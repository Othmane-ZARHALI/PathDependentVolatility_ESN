"""
esn_hyperparam_search.py
======================================

Hyperparameter optimisation pipeline for the ESN A2 stochastic volatility
architecture.


The pipeline simulates the "lnsv" (log-normal stochastic-volatility
fBM) and "bergomi" (rough Bergomi) data-generating processes using the RAW
DAY INDEX t = 1, 2, ..., T as the "time" argument of the fractional Brownian
motion B^H_t, both in the covariance matrix

    Cov(B^H_t, B^H_s) = 1/2 * ( t^(2H) + s^(2H) - |t-s|^(2H) )

and in the variance-matching (Ito) correction of the log-volatility

    log sigma_t = log sigma_0 + lam * B^H_t - 1/2 * lam^2 * t^(2H) .

For H = 0.10 (the "rough" scenarios) this is numerically mild because
t^(2H) grows very slowly (1200^0.2 ~= 4). For H = 0.40 ("fBM_persist")
the exponent 2H = 0.8 is close to 1 (near-Brownian scaling) so t^(2H)
grows almost linearly and reaches 1200^0.8 ~= 291 by the end of a
1200-day path. The subtracted correction term 0.5*lam^2*t^(2H) then
dwarfs the stochastic term (of order lam*t^H ~= lam*1200^0.4 ~= 26),
so log sigma_t drifts to roughly -300 by day 1200 : sigma_t underflows
to ~1e-130 and the simulated path becomes numerically degenerate
(deterministic, non-stochastic, effectively zero). Every downstream
statistic that depends on the variance of returns/vol (H_hat, ACFs,
leverage, Zumbach, kurtosis, ...) is then NaN or garbage, so the score
reference `ref` built for the "fBM_persist" scenario is meaningless.
The CMA-ES/Nelder-Mead calibration therefore has no real signal to
climb for that scenario and its per-scenario 'H' hyperparameter search
collapses to (or near) the lower bound of its admissible domain (0.02)
-- an optimisation artefact, not a genuine best fit.

"""

import math, warnings, time, json, os, copy
import numpy as np
from scipy import optimize, linalg as sp_linalg
warnings.filterwarnings("ignore")

# ============================================================
# 0.  Global constants
# ============================================================

TRADING_DAYS = 252.0
TV_ANN       = 0.20
TV_DAY       = TV_ANN / math.sqrt(TRADING_DAYS)
SIG_MIN      = 0.01  / math.sqrt(TRADING_DAYS)
MU_S         = 0.0

# lam / eta are calibrated on an ANNUALISED scale (consistent with TV_ANN).
# The fBM covariance / Ito correction must use t measured in YEARS, not in
# raw day-count, or the variance-matching correction 0.5*lam^2*t^(2H)
# blows up super-linearly whenever 2H is not small (see module docstring).
TIME_DT      = 1.0 / TRADING_DAYS
# -------------------------------------------------------------------------

# ── DGP parameters (fixed — define the target stylized facts) ──
DGP_SCENARIOS = [
    # (name,  type,          params-dict)
    ("fBM_rough",   "lnsv",     dict(H=0.10, lam=1.0)),
    ("MRW_taylor",  "mrw",      dict(lam2=0.03, T_int=252)),
    ("rBergomi_Z",  "bergomi",  dict(H=0.10, eta=1.0, rho=-0.70)),
    ("fBM_persist", "lnsv",     dict(H=0.40, lam=1.5)),
]

# ── Discrete grid for Nr, Nz ──
NR_GRID = [32, 48, 64, 96]
NZ_GRID = [8, 12, 16, 20]

# ── Inner calibration budget ──
N_DGP_PATHS  = 12    # DGP paths to build score ref
N_CAL_PATHS  = 6     # paths per inner objective evaluation
T_SIM        = 1200  # simulation length (days)
BURN         = 300   # burn-in
T_CAL        = 800   # calibration T (inner loop)
BURN_CAL     = 200   # calibration burn-in

# ── Outer CMA-ES budget ──
MAX_OUTER_EVALS = 300   # total hyperparameter configurations evaluated
POPSIZE         = 12    # CMA-ES population per generation

# ── Inner Nelder-Mead budget (exposed so `fast`/demo runs can shrink it) ──
INNER_NM_MAXITER = 250

# ── Fixed random seed for hyperparameter matrices ──
MATRIX_SEED = 202695547565

# ============================================================
# 1.  Utility — softplus and inverse
# ============================================================

def _sp(x):
    if x >  35: return x
    if x < -35: return math.exp(x)
    return math.log1p(math.exp(x))

def _inv_sp(y):
    y = max(float(y), 1e-15)
    return y if y > 35 else math.log(math.expm1(y))

def _corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() < 1e-14 or b.std() < 1e-14: return np.nan
    return float(np.corrcoef(a, b)[0, 1])

# ============================================================
# 2.  Score function — full 11-term S  (unchanged)
# ============================================================

def compute_statistics(daily_x, daily_var):
    """11 stylized-fact statistics — identical across all pipeline files."""
    x  = np.asarray(daily_x, float)
    v  = np.asarray(daily_var, float)
    sd = np.sqrt(np.maximum(v, 1e-30))
    lv = np.log(v + 1e-30)
    lg = np.arange(1, 21)
    ret_acf = np.array([_corr(x[:-L], x[L:]) for L in lg])
    abs_acf = np.array([_corr(np.abs(x[:-L]), np.abs(x[L:])) for L in lg])
    sq_acf  = np.array([_corr(x[:-L]**2, x[L:]**2) for L in lg])
    vol_acf = np.array([_corr(v[:-L], v[L:]) for L in lg])
    lev     = np.array([_corr(x[:-L], v[L:]) for L in lg])
    xs, ys  = [], []
    for L in [1, 2, 4, 8, 16, 32, 64]:
        if L < len(lv) // 4:
            d = lv[L:] - lv[:-L]; vv = float(np.mean(d * d))
            if vv > 1e-30 and np.isfinite(vv):
                xs.append(math.log(L)); ys.append(math.log(vv))
    H = 0.5 * float(np.polyfit(xs, ys, 1)[0]) if len(xs) >= 3 else np.nan
    zvals = []
    for L in [5, 10, 20]:
        n = len(x) - 2 * L
        if n > 30:
            csx = np.concatenate([[0.], np.cumsum(x)])
            csv = np.concatenate([[0.], np.cumsum(v)])
            idx = L + np.arange(n)
            pR  = csx[idx]-csx[idx-L]; fR = csx[idx+L]-csx[idx]
            pV  = (csv[idx]-csv[idx-L])/L; fV = (csv[idx+L]-csv[idx])/L
            zvals.append(_corr(pR**2, fV) - _corr(pV, fR**2))
    c    = x - x.mean(); var = float(np.var(c))
    kurt = float(np.mean(c**4) / (var**2 + 1e-30) - 3.)
    ann  = math.sqrt(TRADING_DAYS)
    return dict(
        H_hat        = float(H),
        mean_vol_ann = float(sd.mean()) * ann,
        q995_vol_ann = float(np.quantile(sd, 0.995)) * ann,
        max_vol_ann  = float(sd.max()) * ann,
        mean_vol_acf = float(np.nanmean(vol_acf)),
        taylor_gap   = float(np.nanmean(abs_acf) - np.nanmean(sq_acf)),
        taylor_frac  = float(np.nanmean(abs_acf > sq_acf)),
        zumbach      = float(np.nanmean(zvals)),
        leverage     = float(np.nanmean(lev)),
        kurtosis     = kurt,
        max_ret_acf  = float(np.nanmax(np.abs(ret_acf))),
    )

def make_score_ref(dgp_stats):
    _m  = lambda k: float(np.nanmean([s[k] for s in dgp_stats]))
    _sd = lambda k: float(np.nanstd( [s[k] for s in dgp_stats]))
    floors = dict(H_hat=0.01, mean_vol_ann=0.01, q995_vol_ann=0.05,
                  max_vol_ann=0.10, mean_vol_acf=0.005, taylor_gap=0.002,
                  taylor_frac=0.05, zumbach=0.005, leverage=0.005,
                  kurtosis=0.20, max_ret_acf=0.005)
    ref = {}
    for k, fl in floors.items():
        c = _m(k); s = _sd(k)
        ref[k + "_c"] = c
        ref[k + "_s"] = max(s, abs(c) * 0.30, fl)
    ref["stress_mx"] = max(ref["max_vol_ann_c"] * 2.0, 1.5)
    return ref

def score_fn_tent(st, ref):
    """
    ORIGINAL ("tent") score. S >= 0, MINIMISED.  S = sum_k w_k*(1-f_k) + 5*Stress.
      f_k in [0,1] (1 = perfect match, 0 = worst) for two-sided terms,
      g_k in [0,1] similarly for the one-sided (upper-tail) terms.
    S = 0  <=>  a perfect match on every one of the 11 stylised facts and
    no stress event. S = sum_k w_k + 5 = 15.3 is the worst possible value
    (every term maximally wrong AND a stress event fires).

    *** KNOWN LIMITATION (motivates score_fn_smooth below) ***
    f/g are PIECEWISE LINEAR "tents": f(x,c,s)=max(0,1-|x-c|/s). Once a
    term is more than one tolerance-width s away from its target, f=0
    EXACTLY (not just small) -- the local gradient of that term is
    identically zero from that point outward. A local optimiser
    (Nelder-Mead) gets literally no signal to keep improving a term once
    it has "given up" on it, which is a major reason calibration runs
    plateau well above S=0 rather than continuing to descend.
    """
    def f(x, c, s): return max(0., 1. - abs(x - c) / max(s, 1e-8))
    def g(x, c, s): return max(0., 1. - max(0., x - c) / max(s, 1e-8))
    stress = int(st["max_vol_ann"] > ref["stress_mx"]
                 or st["mean_vol_ann"] < 0.05
                 or st["mean_vol_ann"] > 1.50)
    d  = 2.2 * (1. - f(st["H_hat"],        ref["H_hat_c"],       ref["H_hat_s"]))
    d += 1.4 * (1. - f(st["mean_vol_ann"],  ref["mean_vol_ann_c"],ref["mean_vol_ann_s"]))
    d += 0.6 * (1. - g(st["q995_vol_ann"],  ref["q995_vol_ann_c"],ref["q995_vol_ann_s"]))
    d += 0.2 * (1. - g(st["max_vol_ann"],   ref["max_vol_ann_c"], ref["max_vol_ann_s"]))
    d += 0.8 * (1. - f(st["mean_vol_acf"],  ref["mean_vol_acf_c"],ref["mean_vol_acf_s"]))
    d += 1.0 * (1. - f(st["taylor_gap"],    ref["taylor_gap_c"],  ref["taylor_gap_s"]))
    d += 0.8 * (1. - f(st["taylor_frac"],   ref["taylor_frac_c"], ref["taylor_frac_s"]))
    d += 1.0 * (1. - f(st["zumbach"],       ref["zumbach_c"],     ref["zumbach_s"]))
    d += 1.1 * (1. - f(st["leverage"],      ref["leverage_c"],    ref["leverage_s"]))
    d += 0.7 * (1. - f(st["kurtosis"],      ref["kurtosis_c"],    ref["kurtosis_s"]))
    d += 0.5 * (1. - f(st["max_ret_acf"],   ref["max_ret_acf_c"], ref["max_ret_acf_s"]))
    d += 5.0 * stress
    return float(d)


def score_fn_smooth(st, ref):
    """
    "SMOOTHIE" score (new default): identical weights, centres, tolerances
    and Stress definition as score_fn_tent, but the piecewise-linear tent
    f/g are replaced by SMOOTH GAUSSIAN-KERNEL proximity functions:

        f_smooth(x,c,s) = exp( -0.5*((x-c)/s)^2 )                in (0,1]
        g_smooth(x,c,s) = exp( -0.5*(max(0,x-c)/s)^2 )            in (0,1]

    Properties (by design):
      - f_smooth(c,c,s) = 1 exactly (perfect match), and the derivative
        there is exactly 0 -- a genuine smooth minimum of the deviation
        1-f_smooth, not the tent's V-shaped kink.
      - f_smooth NEVER hits exactly 0: there is always a (small but
        non-zero) gradient pulling a badly-off term back toward its
        target, at any distance. This removes the tent's flat,
        zero-gradient plateau beyond one tolerance width, which is the
        main reason Nelder-Mead calibration under score_fn_tent stalls
        well above S=0.
      - Same bounds: d_k = w_k*(1-f_smooth) in [0, w_k), so S is still
        in [0, 15.3) with S=0 <=> a perfect match and no stress event.
    The Stress indicator is left as a hard trigger (unchanged): it flags
    a rare, qualitatively different failure mode (numerical blow-up /
    collapse) that should remain a sharp, deterrent penalty rather than
    a smoothly-traded-off one.
    """
    def f(x, c, s): return math.exp(-0.5 * ((x - c) / max(s, 1e-8)) ** 2)
    def g(x, c, s): return math.exp(-0.5 * (max(0., x - c) / max(s, 1e-8)) ** 2)
    stress = int(st["max_vol_ann"] > ref["stress_mx"]
                 or st["mean_vol_ann"] < 0.05
                 or st["mean_vol_ann"] > 1.50)
    d  = 2.2 * (1. - f(st["H_hat"],        ref["H_hat_c"],       ref["H_hat_s"]))
    d += 1.4 * (1. - f(st["mean_vol_ann"],  ref["mean_vol_ann_c"],ref["mean_vol_ann_s"]))
    d += 0.6 * (1. - g(st["q995_vol_ann"],  ref["q995_vol_ann_c"],ref["q995_vol_ann_s"]))
    d += 0.2 * (1. - g(st["max_vol_ann"],   ref["max_vol_ann_c"], ref["max_vol_ann_s"]))
    d += 0.8 * (1. - f(st["mean_vol_acf"],  ref["mean_vol_acf_c"],ref["mean_vol_acf_s"]))
    d += 1.0 * (1. - f(st["taylor_gap"],    ref["taylor_gap_c"],  ref["taylor_gap_s"]))
    d += 0.8 * (1. - f(st["taylor_frac"],   ref["taylor_frac_c"], ref["taylor_frac_s"]))
    d += 1.0 * (1. - f(st["zumbach"],       ref["zumbach_c"],     ref["zumbach_s"]))
    d += 1.1 * (1. - f(st["leverage"],      ref["leverage_c"],    ref["leverage_s"]))
    d += 0.7 * (1. - f(st["kurtosis"],      ref["kurtosis_c"],    ref["kurtosis_s"]))
    d += 0.5 * (1. - f(st["max_ret_acf"],   ref["max_ret_acf_c"], ref["max_ret_acf_s"]))
    d += 5.0 * stress
    return float(d)


# score_fn is the function used everywhere downstream (inner/outer
# calibration, plotting); switched to the smooth version by default.
score_fn = score_fn_smooth

SCORE_MAX = 2.2+1.4+0.6+0.2+0.8+1.0+0.8+1.0+1.1+0.7+0.5 + 5.0   # = 15.3, worst possible S
SCORE_TERM_WEIGHTS = [
    ("H_hat", 2.2), ("mean_vol_ann", 1.4), ("q995_vol_ann", 0.6),
    ("max_vol_ann", 0.2), ("mean_vol_acf", 0.8), ("taylor_gap", 1.0),
    ("taylor_frac", 0.8), ("zumbach", 1.0), ("leverage", 1.1),
    ("kurtosis", 0.7), ("max_ret_acf", 0.5),
]

def score_breakdown_tent(st, ref):
    """Per-term deviation d_k under the original tent kernel (diagnostics)."""
    def f(x, c, s): return max(0., 1. - abs(x - c) / max(s, 1e-8))
    def g(x, c, s): return max(0., 1. - max(0., x - c) / max(s, 1e-8))
    one_sided = {"q995_vol_ann", "max_vol_ann"}
    out = {}
    for k, w in SCORE_TERM_WEIGHTS:
        fn = g if k in one_sided else f
        out[k] = w * (1. - fn(st[k], ref[k + "_c"], ref[k + "_s"]))
    return out

def score_breakdown_smooth(st, ref):
    """Per-term deviation d_k under the smooth Gaussian-kernel (diagnostics)."""
    def f(x, c, s): return math.exp(-0.5 * ((x - c) / max(s, 1e-8)) ** 2)
    def g(x, c, s): return math.exp(-0.5 * (max(0., x - c) / max(s, 1e-8)) ** 2)
    one_sided = {"q995_vol_ann", "max_vol_ann"}
    out = {}
    for k, w in SCORE_TERM_WEIGHTS:
        fn = g if k in one_sided else f
        out[k] = w * (1. - fn(st[k], ref[k + "_c"], ref[k + "_s"]))
    return out

# score_breakdown is the function used everywhere downstream (plotting);
# switched to the smooth version by default, consistent with score_fn.
score_breakdown = score_breakdown_smooth

# ============================================================
# 3.  DGP generators   *** FIXED: annualised time argument ***
# ============================================================

_CHOL_CACHE_FBM = {}
_CHOL_CACHE_MRW = {}

def _get_fbm_chol(T, H, dt=TIME_DT):
    """
    Cholesky factor of the fBM covariance matrix, with t measured in the
    SAME units as the annualised drift/vol-of-vol parameters (years by
    default: dt = 1/TRADING_DAYS). This is the fix: the original code
    called this with the implicit default dt=1.0 (raw day count), which
    is inconsistent with lam/eta being O(1) annualised parameters and
    causes catastrophic non-stationary blow-up of the correction term
    for H not close to 0 (see module docstring).
    """
    key = (T, H, dt)
    if key not in _CHOL_CACHE_FBM:
        t = np.arange(1, T + 1, dtype=float) * dt
        C = 0.5*(t[:,None]**(2*H) + t[None,:]**(2*H)
                 - np.abs(t[:,None]-t[None,:])**(2*H))
        C += np.eye(T) * 1e-10
        try:    L = sp_linalg.cholesky(C, lower=True)
        except: L = np.linalg.cholesky(C + np.eye(T) * 1e-8)
        _CHOL_CACHE_FBM[key] = L
    return _CHOL_CACHE_FBM[key]

def _get_mrw_chol(T, lam2, T_int):
    # MRW's kernel Cov = lam2 * log(T_int/lag) (lag in days, capped at 0
    # beyond T_int) is already stationary and bounded by construction
    # (the "constant" correction term lam2*log(T_int) does not grow with
    # the simulation horizon), so no fix is needed here.
    key = (T, lam2, T_int)
    if key not in _CHOL_CACHE_MRW:
        t = np.arange(T, dtype=float)
        lag = np.abs(t[:,None] - t[None,:])
        with np.errstate(divide='ignore', invalid='ignore'):
            C = np.where(lag == 0,
                         lam2 * np.log(max(T_int, 1.)),
                         np.maximum(0., lam2 * np.log(T_int / np.maximum(lag, 1.))))
        C = np.where(lag >= T_int, 0., C)
        eigs, V = np.linalg.eigh(C)
        C_psd = V @ np.diag(np.maximum(eigs, 1e-10)) @ V.T
        L = np.linalg.cholesky(C_psd + np.eye(T) * 1e-10)
        _CHOL_CACHE_MRW[key] = L
    return _CHOL_CACHE_MRW[key]

def sim_dgp(scenario_type, params, n_paths, T, burn, seed_base=0):
    """
    Simulate one DGP scenario.  Returns list of stat dicts.

    scenario_type in {'lnsv', 'mrw', 'bergomi'}
    """
    rng_base = np.random.default_rng(seed_base)
    st_all   = []

    if scenario_type == "lnsv":
        H, lam = params["H"], params["lam"]
        L = _get_fbm_chol(T, H, dt=TIME_DT)
        t_   = np.arange(1, T + 1, dtype=float) * TIME_DT      # *** FIX ***
        for p in range(n_paths):
            rng  = np.random.default_rng(int(rng_base.integers(1 << 31)))
            BH   = L @ rng.standard_normal(T)
            lsig = math.log(TV_DAY) + lam * BH - 0.5 * lam**2 * t_**(2*H)
            sig  = np.exp(lsig); v = sig**2
            x    = sig * rng.standard_normal(T)
            st   = compute_statistics(x[burn:], v[burn:])
            st_all.append(st)

    elif scenario_type == "mrw":
        lam2, T_int = params["lam2"], params["T_int"]
        var_omega   = lam2 * math.log(max(T_int, 1.))
        L = _get_mrw_chol(T, lam2, T_int)
        for p in range(n_paths):
            rng   = np.random.default_rng(int(rng_base.integers(1 << 31)))
            omega = L @ rng.standard_normal(T)
            omega -= var_omega / 2.0
            sig   = TV_DAY * np.exp(omega); v = sig**2
            x     = sig * rng.standard_normal(T)
            st    = compute_statistics(x[burn:], v[burn:])
            st_all.append(st)

    elif scenario_type == "bergomi":
        H, eta, rho = params["H"], params["eta"], params["rho"]
        xi0         = TV_DAY**2
        L           = _get_fbm_chol(T, H, dt=TIME_DT)
        t_          = np.arange(1, T + 1, dtype=float) * TIME_DT   # *** FIX ***
        sqrt1m      = math.sqrt(max(1. - rho**2, 0.))
        for p in range(n_paths):
            rng  = np.random.default_rng(int(rng_base.integers(1 << 31)))
            eps1 = rng.standard_normal(T)
            eps2 = rng.standard_normal(T)
            dW2  = rho * eps1 + sqrt1m * eps2
            WH   = L @ dW2
            lV   = math.log(xi0) + eta * WH - 0.5 * eta**2 * t_**(2*H)
            V    = np.exp(lV)
            X    = -0.5 * V + np.sqrt(V) * eps1   # dt=1 (daily return step)
            st   = compute_statistics(X[burn:], V[burn:])
            st_all.append(st)

    else:
        raise ValueError(f"Unknown DGP type: {scenario_type}")

    return st_all

# ============================================================
# 4.  ESN simulator — parameterised by (arch, inner_params)  (unchanged)
# ============================================================

def build_esn_matrices(arch, inner_params):
    """
    Build reservoir matrices given:
      arch         — hyperparameter dict (Nr, Nz, zz_scale, sign_prob_neg, ...;
                     7 shared hyperparameters)
      inner_params — estimated (per-scenario CALIBRATED) parameters dict
                     (H, lam_lo, lam_hi, az_lo, az_hi, rough_scale,
                     zr_lo, zr_hi, m1, m2)

    STRUCTURAL CHANGE (this draft): the quadratic-in-r term of the
    previous draft, kappa_quad/N_r * ||r_t||^2 (i.e. Q = (kappa_quad/N_r)*I,
    a single calibrated scalar), is GENERALISED to
        Q = (1/N_r) * M,   M = m1*I_{N_r} + m2*q q^T,
    so that
        r_t^T Q r_t = (1/N_r) * ( m1*||r_t||^2 + m2*(q^T r_t)^2 ),
    with BOTH m1 and m2 calibrated per DGP scenario (replacing the
    single kappa_quad; net +1 calibrated dimension). M is SYMMETRIC by
    construction for any real (m1,m2) -- I and qq^T are each symmetric,
    and a real linear combination of symmetric matrices is symmetric --
    with NO positive-definiteness constraint imposed (m1, m2 may each be
    positive or negative; M's eigenvalues are m1, with multiplicity
    N_r-1 on q's orthogonal complement, and m1+m2*||q||^2, with
    multiplicity 1 along q, and either or both may be negative). This is
    a genuine matrix (not restricted to a multiple of the identity): the
    m1*I term is the previous draft's isotropic reservoir-energy
    channel, and the NEW m2*qq^T term is a second, independent channel
    that responds specifically to the SQUARED rough factor (q^T r_t)^2
    -- the same projection direction that already drives the model's
    linear rough term (rough_scale*q^T r_t), but now allowed to also
    enter quadratically, with its own independently-calibrated
    coefficient.

    Tractability note (why this parameterisation, not a fully free M):
    a fully free N_r x N_r symmetric M would add N_r(N_r+1)/2 = 2080
    calibrated dimensions at N_r=64 -- far beyond what a per-scenario
    Nelder-Mead multi-start can navigate, and likely to overfit the
    11-term score against a handful of DGP scenarios. M = m1*I + m2*qq^T
    is a genuine, non-trivial symmetric generalisation of "a multiple of
    the identity" (the previous draft's Q) -- diagonal in the eigenbasis
    {q, q's orthogonal complement} -- while adding only one new
    calibrated scalar. See the protocol's discussion of further,
    still-more-general choices (e.g. a small basis of q's at several
    fixed roughness exponents, each with its own calibrated coefficient)
    if this 2-term form proves limiting.

    The z-bank generalisation of the previous draft (per-mode profile
    z_readout_j = geomspace(zr_lo, zr_hi, N_z)[j], calibrated via
    zr_lo, zr_hi) is UNCHANGED here.

    NOTE: this "H" is the ESN's OWN internal memory-kernel exponent
    (shapes q = lambdas^(0.5-H), the reservoir's Volterra weighting),
    calibrated per-scenario. It is a different object from the DGP's
    fBM roughness parameter; the calibration loop searches over it so
    that the resulting SIMULATED H_hat (measured the same way as for
    the DGP, via the log-vol structure function) matches the DGP's
    target H_hat_c. It is not fixed/bugged; the bug was entirely in the
    DGP reference statistics used as the calibration target.
    """
    n_r  = int(arch["Nr"])
    n_z  = int(arch["Nz"])
    H    = float(inner_params["H"])
    lam_lo = float(inner_params.get("lam_lo", 1/3500))
    lam_hi = float(inner_params.get("lam_hi", 2.0))
    az_lo  = float(inner_params.get("az_lo",  1/280))
    az_hi  = float(inner_params.get("az_hi",  1/7))
    rough_scale = float(inner_params.get("rough_scale", 0.40))
    zr_lo  = float(inner_params.get("zr_lo",  0.03))
    zr_hi  = float(inner_params.get("zr_hi",  0.07))
    # m1, m2 are the two calibrated coefficients of the SYMMETRIC matrix
    # M = m1*I + m2*qq^T (no positive-definiteness constraint; each may
    # be positive or negative).
    m1 = float(inner_params.get("m1", 0.0))
    m2 = float(inner_params.get("m2", 0.0))

    # ── Rough bank ────────────────────────────────────────────
    lam_lo = max(lam_lo, 1e-5); lam_hi = max(lam_hi, lam_lo * 2)
    az_lo  = max(az_lo,  1e-5); az_hi  = max(az_hi,  az_lo  * 2)
    lambdas = np.geomspace(lam_lo, lam_hi, n_r)
    q = lambdas ** (0.5 - H)
    q = q / (np.linalg.norm(q) + 1e-15)
    b = np.sqrt(2. * lambdas)
    C = (b[:, None] * b[None, :]) / (lambdas[:, None] + lambdas[None, :])
    q = q / math.sqrt(float(q @ C @ q) + 1e-15)

    # ── z-bank rates and per-mode readout profile ────────────
    az = np.geomspace(az_lo, az_hi, n_z)
    zr_lo_c = max(zr_lo, 1e-4); zr_hi_c = max(zr_hi, zr_lo_c * 1.05)
    z_readout_profile = np.geomspace(zr_lo_c, zr_hi_c, n_z)  # per-mode w_{j,z} source

    # ── Fixed random matrices (drawn once from MATRIX_SEED) ──
    rng    = np.random.default_rng(MATRIX_SEED)
    zz     = rng.uniform(-arch["zz_scale"], arch["zz_scale"], n_z)
    sz     = -rng.choice([-1., 1.], n_z,
                          p=[arch["sign_prob_neg"],
                             1. - arch["sign_prob_neg"]])
    fi = np.array([min(n_r-1, n_r//2+int((n_r//2-1)*j/max(n_z-1,1)))
                   for j in range(n_z)], dtype=int)
    si = np.array([int((n_r//2-1)*(n_z-1-j)/max(n_z-1,1))
                   for j in range(n_z)], dtype=int)

    # ── Bias b0 ───────────────────────────────────────────────
    b0     = _inv_sp(math.sqrt(max(TV_DAY**2 - SIG_MIN**2, 1e-15)))
    kappa0 = (arch["rough_orientation"] * rough_scale
              * float(q @ np.sqrt(2. * lambdas)))

    return dict(lambdas=lambdas, q=q, az=az, zz=zz, sz=sz,
                fi=fi, si=si, b0=b0, kappa0=kappa0, n_r=n_r, n_z=n_z,
                rough_scale=rough_scale,
                z_readout_profile=z_readout_profile,
                zr_lo=zr_lo_c, zr_hi=zr_hi_c,
                m1=m1, m2=m2)


def _sim_esn_with_params(seed, T, arch, matrices, b0_delta=0., scale=1., dt=1.0):
    """
    Simulate ESN given pre-built matrices and calibrated offsets.
    This is the inner-loop simulator used during calibration.
    rough_scale, the per-mode z_readout_profile, and kappa_quad are all
    read from `matrices` (calibrated per scenario by build_esn_matrices
    from inner_params), not from `arch`.
    """
    n_r  = matrices["n_r"]; n_z = matrices["n_z"]
    rng  = np.random.default_rng(int(seed))
    spd  = int(round(1. / dt)); n_st = T * spd; sdt = math.sqrt(dt)
    lam  = matrices["lambdas"]
    al   = np.exp(-lam * dt); cl = np.sqrt(np.maximum(1. - al**2, 1e-14))
    az   = matrices["az"]
    azd  = np.exp(-az * dt); om = 1. - azd
    r    = np.zeros(n_r); z = np.zeros(n_z)
    dx   = np.zeros(T);   dv = np.zeros(T)
    b0   = matrices["b0"] + b0_delta
    # per-mode readout weights w_{j,z} = z_readout_j / sqrt(N_z)
    wz_vec = matrices["z_readout_profile"] / math.sqrt(n_z)
    rc   = arch["rough_orientation"] * matrices["rough_scale"]
    m1   = matrices["m1"] / n_r; m2 = matrices["m2"] / n_r  # Q=(1/N_r)*(m1*I+m2*qq^T)
    qvec = matrices["q"]
    fi   = matrices["fi"]; si = matrices["si"]
    zz   = matrices["zz"]; sz = matrices["sz"]

    for step in range(n_st):
        eps = rng.normal()
        qr  = float(qvec @ r)
        quad_term = m1 * float(r @ r) + m2 * qr * qr   # r_t^T Q r_t, Q=(1/N_r)(m1*I+m2*qq^T)
        eta = (b0 + rc * qr
               + float(wz_vec @ z) + quad_term)
        sp  = _sp(eta); sig = math.sqrt(SIG_MIN**2 + sp**2); var = sig * sig
        day = step // spd
        dx[day] += (MU_S - 0.5 * var) * dt + sig * sdt * eps
        dv[day] += var * dt
        r  = al * r + cl * eps
        m  = max(0., 1. - arch["gamma_norm"] * np.linalg.norm(r) / math.sqrt(n_r))
        zo = z.copy()
        for j in range(n_z):
            p1   = r[fi[j]]; p2 = r[si[j]]
            ev   = 0.7*p1*p1 + 0.3*p2*p2 - 1.; lin = 0.7*p1 + 0.3*p2
            jm   = max(0, j-1); jp = min(n_z-1, j+1)
            lc   = 0.5 * arch["local_z_strength"] * (zo[jm] + zo[jp])
            u    = (sz[j]*p1*p2 + arch["even_strength"]*ev
                    - arch["linear_strength"]*lin + zz[j]*zo[j] + lc)
            z[j] = azd[j]*zo[j] + om[j]*m*arch["z_strength"]*math.tanh(u)

    return dx * scale, dv

# ============================================================
# 5.  Inner calibration
#     *** FIX 2: multi-start Nelder-Mead ***
#     A single fixed starting point x0 = (H=0.10, "fast" reservoir) biases
#     the search toward rough (small-H, fast-mode) solutions. For a
#     genuinely persistent target (fBM_persist, H_hat_c ~ 0.40) the true
#     optimum sits in a completely different, "slow reservoir" region of
#     parameter space (small lam_hi = no fast/day-scale reservoir modes,
#     see protocol Sec. 5.2) that a single Nelder-Mead run initialised in
#     the "rough" region essentially never finds -- it instead reports
#     whatever its local basin's boundary happens to be, which is exactly
#     the symptom originally observed (H pinned at 0.02). Running a small
#     number of geometrically different starts and keeping the best
#     removes this local-optimum/boundary-lock artefact for every
#     scenario, not just fBM_persist.
# ============================================================

INNER_NM_STARTS = [
    # (H0, lam_lo0, lam_hi0, az_lo0, az_hi0, rough_scale0, zr_lo0, zr_hi0, m1_0, m2_0)
    # m1, m2 are the two calibrated coefficients of the SYMMETRIC matrix
    # M = m1*I + m2*qq^T (no positive-definiteness constraint -- each
    # may be positive, negative, or zero).
    # -- "rough / fast reservoir"
    (0.10, 1/3500, 2.0,   1/280, 1/7,  0.40, 0.03, 0.07, 0.0, 0.0),
    # "persistent / slow reservoir": no fast day-scale modes at all
    (0.40, 1e-4,   0.08,  1/280, 1/7,  0.40, 0.03, 0.07, 0.0, 0.0),
    # intermediate
    (0.25, 1e-3,   0.5,   1/280, 1/7,  0.40, 0.03, 0.07, 0.0, 0.0),
    # *** wide-band / multiscale reservoir ***
    # Spans fast (day-scale) AND slow (multi-year) OU modes simultaneously
    # within the SAME reservoir, so a single architecture can express both
    # a rough short-lag component and a persistent long-lag component (the
    # structure-function H_hat estimator only sees lags 1-64 days, but the
    # underlying kernel is free to keep slower modes too) -- this is what
    # lets one calibrated reservoir be flexible enough to fit BOTH the
    # rough (H~0.10) and persistent (H~0.40) DGP families, not just
    # whichever extreme the start happens to land near.
    (0.20, 1e-5,   5.0,   1/400, 1/7,  0.40, 0.03, 0.07, 0.0, 0.0),
    # long-memory extreme (very slow upper rate too): covers targets even
    # more persistent than fBM_persist, for robustness.
    (0.35, 1e-5,   0.02,  1/500, 1/10, 0.40, 0.03, 0.07, 0.0, 0.0),
    # strong amplitude (rough_scale, z-readout profile)
    (0.15, 1e-4,   1.0,   1/300, 1/8,  0.90, 0.10, 0.30, 0.0, 0.0),
    # weak amplitude
    (0.20, 1e-4,   1.0,   1/300, 1/8,  0.15, 0.02, 0.03, 0.0, 0.0),
    # *** ISOTROPIC quadratic term, positive (M~+m1*I) ***
    # probes a regime where slow/cross-scale z-modes (large j, per
    # Eq. 9-10's fast/slow tap separation) matter MORE than fast ones,
    # combined with a genuine, positive, ISOTROPIC energy-driven
    # volatility boost -- the special case m2=0 recovers the previous
    # (scalar kappa_quad) draft's model exactly.
    (0.20, 1e-4,   1.0,   1/300, 1/8,  0.40, 0.02, 0.20, 0.8, 0.0),
    # *** ISOTROPIC quadratic term, negative ***
    (0.20, 1e-4,   1.0,   1/300, 1/8,  0.40, 0.20, 0.02, -0.8, 0.0),
    # *** q-ALIGNED quadratic term (m2), isotropic term off ***
    # probes the genuinely NEW channel this draft adds: a quadratic
    # response specifically to the squared rough factor (q^T r_t)^2,
    # independent of the isotropic reservoir-energy channel (m1=0).
    (0.20, 1e-4,   1.0,   1/300, 1/8,  0.40, 0.03, 0.07, 0.0, 0.8),
    (0.20, 1e-4,   1.0,   1/300, 1/8,  0.40, 0.03, 0.07, 0.0, -0.8),
]

def inner_calibrate(arch, ref, seed_offset=0):
    """
    Given arch (7 shared hyperparameters + Nr, Nz) and ref (DGP score
    reference), calibrate:
      p = [H, log(lam_lo), log(lam_hi), log(az_lo), log(az_hi),
           b0_delta, log(scale), log(rough_scale),
           log(zr_lo), log(zr_hi), m1, m2]

    STRUCTURAL CHANGE (this draft): the quadratic-in-r term is
    generalised from a single scalar kappa_quad (Q=kappa_quad/N_r * I)
    to a genuine 2-term MATRIX Q=(1/N_r)*M, M=m1*I+m2*qq^T -- i.e. M is
    SYMMETRIC by construction for any real (m1,m2) (no
    positive-definiteness constraint is imposed; see
    build_esn_matrices), with BOTH m1 and m2 calibrated. The calibrated
    inner-loop vector grows from 11 to 12 dimensions. Bounds:
    ROUGH_SCALE_BOUNDS, ZR_LO_BOUNDS, ZR_HI_BOUNDS, M1_BOUNDS,
    M2_BOUNDS (both symmetric ranges, m1/m2 may be positive or
    negative).

    Objective: MINIMISE mean score S (>=0, 0=perfect) over N_CAL_PATHS
    paths of T_CAL days, via multi-start Nelder-Mead (see module note
    above). Optimiser: Nelder-Mead (no gradient needed, 12-dimensional).

    Returns (best_inner_params, best_score)  [best_score = lowest S found].
    """
    rs_lo, rs_hi   = ROUGH_SCALE_BOUNDS
    zrl_lo, zrl_hi = ZR_LO_BOUNDS
    zrh_lo, zrh_hi = ZR_HI_BOUNDS
    m1_lo, m1_hi   = M1_BOUNDS
    m2_lo, m2_hi   = M2_BOUNDS

    def decode(p):
        H      = float(np.clip(p[0], 0.02, 0.49))
        # widened lower bound (was -9): allow much longer-memory slow modes
        # (down to lam_lo ~ 1/(40*252) ~ a 40-year mean-reversion time),
        # needed so the SAME reservoir can hold a persistent component
        # alongside a rough one (wide-band/multiscale start, see
        # INNER_NM_STARTS).
        lam_lo = float(np.exp(np.clip(p[1], -11, 1)))
        # widened lower bound (was -3): allow a fully "slow" reservoir
        # with no fast/day-scale modes, needed to reach persistent targets;
        # upper bound slightly raised (was 2) to keep headroom for rough
        # targets even when Nr is large (more modes -> denser geomspace).
        lam_hi = float(np.exp(np.clip(p[2], -6.0, 2.3)))
        az_lo  = float(np.exp(np.clip(p[3], -9, 0)))
        az_hi  = float(np.exp(np.clip(p[4], -5, 0)))
        b0d    = float(p[5])
        scale  = float(np.exp(np.clip(p[6], -1, 1)))
        rough_scale = float(np.clip(math.exp(p[7]), rs_lo, rs_hi))
        zr_lo  = float(np.clip(math.exp(p[8]), zrl_lo, zrl_hi))
        zr_hi  = float(np.clip(math.exp(p[9]), zrh_lo, zrh_hi))
        m1     = float(np.clip(p[10], m1_lo, m1_hi))
        m2     = float(np.clip(p[11], m2_lo, m2_hi))
        lam_hi = max(lam_hi, lam_lo * 1.5)
        az_hi  = max(az_hi,  az_lo  * 1.5)
        return dict(H=H, lam_lo=lam_lo, lam_hi=lam_hi,
                    az_lo=az_lo, az_hi=az_hi,
                    rough_scale=rough_scale,
                    zr_lo=zr_lo, zr_hi=zr_hi, m1=m1, m2=m2)

    def obj(p):
        try:
            ip  = decode(p)
            b0d = float(p[5])
            sc  = float(np.exp(np.clip(p[6], -1, 1)))
            mat = build_esn_matrices(arch, ip)
            if mat["kappa0"] >= 0: return SCORE_MAX   # admissibility check
            scores = []
            for q in range(N_CAL_PATHS):
                x, v = _sim_esn_with_params(
                    1000 + seed_offset + q, T_CAL, arch, mat,
                    b0_delta=b0d, scale=sc)
                scores.append(score_fn(compute_statistics(x[BURN_CAL:], v[BURN_CAL:]), ref))
            return float(np.nanmean(scores))
        except Exception:
            return SCORE_MAX

    best_res, best_val = None, np.inf
    n_starts = len(INNER_NM_STARTS)
    maxiter_per_start = max(30, INNER_NM_MAXITER // n_starts)
    for si, (H0, ll0, lh0, al0, ah0, rs0, zrl0, zrh0, m10, m20) in enumerate(INNER_NM_STARTS):
        x0 = [H0, math.log(ll0), math.log(lh0), math.log(al0), math.log(ah0),
              0.0, 0.0, math.log(rs0), math.log(zrl0), math.log(zrh0),
              m10, m20]
        res = optimize.minimize(obj, x0, method="Nelder-Mead",
                                options={"maxiter": maxiter_per_start,
                                         "xatol": 0.02, "fatol": 0.02})
        if res.fun < best_val:
            best_val = res.fun; best_res = res

    ip  = decode(best_res.x)
    b0d = float(best_res.x[5])
    sc  = float(np.exp(np.clip(best_res.x[6], -1, 1)))
    ip["b0_delta"] = b0d; ip["scale"] = sc
    return ip, float(best_res.fun)

# ============================================================
# 6.  Objective function — given arch, evaluate across all DGPs
# ============================================================

_dgp_cache = {}   # cache DGP paths to avoid regenerating per arch evaluation

def get_dgp_paths(scenario_name, scenario_type, params):
    """Return cached DGP stat list for a scenario."""
    if scenario_name not in _dgp_cache:
        print(f"    [DGP] Generating {scenario_name} ...", flush=True)
        sts = sim_dgp(scenario_type, params,
                      N_DGP_PATHS, T_SIM, BURN,
                      seed_base=abs(hash(scenario_name)) % (1 << 30))
        _dgp_cache[scenario_name] = sts
    return _dgp_cache[scenario_name]


def evaluate_arch(arch, eval_id=0, verbose=False, return_detail=False):
    """
    Evaluate an architectural hyperparameter configuration across all DGP
    scenarios.

    For each scenario:
      1. Get (or generate) DGP paths and build score reference.
      2. Run inner calibration to find optimal (H, lambda-bounds, az-bounds,
         b0_delta, scale) for this arch.
      3. Evaluate mean score on N_CAL_PATHS final paths.

    Returns: mean S across all scenarios (S >= 0, LOWER = better; 0 = perfect).
    If return_detail=True, also returns a list of per-scenario dicts with
    name/score/inner_params/achieved H_hat (used by the diagnostic plots).
    """
    # Admissibility pre-check: rough_orientation must be -1 for leverage<0
    if arch.get("rough_orientation", -1.) > 0:
        return (SCORE_MAX, []) if return_detail else SCORE_MAX

    scenario_scores = []
    detail = []

    for sc_name, sc_type, sc_params in DGP_SCENARIOS:
        try:
            dgp_sts = get_dgp_paths(sc_name, sc_type, sc_params)
            ref     = make_score_ref(dgp_sts)
            inner_params, inner_score = inner_calibrate(arch, ref, seed_offset=eval_id*100)

            # Final evaluation on fresh paths
            mat = build_esn_matrices(arch, inner_params)
            final_scores = []
            final_hhats  = []
            for q in range(N_CAL_PATHS):
                x, v = _sim_esn_with_params(
                    5000 + eval_id * 100 + q, T_CAL, arch, mat,
                    b0_delta=inner_params.get("b0_delta", 0.),
                    scale=inner_params.get("scale", 1.))
                st = compute_statistics(x[BURN_CAL:], v[BURN_CAL:])
                final_scores.append(score_fn(st, ref))
                final_hhats.append(st["H_hat"])

            sc_score = float(np.nanmean(final_scores))
            scenario_scores.append(sc_score)
            detail.append(dict(name=sc_name, score=sc_score,
                                H_esn=inner_params["H"],
                                H_hat_est=float(np.nanmean(final_hhats)),
                                H_hat_true=ref["H_hat_c"],
                                lam_lo=inner_params["lam_lo"],
                                lam_hi=inner_params["lam_hi"],
                                az_lo=inner_params["az_lo"],
                                az_hi=inner_params["az_hi"],
                                rough_scale=inner_params["rough_scale"],
                                zr_lo=inner_params["zr_lo"],
                                zr_hi=inner_params["zr_hi"],
                                m1=inner_params["m1"], m2=inner_params["m2"],
                                inner_params=inner_params,
                                ref=ref))
            if verbose:
                print(f"      {sc_name:20s}: inner={inner_score:.3f}  "
                      f"final={sc_score:.3f}  "
                      f"H={inner_params['H']:.3f}  "
                      f"lam=[{inner_params['lam_lo']:.2e},{inner_params['lam_hi']:.2f}]",
                      flush=True)
        except Exception as e:
            if verbose: print(f"      {sc_name}: ERROR {e}", flush=True)
            scenario_scores.append(SCORE_MAX)
            detail.append(dict(name=sc_name, score=SCORE_MAX, H_esn=np.nan,
                                H_hat_est=np.nan, H_hat_true=np.nan,
                                lam_lo=np.nan, lam_hi=np.nan,
                                az_lo=np.nan, az_hi=np.nan,
                                rough_scale=np.nan, zr_lo=np.nan,
                                zr_hi=np.nan, m1=np.nan, m2=np.nan,
                                inner_params=None, ref=None))

    mean_s = float(np.nanmean(scenario_scores))
    return (mean_s, detail) if return_detail else mean_s

# ============================================================
# 7.  Hyperparameter space — encoding / decoding  (unchanged)
# ============================================================

HP_SPEC = [
    # name                 default   lo      hi      transform
    # STRUCTURAL CHANGE (this draft): rough_scale and z_readout REMOVED
    # from the outer (shared-architecture) search -- they are now
    # per-scenario CALIBRATED parameters (see inner_calibrate), not
    # hyperparameters shared across all 4 DGPs. Rationale: rough_scale
    # sets the amplitude of the (unit-variance) rough factor q^T r_t in
    # the pre-activation, and z_readout sets the z-bank's overall
    # contribution -- both are AMPLITUDE parameters that directly
    # determine how strongly the ESN can express variance-clustering
    # and higher-order (Taylor/Zumbach) effects. A single amplitude
    # shared across four DGPs with very different target stylised
    # facts (H=0.10 vs. H=0.40, different vol-of-vol/kurtosis regimes)
    # is exactly the kind of one-size-fits-all constraint that showed
    # up as a genuine bottleneck in the companion Student-t protocol:
    # promoting the analogous parameter there (rough_scale) from fixed
    # architecture to calibrated let the model reach a much wider,
    # DGP-appropriate range of stylised-fact values. The remaining 7
    # hyperparameters (all "shape", not "amplitude", parameters of the
    # nonlinear z-bank) stay in the outer search.
    ("z_strength",         0.34,     0.05,   2.0,    "log"),
    ("even_strength",      1.50,     0.10,   5.0,    "log"),
    ("linear_strength",    0.25,     0.01,   2.0,    "log"),
    ("gamma_norm",         1.00,     0.10,   3.0,    "log"),
    ("local_z_strength",   0.03,     0.001,  0.5,    "log"),
    ("zz_scale",           0.08,     0.01,   0.5,    "log"),
    ("sign_prob_neg",      0.22,     0.05,   0.50,   "linear"),
]

# Calibrated (per-scenario) bounds for rough_scale, z_readout -- same
# admissible ranges as the OLD outer-search HP_SPEC rows they replace,
# now applied inside inner_calibrate instead of the CMA-ES outer loop.
ROUGH_SCALE_BOUNDS = (0.05, 2.0)
# Z_READOUT_BOUNDS replaced by a LO/HI profile pair (geomspace across
# the N_z z-modes -- see build_esn_matrices) plus a new quadratic-term
# MATRIX Q=(1/N_r)*M, M SYMMETRIC by construction: M = m1*I + m2*qq^T,
# a real linear combination of two symmetric matrices (I and qq^T) is
# always symmetric, for ANY real m1, m2 -- no positive-definiteness
# constraint is imposed. Bounds for zr_lo/zr_hi match the OLD single
# Z_READOUT_BOUNDS range; M1_BOUNDS/M2_BOUNDS are symmetric O(1) ranges
# (same range as the old signed KAPPA_QUAD_BOUNDS), allowing each
# coefficient to be positive, negative, or (near) zero.
ZR_LO_BOUNDS   = (0.01, 0.5)
ZR_HI_BOUNDS   = (0.01, 0.5)
M1_BOUNDS = (-2.0, 2.0)
M2_BOUNDS = (-2.0, 2.0)

def encode_hp(arch):
    x = []
    for name, default, lo, hi, tr in HP_SPEC:
        v = arch.get(name, default)
        if tr == "log":
            x.append(math.log(max(v, 1e-9)))
        else:
            x.append(float(v))
    return np.array(x)

def decode_hp(x, Nr, Nz):
    arch = {"rough_orientation": -1.0, "Nr": Nr, "Nz": Nz,
            "matrix_seed": MATRIX_SEED}
    for i, (name, default, lo, hi, tr) in enumerate(HP_SPEC):
        if tr == "log":
            v = math.exp(x[i])
        else:
            v = float(x[i])
        arch[name] = float(np.clip(v, lo, hi))
    return arch

def default_arch(Nr=64, Nz=12):
    """Return the A2-981003 default hyperparameters (7 shared + Nr, Nz;
    rough_scale/z_readout are no longer part of arch -- see HP_SPEC)."""
    return {
        "rough_orientation": -1.0, "Nr": Nr, "Nz": Nz,
        "matrix_seed": MATRIX_SEED,
        "z_strength": 0.34,
        "even_strength": 1.50, "linear_strength": 0.25, "gamma_norm": 1.00,
        "local_z_strength": 0.03, "zz_scale": 0.08, "sign_prob_neg": 0.22,
    }

# ── Reported optimal architecture ───────────────────────────────────────
# The 9 shared hyperparameters found by the outer CMA-ES search (originally
# reported at the Nr=64, Nz=8 cell under the tent score), together with
# Nr, Nz re-derived from the full (Nr,Nz) grid-matrix sweep under the
# SMOOTH score (score_fn_smooth, see Sec. "Smooth score" of the protocol):
# Nr=64, Nz=8 (S=3.868) and Nr=64, Nz=12 (S=3.869) are statistically tied
# for best across the full 16-cell grid, both far below the tent score's
# best of S=7.099 -- direct evidence the smooth kernel lets calibration
# get substantially closer to S=0 with the SAME architecture/budget.
OPTIMAL_ARCH = {
    "rough_orientation": -1.0, "Nr": 64, "Nz": 8,
    "matrix_seed": MATRIX_SEED,
    "z_strength": 0.19288535085059877,
    "even_strength": 2.2291053209265703,
    "linear_strength": 0.1930215632408165,
    "gamma_norm": 0.9651266728165394,
    "local_z_strength": 0.04759055003377607,
    "zz_scale": 0.039730814766135034,
    "sign_prob_neg": 0.22189311131648579,
}
# rough_scale and z_readout are NO LONGER part of the shared architecture
# (see HP_SPEC) -- they are calibrated per-scenario by inner_calibrate,
# alongside H, lam_lo, lam_hi, az_lo, az_hi, b0_delta, scale. Their
# per-scenario calibrated values are reported in
# SMOOTH_DUAL_REGIME_RESULTS below and in the protocol's new section.

# ── Reference dual-regime fit quality at OPTIMAL_ARCH under the smooth
# score (fresh held-out re-evaluation; see protocol V5, Sec. 8.5-8.6).
# H_hat_est/H_hat_true are the realised vs. target structure-function
# Hurst exponents (Remark on the two different "H" objects, Sec. 1.2);
# lam_lo/lam_hi are the per-scenario calibrated OU-rate window. Kept here
# purely for documentation/reproducibility -- not read by any function.
SMOOTH_DUAL_REGIME_RESULTS = {
    "fBM_rough":   dict(H_hat_est=0.091, H_hat_true=0.080, delta_H=0.011,
                        lam_lo=2.2e-4, lam_hi=2.74),
    "MRW_taylor":  dict(H_hat_est=0.288, H_hat_true=0.297, delta_H=0.008,
                        lam_lo=3.4e-5, lam_hi=0.33),
    "rBergomi_Z":  dict(H_hat_est=0.117, H_hat_true=0.107, delta_H=0.010,
                        lam_lo=1.7e-5, lam_hi=2.43),
    "fBM_persist": dict(H_hat_est=0.295, H_hat_true=0.347, delta_H=0.052,
                        lam_lo=1.4e-2, lam_hi=0.25),
}
# For comparison, the same quantities under the ORIGINAL tent score at
# the tent-optimal architecture (Nr=64, Nz=20): fBM_rough delta_H=0.065,
# MRW_taylor delta_H=0.048, rBergomi_Z delta_H=0.208, fBM_persist
# delta_H=0.053. The smooth kernel tightens delta_H by ~5-20x for three
# of the four scenarios (most dramatically rBergomi_Z), and leaves
# fBM_persist essentially unchanged (it was already the least
# gradient-starved of the four under the tent kernel).

# ============================================================
# 8.  CMA-ES optimiser (self-contained, no external dependency) (unchanged)
# ============================================================

class CMAESOptimiser:
    def __init__(self, x0, sigma0=0.3, popsize=12, seed=42):
        self.dim     = len(x0)
        self.mean    = np.array(x0, dtype=float)
        self.sigma   = sigma0
        self.popsize = popsize
        self.mu      = popsize // 2
        self.rng     = np.random.default_rng(seed)

        raw_w = math.log(self.mu + 0.5) - np.log(np.arange(1, self.mu + 1))
        self.w    = raw_w / raw_w.sum()
        self.mueff = 1. / (self.w**2).sum()

        n = self.dim
        self.cs   = (self.mueff + 2.) / (n + self.mueff + 5.)
        self.ds   = (1. + 2.*max(0., math.sqrt((self.mueff-1.)/(n+1.))-1.)
                     + self.cs)
        self.chiN = math.sqrt(n) * (1. - 1./(4.*n) + 1./(21.*n**2))

        self.cc   = (4. + self.mueff/n) / (n + 4. + 2.*self.mueff/n)
        self.c1   = 2. / ((n + 1.3)**2 + self.mueff)
        self.cmu  = min(1. - self.c1,
                        2.*(self.mueff - 2. + 1./self.mueff)
                        / ((n+2.)**2 + self.mueff))

        self.ps   = np.zeros(n)
        self.pc   = np.zeros(n)
        self.B    = np.eye(n)
        self.D    = np.ones(n)
        self.C    = np.eye(n)
        self.invsqrtC = np.eye(n)
        self.eigeneval = 0
        self.counteval = 0
        self.generation = 0

    def ask(self):
        n  = self.dim
        BD = self.B * self.D[None, :]
        samples = [self.mean + self.sigma * BD @ self.rng.standard_normal(n)
                   for _ in range(self.popsize)]
        return samples

    def tell(self, solutions, fitnesses):
        n      = self.dim
        mu     = self.mu
        ranked = sorted(zip(fitnesses, range(self.popsize)))
        best_idx = [ranked[i][1] for i in range(mu)]

        x_old   = self.mean.copy()
        x_new   = sum(self.w[i] * solutions[best_idx[i]] for i in range(mu))
        self.mean = x_new

        invsqrtC_x = self.invsqrtC @ ((x_new - x_old) / self.sigma)
        self.ps = ((1. - self.cs) * self.ps
                   + math.sqrt(self.cs*(2.-self.cs)*self.mueff) * invsqrtC_x)
        hs = (np.linalg.norm(self.ps) /
              math.sqrt(1.-(1.-self.cs)**(2.*(self.counteval+1)/self.popsize))
              / self.chiN < 1.4 + 2./(n+1))
        self.sigma *= math.exp((self.cs / self.ds)
                               * (np.linalg.norm(self.ps) / self.chiN - 1.))
        self.sigma  = float(np.clip(self.sigma, 1e-5, 5.))

        self.pc = ((1. - self.cc) * self.pc
                   + hs * math.sqrt(self.cc*(2.-self.cc)*self.mueff)
                   * (x_new - x_old) / self.sigma)
        artmp = [(solutions[best_idx[i]] - x_old) / self.sigma for i in range(mu)]
        self.C = ((1. - self.c1 - self.cmu) * self.C
                  + self.c1 * (np.outer(self.pc, self.pc) + (1-hs)*self.cc*(2.-self.cc)*self.C)
                  + self.cmu * sum(self.w[i] * np.outer(artmp[i], artmp[i])
                                   for i in range(mu)))

        self.counteval += self.popsize; self.generation += 1
        if self.counteval - self.eigeneval > self.popsize / (self.c1 + self.cmu) / n / 10:
            self.eigeneval = self.counteval
            self.C = np.triu(self.C) + np.triu(self.C, 1).T
            try:
                self.D, self.B = np.linalg.eigh(self.C)
                self.D = np.sqrt(np.maximum(self.D, 1e-20))
                self.invsqrtC = self.B @ np.diag(1./self.D) @ self.B.T
            except np.linalg.LinAlgError:
                self.D = np.ones(n); self.B = np.eye(n)
                self.invsqrtC = np.eye(n)

    @property
    def best_sigma(self): return self.sigma

# ============================================================
# 9.  Main pipeline  (unchanged, other than the DGP fix upstream)
# ============================================================

def run(max_evals=MAX_OUTER_EVALS, popsize=POPSIZE,
        nr_grid=None, nz_grid=None,
        save_prefix="esn_hyperparam_search",
        verbose=True, fast=False):
    """
    CMA-ES MINIMISES S >= 0 directly (S=0 is a perfect match to all 4 DGPs,
    lower is always better -- no sign flips anywhere in this function).
    """
    if fast:
        global N_DGP_PATHS, N_CAL_PATHS, T_SIM, BURN, T_CAL, BURN_CAL, INNER_NM_MAXITER
        N_DGP_PATHS = 4; N_CAL_PATHS = 3; T_SIM = 600; BURN = 150
        T_CAL = 400; BURN_CAL = 100; max_evals = 20; popsize = 6
        INNER_NM_MAXITER = 60

    if nr_grid is None: nr_grid = NR_GRID if not fast else [32, 64]
    if nz_grid is None: nz_grid = NZ_GRID if not fast else [8, 12]

    n_cells   = len(nr_grid) * len(nz_grid)
    evals_cell = max(popsize, max_evals // n_cells)

    history  = []
    best_score = np.inf
    best_arch  = None
    best_detail = None
    cell_best_table = {}   # (Nr,Nz) -> best S in that cell

    print(f"\n{'='*65}")
    print(f"ESN Hyperparameter Search  --  CMA-ES MINIMISES S >= 0 directly")
    print(f"DGP scenarios   : {[s[0] for s in DGP_SCENARIOS]}")
    print(f"(Nr, Nz) grid   : {[(nr,nz) for nr in nr_grid for nz in nz_grid]}")
    print(f"Evals per cell  : {evals_cell}  (popsize={popsize})")
    print(f"Total budget    : ~{evals_cell * n_cells} evaluations")
    print(f"Score range     : S=0 (perfect) .. S={SCORE_MAX:.1f} (worst)")
    print(f"{'='*65}\n")

    print("Pre-generating DGP paths ...")
    for sc_name, sc_type, sc_params in DGP_SCENARIOS:
        get_dgp_paths(sc_name, sc_type, sc_params)
    print("DGP paths ready.\n")
    for sc_name, _, _ in DGP_SCENARIOS:
        ref = make_score_ref(_dgp_cache[sc_name])
        print(f"  [ref] {sc_name:14s} H_hat_c={ref['H_hat_c']:.4f}  "
              f"mean_vol_ann_c={ref['mean_vol_ann_c']:.4f}")
    print()

    eval_counter = [0]
    t0 = time.time()

    for Nr in nr_grid:
        for Nz in nz_grid:
            cell_label = f"Nr={Nr} Nz={Nz}"
            print(f"\n{'─'*55}  {cell_label}")

            arch0  = default_arch(Nr=Nr, Nz=Nz)
            x0     = encode_hp(arch0)
            cmaes  = CMAESOptimiser(x0, sigma0=0.35, popsize=popsize, seed=Nr*100+Nz)
            cell_best = np.inf

            n_gens = max(1, evals_cell // popsize)
            for gen in range(n_gens):
                candidates = cmaes.ask()
                fits = []
                for x_cand in candidates:
                    arch = decode_hp(x_cand, Nr, Nz)
                    eid  = eval_counter[0]; eval_counter[0] += 1
                    sc, detail = evaluate_arch(arch, eval_id=eid, verbose=False, return_detail=True)
                    fits.append(sc)   # CMA-ES minimises S directly -- no sign flip

                    record = {"eval_id": eid, "Nr": Nr, "Nz": Nz,
                              "score": sc, "arch": {k: v for k, v in arch.items()
                                                     if k != "matrix_seed"},
                              "detail": [{k: v for k, v in d.items() if k != "ref"}
                                         for d in detail]}
                    history.append(record)

                    if sc < cell_best:
                        cell_best = sc
                        if verbose:
                            print(f"  gen {gen:3d}  eval {eid:4d}  "
                                  f"S={sc:.4f}  {cell_label}  *** new cell best",
                                  flush=True)
                    if sc < best_score:
                        best_score = sc
                        best_arch  = copy.deepcopy(arch)
                        best_detail = detail
                        if verbose:
                            print(f"  gen {gen:3d}  eval {eid:4d}  "
                                  f"S={sc:.4f}  *** GLOBAL BEST ***",
                                  flush=True)
                    elif verbose and eid % popsize == 0:
                        elapsed = time.time() - t0
                        print(f"  gen {gen:3d}  eval {eid:4d}  "
                              f"S={sc:.4f}  "
                              f"cell_best={cell_best:.4f}  "
                              f"global_best={best_score:.4f}  "
                              f"[{elapsed:.0f}s]",
                              flush=True)

                cmaes.tell(candidates, fits)

            cell_best_table[(Nr, Nz)] = cell_best

            # Incremental checkpoint (robust to interruption on wide grids)
            try:
                ckpt = {
                    "cell_best": {f"Nr={k[0]}_Nz={k[1]}": v for k, v in cell_best_table.items()},
                    "best_score": best_score,
                    "best_arch": best_arch,
                    "n_evals_so_far": eval_counter[0],
                    "elapsed_s_so_far": time.time() - t0,
                }
                with open(f"{save_prefix}_checkpoint.json", "w") as f:
                    json.dump(ckpt, f, indent=2)
            except Exception as e:
                print(f"Checkpoint save failed: {e}")

    elapsed = time.time() - t0
    print(f"\n{'='*65}")
    print(f"Search complete in {elapsed:.0f}s  ({eval_counter[0]} evaluations)")
    print(f"Best global score: S* = {best_score:.4f}   (0=perfect, {SCORE_MAX:.1f}=worst)")
    print(f"Best architecture:")
    for k, v in best_arch.items():
        if k != "matrix_seed":
            print(f"  {k:<22} = {v}")
    print(f"{'='*65}\n")

    print("Final evaluation of best architecture (verbose):")
    final_score, final_detail = evaluate_arch(best_arch, eval_id=9999, verbose=True, return_detail=True)
    print(f"Final mean score: {final_score:.4f}")
    # Prefer the freshly recomputed detail (same seed convention as the
    # printed per-scenario breakdown) for the report / figures.
    best_detail = final_detail

    results = {
        "best_score":  best_score,
        "final_score": final_score,
        "best_arch":   best_arch,
        "n_evals":     eval_counter[0],
        "elapsed_s":   elapsed,
        "cell_best":   {f"Nr={k[0]}_Nz={k[1]}": v for k, v in cell_best_table.items()},
        "per_scenario": [{k: v for k, v in d.items() if k not in ("ref", "inner_params")}
                          for d in best_detail],
        "history":     history,
    }
    json_path = f"{save_prefix}_results.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Results saved -> {json_path}")

    best_path = f"{save_prefix}_best.json"
    with open(best_path, "w") as f:
        json.dump({"best_arch": best_arch, "best_score": best_score,
                   "final_score": final_score,
                   "per_scenario": results["per_scenario"]}, f, indent=2)
    print(f"Best arch saved -> {best_path}")

    try:
        make_figure_overview(history, best_arch, best_score, best_detail,
                              nr_grid, nz_grid, cell_best_table,
                              save_path=f"{save_prefix}_overview.png")
        print(f"Overview figure saved -> {save_prefix}_overview.png")
    except Exception as e:
        print(f"Overview figure failed: {e}")

    try:
        make_figure_deviation(history, best_arch, best_detail,
                              nr_grid, nz_grid,
                              save_path=f"{save_prefix}_deviation.png")
        print(f"Deviation figure saved -> {save_prefix}_deviation.png")
    except Exception as e:
        print(f"Deviation figure failed: {e}")

    return best_arch, best_score, history


# ============================================================
# 9b.  Diagnostic figures
# ============================================================

def make_figure_overview(history, best_arch, best_score, best_detail,
                          nr_grid, nz_grid, cell_best_table, save_path):
    """
    6-panel overview:
      1. Convergence (S minimised directly, running MIN)
      2. Score distribution vs Nr (boxplot)
      3. Best S per (Nr,Nz) cell (heatmap, optimal cell highlighted)
      4. Per-DGP score at the optimal architecture
      5. Hyperparameter shift from the A2-981003 default (ratio)
      6. Estimated H (per-DGP calibrated) vs the DGP's true H
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    eids   = [h["eval_id"] for h in history]
    scores = [h["score"]   for h in history]
    running_min = np.minimum.accumulate(scores)
    nr_used = sorted(set(h["Nr"] for h in history))
    nz_used = sorted(set(h["Nz"] for h in history))

    fig = plt.figure(figsize=(19, 11))
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.28, top=0.85, bottom=0.06)

    fig.suptitle(
        "ESN Hyperparameter Search — CMA-ES minimising $S\\geq 0$ directly\n"
        f"Score: $S=\\sum_k w_k(1-f_k)+5\\cdot\\mathrm{{Stress}}\\in[0,{SCORE_MAX:.1f}]$"
        "  |  $S=0$: perfect match  |  lower = better\n"
        f"Optimal: $N_r^*={best_arch['Nr']}$, $N_z^*={best_arch['Nz']}$  |  "
        f"$S^*={best_score:.3f}$  |  {len(history)} evals\n"
        "fBM-persist DGP fixed: annualised (stationary) time argument -- no volatility collapse",
        fontsize=12.5, fontweight="bold")

    # --- Panel 1: convergence -------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(eids, scores, s=10, alpha=0.35, color="#5b84b1", label="$S$ per eval")
    ax.plot(eids, running_min, color="#c0392b", lw=1.8, label="Running min $S$")
    ax.axhline(best_score, color="#1D9E75", lw=1.2, ls="--", label=f"$S^*={best_score:.3f}$")
    ax.axhline(0, color="gray", lw=0.8, ls=":", label="$S=0$ ideal")
    ax.axhline(SCORE_MAX, color="lightgray", lw=0.8, ls=":", label=f"$S={SCORE_MAX:.1f}$ worst")
    ax.set_xlabel("Evaluation"); ax.set_ylabel("Score $S\\geq0$")
    ax.set_title("① Convergence\nCMA-ES minimises $S$ directly")
    ax.legend(fontsize=7.5, loc="upper right"); ax.grid(True, alpha=0.2)

    # --- Panel 2: score distribution vs Nr -------------------------------
    ax = fig.add_subplot(gs[0, 1])
    box_data = [[h["score"] for h in history if h["Nr"] == nr] for nr in nr_used]
    bp = ax.boxplot(box_data, labels=[f"$N_r={nr}$" for nr in nr_used], patch_artist=True)
    palette = ["#5b84b1", "#e08e6d", "#59b797", "#c19be0"]
    for patch, col in zip(bp["boxes"], palette):
        patch.set_facecolor(col); patch.set_alpha(0.55)
    ax.axhline(best_score, color="#c0392b", lw=1.0, ls="--", label=f"$S^*={best_score:.2f}$")
    ax.set_ylabel("$S\\geq0$ (lower=better)")
    ax.set_title("② Score distribution vs $N_r$\n(lower box = better)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)

    # --- Panel 3: best S per (Nr,Nz) cell heatmap ------------------------
    ax = fig.add_subplot(gs[0, 2])
    grid = np.full((len(nr_used), len(nz_used)), np.nan)
    for (nr, nz), v in cell_best_table.items():
        if nr in nr_used and nz in nz_used:
            grid[nr_used.index(nr), nz_used.index(nz)] = v
    im = ax.imshow(grid, cmap="RdYlGn_r", aspect="auto")
    for i in range(len(nr_used)):
        for j in range(len(nz_used)):
            if np.isfinite(grid[i, j]):
                is_best = (nr_used[i] == best_arch["Nr"] and nz_used[j] == best_arch["Nz"])
                ax.text(j, i, f"{grid[i,j]:.2f}", ha="center", va="center",
                        fontsize=11, fontweight="bold",
                        color="black")
                if is_best:
                    rect = plt.Rectangle((j-0.5, i-0.5), 1, 1, fill=False,
                                          edgecolor="gold", lw=3.5)
                    ax.add_patch(rect)
    ax.set_xticks(range(len(nz_used))); ax.set_xticklabels([f"$N_z={nz}$" for nz in nz_used])
    ax.set_yticks(range(len(nr_used))); ax.set_yticklabels([f"$N_r={nr}$" for nr in nr_used])
    plt.colorbar(im, ax=ax, label="Best $S$ (lower=better)", shrink=0.85)
    ax.set_title("③ Best $S$ per $(N_r,N_z)$ cell\n(gold border = optimal)")

    # --- Panel 4: per-DGP score at optimal arch --------------------------
    ax = fig.add_subplot(gs[1, 0])
    names  = [d["name"] for d in best_detail]
    dscores = [d["score"] for d in best_detail]
    colors  = ["#5b84b1", "#e08e6d", "#59b797", "#9a9a9a"]
    bars = ax.bar(names, dscores, color=colors[:len(names)])
    ax.axhline(np.mean(dscores), color="#c0392b", ls="--", lw=1.2,
               label=f"Mean $S$={np.mean(dscores):.3f}")
    ax.axhline(0, color="gray", ls=":", lw=0.8, label="$S=0$ (ideal)")
    ax.axhline(SCORE_MAX, color="lightgray", ls=":", lw=0.8, label=f"$S={SCORE_MAX:.1f}$ (worst)")
    for b, v in zip(bars, dscores):
        ax.text(b.get_x()+b.get_width()/2, v+0.05, f"{v:.3f}", ha="center", fontsize=9, fontweight="bold")
    ax.set_ylabel("$S\\geq0$ (lower=better)")
    ax.set_title("④ Per-DGP score — optimal arch\n(final held-out paths)")
    ax.legend(fontsize=7.5); ax.grid(True, alpha=0.2, axis="y")
    plt.setp(ax.get_xticklabels(), rotation=0, fontsize=8)

    # --- Panel 5: hyperparameter shift from default ----------------------
    ax = fig.add_subplot(gs[1, 1])
    arch0 = default_arch(Nr=best_arch["Nr"], Nz=best_arch["Nz"])
    hp_names = [spec[0] for spec in HP_SPEC]
    ratios = [best_arch[n] / arch0[n] for n in hp_names]
    bar_colors = ["#c0392b" if r > 1.15 else ("#3a6fb5" if r < 0.85 else "#9a9a9a") for r in ratios]
    ax.barh(hp_names, ratios, color=bar_colors)
    ax.axvline(1.0, color="black", ls="--", lw=1.0, label="Default A2-981003")
    for i, r in enumerate(ratios):
        ax.text(r + 0.02, i, f"{best_arch[hp_names[i]]:.4f}", va="center", fontsize=7.5)
    ax.set_xlabel("Optimal / Default ratio")
    ax.set_title("⑤ Hyperparameter shift from A2-981003\n(red ↑>15%, blue ↓>15%, grey ≈unchanged)")
    ax.legend(fontsize=8)

    # --- Panel 6: estimated H vs true H per DGP --------------------------
    ax = fig.add_subplot(gs[1, 2])
    x = np.arange(len(best_detail))
    h_est  = [d["H_hat_est"]  for d in best_detail]
    h_true = [d["H_hat_true"] for d in best_detail]
    bars = ax.bar(x, h_est, color=colors[:len(x)], label="Est. $\\hat H$ (achieved)")
    ax.scatter(x, h_true, color="#c0392b", marker="D", s=70, zorder=5, label="True $\\hat H_c$ (DGP target)")
    for xi, (e, t) in enumerate(zip(h_est, h_true)):
        ax.text(xi, e + 0.012, f"{e:.3f}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d['name']}" for d in best_detail], fontsize=8)
    ax.set_ylabel("Hurst / roughness exponent $\\hat H$")
    ax.set_title("⑥ Achieved $\\hat H$ vs DGP target $\\hat H_c$\n(per-scenario, calibrated, multi-start NM)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, axis="y")

    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def make_figure_deviation(history, best_arch, best_detail, nr_grid, nz_grid, save_path):
    """
    2-panel deviation figure:
      1. Per-(Nr,Nz)-cell score trajectory (raw + running min)
      2. 11-term score deviation breakdown at the optimal architecture,
         averaged across all 4 DGP scenarios and their calibration paths.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(15.5, 5.6))
    gs = fig.add_gridspec(1, 2, top=0.78, wspace=0.28)
    fig.suptitle(
        "Score $S\\geq 0$ per cell  |  11-term deviation breakdown\n"
        "$S=\\sum_k w_k(1-f_k)+5\\cdot\\mathrm{Stress}$  |  $d_k=0$ = perfect, $d_k=w_k$ = worst",
        fontsize=12, fontweight="bold")

    # --- Panel 1: per-cell trajectories -----------------------------------
    ax = fig.add_subplot(gs[0, 0])
    cmap = plt.get_cmap("tab10")
    cells = sorted(set((h["Nr"], h["Nz"]) for h in history))
    for ci, (nr, nz) in enumerate(cells):
        sub = [h for h in history if h["Nr"] == nr and h["Nz"] == nz]
        sub = sorted(sub, key=lambda h: h["eval_id"])
        eids = [h["eval_id"] for h in sub]
        scs  = [h["score"]   for h in sub]
        running_min = np.minimum.accumulate(scs)
        color = cmap(ci % 10)
        ax.scatter(eids, scs, s=14, alpha=0.55, color=color)
        ax.plot(eids, running_min, color=color, lw=2.0, label=f"$N_r={nr}, N_z={nz}$")
    best_score = min(h["score"] for h in history)
    ax.axhline(best_score, color="#c0392b", ls="--", lw=1.2, label=f"$S^*={best_score:.3f}$")
    ax.axhline(0, color="gray", ls=":", lw=0.8, label="$S=0$ ideal")
    ax.set_xlabel("Evaluation"); ax.set_ylabel("$S\\geq0$ (lower=better)")
    ax.set_title("Score trajectory per $(N_r,N_z)$ cell\n(dots=raw, line=running min $S$)")
    ax.legend(fontsize=7, ncol=1); ax.grid(True, alpha=0.2)

    # --- Panel 2: 11-term deviation breakdown at the optimum --------------
    ax = fig.add_subplot(gs[0, 1])
    accum = {k: [] for k, _ in SCORE_TERM_WEIGHTS}
    for d in best_detail:
        if d.get("inner_params") is None or d.get("ref") is None:
            continue
        mat = build_esn_matrices(best_arch, d["inner_params"])
        for q in range(N_CAL_PATHS):
            x, v = _sim_esn_with_params(
                7000 + q, T_CAL, best_arch, mat,
                b0_delta=d["inner_params"].get("b0_delta", 0.),
                scale=d["inner_params"].get("scale", 1.))
            st = compute_statistics(x[BURN_CAL:], v[BURN_CAL:])
            dev = score_breakdown(st, d["ref"])
            for k in accum: accum[k].append(dev[k])

    labels_map = {"H_hat": "$\\hat H$", "mean_vol_ann": "$\\bar\\sigma$",
                  "q995_vol_ann": "$q_{995}$", "max_vol_ann": "$\\sigma_{\\max}$",
                  "mean_vol_acf": "$\\bar V_{ACF}$", "taylor_gap": "$F_T$",
                  "taylor_frac": "$G_T$", "zumbach": "$Z$", "leverage": "$L$",
                  "kurtosis": "$K$", "max_ret_acf": "$A$"}
    weights = dict(SCORE_TERM_WEIGHTS)
    order = ["max_ret_acf","kurtosis","leverage","zumbach","taylor_frac","taylor_gap",
             "mean_vol_acf","max_vol_ann","q995_vol_ann","mean_vol_ann","H_hat"]
    means = [float(np.mean(accum[k])) if accum[k] else 0.0 for k in order]
    wts   = [weights[k] for k in order]
    ylabels = [f"{wts[i]:.1f} (1-{labels_map[order[i]].replace('$','')})" for i in range(len(order))]
    bar_colors = []
    for m, w in zip(means, wts):
        if m < 0.4*w: bar_colors.append("#59b797")
        elif m < 0.7*w: bar_colors.append("#e08e6d")
        else: bar_colors.append("#c0392b")
    ax.barh(range(len(order)), means, color=bar_colors)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(ylabels, fontsize=9)
    for i, (m, w) in enumerate(zip(means, wts)):
        ax.text(m + 0.02, i, f"{m:.3f}/{w:.1f}", va="center", fontsize=8)
    ax.axvline(0, color="gray", ls=":", lw=1.0, label="$d_k=0$ (perfect term)")
    ax.set_xlabel("Deviation $d_k=w_k(1-f_k)$  (0=perfect, $w_k$=worst)")
    ax.set_title(f"11-term deviation at optimal arch\n"
                 f"(green $d_k<0.4w_k$, orange $<0.7w_k$, red $\\geq0.7w_k$)\n"
                 f"Mean over {len(best_detail)} DGPs × {N_CAL_PATHS} paths  |  "
                 f"$\\sum d_k$={sum(means):.3f}")
    ax.legend(fontsize=8, loc="lower right")

    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# 9c.  Full-grid score matrix for a FIXED hyperparameter set
#      (varies only Nr, Nz; re-calibrates H_ESN/lambda per cell as usual)
# ============================================================

def evaluate_fixed_arch_over_grid(hp, nr_grid=None, nz_grid=None,
                                   eval_id=0, verbose=True,
                                   checkpoint_path=None):
    """
    Evaluate ONE fixed hyperparameter set (the 9 shared architecture
    hyperparameters of HP_SPEC, e.g. the outer search's reported optimum)
    at every cell of a (Nr, Nz) grid, sweeping ONLY the reservoir sizes.

    This isolates the effect of reservoir size from the effect of the
    other 9 hyperparameters (which are held fixed at whatever `hp`
    specifies), unlike `run()` which re-optimises all 9 within every
    cell. Each cell still runs the usual per-DGP-scenario multi-start
    inner calibration (Sec. 5 / INNER_NM_STARTS) -- only Nr, Nz change.

    Args:
        hp:  dict with the 9 shared hyperparameters, e.g.
             {"rough_scale":..., "z_strength":..., "z_readout":...,
              "even_strength":..., "linear_strength":..., "gamma_norm":...,
              "local_z_strength":..., "zz_scale":..., "sign_prob_neg":...}
             (rough_orientation is always fixed to -1.0; matrix_seed is
             always MATRIX_SEED -- neither needs to be supplied).
        nr_grid, nz_grid: grid to sweep (defaults to the full production
             grid NR_GRID x NZ_GRID).
        checkpoint_path: if given, partial results are written to this
             JSON path after every cell (robust to interruption on wide
             grids).

    Returns: dict {"Nr={nr}_Nz={nz}": score, ...}
    """
    if nr_grid is None: nr_grid = NR_GRID
    if nz_grid is None: nz_grid = NZ_GRID

    for sc_name, sc_type, sc_params in DGP_SCENARIOS:
        get_dgp_paths(sc_name, sc_type, sc_params)

    results = {}
    t0 = time.time()
    for Nr in nr_grid:
        for Nz in nz_grid:
            arch = {"rough_orientation": -1.0, "Nr": Nr, "Nz": Nz,
                    "matrix_seed": MATRIX_SEED, **hp}
            t1 = time.time()
            sc = evaluate_arch(arch, eval_id=eval_id, verbose=False,
                                return_detail=False)
            key = f"Nr={Nr}_Nz={Nz}"
            results[key] = sc
            if verbose:
                print(f"  Nr={Nr:3d} Nz={Nz:3d}  S={sc:.4f}   "
                      f"({time.time()-t1:.1f}s)", flush=True)
            if checkpoint_path:
                with open(checkpoint_path, "w") as f:
                    json.dump({"results": results,
                               "elapsed_s": time.time() - t0}, f, indent=2)
    return results


def make_figure_grid_matrix(results, nr_grid, nz_grid, hp, save_path):
    """
    Heatmap of `results` (as returned by evaluate_fixed_arch_over_grid)
    over the full (Nr, Nz) grid: one cell per (Nr, Nz), annotated with
    its score S, gold border on the best (lowest-S) cell.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    grid = np.full((len(nr_grid), len(nz_grid)), np.nan)
    for i, Nr in enumerate(nr_grid):
        for j, Nz in enumerate(nz_grid):
            grid[i, j] = results[f"Nr={Nr}_Nz={Nz}"]

    best_idx = np.unravel_index(np.nanargmin(grid), grid.shape)
    best_nr, best_nz = nr_grid[best_idx[0]], nz_grid[best_idx[1]]
    best_val = grid[best_idx]

    fig, ax = plt.subplots(figsize=(8.2, 6.4))
    im = ax.imshow(grid, cmap="RdYlGn_r", aspect="auto",
                   vmin=np.nanmin(grid), vmax=np.nanmax(grid))

    for i in range(len(nr_grid)):
        for j in range(len(nz_grid)):
            val = grid[i, j]
            if not np.isfinite(val):
                continue
            is_best = (i == best_idx[0] and j == best_idx[1])
            ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                    fontsize=12.5, fontweight="bold", color="black")
            if is_best:
                rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                      edgecolor="gold", lw=4)
                ax.add_patch(rect)

    ax.set_xticks(range(len(nz_grid)))
    ax.set_xticklabels([f"$N_z={nz}$" for nz in nz_grid], fontsize=11)
    ax.set_yticks(range(len(nr_grid)))
    ax.set_yticklabels([f"$N_r={nr}$" for nr in nr_grid], fontsize=11)
    cbar = plt.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("Score $S\\geq0$  (lower = better)", fontsize=11)

    hp_str = ", ".join(f"{k.replace('_',' ')}={v:.4f}"
                        for k, v in list(hp.items())[:2])
    ax.set_title(
        "Score $S$ over the full $(N_r,N_z)$ grid\n"
        f"hyperparameters fixed at the reported optimum ({hp_str}, ...)\n"
        f"best cell: $N_r={best_nr}$, $N_z={best_nz}$  ($S^*={best_val:.4f}$)",
        fontsize=12)
    plt.tight_layout()
    fig.savefig(save_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return best_nr, best_nz, best_val


def run_grid_matrix(hp=None, nr_grid=None, nz_grid=None,
                     save_prefix="esn_hyperparam_search_gridmatrix",
                     fast=False):
    """
    Convenience wrapper: evaluate a fixed hyperparameter set over the
    full (Nr, Nz) grid, save the raw scores to JSON, and save the
    heatmap figure. Returns (results_dict, (best_nr, best_nz, best_val)).

    If `hp` is omitted, defaults to the 9 shared hyperparameters of
    OPTIMAL_ARCH. Sweeping the full grid with these fixed hyperparameters
    finds Nr=64, Nz=20 as the best cell (S*=7.099), not the Nr=64, Nz=8
    cell originally reported from the smaller demonstration search --
    OPTIMAL_ARCH has been updated accordingly (see its definition, Sec.
    6.2/6.3 of the protocol).
    """
    if fast:
        global N_DGP_PATHS, N_CAL_PATHS, T_SIM, BURN, T_CAL, BURN_CAL, INNER_NM_MAXITER
        N_DGP_PATHS = 2; N_CAL_PATHS = 2; T_SIM = 200; BURN = 50
        T_CAL = 120; BURN_CAL = 30; INNER_NM_MAXITER = 20

    if hp is None:
        hp = {k: v for k, v in OPTIMAL_ARCH.items()
              if k not in ("rough_orientation", "Nr", "Nz", "matrix_seed")}
    if nr_grid is None: nr_grid = NR_GRID
    if nz_grid is None: nz_grid = NZ_GRID

    ckpt_path = f"{save_prefix}_checkpoint.json"
    results = evaluate_fixed_arch_over_grid(
        hp, nr_grid=nr_grid, nz_grid=nz_grid, verbose=True,
        checkpoint_path=ckpt_path)

    with open(f"{save_prefix}_results.json", "w") as f:
        json.dump({"hp": hp, "nr_grid": nr_grid, "nz_grid": nz_grid,
                   "results": results}, f, indent=2)

    best_nr, best_nz, best_val = make_figure_grid_matrix(
        results, nr_grid, nz_grid, hp, save_path=f"{save_prefix}_heatmap.png")
    print(f"Grid-matrix heatmap saved -> {save_prefix}_heatmap.png")
    print(f"Best cell: Nr={best_nr}, Nz={best_nz}  (S*={best_val:.4f})")
    return results, (best_nr, best_nz, best_val)


# ============================================================
# 9d.  Full-grid (all 16 cells) versions of the overview/deviation figures
# ============================================================

def make_figure_overview_gridmatrix(grid_results, nr_grid, nz_grid,
                                     best_arch, best_detail, save_path):
    """
    Same 6-panel layout as `make_figure_overview`, but panels 1-3 are
    built from the COMPLETE (Nr,Nz) grid-matrix sweep (all 16 cells of
    the production grid) instead of a partial CMA-ES search history:
      1. Score S at every grid cell, in (Nr,Nz) sweep order, + running min
      2. Score distribution vs Nr (boxplot, ALL Nz values per Nr)
      3. Best S per (Nr,Nz) cell -- the full 4x4 heatmap
      4. Per-DGP score at the optimal architecture (held-out paths)
      5. Hyperparameter shift from the A2-981003 default
      6. Estimated H (per-DGP calibrated) vs the DGP's true H
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cells = [(nr, nz) for nr in nr_grid for nz in nz_grid]
    scores = [grid_results[f"Nr={nr}_Nz={nz}"] for nr, nz in cells]
    eids = list(range(len(cells)))
    running_min = np.minimum.accumulate(scores)
    best_i = int(np.argmin(scores))
    best_score = scores[best_i]

    fig = plt.figure(figsize=(19, 11))
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.28, top=0.85, bottom=0.06)

    fig.suptitle(
        "ESN Hyperparameter Search — full $(N_r,N_z)$ grid-matrix sweep\n"
        f"Score: $S=\\sum_k w_k(1-f_k)+5\\cdot\\mathrm{{Stress}}\\in[0,{SCORE_MAX:.1f}]$"
        "  |  $S=0$: perfect match  |  lower = better\n"
        f"Optimal: $N_r^*={best_arch['Nr']}$, $N_z^*={best_arch['Nz']}$  |  "
        f"$S^*={best_score:.3f}$  |  {len(cells)} cells (complete production grid)\n"
        "fBM-persist DGP fixed: annualised (stationary) time argument -- no volatility collapse",
        fontsize=12.5, fontweight="bold")

    # --- Panel 1: score across the full sweep -----------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(eids, scores, s=22, alpha=0.6, color="#5b84b1", label="$S$ per cell")
    ax.plot(eids, running_min, color="#c0392b", lw=1.8, label="Running min $S$")
    ax.axhline(best_score, color="#1D9E75", lw=1.2, ls="--", label=f"$S^*={best_score:.3f}$")
    ax.axhline(0, color="gray", lw=0.8, ls=":", label="$S=0$ ideal")
    ax.set_xlabel("Grid cell (sweep order: $N_r$ outer, $N_z$ inner)")
    ax.set_ylabel("Score $S\\geq0$")
    ax.set_title("① Score across the full grid sweep\n(all 16 cells, fixed hyperparameters)")
    ax.legend(fontsize=7.5, loc="upper right"); ax.grid(True, alpha=0.2)

    # --- Panel 2: score distribution vs Nr (all Nz per Nr) -----------------
    ax = fig.add_subplot(gs[0, 1])
    box_data = [[grid_results[f"Nr={nr}_Nz={nz}"] for nz in nz_grid] for nr in nr_grid]
    bp = ax.boxplot(box_data, labels=[f"$N_r={nr}$" for nr in nr_grid], patch_artist=True)
    palette = ["#5b84b1", "#e08e6d", "#59b797", "#c19be0"]
    for patch, col in zip(bp["boxes"], palette):
        patch.set_facecolor(col); patch.set_alpha(0.55)
    ax.axhline(best_score, color="#c0392b", lw=1.0, ls="--", label=f"$S^*={best_score:.2f}$")
    ax.set_ylabel("$S\\geq0$ (lower=better)")
    ax.set_title(f"② Score distribution vs $N_r$\n(all {len(nz_grid)} $N_z$ values per box)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)

    # --- Panel 3: best S per (Nr,Nz) cell -- FULL heatmap ------------------
    ax = fig.add_subplot(gs[0, 2])
    grid = np.array([[grid_results[f"Nr={nr}_Nz={nz}"] for nz in nz_grid]
                      for nr in nr_grid])
    im = ax.imshow(grid, cmap="RdYlGn_r", aspect="auto")
    for i in range(len(nr_grid)):
        for j in range(len(nz_grid)):
            is_best = (nr_grid[i] == best_arch["Nr"] and nz_grid[j] == best_arch["Nz"])
            ax.text(j, i, f"{grid[i,j]:.2f}", ha="center", va="center",
                    fontsize=11, fontweight="bold", color="black")
            if is_best:
                rect = plt.Rectangle((j-0.5, i-0.5), 1, 1, fill=False,
                                      edgecolor="gold", lw=3.5)
                ax.add_patch(rect)
    ax.set_xticks(range(len(nz_grid))); ax.set_xticklabels([f"$N_z={nz}$" for nz in nz_grid])
    ax.set_yticks(range(len(nr_grid))); ax.set_yticklabels([f"$N_r={nr}$" for nr in nr_grid])
    plt.colorbar(im, ax=ax, label="$S$ (lower=better)", shrink=0.85)
    ax.set_title(f"③ $S$ per $(N_r,N_z)$ cell -- complete grid\n"
                 f"({len(nr_grid)}$\\times${len(nz_grid)}={len(cells)} cells; gold = optimal)")

    # --- Panel 4: per-DGP score at optimal arch ----------------------------
    ax = fig.add_subplot(gs[1, 0])
    names  = [d["name"] for d in best_detail]
    dscores = [d["score"] for d in best_detail]
    colors  = ["#5b84b1", "#e08e6d", "#59b797", "#9a9a9a"]
    bars = ax.bar(names, dscores, color=colors[:len(names)])
    ax.axhline(np.mean(dscores), color="#c0392b", ls="--", lw=1.2,
               label=f"Mean $S$={np.mean(dscores):.3f}")
    ax.axhline(0, color="gray", ls=":", lw=0.8, label="$S=0$ (ideal)")
    ax.axhline(SCORE_MAX, color="lightgray", ls=":", lw=0.8, label=f"$S={SCORE_MAX:.1f}$ (worst)")
    for b, v in zip(bars, dscores):
        ax.text(b.get_x()+b.get_width()/2, v+0.05, f"{v:.3f}", ha="center", fontsize=9, fontweight="bold")
    ax.set_ylabel("$S\\geq0$ (lower=better)")
    ax.set_title("④ Per-DGP score — optimal arch\n(final held-out paths)")
    ax.legend(fontsize=7.5); ax.grid(True, alpha=0.2, axis="y")
    plt.setp(ax.get_xticklabels(), rotation=0, fontsize=8)

    # --- Panel 5: hyperparameter shift from default ------------------------
    ax = fig.add_subplot(gs[1, 1])
    arch0 = default_arch(Nr=best_arch["Nr"], Nz=best_arch["Nz"])
    hp_names = [spec[0] for spec in HP_SPEC]
    ratios = [best_arch[n] / arch0[n] for n in hp_names]
    bar_colors = ["#c0392b" if r > 1.15 else ("#3a6fb5" if r < 0.85 else "#9a9a9a") for r in ratios]
    ax.barh(hp_names, ratios, color=bar_colors)
    ax.axvline(1.0, color="black", ls="--", lw=1.0, label="Default A2-981003")
    for i, r in enumerate(ratios):
        ax.text(r + 0.02, i, f"{best_arch[hp_names[i]]:.4f}", va="center", fontsize=7.5)
    ax.set_xlabel("Optimal / Default ratio")
    ax.set_title("⑤ Hyperparameter shift from A2-981003\n(red ↑>15%, blue ↓>15%, grey ≈unchanged)")
    ax.legend(fontsize=8)

    # --- Panel 6: estimated H vs true H per DGP ----------------------------
    ax = fig.add_subplot(gs[1, 2])
    x = np.arange(len(best_detail))
    h_est  = [d["H_hat_est"]  for d in best_detail]
    h_true = [d["H_hat_true"] for d in best_detail]
    bars = ax.bar(x, h_est, color=colors[:len(x)], label="Est. $\\hat H$ (achieved)")
    ax.scatter(x, h_true, color="#c0392b", marker="D", s=70, zorder=5, label="True $\\hat H_c$ (DGP target)")
    for xi, (e, t) in enumerate(zip(h_est, h_true)):
        ax.text(xi, e + 0.012, f"{e:.3f}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d['name']}" for d in best_detail], fontsize=8)
    ax.set_ylabel("Hurst / roughness exponent $\\hat H$")
    ax.set_title("⑥ Achieved $\\hat H$ vs DGP target $\\hat H_c$\n(per-scenario, calibrated, multi-start NM)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, axis="y")

    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def make_figure_deviation_gridmatrix(grid_results, nr_grid, nz_grid,
                                      best_arch, best_detail, save_path):
    """
    Same 2-panel layout as `make_figure_deviation`, but the left panel
    shows S at EVERY cell of the complete (Nr,Nz) grid (bar chart grouped
    by Nr, best cell starred) instead of a partial CMA-ES search history.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(15.5, 5.6))
    gs = fig.add_gridspec(1, 2, top=0.78, wspace=0.28)
    fig.suptitle(
        "Score $S\\geq 0$ over the complete $(N_r,N_z)$ grid  |  11-term deviation breakdown\n"
        "$S=\\sum_k w_k(1-f_k)+5\\cdot\\mathrm{Stress}$  |  $d_k=0$ = perfect, $d_k=w_k$ = worst",
        fontsize=12, fontweight="bold")

    # --- Panel 1: S at every cell, grouped by Nr --------------------------
    ax = fig.add_subplot(gs[0, 0])
    cmap = plt.get_cmap("tab10")
    best_score = min(grid_results.values())
    x0 = 0
    xticks, xlabels = [], []
    for ci, nr in enumerate(nr_grid):
        color = cmap(ci % 10)
        vals = [grid_results[f"Nr={nr}_Nz={nz}"] for nz in nz_grid]
        xs = [x0 + k for k in range(len(nz_grid))]
        bars = ax.bar(xs, vals, color=color, alpha=0.85, label=f"$N_r={nr}$")
        for xi, v, nz in zip(xs, vals, nz_grid):
            if abs(v - best_score) < 1e-9:
                ax.scatter([xi], [v + 0.03], marker="*", s=140, color="gold",
                           edgecolor="black", zorder=5)
        xticks += xs
        xlabels += [f"$N_z${nz}" for nz in nz_grid]
        x0 += len(nz_grid) + 1
    ax.axhline(best_score, color="#c0392b", ls="--", lw=1.2, label=f"$S^*={best_score:.3f}$")
    ax.axhline(0, color="gray", ls=":", lw=0.8, label="$S=0$ ideal")
    ax.set_xticks(xticks); ax.set_xticklabels(xlabels, fontsize=6.5, rotation=90)
    ax.set_ylabel("$S\\geq0$ (lower=better)")
    ax.set_title("Score $S$ at every cell of the complete grid\n"
                 "(grouped by $N_r$; gold star = optimal cell)")
    ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.2, axis="y")

    # --- Panel 2: 11-term deviation breakdown at the optimum --------------
    ax = fig.add_subplot(gs[0, 1])
    accum = {k: [] for k, _ in SCORE_TERM_WEIGHTS}
    for d in best_detail:
        if d.get("inner_params") is None or d.get("ref") is None:
            continue
        mat = build_esn_matrices(best_arch, d["inner_params"])
        for q in range(N_CAL_PATHS):
            x, v = _sim_esn_with_params(
                7000 + q, T_CAL, best_arch, mat,
                b0_delta=d["inner_params"].get("b0_delta", 0.),
                scale=d["inner_params"].get("scale", 1.))
            st = compute_statistics(x[BURN_CAL:], v[BURN_CAL:])
            dev = score_breakdown(st, d["ref"])
            for k in accum: accum[k].append(dev[k])

    labels_map = {"H_hat": "$\\hat H$", "mean_vol_ann": "$\\bar\\sigma$",
                  "q995_vol_ann": "$q_{995}$", "max_vol_ann": "$\\sigma_{\\max}$",
                  "mean_vol_acf": "$\\bar V_{ACF}$", "taylor_gap": "$F_T$",
                  "taylor_frac": "$G_T$", "zumbach": "$Z$", "leverage": "$L$",
                  "kurtosis": "$K$", "max_ret_acf": "$A$"}
    weights = dict(SCORE_TERM_WEIGHTS)
    order = ["max_ret_acf","kurtosis","leverage","zumbach","taylor_frac","taylor_gap",
             "mean_vol_acf","max_vol_ann","q995_vol_ann","mean_vol_ann","H_hat"]
    means = [float(np.mean(accum[k])) if accum[k] else 0.0 for k in order]
    wts   = [weights[k] for k in order]
    ylabels = [f"{wts[i]:.1f} (1-{labels_map[order[i]].replace('$','')})" for i in range(len(order))]
    bar_colors = []
    for m, w in zip(means, wts):
        if m < 0.4*w: bar_colors.append("#59b797")
        elif m < 0.7*w: bar_colors.append("#e08e6d")
        else: bar_colors.append("#c0392b")
    ax.barh(range(len(order)), means, color=bar_colors)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(ylabels, fontsize=9)
    for i, (m, w) in enumerate(zip(means, wts)):
        ax.text(m + 0.02, i, f"{m:.3f}/{w:.1f}", va="center", fontsize=8)
    ax.axvline(0, color="gray", ls=":", lw=1.0, label="$d_k=0$ (perfect term)")
    ax.set_xlabel("Deviation $d_k=w_k(1-f_k)$  (0=perfect, $w_k$=worst)")
    ax.set_title(f"11-term deviation at optimal arch\n"
                 f"(green $d_k<0.4w_k$, orange $<0.7w_k$, red $\\geq0.7w_k$)\n"
                 f"Mean over {len(best_detail)} DGPs × {N_CAL_PATHS} paths  |  "
                 f"$\\sum d_k$={sum(means):.3f}")
    ax.legend(fontsize=8, loc="lower right")

    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run_full_grid_report(hp=None, nr_grid=None, nz_grid=None,
                          save_prefix="esn_hyperparam_search_fullgrid",
                          fast=False):
    """
    End-to-end convenience wrapper producing the complete-grid analogues
    of both the overview and deviation figures (Sec. 6 of the protocol):
      1. Sweep the fixed hyperparameter set over the full (Nr,Nz) grid
         (`evaluate_fixed_arch_over_grid`).
      2. Re-evaluate the resulting best cell with `return_detail=True`
         to get the per-scenario detail (H_hat, ref, inner_params, ...)
         needed for panels 4-6 / the deviation breakdown.
      3. Save the grid-matrix heatmap, the full-grid overview figure, and
         the full-grid deviation figure.
    """
    if fast:
        global N_DGP_PATHS, N_CAL_PATHS, T_SIM, BURN, T_CAL, BURN_CAL, INNER_NM_MAXITER
        N_DGP_PATHS = 2; N_CAL_PATHS = 2; T_SIM = 200; BURN = 50
        T_CAL = 120; BURN_CAL = 30; INNER_NM_MAXITER = 20

    if hp is None:
        hp = {k: v for k, v in OPTIMAL_ARCH.items()
              if k not in ("rough_orientation", "Nr", "Nz", "matrix_seed")}
    if nr_grid is None: nr_grid = NR_GRID
    if nz_grid is None: nz_grid = NZ_GRID

    results, (best_nr, best_nz, best_val) = run_grid_matrix(
        hp=hp, nr_grid=nr_grid, nz_grid=nz_grid, save_prefix=save_prefix, fast=False)

    best_arch = {"rough_orientation": -1.0, "Nr": best_nr, "Nz": best_nz,
                 "matrix_seed": MATRIX_SEED, **hp}
    _, best_detail = evaluate_arch(best_arch, eval_id=999, verbose=True, return_detail=True)

    make_figure_overview_gridmatrix(
        results, nr_grid, nz_grid, best_arch, best_detail,
        save_path=f"{save_prefix}_overview.png")
    make_figure_deviation_gridmatrix(
        results, nr_grid, nz_grid, best_arch, best_detail,
        save_path=f"{save_prefix}_deviation.png")
    print(f"Full-grid overview  -> {save_prefix}_overview.png")
    print(f"Full-grid deviation -> {save_prefix}_deviation.png")
    return results, best_arch, best_detail


if __name__ == "__main__":
    import argparse
    pa = argparse.ArgumentParser(
        description="ESN hyperparameter search via CMA-ES (fBM_persist fix applied)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pa.add_argument("--max_evals", type=int,   default=MAX_OUTER_EVALS)
    pa.add_argument("--popsize",   type=int,   default=POPSIZE)
    pa.add_argument("--nr",   nargs="+", type=int, default=NR_GRID)
    pa.add_argument("--nz",   nargs="+", type=int, default=NZ_GRID)
    pa.add_argument("--out",  type=str,  default="esn_hyperparam_search")
    pa.add_argument("--fast", action="store_true",
                    help="Quick test run (reduced budget)")
    pa.add_argument("--grid_matrix_only", action="store_true",
                    help="Skip the CMA-ES search; instead evaluate the "
                         "reported-optimal (or --out-prefixed best.json) "
                         "hyperparameter set over the full (Nr,Nz) grid "
                         "and save the score-matrix heatmap (Sec. 6.2/6.3 "
                         "of the protocol).")
    pa.add_argument("--full_grid_report", action="store_true",
                    help="Like --grid_matrix_only, but also regenerates "
                         "the 6-panel overview and 2-panel deviation "
                         "figures using the COMPLETE (Nr,Nz) grid instead "
                         "of a partial CMA-ES search history.")
    args = pa.parse_args()
    if args.full_grid_report:
        run_full_grid_report(nr_grid=args.nr, nz_grid=args.nz,
                              save_prefix=f"{args.out}_fullgrid", fast=args.fast)
    elif args.grid_matrix_only:
        hp = None
        best_json = f"{args.out}_best.json"
        if os.path.exists(best_json):
            with open(best_json) as f:
                saved = json.load(f)
            arch = saved["best_arch"]
            hp = {k: arch[k] for k, *_ in HP_SPEC}
            print(f"Using hyperparameters from {best_json}")
        run_grid_matrix(hp=hp, nr_grid=args.nr, nz_grid=args.nz,
                        save_prefix=f"{args.out}_gridmatrix", fast=args.fast)
    else:
        run(max_evals=args.max_evals, popsize=args.popsize,
            nr_grid=args.nr, nz_grid=args.nz,
            save_prefix=args.out, verbose=True, fast=args.fast)
