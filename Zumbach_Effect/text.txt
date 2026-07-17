"""
student_t_fat_tail_v3.py
========================

DATA-GENERATING PROCESS
------------------------
    r_t = sigma * eps_t,   eps_t ~ t_nu   (i.i.d. location-scale Student-t)
    sigma = 0.20 / sqrt(252)              (daily vol = 20% annualised)

True tail index alpha = nu (analytically known).
Grid:  nu in {3, 5, 10, 30}.

THREE MODELS — all use GAUSSIAN Brownian dynamics (their correct specification)
-------------------------------------------------------------------------------
The fat tails in each model must EMERGE from the stochastic volatility dynamics,
NOT from injecting Student-t noise.  The DGP has i.i.d. fat tails; the SV models
have Gaussian innovations but stochastic vol that generates (mild) excess kurtosis.
The calibration score measures how closely each model's Gaussian-SV distribution
matches the Student-t DGP statistics.  The residual gap IS the finding.

1. ESN A2-981003
   Exact replica of the notebook `desn_A2_lowdim_selected_981003_minimal_colab.ipynb`.
   Return:     dX_t = (mu - V_t/2) dt + sigma_t sqrt(dt) * eps_t,  eps_t ~ N(0,1)
   Reservoir:  r_{i,t} = alpha_i r_{i,t-1} + c_i * eps_t  (same Gaussian eps_t)
   Volatility: sigma_t = sqrt(sigma_min^2 + softplus(eta_t)^2)
   This is the EXACT notebook simulate_path loop.

2. Quadratic Rough Heston (QRH)
   Bourgey & Gatheral (2026), SSRN 5239929.
   dS_t / S_t = -sqrt(V_t) dW_t,   V_t = Y_t^2 + c   (pure Brownian, c > 0)
   dY_t = kappa(T-t) sqrt(V_t) dW_t   (SAME Brownian as spot — perfect leverage -1)
   Gamma kernel: kappa(tau) = nu_vol/Gamma(alpha) * tau^{alpha-1} * exp(-lam*tau)
   All innovations are N(0,1).  Euler-Volterra, ring-buffer 252 steps.

3. PDV Guyon–Lekeufack (GL)
   Guyon & Lekeufack (2023), SSRN 4174589; impl. follows arXiv:2507.09412.
   dS_t / S_t = sqrt(V_t) dW_t   (Brownian)
   V_t = (1-beta)*V_bar + beta*[alpha1*F_r(t) + (1-alpha1)*F_v(t)]
   F_t(gamma) = sum_k k^{gamma-0.5} r_{t-k}^2/dt / sum_k k^{gamma-0.5}
   All innovations are N(0,1).

CALIBRATION
-----------
Data-adaptive score S (Eq. 9 in the problem statement):
    S = 2.2 f(H_hat,   c_H,   s_H)
      + 1.4 f(vol,     c_vol, s_vol)
      + 0.6 g(q995,    c_q,   s_q)
      + 0.2 g(maxV,    c_m,   s_m)
      + 0.8 f(V_ACF,   c_V,   s_V)
      + 1.0 f(G_T,     c_GT,  s_GT)
      + 0.8 f(F_T,     c_FT,  s_FT)
      + 1.0 f(Z,       c_Z,   s_Z)
      + 1.1 f(L,       c_L,   s_L)
      + 0.7 f(K,       c_K,   s_K)
      + 0.5 f(A,       c_A,   s_A)
      - 5.0 * Stress

where f(x,c,s) = max(0, 1-|x-c|/s) and g(x,c,s) = max(0, 1-(x-c)_+/s).
Each centre c and tolerance s is the cross-path mean and std of the DGP
observations for the given nu — strictly model-vs-data for all eleven terms.

FIGURES
-------
Figure 1  (returns): histogram (log-y), Hill plot alpha(k), kurtosis bar,
                     QQ-plot vs Gaussian — 4 rows x 4 columns (one per nu).
Figure 2  (summary): Hill alpha vs nu, excess kurtosis vs nu, score S vs nu.

OUTPUT FILES
------------
student_t_v3_returns.png
student_t_v3_summary.png
student_t_v3_protocol.tex
student_t_v3_protocol.pdf   (compiled if pdflatex available)
"""

# ============================================================================
# 0.  Imports and global constants
# ============================================================================

import math, warnings, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats, optimize
from scipy.special import gamma as Gamma

warnings.filterwarnings("ignore")

# ── physical constants (notebook-exact) ─────────────────────────────────────
TRADING_DAYS    = 252.0
TARGET_ANN_VOL  = 0.20
TARGET_DAY_VOL  = TARGET_ANN_VOL / math.sqrt(TRADING_DAYS)   # sigma in notebook
SIGMA_MIN       = 0.01 / math.sqrt(TRADING_DAYS)             # floor (notebook)
MU_S            = 0.0

# ── architecture 981003 (notebook-exact) ────────────────────────────────────
ARCH = dict(
    architecture_id=981003,
    matrix_seed=202695547565,
    n_r=64,  n_z=12,
    H_target=0.08,
    rough_scale=0.40,   z_strength=0.34,   z_readout=0.05,
    even_strength=1.50, linear_strength=0.25, gamma_norm=1.00,
    local_z_strength=0.03, k_sparse=1, zz_scale=0.08,
    sign_prob_neg=0.22, rough_orientation=-1.0,
)

# ── experiment grid ──────────────────────────────────────────────────────────
NU_GRID = [3, 5, 10, 30]

COLORS = {
    "DGP $t_\\nu$":  "#888780",
    "ESN A2-981003": "#1a4d80",
    "QRH":           "#D85A30",
    "PDV-GL":        "#1D9E75",
}
GRAY  = "#888780"
CORAL = "#c0392b"
MODEL_NAMES = ["ESN A2-981003", "QRH", "PDV-GL"]

NU_LABELS = {
    3:  "$\\nu=3$\n(very heavy)",
    5:  "$\\nu=5$\n(equity-like)",
    10: "$\\nu=10$\n(mild tails)",
    30: "$\\nu=30$\n(near-Gaussian)",
}

# Do not display extremely small densities
DENSITY_THRESHOLD = 1e-4  # try 1e-5 or 1e-6 if you want more tails

# ============================================================================
# 1.  Shared utilities: softplus, statistics, score
# ============================================================================

def _inv_sp(y: float) -> float:
    y = max(float(y), 1e-15)
    return y if y > 35 else math.log(math.expm1(y))

def _sp(x: float) -> float:
    """Scalar softplus (notebook-exact)."""
    if x > 35:   return float(x)
    if x < -35:  return math.exp(float(x))
    return math.log1p(math.exp(float(x)))

def _sp_arr(x):
    """Vectorised softplus."""
    x = np.asarray(x, float)
    return np.where(x > 35, x, np.where(x < -35, np.exp(x), np.log1p(np.exp(x))))

def _corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() < 1e-14 or b.std() < 1e-14: return np.nan
    return float(np.corrcoef(a, b)[0, 1])

def compute_statistics(daily_x, daily_var, burn=0):
    """
    Compute all eleven stylized-fact statistics from one path.
    Mirrors the notebook `statistics()` function exactly.

    Parameters
    ----------
    daily_x   : (T,) daily log-returns (post-Ito correction already included)
    daily_var : (T,) daily integrated variance
    burn      : burn-in days to skip (already sliced out before this call)
    """
    x = np.asarray(daily_x, float)
    v = np.asarray(daily_var, float)
    sigma_daily = np.sqrt(np.maximum(v, 1e-30))
    lv          = np.log(v + 1e-30)
    lags        = np.arange(1, 21)

    ret_acf = np.array([_corr(x[:-L], x[L:]) for L in lags])
    abs_acf = np.array([_corr(np.abs(x[:-L]), np.abs(x[L:])) for L in lags])
    sq_acf  = np.array([_corr(x[:-L]**2, x[L:]**2) for L in lags])
    vol_acf = np.array([_corr(v[:-L], v[L:]) for L in lags])
    lev     = np.array([_corr(x[:-L], v[L:]) for L in lags])

    # Hurst — dyadic variogram (notebook-exact)
    xs, ys = [], []
    for L in [1, 2, 4, 8, 16, 32, 64]:
        if L < len(lv) // 4:
            d = lv[L:] - lv[:-L]; vv = float(np.mean(d*d))
            if vv > 1e-30 and np.isfinite(vv):
                xs.append(math.log(L)); ys.append(math.log(vv))
    H = 0.5 * float(np.polyfit(xs, ys, 1)[0]) if len(xs) >= 3 else np.nan

    # Zumbach — windowed integrated variance (notebook-exact)
    zvals = []
    for L in [5, 10, 20]:
        n = len(x) - 2*L
        if n > 30:
            csx = np.concatenate([[0.], np.cumsum(x)])
            csv = np.concatenate([[0.], np.cumsum(v)])
            idx = L + np.arange(n)
            pR = csx[idx] - csx[idx-L]; fR = csx[idx+L] - csx[idx]
            pV = (csv[idx] - csv[idx-L]) / L
            fV = (csv[idx+L] - csv[idx]) / L
            zvals.append(_corr(pR**2, fV) - _corr(pV, fR**2))

    c   = x - x.mean(); var = float(np.var(c))
    kurt = float(np.mean(c**4) / (var**2 + 1e-30) - 3.)
    ann  = math.sqrt(TRADING_DAYS)

    return dict(
        H_hat        = float(H),
        mean_vol_ann = float(sigma_daily.mean()) * ann,
        q995_vol_ann = float(np.quantile(sigma_daily, 0.995)) * ann,
        max_vol_ann  = float(sigma_daily.max()) * ann,
        mean_vol_acf = float(np.nanmean(vol_acf)),
        taylor_gap   = float(np.nanmean(abs_acf) - np.nanmean(sq_acf)),
        taylor_frac  = float(np.nanmean(abs_acf > sq_acf)),
        zumbach      = float(np.nanmean(zvals)),
        leverage     = float(np.nanmean(lev)),
        kurtosis     = kurt,
        max_ret_acf  = float(np.nanmax(np.abs(ret_acf))),
        logret       = x,       # kept for pooling
    )


def make_score_ref(dgp_stats: list) -> dict:
    """
    Build data-adaptive calibration centres and tolerances from DGP observations.
    For every statistic k:
        c_k = cross-path mean of DGP k
        s_k = max( cross-path std of DGP k,  0.30 * |c_k|,  floor )
    """
    def _m(k):  return float(np.nanmean([s[k] for s in dgp_stats]))
    def _sd(k): return float(np.nanstd( [s[k] for s in dgp_stats]))

    ref = {}
    floors = dict(
        H_hat=0.01, mean_vol_ann=0.01, q995_vol_ann=0.05, max_vol_ann=0.10,
        mean_vol_acf=0.005, taylor_gap=0.002, taylor_frac=0.05,
        zumbach=0.005, leverage=0.005, kurtosis=0.20, max_ret_acf=0.005)
    for k, fl in floors.items():
        c = _m(k); s = _sd(k)
        ref[k + "_c"] = c
        ref[k + "_s"] = max(s, abs(c) * 0.30, fl)
    ref["stress_mx"] = max(ref["max_vol_ann_c"] * 2.0, 1.5)
    return ref


def score_fn(st: dict, ref: dict) -> float:
    """
    Data-adaptive stylized-fact score S (Eq. 9).
    f(x,c,s) = max(0, 1 - |x-c|/s)   symmetric proximity
    g(x,c,s) = max(0, 1 - (x-c)_+/s) one-sided (penalises excess only)
    """
    def f(x, c, s): return max(0., 1. - abs(x - c) / max(s, 1e-8))
    def g(x, c, s): return max(0., 1. - max(0., x - c) / max(s, 1e-8))

    H    = st['H_hat'];          vol  = st['mean_vol_ann']
    q995 = st['q995_vol_ann'];   mx   = st['max_vol_ann']
    V    = st['mean_vol_acf'];   GT   = st['taylor_gap']
    FT   = st['taylor_frac'];    Z    = st['zumbach']
    L    = st['leverage'];       K    = st['kurtosis']
    A    = st['max_ret_acf']

    stress = int(mx > ref['stress_mx'] or vol < 0.05 or vol > 1.50)

    s  = 2.2 * f(H,   ref['H_hat_c'],        ref['H_hat_s'])
    s += 1.4 * f(vol, ref['mean_vol_ann_c'],  ref['mean_vol_ann_s'])
    s += 0.6 * g(q995,ref['q995_vol_ann_c'],  ref['q995_vol_ann_s'])
    s += 0.2 * g(mx,  ref['max_vol_ann_c'],   ref['max_vol_ann_s'])
    s += 0.8 * f(V,   ref['mean_vol_acf_c'],  ref['mean_vol_acf_s'])
    s += 1.0 * f(GT,  ref['taylor_gap_c'],    ref['taylor_gap_s'])
    s += 0.8 * f(FT,  ref['taylor_frac_c'],   ref['taylor_frac_s'])
    s += 1.0 * f(Z,   ref['zumbach_c'],        ref['zumbach_s'])
    s += 1.1 * f(L,   ref['leverage_c'],       ref['leverage_s'])
    s += 0.7 * f(K,   ref['kurtosis_c'],       ref['kurtosis_s'])
    s += 0.5 * f(A,   ref['max_ret_acf_c'],    ref['max_ret_acf_s'])
    s -= 5.0 * stress
    return float(s)


def hill(x, k=None):
    """Hill tail-index estimator."""
    sx = np.sort(np.abs(x))[::-1]; n = len(sx)
    k  = min(k or max(30, n // 10), n - 1)
    lr = np.log(sx[:k]) - np.log(sx[k]); lr = lr[lr > 0]
    return 1. / np.mean(lr) if len(lr) > 0 else np.nan

def hill_curve(x, n_pts=40):
    n = len(x); sx = np.sort(np.abs(x))[::-1]
    ks = np.unique(np.round(np.exp(
        np.linspace(np.log(5), np.log(max(n // 5, 6)), n_pts))).astype(int))
    out = []
    for k in ks:
        k = min(k, n - 1); lr = np.log(sx[:k]) - np.log(sx[k]); lr = lr[lr > 0]
        out.append(1. / np.mean(lr) if len(lr) > 0 else np.nan)
    return ks, np.array(out)


# ============================================================================
# 2.  DGP — i.i.d. Student-t
# ============================================================================

def dgp_student_t(nu: float, n_paths: int, T: int, burn: int,
                  sigma: float = None, seed_base: int = 42):
    """
    r_t = sigma * eps_t,  eps_t ~ t_nu.
    The t_nu variate is scaled so that Var(r_t) = sigma^2  (unit-variance t_nu
    requires dividing by sqrt(nu/(nu-2)) for nu > 2).
    Returns: list of stat dicts, pooled dict {logret: array}.
    """
    if sigma is None: sigma = TARGET_DAY_VOL
    # scale factor so that Var(eps_t) = 1
    if nu > 2:
        t_scale = math.sqrt((nu - 2.) / nu)
    else:
        t_scale = 0.5   # nu <= 2: infinite variance; use small scale
    v_const = sigma**2  # deterministic daily variance

    rng_base = np.random.default_rng(seed_base)
    st_all   = []; lr_all = []

    for p in range(n_paths):
        rng = np.random.default_rng(int(rng_base.integers(1 << 31)))
        eps = rng.standard_t(nu, size=T) * sigma * t_scale
        # Tiny noise on variance so ACF computations don't degenerate
        v   = np.full(T, v_const) + rng.normal(0, v_const * 0.005, T)**2
        st  = compute_statistics(eps[burn:], v[burn:])
        st_all.append(st)
        lr_all.append(eps[burn:])

    pool = dict(logret=np.concatenate(lr_all))
    return st_all, pool


# ============================================================================
# 3.  Model A — ESN A2-981003 (notebook-exact)
# ============================================================================

def _build_esn_params(arch=ARCH):
    """Build reservoir matrices — exact replica of notebook `build_architecture`."""
    rng    = np.random.default_rng(int(arch['matrix_seed']))
    n_r, n_z = arch['n_r'], arch['n_z']

    lambdas = np.geomspace(1/3500, 2.0, int(n_r)).astype(float)
    q = lambdas ** (0.5 - arch['H_target'])
    q = q / (np.linalg.norm(q) + 1e-15)
    b = np.sqrt(2.0 * lambdas)
    C = (b[:, None] * b[None, :]) / (lambdas[:, None] + lambdas[None, :])
    q = q / math.sqrt(float(q @ C @ q) + 1e-15)

    az_rates = np.geomspace(1/280, 1/7, n_z).astype(float)
    zz_diag  = rng.uniform(-arch['zz_scale'], arch['zz_scale'], n_z)
    sign_z   = -rng.choice([-1., 1.], n_z,
                            p=[arch['sign_prob_neg'], 1. - arch['sign_prob_neg']])
    fast_idx = np.array([min(n_r-1, n_r//2 + int((n_r//2-1)*j/max(n_z-1,1)))
                          for j in range(n_z)], dtype=int)
    slow_idx = np.array([int((n_r//2-1)*(n_z-1-j)/max(n_z-1,1))
                          for j in range(n_z)], dtype=int)

    b0     = _inv_sp(math.sqrt(max(TARGET_DAY_VOL**2 - SIGMA_MIN**2, 1e-15)))
    kappa0 = arch['rough_orientation'] * arch['rough_scale'] * float(q @ np.sqrt(2.*lambdas))
    return dict(lambdas=lambdas, q=q, az_rates=az_rates, zz_diag=zz_diag,
                sign_z=sign_z, fast_idx=fast_idx, slow_idx=slow_idx,
                b0=b0, kappa0=kappa0)


_ESN_PARAMS = None   # initialised once in run()


def _simulate_esn(seed: int, T: int,
                  dt: float = 1.0, arch=ARCH, params=None):
    """
    ESN A2-981003 — EXACT replica of notebook simulate_path().

    Single Gaussian innovation eps_t drives BOTH the Ito log-return AND the
    reservoir.  This is the notebook's design: the rough bank is driven by the
    normalised martingale innovation eps_t, not a separate draw.

    dX_t = (mu - V_t/2) dt + sigma_t sqrt(dt) eps_t,   eps_t ~ N(0,1)
    r_{i,t} = alpha_i r_{i,t-1} + c_i eps_t              (same eps_t)

    Returns daily_x (T,), daily_var (T,).
    """
    if params is None: params = _ESN_PARAMS
    rng   = np.random.default_rng(int(seed))
    n_r, n_z = arch['n_r'], arch['n_z']

    spd     = int(round(1.0 / dt))
    n_steps = T * spd
    sqrt_dt = math.sqrt(dt)

    alpha_r = np.exp(-params['lambdas'] * dt)
    c_r     = np.sqrt(np.maximum(1. - alpha_r**2, 1e-14))
    az      = np.exp(-params['az_rates'] * dt)
    omaz    = 1. - az

    r = np.zeros(n_r); z = np.zeros(n_z)
    daily_x   = np.zeros(T)
    daily_var = np.zeros(T)
    wz  = arch['z_readout'] / math.sqrt(n_z)
    rc  = arch['rough_orientation'] * arch['rough_scale']

    for step in range(n_steps):
        eps = rng.normal()     # single N(0,1): drives BOTH return and reservoir

        eta = params['b0'] + rc * float(params['q'] @ r) + wz * float(np.sum(z))
        sp  = _sp(eta)
        sig = math.sqrt(SIGMA_MIN**2 + sp**2)
        var = sig * sig

        day = step // spd
        daily_x[day]   += (MU_S - 0.5 * var) * dt + sig * sqrt_dt * eps
        daily_var[day] += var * dt

        # Diagonal innovation-normalised rough update (notebook-exact)
        r = alpha_r * r + c_r * eps
        m = max(0., 1. - arch['gamma_norm'] * np.linalg.norm(r) / math.sqrt(n_r))

        z_old = z.copy()
        for j in range(n_z):
            p1 = r[params['fast_idx'][j]]
            p2 = r[params['slow_idx'][j]]
            ev   = 0.7*p1*p1 + 0.3*p2*p2 - 1.
            lin  = 0.7*p1 + 0.3*p2
            jm   = max(0, j-1); jp = min(n_z-1, j+1)
            lc   = 0.5 * arch['local_z_strength'] * (z_old[jm] + z_old[jp])
            u = (params['sign_z'][j]*p1*p2 + arch['even_strength']*ev
                 - arch['linear_strength']*lin
                 + params['zz_diag'][j]*z_old[j] + lc)
            z[j] = az[j]*z_old[j] + omaz[j]*m*arch['z_strength']*math.tanh(u)

    return daily_x, daily_var


def _esn_loss(b0_delta: float, scale: float,
              ref: dict, n_cal: int, T_cal: int, burn_cal: int, dt: float):
    """
    Calibration objective for ESN (Gaussian dynamics, notebook-exact).
    b0_delta: readout bias offset (shifts mean vol level).
    scale:    post-hoc vol rescaling (linear adjustment).
    """
    global _ESN_PARAMS
    scale = max(scale, 0.05)
    old_b0 = _ESN_PARAMS['b0']
    _ESN_PARAMS['b0'] = old_b0 + b0_delta
    scores = []
    for q in range(n_cal):
        try:
            x, v = _simulate_esn(3000+q, T_cal, dt)
            x    = x * scale
            st   = compute_statistics(x[burn_cal:], v[burn_cal:])
            scores.append(score_fn(st, ref))
        except Exception:
            scores.append(-10.)
    _ESN_PARAMS['b0'] = old_b0
    return -float(np.nanmean(scores))


def calibrate_esn(nu: float, ref: dict,
                  n_cal=4, T_cal=800, burn_cal=200, dt=1.0) -> dict:
    """Calibrate ESN to DGP statistics.  nu passed for API consistency only."""
    print(f"  ESN calibration  nu={nu}  (Gaussian SV) ...", flush=True)
    def obj(p): return _esn_loss(p[0], p[1], ref, n_cal, T_cal, burn_cal, dt)
    res = optimize.minimize(obj, [0.0, 1.0], method="Nelder-Mead",
                            options={"maxiter":120,"xatol":0.05,"fatol":0.05})
    b0d, sc = res.x; sc = max(float(sc), 0.05)
    print(f"    b0_delta={b0d:.4f} scale={sc:.4f} score={-res.fun:.3f}")
    return dict(b0_delta=float(b0d), scale=sc, nu=nu)


def simulate_esn_paths(cal: dict, n_paths: int, T: int, burn: int,
                       dt=1.0) -> list:
    """Simulate ESN paths with Gaussian innovations (notebook-exact)."""
    global _ESN_PARAMS
    old_b0 = _ESN_PARAMS['b0']
    _ESN_PARAMS['b0'] += cal['b0_delta']
    st_all = []
    for p in range(n_paths):
        x, v = _simulate_esn(9000+p, T, dt)
        x = x * cal['scale']
        st_all.append(compute_statistics(x[burn:], v[burn:]))
    _ESN_PARAMS['b0'] = old_b0
    return st_all


# ============================================================================
# 4.  Model B — Quadratic Rough Heston (QRH)
#     Bourgey & Gatheral (2026), SSRN 5239929
#     Euler–Volterra with ring-buffer
# ============================================================================

def _rl_weights(n: int, al: float, dt: float) -> np.ndarray:
    """
    Riemann–Liouville quadrature weights for the Euler–Volterra scheme.
    w_j = dt^alpha / (Gamma(alpha) * alpha) * [(j+1)^alpha - j^alpha]
    (Bourgey-Gatheral §3, discretisation of the Volterra integral).
    """
    j = np.arange(n, dtype=float)
    return dt**al / (Gamma(al) * al) * ((j+1)**al - j**al)


def _simulate_qrh(seed: int, T: int,
                  H: float, nu_vol: float, lam: float, c_frac: float,
                  dt: float = 1.0, window: int = 252):
    """
    QRH — Bourgey & Gatheral (2026), SSRN 5239929.

    EXACT model dynamics (§2):
        dS_t / S_t = -sqrt(V_t) dW_t        (spot; SAME Brownian as vol)
        V_t = Y_t^2 + c,   c > 0
        dY_t(u) = kappa(u-t) sqrt(V_t) dW_t

    ALL innovations are N(0,1) Brownian increments.
    The Euler-Volterra scheme (§3) uses one Gaussian per step.

    Gamma kernel:
        kappa(tau) = nu_vol/Gamma(alpha) * tau^{alpha-1} * exp(-lam*tau)
        alpha = H + 1/2

    RL quadrature weights:
        w_j = dt^alpha / (Gamma(alpha)*alpha) * [(j+1)^alpha - j^alpha]

    Ring buffer of `window` steps prevents unbounded Volterra accumulation.
    """
    al     = H + 0.5
    nu_hat = nu_vol * math.sqrt(Gamma(2.*H)) / Gamma(al)
    xi0    = TARGET_DAY_VOL**2
    c      = c_frac * xi0
    Y0     = math.sqrt(max(xi0 - c, 0.))

    spd     = int(round(1.0 / dt))
    n_steps = T * spd
    sqrt_dt = math.sqrt(dt)

    rng       = np.random.default_rng(int(seed))
    daily_x   = np.zeros(T)
    daily_var = np.zeros(T)

    buf = np.zeros(window, dtype=float)  # ring buffer for Volterra integrand
    bi  = 0; nf = 0

    for step in range(n_steps):
        nw   = min(nf + 1, window)
        w    = _rl_weights(nw, al, dt)
        idxs = [(bi - 1 - k) % window for k in range(nw)]
        Y    = Y0 + nu_hat * float(np.dot(w, buf[idxs]))
        V    = max(min(Y*Y + c, xi0*50.), 1e-12)
        sig  = math.sqrt(V)

        # Single N(0,1) Brownian increment — drives BOTH spot and Y
        eps  = rng.normal()

        day = step // spd
        # dS_t/S_t = -sqrt(V_t) dW_t  (paper eq.)
        # Euler on S: S_{n+1} = S_n * (1 - sqrt(V_n*dt)*eps)
        # Daily log-return: ln(S_{n+1}/S_n) ≈ -sqrt(V_n*dt)*eps  (no Ito drift)
        daily_x[day]   -= sig * sqrt_dt * eps
        daily_var[day] += V * dt

        # Volterra integrand stored for RL convolution
        buf[bi % window] = sig * eps
        bi += 1; nf = min(nf + 1, window)

    return daily_x, daily_var


def _qrh_loss(params, ref, n_cal, T_cal, burn_cal, dt):
    H, nu_vol, lam, c_frac = params
    if not (0.01<H<0.49 and 0.01<nu_vol<3. and 0.1<lam<20. and 0.01<c_frac<0.99):
        return 20.
    scores = []
    for q in range(n_cal):
        try:
            x, v = _simulate_qrh(500+q, T_cal, H, nu_vol, lam, c_frac, dt)
            scores.append(score_fn(compute_statistics(x[burn_cal:], v[burn_cal:]), ref))
        except Exception:
            scores.append(-10.)
    return -float(np.nanmean(scores))


def calibrate_qrh(nu: float, ref: dict,
                  n_cal=4, T_cal=800, burn_cal=200, dt=1.0) -> dict:
    """Calibrate QRH (Gaussian SV) to DGP statistics."""
    print(f"  QRH calibration  nu={nu}  (Gaussian SV, dS/S=-sqrt(V)dW) ...", flush=True)
    tH = ref['H_hat_c']
    p0 = [max(0.02, min(tH, 0.48)), 0.30, 2.0, 0.50]
    def obj(p): return _qrh_loss(p, ref, n_cal, T_cal, burn_cal, dt)
    res = optimize.minimize(obj, p0, method="Nelder-Mead",
                            options={"maxiter":300,"xatol":0.01,"fatol":0.02})
    H_o, nv_o, lk_o, cf_o = res.x
    H_o  = float(np.clip(H_o,  0.01, 0.48))
    nv_o = float(np.clip(nv_o, 0.01, 3.))
    lk_o = float(np.clip(lk_o, 0.1,  20.))
    cf_o = float(np.clip(cf_o, 0.01, 0.99))
    print(f"    H={H_o:.4f} nu_vol={nv_o:.4f} lam={lk_o:.4f} c={cf_o:.4f} "
          f"score={-res.fun:.3f}")
    return dict(H=H_o, nu_vol=nv_o, lam=lk_o, c_frac=cf_o)


def simulate_qrh_paths(cal: dict, n_paths: int, T: int, burn: int,
                       dt=1.0) -> list:
    """Simulate QRH paths with Gaussian innovations."""
    return [compute_statistics(
                *[arr[burn:] for arr in _simulate_qrh(
                    7000+p, T,
                    cal['H'], cal['nu_vol'], cal['lam'], cal['c_frac'], dt)])
            for p in range(n_paths)]


# ============================================================================
# 5.  Model C — PDV Guyon–Lekeufack (GL)
#     Guyon & Lekeufack (2023), SSRN 4174589; arXiv:2507.09412
#
#  EXACT MODEL (from the paper):
#     dS_t = sigma_t dW_t                    (arithmetic Brownian, not geometric)
#     sigma_t = sigma(R_{1,t}, R_{2,t})
#             = beta0 + beta1 * R_{1,t} + beta2 * sqrt(R_{2,t})
#
#     R_{1,t} = int_{-inf}^{t} K1(t-u) sigma_u dW_u   (kernel-weighted vol-of-vol)
#     R_{2,t} = int_{-inf}^{t} K2(t-u) sigma_u^2 du   (kernel-weighted variance)
#
#     K1(tau) = (1/theta1) * exp(-tau/theta1)    (exponential kernel)
#     K2(tau) = (1/theta2) * exp(-tau/theta2)
#
#  Discretisation (arXiv:2507.09412):
#     R_{1,n+1} = exp(-dt/theta1) * R_{1,n} + sigma_n * sqrt(dt) * eps_n
#     R_{2,n+1} = exp(-dt/theta2) * R_{2,n} + sigma_n^2 * dt
#     r_n       = sigma_n * sqrt(dt) * eps_n    (return, NO Ito correction)
#     sigma_n   = max( beta0 + beta1*R_{1,n} + beta2*sqrt(R_{2,n}), sigma_min )
# ============================================================================

def _simulate_gl(seed: int, T: int,
                 beta0: float, beta1: float, beta2: float,
                 theta1: float, theta2: float,
                 dt: float = 1.0):
    """
    PDV Guyon-Lekeufack (2023) — exact paper specification.

    MODEL:
        dS_t = sigma_t dW_t    (arithmetic Brownian — daily return = sigma_t * sqrt(dt) * eps)
        sigma_t = beta0 + beta1 * R1_t + beta2 * sqrt(R2_t)
        R1_t = int K1(t-u) sigma_u dWu    (exponential kernel, decay theta1)
        R2_t = int K2(t-u) sigma_u^2 du   (exponential kernel, decay theta2)

    Euler discretisation:
        R1_{n+1} = a1 * R1_n + sigma_n * sqrt(dt) * eps_n
        R2_{n+1} = a2 * R2_n + sigma_n^2 * dt
        r_n      = sigma_n * sqrt(dt) * eps_n          (daily return, NO Ito term)

    where a1 = exp(-dt/theta1), a2 = exp(-dt/theta2).

    Parameters
    ----------
    beta0  : floor vol (ensures sigma > 0; ~ TARGET_DAY_VOL initially)
    beta1  : loading on past signed-vol integral R1
    beta2  : loading on sqrt of past variance integral R2
    theta1 : mean-reversion time of R1 (days)
    theta2 : mean-reversion time of R2 (days)
    """
    spd     = int(round(1.0 / dt))
    n_steps = T * spd
    sqrt_dt = math.sqrt(dt)

    a1 = math.exp(-dt / max(theta1, 1e-6))
    # a2 is computed in stationary init block below

    rng       = np.random.default_rng(int(seed))
    daily_x   = np.zeros(T)   # daily return  dS/S  (arithmetic)
    daily_var = np.zeros(T)   # daily quadratic variation

    # Stationary R2: self-consistent solution sigma0 = beta0 + beta2*sqrt(R2_stat)
    # R2_stat = sigma0^2 * dt / (1-a2) => sigma0 = beta0 / (1 - beta2*sqrt(dt/(1-a2)))
    a2    = math.exp(-dt / max(theta2, 1e-6))
    denom = max(1. - beta2 * math.sqrt(dt / max(1. - a2, 1e-10)), 0.1)
    sig0  = beta0 / denom
    sig0  = min(sig0, TARGET_DAY_VOL * 5.)   # safety cap
    R1 = 0.
    R2 = sig0**2 * dt / max(1. - a2, 1e-10)

    for step in range(n_steps):
        # Volatility function: sigma(R1, R2) — clipped for stability
        sig = max(beta0 + beta1 * R1 + beta2 * math.sqrt(max(R2, 0.)), 1e-6)
        sig = min(sig, TARGET_DAY_VOL * 10.)  # hard cap at 10x target vol

        eps = rng.normal()
        r_t = sig * sqrt_dt * eps   # dS_t = sigma_t dW_t  =>  r_t = sig * sqrt(dt) * eps

        day = step // spd
        daily_x[day]   += r_t       # arithmetic return (no 1/2 Ito term)
        daily_var[day] += sig**2 * dt

        # OU updates for the two kernel factors
        R1 = a1 * R1 + sig * sqrt_dt * eps   # R1: weighted past signed vol shocks
        R2 = a2 * R2 + sig**2 * dt           # R2: weighted past variance

    return daily_x, daily_var


def _gl_loss(params, ref, n_cal, T_cal, burn_cal, dt):
    beta0, beta1, beta2, theta1, theta2 = params
    # Stability: require beta2 * sqrt(theta2) < 0.85 to keep sigma mean-reverting
    if not (1e-5 < beta0 < 0.05 and abs(beta1) < 1. and 0. <= beta2 < 1.
            and 1. <= theta1 < 100. and 5. <= theta2 < 500.
            and beta2 * math.sqrt(theta2) < 0.85):
        return 20.
    scores = []
    for q in range(n_cal):
        try:
            x, v = _simulate_gl(600+q, T_cal, beta0, beta1, beta2,
                                  theta1, theta2, dt)
            scores.append(score_fn(compute_statistics(x[burn_cal:], v[burn_cal:]), ref))
        except Exception:
            scores.append(-10.)
    return -float(np.nanmean(scores))


def calibrate_gl(nu: float, ref: dict,
                 n_cal=4, T_cal=800, burn_cal=200, dt=1.0) -> dict:
    """
    Calibrate GL to DGP statistics.
    Parameters: beta0 (floor vol), beta1 (R1 loading), beta2 (R2 loading),
                theta1 (R1 decay, days), theta2 (R2 decay, days).
    """
    print(f"  PDV-GL calibration  nu={nu}  (dS=sigma dW, exact GL) ...", flush=True)
    # Initial guess: beta0 = 80% of target vol, small beta2, moderate memory scales
    v0 = TARGET_DAY_VOL
    p0 = [v0 * 0.8, 0.0, 0.05, 5., 50.]
    def obj(p): return _gl_loss(p, ref, n_cal, T_cal, burn_cal, dt)
    res = optimize.minimize(obj, p0, method="Nelder-Mead",
                            options={"maxiter":500,"xatol":1e-3,"fatol":0.02})
    b0, b1, b2, th1, th2 = res.x
    b0  = float(np.clip(b0,  1e-6, 0.5))
    b1  = float(np.clip(b1,  -5.,  5.))
    b2  = float(np.clip(b2,  0.,   5.))
    th1 = float(np.clip(th1, 1.,   500.))
    th2 = float(np.clip(th2, 1.,   2000.))
    print(f"    beta0={b0:.5f} beta1={b1:.4f} beta2={b2:.4f} "
          f"theta1={th1:.1f} theta2={th2:.1f}  score={-res.fun:.3f}")
    return dict(beta0=b0, beta1=b1, beta2=b2, theta1=th1, theta2=th2)


def simulate_gl_paths(cal: dict, n_paths: int, T: int, burn: int,
                      dt=1.0) -> list:
    """Simulate GL paths using exact paper dynamics."""
    return [compute_statistics(
                *[arr[burn:] for arr in _simulate_gl(
                    8000+p, T,
                    cal['beta0'], cal['beta1'], cal['beta2'],
                    cal['theta1'], cal['theta2'], dt)])
            for p in range(n_paths)]


# ============================================================================
# 6.  Figures
# ============================================================================

def _pool_ret(st_list):
    return np.concatenate([s['logret'] for s in st_list])


def plot_returns_grid(results: dict, nu_grid: list,
                      save_path: str = None) -> plt.Figure:
    """
    Figure 1 — fat-tails of daily log-returns.
    4 rows x len(nu_grid) columns.
    Rows: histogram | Hill plot | kurtosis bar | QQ vs Gaussian.
    """
    n_nu  = len(nu_grid)
    fig, axes = plt.subplots(4, n_nu, figsize=(5*n_nu, 16))
    if n_nu == 1: axes = axes[:, None]
    fig.suptitle(
        "Fat-Tailedness of Daily Log-Returns\n"
        "DGP = i.i.d. Student-$t_\\nu$  |  "
        "True tail index $\\alpha = \\nu$  |  "
        "Three models calibrated by data-adaptive score $S$",
        fontsize=12, fontweight="bold")

    for ci, nu in enumerate(nu_grid):
        data    = results[nu]
        dgp_r   = data['dgp_pool']['logret']
        all_r   = np.concatenate([dgp_r] +
                                  [_pool_ret(data[m]) for m in MODEL_NAMES])
        lo, hi  = np.percentile(all_r, 0.3), np.percentile(all_r, 99.7)
        axes[0, ci].set_title(NU_LABELS[nu], fontsize=11, fontweight="bold")

        # ── Row 0: histogram (log-y) ────────────────────────────────
        ax = axes[0, ci]
        ax.hist(dgp_r, bins=120, density=True, histtype="step",
                lw=2.2, color=GRAY, label="DGP $t_\\nu$", zorder=5)
        xr = np.linspace(lo, hi, 400)
        ax.plot(xr, stats.norm.pdf(xr, 0, dgp_r.std()),
                color=GRAY, lw=1., ls=":", label="Gaussian ref")
        for m in MODEL_NAMES:
            ax.hist(_pool_ret(data[m]), bins=120, density=True,
                    histtype="step", lw=1.5, color=COLORS[m], label=m)

        ax.set_yscale("log")
        ax.set_ylim(bottom=DENSITY_THRESHOLD)
        ax.set_yscale("log"); ax.set_xlim(lo, hi)
        ax.set_xlabel("$r_t$", fontsize=9)
        if ci == 0: ax.set_ylabel("Density (log)", fontsize=9)
        ax.legend(fontsize=7); ax.grid(True, alpha=0.18)

        # ── Row 1: Hill plot ─────────────────────────────────────────
        ax = axes[1, ci]
        ks, al = hill_curve(dgp_r)
        ax.plot(ks, al, color=GRAY, lw=2.2, label="DGP $t_\\nu$")
        for m in MODEL_NAMES:
            ks2, al2 = hill_curve(_pool_ret(data[m]))
            ax.plot(ks2, al2, color=COLORS[m], lw=1.5, label=m)
        ax.axhline(nu, color=CORAL, lw=1.4, ls="--",
                   label=f"True $\\alpha={nu}$")
        ax.axhline(4,  color=GRAY,  lw=0.8, ls=":", alpha=0.7)
        ax.set_ylim(0, min(nu*3, 40))
        ax.set_xlabel("$k$", fontsize=9)
        if ci == 0: ax.set_ylabel(r"$\hat\alpha$ Hill index", fontsize=9)
        ax.legend(fontsize=7); ax.grid(True, alpha=0.18)

        # ── Row 2: kurtosis bar ──────────────────────────────────────
        ax = axes[2, ci]
        labels = ["DGP $t_\\nu$"] + MODEL_NAMES
        bar_cols = [GRAY] + [COLORS[m] for m in MODEL_NAMES]
        dgp_k   = float(np.nanmean([s['kurtosis'] for s in data['dgp_stats']]))
        kvals   = [dgp_k] + [float(np.nanmean([s['kurtosis'] for s in data[m]]))
                              for m in MODEL_NAMES]
        bars = ax.bar(labels, kvals, color=bar_cols, alpha=0.82, width=0.6)
        if nu > 4:
            theo_k = 6. / (nu - 4.)
            ax.axhline(theo_k, color=CORAL, lw=1.4, ls="--",
                       label=f"Theoretical $6/(\\nu-4)={theo_k:.2f}$")
        ax.axhline(0, color=GRAY, lw=0.6)
        ax.set_ylabel("Excess kurtosis", fontsize=9)
        ax.tick_params(axis='x', labelsize=7, rotation=20)
        ax.legend(fontsize=7); ax.grid(True, axis="y", alpha=0.18)
        for bar, v in zip(bars, kvals):
            ax.text(bar.get_x()+bar.get_width()/2, max(v,0)+0.1,
                    f"{v:.1f}", ha="center", fontsize=7)

        # ── Row 3: QQ-plot vs Gaussian ───────────────────────────────
        ax = axes[3, ci]
        for r_arr, lbl, col, lw_ in (
                [(dgp_r, "DGP $t_\\nu$", GRAY, 2.2)] +
                [(_pool_ret(data[m]), m, COLORS[m], 1.5)
                 for m in MODEL_NAMES]):
            z  = (r_arr - r_arr.mean()) / (r_arr.std() + 1e-14)
            qs = np.linspace(0.01, 0.99, 200)
            ax.plot(stats.norm.ppf(qs), np.quantile(z, qs),
                    color=col, lw=lw_, label=lbl)
        ax.plot([-4, 4], [-4, 4], color=GRAY, lw=0.8, ls=":",
                label="Gaussian $y=x$")
        ax.set_xlabel("Gaussian theoretical quantile", fontsize=9)
        if ci == 0: ax.set_ylabel("Empirical quantile", fontsize=9)
        ax.set_title("QQ vs Gaussian", fontsize=9)
        ax.legend(fontsize=7); ax.grid(True, alpha=0.18)

    for ax in axes.ravel(): ax.tick_params(labelsize=8)
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved -> {save_path}")
    return fig


def plot_summary(results: dict, nu_grid: list,
                 save_path: str = None) -> plt.Figure:
    """
    Figure 2 — summary across all nu values.
    3 rows: Hill alpha | excess kurtosis | score S.
    """
    fig, axes = plt.subplots(3, 1, figsize=(9, 12))
    fig.suptitle(
        "Summary: Tail Properties vs Degrees of Freedom $\\nu$\n"
        "DGP = i.i.d. Student-$t_\\nu$  |  "
        "Models calibrated by data-adaptive score $S$",
        fontsize=12, fontweight="bold")
    marks = ["o", "s", "^", "D"]

    # Row 0: Hill alpha of returns vs nu
    ax = axes[0]
    dgp_hill = [hill(_pool_ret(results[nu]['dgp_stats'])) for nu in nu_grid]
    ax.plot(nu_grid, nu_grid,   color=GRAY,  lw=2.,  marker="o", ms=8,
            label="DGP $t_\\nu$ (true $\\alpha=\\nu$)")
    ax.plot(nu_grid, dgp_hill,  color=GRAY,  lw=1.4, marker="o", ms=6,
            ls="--", label="DGP $t_\\nu$ Hill $\\hat\\alpha$")
    for m, mk in zip(MODEL_NAMES, marks[1:]):
        vals = [hill(_pool_ret(results[nu][m])) for nu in nu_grid]
        ax.plot(nu_grid, vals, color=COLORS[m], lw=1.8, marker=mk, ms=7, label=m)
    ax.set_xticks(nu_grid); ax.set_xlabel("$\\nu$", fontsize=10)
    ax.set_ylabel(r"Hill $\hat\alpha$ of $|r_t|$ at $k=n/10$", fontsize=10)
    ax.set_title(r"Tail index vs $\nu$  (closer to $\alpha=\nu$ = better)", fontsize=11)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.22)

    # Row 1: excess kurtosis vs nu
    ax = axes[1]
    theo_k = [6./(nu-4.) if nu > 4 else float('nan') for nu in nu_grid]
    ax.plot(nu_grid, theo_k, color=GRAY, lw=2., marker="o", ms=8,
            label="Theoretical $6/(\\nu-4)$ [DGP]")
    for m, mk in zip(MODEL_NAMES, marks[1:]):
        vals = [float(np.nanmean([s['kurtosis'] for s in results[nu][m]]))
                for nu in nu_grid]
        ax.plot(nu_grid, vals, color=COLORS[m], lw=1.8, marker=mk, ms=7, label=m)
    ax.set_xticks(nu_grid); ax.set_xlabel("$\\nu$", fontsize=10)
    ax.set_ylabel("Excess kurtosis of $r_t$", fontsize=10)
    ax.set_title("Excess kurtosis vs $\\nu$", fontsize=11)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.22)

    # Row 2: score S vs nu
    ax = axes[2]
    for m, mk in zip(MODEL_NAMES, marks[1:]):
        vals = [float(np.nanmean([score_fn(s, results[nu]['ref']) for s in results[nu][m]]))
                for nu in nu_grid]
        ax.plot(nu_grid, vals, color=COLORS[m], lw=1.8, marker=mk, ms=7, label=m)
    dgp_sc = [float(np.nanmean([score_fn(s, results[nu]['ref'])
                                 for s in results[nu]['dgp_stats']]))
              for nu in nu_grid]
    ax.plot(nu_grid, dgp_sc, color=GRAY, lw=2., marker="o", ms=8,
            ls="--", label="DGP $t_\\nu$ (self-score)")
    ax.set_xticks(nu_grid); ax.set_xlabel("$\\nu$", fontsize=10)
    ax.set_ylabel("Mean score $S$ (model vs DGP)", fontsize=10)
    ax.set_title(r"Score $S$ vs $\nu$  (higher = closer to DGP)", fontsize=11)
    ax.axhline(0, color=GRAY, lw=0.6)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.22)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved -> {save_path}")
    return fig


# ============================================================================
# 7.  LaTeX protocol
# ============================================================================

LATEX = r"""
\documentclass[11pt,a4paper]{article}
\usepackage{amsmath,amssymb,amsthm,amsfonts}
\usepackage{geometry}\geometry{left=2.3cm,right=2.3cm,top=2.5cm,bottom=2.5cm}
\usepackage{setspace}\setstretch{1.20}
\usepackage{hyperref}\hypersetup{colorlinks=true,linkcolor=black,citecolor=black}
\usepackage[T1]{fontenc}
\usepackage{parskip,booktabs,array,xcolor,titlesec,enumitem}
\titleformat{\section}{\large\bfseries}{}{0em}{}[\titlerule]
\titlespacing{\section}{0pt}{14pt}{6pt}
\titleformat{\subsection}{\normalsize\bfseries}{\thesubsection\;}{0em}{}
\newtheorem{remark}{Remark}

\begin{document}

\begin{center}
{\LARGE\bfseries Numerical Protocol}\\[0.3cm]
{\large Fat-Tailedness of Returns under the\\
i.i.d.\ Student-$t$ Data-Generating Process}\\[0.2cm]
{\normalsize Othmane Zarhali --- Paris Dauphine / CNRS}
\end{center}
\vspace{0.15cm}\noindent\rule{\textwidth}{0.8pt}

\begin{abstract}
We evaluate how well three stochastic volatility models --- ESN A2-981003,
Quadratic Rough Heston (QRH), and the path-dependent volatility model of
Guyon--Lekeufack (PDV-GL) --- can reproduce the fat-tail properties of daily
log-returns under an i.i.d.\ location-scale Student-$t$ data-generating process
(DGP) with analytically known tail index $\alpha=\nu\in\{3,5,10,30\}$.

Crucially, all three models use \textbf{Gaussian Brownian innovations} --- this
is their correct mathematical specification (ESN: $\varepsilon_t\sim\mathcal{N}(0,1)$;
QRH: $dS/S_t=-\sqrt{V_t}\,dW_t$; GL: $dS/S_t=\sqrt{V_t}\,dW_t$).
Fat tails in the model returns must \emph{emerge} from the stochastic volatility
dynamics (variance clustering, rough memory, path-dependence) rather than from
heavy-tailed noise.  The residual gap between the model's emergent kurtosis and
the true Student-$t$ tail index is the central finding of this protocol.
\end{abstract}

\tableofcontents\vspace{0.4cm}

%----------------------------------------------------------------------
\section{Data-Generating Process}
%----------------------------------------------------------------------

The DGP is the i.i.d.\ location-scale Student-$t$ model:
\begin{equation}
  r_t \;=\; \sigma\,\varepsilon_t, \qquad
  \varepsilon_t \;\sim\; t_\nu, \qquad
  \sigma \;=\; \frac{0.20}{\sqrt{252}},
  \label{eq:dgp}
\end{equation}
where $t_\nu$ denotes the standardised Student-$t$ distribution.
The scale factor is $\sigma = 0.20/\sqrt{252}$ (daily vol $=20\%$ annual).
For $\nu>2$ the innovations are further rescaled by $\sqrt{(\nu-2)/\nu}$
so that $\mathrm{Var}(r_t) = \sigma^2$ exactly.

\paragraph{Analytical tail index.}  By the regular variation of $t_\nu$:
\begin{equation}
  P(|r_t|>u) \;\sim\; C\,u^{-\nu} \quad\text{as }u\to\infty,
  \qquad \alpha \;=\; \nu.
  \label{eq:tail}
\end{equation}

\paragraph{Excess kurtosis.}
For $\nu>4$: $\kappa=6/(\nu-4)$; infinite for $\nu\le4$.

\paragraph{Parameter grid.}

\begin{center}
\renewcommand{\arraystretch}{1.4}
\begin{tabular}{cccc}\toprule
$\nu$ & Tail behaviour & Excess kurtosis & Context \\\midrule
3  & Very heavy ($\alpha=3$) & $\infty$ & Stress scenario \\
5  & Heavy ($\alpha=5$)      & 6.00     & Typical equity  \\
10 & Mild ($\alpha=10$)      & 1.00     & Diversified portfolio \\
30 & Near-Gaussian           & 0.24     & Large-cap index \\
\bottomrule\end{tabular}
\end{center}

%----------------------------------------------------------------------
\section{Models}
%----------------------------------------------------------------------

\subsection{ESN A2-981003}

The architecture is the low-dimensional A2 candidate identified in
Zarhali (2026). It is a block-triangular continuous-time Deep ESN / CDE
with latent dimension $N_r+N_z=64+12=76$.

\paragraph{Rough memory bank.}
$r_t\in\mathbb{R}^{64}$ with diagonal OU spectrum:
\begin{equation}
  r_{i,t+\Delta t} \;=\; e^{-\lambda_i\Delta t}\,r_{i,t}
  + \sqrt{1-e^{-2\lambda_i\Delta t}}\;\varepsilon_t^G,
  \qquad \varepsilon_t^G\sim\mathcal{N}(0,1),
  \quad \lambda_i\in\mathrm{Geomspace}(1/3500,\,2),
  \label{eq:rough_bank}
\end{equation}
so each node is an OU process driven by the \emph{same} scalar Gaussian
innovation $\varepsilon_t^G$ (no cross-talk, innovation-normalised).

\paragraph{Interaction bank.}
$z_t\in\mathbb{R}^{12}$ with $a_j\in\mathrm{Geomspace}(1/280,1/7)$:
\begin{equation}
  z_{j,t+\Delta t}
  = e^{-a_j\Delta t}\,z_{j,t}
  + (1-e^{-a_j\Delta t})\,m_t\cdot 0.34\,\tanh(u_{j,t}),
  \label{eq:z_update}
\end{equation}
\begin{equation}
  u_{j,t} = s_j\,p_{1,j}p_{2,j}
    + 1.50\underbrace{(0.7p_{1,j}^2+0.3p_{2,j}^2-1)}_{\text{weighted even response}}
    - 0.25(0.7p_{1,j}+0.3p_{2,j})
    + d_j\,z_{j,t},
\end{equation}
where $p_{1,j}=r_{f_j,t}$ (fast projection), $p_{2,j}=r_{s_j,t}$ (slow
projection), $m_t=\max(0,1-\|r_t\|/\sqrt{N_r})$ (damping).

\paragraph{Volatility readout.}
\begin{equation}
  \eta_t = b_0 - 0.40\,q^\top r_t + \frac{0.05}{\sqrt{12}}\sum_{j=1}^{12}z_{j,t},
  \qquad
  \sigma_t = \sqrt{\sigma_{\min}^2 + \mathrm{softplus}(\eta_t)^2}.
  \label{eq:esn_readout}
\end{equation}

\paragraph{Return (Gaussian Brownian innovation — notebook-exact).}
\begin{equation}
  \Delta X_t = \Bigl(\mu-\tfrac12\sigma_t^2\Bigr)\Delta t
              + \sigma_t\sqrt{\Delta t}\,\varepsilon_t,
  \qquad \varepsilon_t \;\sim\; \mathcal{N}(0,1).
  \label{eq:esn_return}
\end{equation}
The \emph{same} $\varepsilon_t$ drives both the Ito log-return and the rough
bank update $r_{i,t}\leftarrow\alpha_i r_{i,t}+c_i\varepsilon_t$ (notebook-exact).
The fat-tail properties of the model emerge entirely from the stochastic
volatility dynamics, not from the innovation distribution.

\paragraph{Admissibility.}
The instantaneous leverage coefficient $\kappa_0=-0.472416<0$ is verified
analytically (Section 4 of the architecture report).

\subsection{Quadratic Rough Heston (QRH)}

Bourgey \& Gatheral \cite{bourgey2026}, SSRN 5239929.

\paragraph{Model.}
\begin{equation}
  V_t \;=\; Y_t^2 + c, \qquad c>0,
  \label{eq:qrh_var}
\end{equation}
\begin{equation}
  \mathrm{d}Y_t(u) \;=\; \kappa(u-t)\sqrt{V_t}\,\mathrm{d}W_t^{\rm vol},
  \label{eq:qrh_dy}
\end{equation}
with the gamma kernel:
\begin{equation}
  \kappa(\tau)
  = \frac{\nu_{\rm vol}}{\Gamma(\alpha)}\,\tau^{\alpha-1}\,e^{-\lambda\tau},
  \qquad \alpha = H+\tfrac12,
  \label{eq:qrh_kernel}
\end{equation}
and calibrated parameters $(H,\nu_{\rm vol},\lambda,c)$.

\paragraph{Simulation.}
Euler--Volterra discretisation (Bourgey--Gatheral, \S3) with
Riemann--Liouville quadrature weights:
\begin{equation}
  w_j = \frac{\Delta t^\alpha}{\Gamma(\alpha)\,\alpha}
        \bigl[(j+1)^\alpha - j^\alpha\bigr],
  \qquad j=0,\ldots,N-1.
  \label{eq:rl_weights}
\end{equation}
A ring buffer of 252 steps (one trading year at $\Delta t=1$) prevents
unbounded accumulation of the Volterra sum while preserving the one-year
fractional memory.

\paragraph{Simulation of returns.}
The model SDE is $\mathrm{d}S_t/S_t = -\sqrt{V_t}\,\mathrm{d}W_t$.
In the Euler scheme the daily return is:
\begin{equation}
  r_n \;=\; \frac{S_{n+1}-S_n}{S_n}
       \;\approx\; -\sqrt{V_n\,\Delta t}\,\varepsilon_n,
  \qquad \varepsilon_n\sim\mathcal{N}(0,1).
  \label{eq:qrh_return}
\end{equation}
No Itô correction is added — the model works with the price $S_t$ directly.
Excess kurtosis emerges from variance clustering ($V_t = Y_t^2+c$ amplified by
the rough kernel).

\subsection{PDV Guyon--Lekeufack (GL)}

Guyon \& Lekeufack \cite{guyon2023}, SSRN 4174589;
implementation follows \cite{bourgey2025}, arXiv:2507.09412.

\paragraph{Model.}
\begin{equation}
  V_t = (1-\beta)\bar V
        + \beta\bigl[\alpha_1 F_t(\alpha_r) + (1-\alpha_1)F_t(\alpha_v)\bigr],
  \label{eq:gl_var}
\end{equation}
\begin{equation}
  F_t(\gamma)
  = \frac{\displaystyle\sum_{k=1}^{t} k^{\gamma-1/2}\,r_{t-k}^2/\Delta t}
         {\displaystyle\sum_{k=1}^{t} k^{\gamma-1/2}},
  \label{eq:gl_factor}
\end{equation}
where $\alpha_r$ (short-memory, rough factor) and $\alpha_v$ (long-memory
factor) are calibrated exponents.
The factor $F_t(\gamma)$ is computed online from all available history at each step.

\paragraph{Simulation of returns.}
The model SDE is $\mathrm{d}S_t/S_t = -\sqrt{V_t}\,\mathrm{d}W_t$.
In the Euler scheme the daily return is:
\begin{equation}
  r_n \;=\; \frac{S_{n+1}-S_n}{S_n}
       \;\approx\; -\sqrt{V_n\,\Delta t}\,\varepsilon_n,
  \qquad \varepsilon_n\sim\mathcal{N}(0,1).
  \label{eq:qrh_return}
\end{equation}
No Itô correction is added — the model works with the price $S_t$ directly.
Excess kurtosis emerges from variance clustering ($V_t = Y_t^2+c$ amplified by
the rough kernel).

%----------------------------------------------------------------------
\section{Data-Adaptive Calibration Score}
%----------------------------------------------------------------------

\subsection{Reference construction}

Given $n_{\rm dgp}$ DGP paths, each producing a set of statistics
$\{k^{(p)}\}_{p=1}^{n_{\rm dgp}}$, define for each statistic $k$:
\begin{equation}
  c_k^{\rm DGP} = \frac{1}{n_{\rm dgp}}\sum_p k^{(p)}, \qquad
  s_k^{\rm DGP} = \max\!\Bigl(\mathrm{std}_p(k^{(p)}),\;
                          0.30\,|c_k^{\rm DGP}|,\;
                          s_k^{\min}\Bigr),
  \label{eq:ref}
\end{equation}
where $s_k^{\min}>0$ is a small floor (preventing infinite sensitivity
when the DGP statistic is nearly constant).

\subsection{Score function}

The data-adaptive score is:
\begin{equation}
\begin{aligned}
S &= 2.2\,f\!\bigl(\hat H,\,c_H,s_H\bigr)
   + 1.4\,f\!\bigl(\bar\sigma,\,c_\sigma,s_\sigma\bigr)
   + 0.6\,g\!\bigl(q_{995},\,c_q,s_q\bigr)
   + 0.2\,g\!\bigl(\sigma_{\max},\,c_m,s_m\bigr) \\
  &+ 0.8\,f\!\bigl(\bar V_{\rm ACF},\,c_V,s_V\bigr)
   + 1.0\,f\!\bigl(G_T,\,c_{GT},s_{GT}\bigr)
   + 0.8\,f\!\bigl(F_T,\,c_{FT},s_{FT}\bigr) \\
  &+ 1.0\,f\!\bigl(Z,\,c_Z,s_Z\bigr)
   + 1.1\,f\!\bigl(L,\,c_L,s_L\bigr)
   + 0.7\,f\!\bigl(K,\,c_K,s_K\bigr)
   + 0.5\,f\!\bigl(A,\,c_A,s_A\bigr)
   - 5\cdot\mathrm{Stress},
\end{aligned}
\label{eq:score}
\end{equation}
where
$f(x,c,s) = \max(0,1-|x-c|/s)$ (symmetric proximity) and
$g(x,c,s) = \max(0,1-(x-c)_+/s)$ (one-sided, penalises excess).
All centres $c_k$ and tolerances $s_k$ are the DGP-observed values from
\eqref{eq:ref}, so every term measures \emph{model vs data}.

\paragraph{Statistics in $S$.}
$\hat H$ = Hurst exponent (dyadic variogram on $\log V_t$);
$\bar\sigma$ = mean annualised volatility;
$q_{995}$ = 99.5th quantile of annualised vol;
$\sigma_{\max}$ = path-maximum annualised vol;
$\bar V_{\rm ACF}$ = mean vol ACF at lags 1--20;
$G_T$ = Taylor gap ($\bar\rho_{|r|}-\bar\rho_{r^2}$);
$F_T$ = Taylor fraction;
$Z$ = Zumbach statistic (windowed, lags 5/10/20);
$L$ = dynamic leverage ($\mathrm{Corr}(r_t,V_{t+L})$, lags 1--20);
$K$ = excess kurtosis of $r_t$;
$A$ = maximum absolute raw-return ACF;
Stress $=\mathbf{1}[\sigma_{\max}>2c_m\;\text{or}\;\bar\sigma<5\%]$.

\subsection{Optimisation}

Calibration minimises $-\bar S$ (mean over $n_{\rm cal}$ paths) via
Nelder--Mead (\texttt{scipy.optimize.minimize}).

\begin{center}
\renewcommand{\arraystretch}{1.4}
\begin{tabular}{lll}\toprule
Model & Parameters calibrated & Remark \\\midrule
ESN A2-981003 & $\Delta b_0$, scale $s$ &
  Architecture fixed; bias offset + vol rescaling only \\
QRH           & $H$, $\nu_{\rm vol}$, $\lambda$, $c_{\rm frac}$ &
  $c = c_{\rm frac}\cdot\sigma_{\rm daily}^2$ \\
PDV-GL        & $\beta$, $\alpha_1$, $\alpha_r$, $\alpha_v$ &
  $\bar V = \sigma_{\rm daily}^2\cdot\Delta t$ \\
\bottomrule\end{tabular}
\end{center}

%----------------------------------------------------------------------
\section{Fat-Tail Evaluation}
%----------------------------------------------------------------------

\subsection{Hill estimator}

For $n$ observations $|r_1|,\ldots,|r_n|$ with order statistics
$|r|_{(1)}\ge\cdots\ge|r|_{(n)}$:
\begin{equation}
  \hat\alpha_k \;=\; \left[\frac{1}{k}\sum_{i=1}^k
    \log\frac{|r|_{(i)}}{|r|_{(k+1)}}\right]^{-1}.
  \label{eq:hill}
\end{equation}
Bar charts use $k=\lfloor n/10\rfloor$; Hill plots show $\hat\alpha_k$ vs $k$.
Under the DGP, $\hat\alpha_k\to\nu$ as $n\to\infty$.

\subsection{Excess kurtosis}

$\hat\kappa = m_4/m_2^2-3$.
Under the DGP: $\mathrm{E}[\hat\kappa]\to 6/(\nu-4)$ for $\nu>4$ (infinite otherwise).

%----------------------------------------------------------------------
\section{Monte Carlo Setup}
%----------------------------------------------------------------------

\begin{center}
\renewcommand{\arraystretch}{1.4}
\begin{tabular}{ll}\toprule
Parameter & Value \\\midrule
DGP paths $n_{\rm dgp}$             & 30 \\
Calibration paths $n_{\rm cal}$     & 5  \\
Calibration trading days $T_{\rm cal}$ & 800 \\
Final simulation paths $n_{\rm sim}$& 30 \\
Final trading days $T_{\rm sim}$    & 2500 \\
Burn-in                             & 500 \\
Time step $\Delta t$                & 1 (daily) \\
$\nu$ grid                          & $\{3,5,10,30\}$ \\
QRH ring-buffer window              & 252 steps \\
\bottomrule\end{tabular}
\end{center}

%----------------------------------------------------------------------
\section{Output Figures}
%----------------------------------------------------------------------

\textbf{Figure 1.}  Four columns (one per $\nu$), four rows:
(i) log-scale histogram of $r_t$ with Gaussian overlay;
(ii) Hill plot $\hat\alpha(k)$ vs $k$ with true $\alpha=\nu$ reference;
(iii) excess kurtosis bar chart with theoretical $6/(\nu-4)$ reference;
(iv) QQ-plot of standardised returns against Gaussian.

\textbf{Figure 2.}  Three rows: Hill $\hat\alpha$ vs $\nu$,
excess kurtosis vs $\nu$, score $S$ vs $\nu$.

%----------------------------------------------------------------------
\begin{thebibliography}{9}
\bibitem{arch_report}
O.\ Zarhali.
\textit{Low-Dimensional A2 Architecture Study}.
Internal report, Paris Dauphine / CNRS, June 2026.

\bibitem{bourgey2026}
F.\ Bourgey, J.\ Gatheral.
\textit{Quadratic Rough Heston: SPX, VIX, and the SSR}.
SSRN:5239929, 2026.

\bibitem{guyon2023}
J.\ Guyon, J.\ Lekeufack.
\textit{Volatility is (mostly) path-dependent}.
Quantitative Finance 23(9):1221--1258, 2023.

\bibitem{bourgey2025}
F.\ Bourgey.
\textit{Almost path-dependent volatility model: calibration and simulation}.
arXiv:2507.09412, 2025.

\bibitem{hill1975}
B.M.\ Hill.
\textit{A simple general approach to inference about the tail of a distribution}.
Annals of Statistics 3(5):1163--1174, 1975.

\bibitem{cont2001}
R.\ Cont.
\textit{Empirical properties of asset returns: stylized facts and statistical issues}.
Quantitative Finance 1(2):223--236, 2001.
\end{thebibliography}

\end{document}
"""


def write_latex(path: str):
    with open(path, "w") as f: f.write(LATEX.lstrip())
    print(f"  Saved -> {path}")


# ============================================================================
# 8.  Main pipeline
# ============================================================================

def run(
    nu_grid   = None,
    n_dgp     = 30,
    n_cal     = 5,
    n_sim     = 30,
    T_dgp     = 2500, burn_dgp = 500,
    T_cal     = 800,  burn_cal = 200,
    T_sim     = 2500, burn_sim = 500,
    dt        = 1.0,
    prefix    = "student_t_v3",
):
    global _ESN_PARAMS
    if nu_grid is None: nu_grid = NU_GRID

    # Build ESN reservoir once
    _ESN_PARAMS = _build_esn_params(ARCH)
    assert _ESN_PARAMS['kappa0'] < 0, \
        f"kappa0={_ESN_PARAMS['kappa0']:.4f} must be negative"
    print(f"ESN kappa0 = {_ESN_PARAMS['kappa0']:.6f}  ✓\n")

    results = {}

    for nu in nu_grid:
        print(f"{'='*60}")
        print(f"  nu = {nu}  (true alpha = {nu})")
        print(f"{'='*60}")

        # ── STEP 1: DGP ─────────────────────────────────────────────
        print(f"  DGP (t_{nu}) ...", flush=True)
        dgp_sts, dgp_pool = dgp_student_t(nu, n_dgp, T_dgp, burn_dgp)
        ref = make_score_ref(dgp_sts)
        print(f"  Score reference: H_c={ref['H_hat_c']:.4f}  "
              f"vol_c={ref['mean_vol_ann_c']*100:.1f}%  "
              f"K_c={ref['kurtosis_c']:.3f}  "
              f"Z_c={ref['zumbach_c']:.4f}")

        # ── STEP 2: Calibrate ────────────────────────────────────────
        esn_cal = calibrate_esn(nu, ref, n_cal, T_cal, burn_cal, dt)
        qrh_cal = calibrate_qrh(nu, ref, n_cal, T_cal, burn_cal, dt)
        gl_cal  = calibrate_gl(nu,  ref, n_cal, T_cal, burn_cal, dt)

        # ── STEP 3: Simulate final paths ─────────────────────────────
        print(f"  Simulating final paths ...", flush=True)
        esn_sts = simulate_esn_paths(esn_cal, n_sim, T_sim, burn_sim, dt)
        qrh_sts = simulate_qrh_paths(qrh_cal, n_sim, T_sim, burn_sim, dt)
        gl_sts  = simulate_gl_paths(gl_cal,  n_sim, T_sim, burn_sim, dt)

        results[nu] = {
            "dgp_stats":     dgp_sts,
            "dgp_pool":      dgp_pool,
            "ref":           ref,
            "ESN A2-981003": esn_sts,
            "QRH":           qrh_sts,
            "PDV-GL":        gl_sts,
        }

        # Summary table
        def _m(sts, k): return float(np.nanmean([s[k] for s in sts]))
        print(f"\n  Results nu={nu}:")
        print(f"  {'Model':<22} {'H':>7} {'vol%':>6} {'Z':>8} "
              f"{'Lev':>7} {'Kurt':>7} {'Score':>7}")
        print("  " + "-"*68)
        for lbl, sts in [("DGP (target)", dgp_sts),
                          ("ESN A2-981003", esn_sts),
                          ("QRH",           qrh_sts),
                          ("PDV-GL",        gl_sts)]:
            sc = float(np.nanmean([score_fn(s, ref) for s in sts]))
            print(f"  {lbl:<22} {_m(sts,'H_hat'):>7.4f} "
                  f"{_m(sts,'mean_vol_ann')*100:>5.1f}% "
                  f"{_m(sts,'zumbach'):>8.4f} "
                  f"{_m(sts,'leverage'):>7.4f} "
                  f"{_m(sts,'kurtosis'):>7.3f} {sc:>7.3f}")
        print()

    # ── Figures ──────────────────────────────────────────────────────
    print(f"\n{'='*60}\nProducing figures ...")
    fig1 = plot_returns_grid(results, nu_grid, f"{prefix}_returns.png")
    plt.close(fig1)
    fig2 = plot_summary(results, nu_grid, f"{prefix}_summary.png")
    plt.close(fig2)

    # ── LaTeX ────────────────────────────────────────────────────────
    tex_path = f"{prefix}_protocol.tex"
    write_latex(tex_path)

    # Compile PDF if pdflatex available
    if os.system("which pdflatex > /dev/null 2>&1") == 0:
        os.system(f"pdflatex -interaction=nonstopmode {tex_path} > /dev/null 2>&1")
        os.system(f"pdflatex -interaction=nonstopmode {tex_path} > /dev/null 2>&1")
        pdf = tex_path.replace(".tex", ".pdf")
        if os.path.exists(pdf):
            print(f"  Compiled -> {pdf}")

    return results


# ============================================================================
# 9.  CLI
# ============================================================================

if __name__ == "__main__":
    import argparse
    pa = argparse.ArgumentParser(
        description="Student-t DGP fat-tail comparison (corrected model simulations)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pa.add_argument("--nu",    nargs="+", type=int, default=[3,5,10,30])
    pa.add_argument("--n_dgp", type=int, default=30)
    pa.add_argument("--n_cal", type=int, default=5)
    pa.add_argument("--n_sim", type=int, default=30)
    pa.add_argument("--T",     type=int, default=2500)
    pa.add_argument("--burn",  type=int, default=500)
    pa.add_argument("--dt",    type=float, default=1.0)
    pa.add_argument("--out",   type=str, default="student_t_v3")
    pa.add_argument("--fast",  action="store_true",
                    help="Quick test: nu=[5,10], n=8, T=800")
    args = pa.parse_args()
    if args.fast:
        args.nu=[5,10]; args.n_dgp=8; args.n_cal=3
        args.n_sim=8; args.T=800; args.burn=200
    run(nu_grid=args.nu,
        n_dgp=args.n_dgp, T_dgp=args.T, burn_dgp=args.burn,
        n_cal=args.n_cal, T_cal=min(args.T,800), burn_cal=min(args.burn,200),
        n_sim=args.n_sim, T_sim=args.T, burn_sim=args.burn,
        dt=args.dt, prefix=args.out)
