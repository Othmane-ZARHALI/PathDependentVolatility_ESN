"""
esn_base.py (reconstructed)
======================================
Core ESN-A2 simulation engine, extracted/reconstructed from
esn_hyperparam_search.py (same MATRIX_SEED, build_esn_matrices,
_sim_esn_with_params, make_score_ref, score_fn_smooth, SCORE_TERM_WEIGHTS
-- verified byte-identical function signatures) after the working
container environment reset and these modules were lost. Validated by
reproducing the exact target H_hat_c values documented in
esn_market_protocol.tex Table 2 for the original 5 assets
(GPC=0.0467, QCOM=0.0529, INTC=0.0558, LUV=0.0459, NI=0.0543).
"""
import math
import numpy as np

TRADING_DAYS = 252.0
TV_ANN       = 0.20
TV_DAY       = TV_ANN / math.sqrt(TRADING_DAYS)
SIG_MIN      = 0.01  / math.sqrt(TRADING_DAYS)
MU_S         = 0.0
TIME_DT      = 1.0 / TRADING_DAYS

MATRIX_SEED = 202695547565

# The 7 shared z-bank hyperparameters: (name, A2-981003 default, lo, hi,
# transform). Bounds match the companion note's own HP_SPEC table
# (§1.1); transform="linear" throughout for this note's own bounded-LM
# (TRF) search, which enforces box bounds natively and does not need
# the log/linear transform some other (e.g. CMA-ES) scripts in this
# family apply.
HP_SPEC = [
    ("z_strength",      0.34,  0.05, 2.0, "linear"),
    ("even_strength",   1.50,  0.10, 5.0, "linear"),
    ("linear_strength", 0.25,  0.01, 2.0, "linear"),
    ("gamma_norm",      1.00,  0.10, 3.0, "linear"),
    ("local_z_strength",0.030, 0.001,0.5, "linear"),
    ("zz_scale",        0.080, 0.01, 0.5, "linear"),
    ("sign_prob_neg",   0.22,  0.05, 0.50,"linear"),
]

# The bounded-LM-found optimal values of the 7 shared hyperparameters
# (Nr=22, Nz=29 alongside these -- see esn_market_protocol.tex Table 1).
OPTIMAL_ARCH = dict(
    z_strength=0.27136029620674146,
    even_strength=3.1842560105665587,
    linear_strength=0.18199009132966565,
    gamma_norm=1.308098702657785,
    local_z_strength=0.06081431502674114,
    zz_scale=0.03120698941082984,
    sign_prob_neg=0.222781582556174,
)

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

SCORE_MAX = 2.2+1.4+0.6+0.2+0.8+1.0+0.8+1.0+1.1+0.7+0.5 + 5.0   # = 15.3
SCORE_TERM_WEIGHTS = [
    ("H_hat", 2.2), ("mean_vol_ann", 1.4), ("q995_vol_ann", 0.6),
    ("max_vol_ann", 0.2), ("mean_vol_acf", 0.8), ("taylor_gap", 1.0),
    ("taylor_frac", 0.8), ("zumbach", 1.0), ("leverage", 1.1),
    ("kurtosis", 0.7), ("max_ret_acf", 0.5),
]

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
