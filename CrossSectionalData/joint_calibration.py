"""
joint_calibration.py
======================================
Calibrates the ESN-A2 architecture against ALL 5 stylised facts this
note family has examined separately -- roughness (moment-based Hurst
exponent), the leverage effect, the Taylor effect, the Zumbach effect,
and excess kurtosis -- SIMULTANEOUSLY, per asset, by building a single
combined residual vector that concatenates each individual note's own
established objective, and letting bounded Levenberg-Marquardt
minimise the sum of squares of all of them jointly.

FREE PARAMETERS: theta = (H, rough_scale, lam_lo, lam_hi, m1), 5
dimensions:
  H           in [0.01, 2.5]
  rough_scale in [0.2, 0.7]
  lam_lo      in [1e-5, 0.1]
  lam_hi      in [0.001, 5.0]
  m1          in [-1.0, 0.5]
m2 and the other 6 entries of p are held fixed at their own grid-
midpoint values (see any of this note family's own "how the non-
calibrated parameters are fixed" sections).

ONE simulated variance signal v_t (= sigma_t^2) drives returns AND
every stylised fact estimated from a given path -- Hurst, leverage,
Taylor, Zumbach, and kurtosis all read the same simulated (x_t, v_t),
1 simulation per replicate.

THE COMBINED RESIDUAL VECTOR, per asset, concatenates 5 blocks, each
normalised by a cross-asset tolerance s_k AND by sqrt(number of terms
in that block), so no single fact numerically dominates just by having
more terms or a larger natural scale:
  - roughness:  1 term,  (Hhat_ESN - Hhat_real)
  - leverage:   40 terms, lags 1-40, 1/(L+1)^2-weighted (strong
                short-lag emphasis)
  - Taylor:     12 terms, lags {1,3,5,10,15,20}, acf(|x|) & acf(x^2)
  - Zumbach:    5 terms,  L in {2,10,20,30,40}, Z(L)
  - kurtosis:   1 term,  (kurtosis_ESN - kurtosis_real)
= 59 residual entries in total per asset. Hurst and Zumbach additionally
carry extra weight multipliers (3x, 4x) on top of the tolerance/count
normalisation.

STARTING POINTS matter more than any weighting scheme tried: a single
shared starting point, or even a coarse stock-vs-index split, leaves
roughness coverage stuck near 0/9 because H has almost no effective
control over Hhat in the wrong region of lam_hi. The reliable recipe is
to invert an H-to-Hhat map (built by a direct sweep at a wide lam_hi,
e.g. 3.5) against each asset's own real Hurst target, giving a precise
per-asset H0 (and, for targets below what lam_hi=3.5 reaches, a wider
lam_hi0=5.0 too).

Requires: esn_base.py (core simulation module) and the 5 real-target
JSONs: bench9_real_hurst.json, bench9_real_leverage_wide.json,
bench9_real_abs_acf_40.json, bench9_real_sq_acf_40.json,
bench9_real_zumbach_profile.json, bench9_real_kurtosis.json.
"""
import sys, json, time
sys.path.insert(0, ".")
import numpy as np
from scipy.optimize import least_squares
import esn_base as E

OPT_ARCH = dict(
    rough_orientation=-1.0, Nr=22, Nz=29, matrix_seed=E.MATRIX_SEED,
    z_strength=0.27136029620674146, even_strength=3.1842560105665587,
    linear_strength=0.18199009132966565, gamma_norm=1.308098702657785,
    local_z_strength=0.06081431502674114, zz_scale=0.03120698941082984,
    sign_prob_neg=0.222781582556174,
)
P_FIXED = dict(az_lo=1/400, az_hi=1/7, zr_lo=0.055, zr_hi=0.25,
               m2=0.0, b0_delta=0.0, scale=1.0)

H_BOUNDS = (0.01, 2.5)
RS_BOUNDS = (0.2, 0.7)
LL_BOUNDS = (1e-5, 0.1)
LH_BOUNDS = (0.001, 5.0)
M1_BOUNDS = (-1.0, 0.5)  # upper widened from 0.1: confirmed by direct sweep that
ALPHA_BOUNDS = (0.02, 1.0)  # feedback-layer smoothing factor; 1.0 recovers the single-signal model
BOUNDS_LO = [H_BOUNDS[0], RS_BOUNDS[0], LL_BOUNDS[0], LH_BOUNDS[0], M1_BOUNDS[0], ALPHA_BOUNDS[0]]
BOUNDS_HI = [H_BOUNDS[1], RS_BOUNDS[1], LL_BOUNDS[1], LH_BOUNDS[1], M1_BOUNDS[1], ALPHA_BOUNDS[1]]

LEV_LAGS = list(range(1, 41))
TAYLOR_LAGS = [1, 3, 5, 10, 15, 20]
ZUMBACH_LAGS = [2, 10, 20, 30, 40]
HURST_LAGS_ESN = np.arange(1, 41)  # ESN paths have no proxy noise: full range
MOMENT_QS = [0.5, 1.0, 1.5, 2.0]

N_PATHS_EVAL = 6
T_EVAL = 4000
MAX_NFEV = 40
T_CONFIRM = 12000
N_REPLICATES_CI = 40
T_REP_CI = 4000
NAMES9 = ["LHX", "BK", "HUM", "TSN", "IEX", "KMB", "SP500", "RUSSELL2000", "NASDAQ100"]


# ---------- per-path statistic extraction (shared across all 5 facts) ----------

def _corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() < 1e-14 or b.std() < 1e-14:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def smooth_curve(arr, window=3):
    """Light 3-point moving average across adjacent entries of a
    lag-indexed statistic (leverage, Taylor acf, Zumbach Z(L)),
    reducing raw sampling noise in each individual point estimate.
    Edges use a shrinking window (no wraparound, no padding with
    zeros) so the smoothed curve's own endpoints are not biased toward
    0. Applied identically to the ESN's own simulated statistics and
    the real target curves, so the comparison itself is unaffected by
    the smoothing's own choice of window -- only the per-point noise
    each side must fit is reduced."""
    arr = np.asarray(arr, float)
    n = len(arr)
    half = window // 2
    out = np.empty(n)
    for i in range(n):
        lo = max(0, i - half); hi = min(n, i + half + 1)
        seg = arr[lo:hi]
        seg = seg[~np.isnan(seg)]
        out[i] = seg.mean() if len(seg) else np.nan
    return out


def leverage_at_lags(x, v, lags):
    raw = np.array([_corr(x[:-L], v[L:]) for L in lags])
    return smooth_curve(raw)


def acf_at_lags(series, lags):
    raw = np.array([_corr(series[:-L], series[L:]) for L in lags])
    return smooth_curve(raw)


def zumbach_Z(x, v, lags):
    n = len(x)
    out = []
    for L in lags:
        m = n - 2 * L
        if m > 30:
            csx = np.concatenate([[0.0], np.cumsum(x)])
            csv = np.concatenate([[0.0], np.cumsum(v)])
            idx = L + np.arange(m)
            pR = csx[idx] - csx[idx - L]
            fR = csx[idx + L] - csx[idx]
            pV = (csv[idx] - csv[idx - L]) / L
            fV = (csv[idx + L] - csv[idx]) / L
            out.append(_corr(pR ** 2, fV) - _corr(pV, fR ** 2))
        else:
            out.append(np.nan)
    return smooth_curve(np.array(out))


def excess_kurtosis(x):
    c = x - x.mean()
    var = np.var(c)
    return float(np.mean(c ** 4) / (var ** 2 + 1e-30) - 3.0)


def estimate_H(log_sigma, lags):
    zetas = []
    for q in MOMENT_QS:
        log_ms = []
        for L in lags:
            d = np.abs(log_sigma[L:] - log_sigma[:-L]) ** q
            log_ms.append(np.log(max(np.mean(d), 1e-300)))
        slope, _ = np.polyfit(np.log(lags), log_ms, 1)
        zetas.append(slope)
    Hhat, _ = np.polyfit(MOMENT_QS, zetas, 1)
    return float(Hhat)


def all_stats_one_path(x, v, v_raw=None):
    """Every one of the 5 stylised facts, extracted from the SAME
    simulated path. Hhat, leverage, and Zumbach are estimated from
    v_raw (the unsmoothed signal) -- an EWMA is itself autocorrelated,
    and was found to inject a spurious dip-and-recovery shape into
    leverage's own corr(x_t,v_{t+L}) when read from a smoothed v.
    Taylor and kurtosis, functions of x alone, benefit from whatever
    persistence alpha's own feedback layer injected into x. Falls back
    to v (single-signal behaviour) if v_raw is not given. Each
    lag-indexed statistic is smoothed via smooth_curve (leverage_
    at_lags/acf_at_lags/zumbach_Z all apply it internally)."""
    v_for_structural = v_raw if v_raw is not None else v
    log_sigma = 0.5 * np.log(np.maximum(v_for_structural, 1e-30))
    return dict(
        Hhat=estimate_H(log_sigma, HURST_LAGS_ESN),
        leverage=leverage_at_lags(x, v_for_structural, LEV_LAGS),
        abs_acf=acf_at_lags(np.abs(x), TAYLOR_LAGS),
        sq_acf=acf_at_lags(x ** 2, TAYLOR_LAGS),
        zumbach=zumbach_Z(x, v_for_structural, ZUMBACH_LAGS),
        kurtosis=excess_kurtosis(x),
    )


# ---------- joint residual and calibration ----------

def build_target_vector(real_hurst, real_lev, real_abs, real_sq, real_zum, real_kurt, asset_name):
    lev_idx = [l - 1 for l in LEV_LAGS]         # real_lev has lag 0..40
    taylor_idx = [l - 1 for l in TAYLOR_LAGS]   # real_abs/sq have lag 1..40 (index 0 = lag1)
    zum_lags_all = real_zum["lags"]
    zum_idx = [zum_lags_all.index(L) for L in ZUMBACH_LAGS]
    real_zum_Z = np.array(real_zum["pR2_fV"][asset_name]) - np.array(real_zum["pV_fR2"][asset_name])
    # Smoothed identically to the ESN's own simulated statistics
    # (leverage_at_lags/acf_at_lags/zumbach_Z, see smooth_curve's own
    # docstring), so the comparison is fair -- neither side is
    # penalised or favoured by the smoothing's own choice of window.
    real_lev_a = smooth_curve(np.array(real_lev)[lev_idx])
    real_abs_a = smooth_curve(np.array(real_abs)[taylor_idx])
    real_sq_a = smooth_curve(np.array(real_sq)[taylor_idx])
    real_zum_a = smooth_curve(real_zum_Z[zum_idx])
    return dict(
        Hhat=real_hurst,
        leverage=real_lev_a,
        abs_acf=real_abs_a,
        sq_acf=real_sq_a,
        zumbach=real_zum_a,
        kurtosis=real_kurt,
    )


# ---------- tolerances s_k, cross-asset spread of each fact's own real value ----------
# (analogous to the main protocol's own per-asset s_k, computed here across the 9
# assets rather than across 1 asset's own 6 time windows, for tractability)
TOLERANCES = dict(
    Hhat=0.0166, leverage=0.0247, abs_acf=0.0768, sq_acf=0.0894,
    zumbach=0.0807, kurtosis=3.3962,
)


def joint_residual(avg_stats, targets, normalize=True):
    """Concatenates all 5 objectives' own residuals into a single
    vector. Each block is divided by its own cross-asset tolerance s_k
    (TOLERANCES) AND by sqrt(number of terms in that block), so that
    no single stylised fact dominates the total sum of squares merely
    by having a larger natural scale (kurtosis) or more residual
    entries -- both would otherwise happen. Leverage's own lag
    weighting uses 1/(L+1)^2 (strong short-lag emphasis); Hurst
    carries its own extra weight multiplier on top of this
    normalisation. Per direct request, Taylor's own 2 curves
    (acf(|x|), acf(x^2)) are each fit separately against their own
    real target (not the difference between them), and Zumbach is fit
    on its own raw level (not a shape-normalised/z-scored version);
    neither carries an extra weight multiplier beyond the standard
    tolerance/block-count normalisation below. All 4 lag-indexed
    statistics (leverage, both Taylor curves, Zumbach) are smoothed
    via smooth_curve before this function ever sees them (applied
    inside leverage_at_lags/acf_at_lags/zumbach_Z for the ESN's own
    simulated statistics, and directly in build_target_vector for the
    real targets)."""
    s = TOLERANCES if normalize else {k: 1.0 for k in TOLERANCES}
    block_n = dict(Hhat=1, leverage=len(LEV_LAGS), abs_acf=len(TAYLOR_LAGS),
                   sq_acf=len(TAYLOR_LAGS), zumbach=len(ZUMBACH_LAGS), kurtosis=1)
    bn = {k: (np.sqrt(v) if normalize else 1.0) for k, v in block_n.items()}
    # Additional per-fact weight multipliers, applied on top of the
    # tolerance/block-count normalisation above. Taylor and Zumbach
    # carry NO extra multiplier (1.0), per direct request.
    fact_boost = dict(Hhat=3.0, leverage=1.0, abs_acf=1.0, sq_acf=1.0, zumbach=1.0, kurtosis=1.0)

    r_hurst = np.array([avg_stats["Hhat"] - targets["Hhat"]]) / (s["Hhat"] * bn["Hhat"] / fact_boost["Hhat"])
    lev_w = 1.0 / (np.array(LEV_LAGS) + 1) ** 2.0  # steepened from **1.5: more small-lag weight.
    r_lev = lev_w * (avg_stats["leverage"] - targets["leverage"]) / (s["leverage"] * bn["leverage"] / fact_boost["leverage"])

    r_abs = (avg_stats["abs_acf"] - targets["abs_acf"]) / (s["abs_acf"] * bn["abs_acf"] / fact_boost["abs_acf"])
    r_sq = (avg_stats["sq_acf"] - targets["sq_acf"]) / (s["sq_acf"] * bn["sq_acf"] / fact_boost["sq_acf"])

    r_zum = (avg_stats["zumbach"] - targets["zumbach"]) / (s["zumbach"] * bn["zumbach"] / fact_boost["zumbach"])

    r_kurt = np.array([avg_stats["kurtosis"] - targets["kurtosis"]]) / (s["kurtosis"] * bn["kurtosis"] / fact_boost["kurtosis"])
    return np.concatenate([r_hurst, r_lev, r_abs, r_sq, r_zum, r_kurt])


def bisect_m1_for_kurtosis(H, rs, ll, lh, target_kurtosis, n_paths=20, T_eval=4000,
                            seed_offset=0, alpha=1.0):
    """Stage-2 refinement: the stage-1 joint search's own finite-
    difference gradient leaves m1 stuck near its own starting value --
    verified NOT because the gradient is flat (a direct sweep shows
    kurtosis moves cleanly across m1's own range at a fixed calibrated
    point) but because the full combined objective's own evaluation-
    to-evaluation noise swamps m1's own marginal contribution to a
    numerical derivative taken across all terms at once. Bound
    bisection (not a joint least-squares gradient step) targeting
    kurtosis specifically, holding H, rough_scale, lam_lo, lam_hi,
    alpha fixed at their own stage-1 values, with heavy, fixed-seed
    averaging (20 paths, T=4000, non-incrementing across brentq's own
    internal calls), resolves this reliably."""
    from scipy.optimize import brentq

    def kurt_at_m1(m1):
        ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh, m1=m1)
        mat = E.build_esn_matrices(OPT_ARCH, ip)
        if mat["kappa0"] >= 0:
            return 1000.0
        ks = []
        for q in range(n_paths):
            seed = 16_600_000 + seed_offset + q
            x_, v_, v_raw_ = E._sim_esn_with_params(seed, T_eval, OPT_ARCH, mat, b0_delta=0.0, scale=1.0, alpha=alpha)
            ks.append(excess_kurtosis(x_))
        return float(np.mean(ks))

    def f(m1):
        return kurt_at_m1(m1) - target_kurtosis

    lo, hi = M1_BOUNDS
    f_lo, f_hi = f(lo), f(hi)
    if f_lo > 0 and f_hi > 0:
        m1_star = lo
    elif f_lo < 0 and f_hi < 0:
        m1_star = hi
    else:
        m1_star = brentq(f, lo, hi, xtol=0.002, maxiter=12)
    return float(m1_star), kurt_at_m1(m1_star)


def calibrate_joint(targets, x0, n_paths=N_PATHS_EVAL, T_eval=T_EVAL,
                     max_nfev=MAX_NFEV, seed_offset=0, cost_log=None):
    """Joint calibration over theta=(H, rough_scale, lam_lo, lam_hi,
    m1, alpha), minimising the combined 53-term residual (Hurst +
    leverage + Taylor + Zumbach + kurtosis) via bounded
    Levenberg-Marquardt. alpha is a feedback-layer EWMA smoothing
    factor splitting the raw variance signal (drives Hurst, leverage,
    Zumbach) from the smoothed signal actually driving returns (hence
    Taylor and kurtosis). x0 = (H0, rough_scale0, lam_lo0, lam_hi0,
    m1_0, alpha0). If cost_log is a list, the sum-of-squares cost of
    every residual evaluation is appended to it, in call order --
    diagnostic use only, no effect on the optimisation itself.

    COMMON RANDOM NUMBERS: seed_base is FIXED across every call to
    resid() within one calibrate_joint() run (not incremented per
    call). scipy's own finite-difference Jacobian estimate probes
    theta+d*e_i for each parameter i; with a fresh random seed at
    every probe, the estimated derivative (r(theta+d*e_i)-r(theta))/d
    was dominated by 2 INDEPENDENT stochastic path realisations' own
    sampling noise, not the true local sensitivity of the residual to
    theta -- confirmed to be the dominant cause of erratic, barely-
    decaying evaluation traces when tracking cost vs. evaluation
    number. Fixing the seed means every probe at a given outer
    iteration reuses the SAME underlying random paths as the base
    point, so most of the sampling noise cancels in the difference,
    leaving the genuine local gradient. The seed is still asset-
    specific (via seed_offset) and still changes between separate
    calibrate_joint() calls, so different assets and different
    calibration runs are not correlated -- only evaluations WITHIN one
    run share randomness."""
    counter = [0]

    def resid(x):
        H, rs, ll, lh, m1, alpha = x
        H = float(np.clip(H, *H_BOUNDS)); rs = float(np.clip(rs, *RS_BOUNDS))
        ll = float(np.clip(ll, *LL_BOUNDS)); lh = float(np.clip(lh, *LH_BOUNDS))
        lh = max(lh, ll * 1.5)
        m1 = float(np.clip(m1, *M1_BOUNDS))
        alpha = float(np.clip(alpha, *ALPHA_BOUNDS))
        counter[0] += 1
        seed_base = 17_500_000 + seed_offset  # FIXED (common random numbers) -- was
                                                # + counter[0]*50, which gave every
                                                # finite-difference probe independent
                                                # noise and swamped the true gradient
        ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh, m1=m1)
        mat = E.build_esn_matrices(OPT_ARCH, ip)
        if mat["kappa0"] >= 0:
            r = np.full(53, 5.0)
            if cost_log is not None:
                cost_log.append(float(np.sum(r ** 2)))
            return r
        path_stats = []
        for q in range(n_paths):
            x_, v_, v_raw_ = E._sim_esn_with_params(seed_base + q, T_eval, OPT_ARCH, mat,
                                             b0_delta=0.0, scale=1.0, alpha=alpha)
            path_stats.append(all_stats_one_path(x_, v_, v_raw_))
        avg = dict(
            Hhat=np.nanmean([s["Hhat"] for s in path_stats]),
            leverage=np.nanmean([s["leverage"] for s in path_stats], axis=0),
            abs_acf=np.nanmean([s["abs_acf"] for s in path_stats], axis=0),
            sq_acf=np.nanmean([s["sq_acf"] for s in path_stats], axis=0),
            zumbach=np.nanmean([s["zumbach"] for s in path_stats], axis=0),
            kurtosis=np.nanmean([s["kurtosis"] for s in path_stats]),
        )
        r = joint_residual(avg, targets)
        r = np.nan_to_num(r, nan=1.0)
        if cost_log is not None:
            cost_log.append(float(np.sum(r ** 2)))
        return r

    res = least_squares(resid, np.array(x0), bounds=(BOUNDS_LO, BOUNDS_HI),
                         method="trf", max_nfev=max_nfev,
                         xtol=1e-8, ftol=1e-8, gtol=1e-8, diff_step=0.15)
    H, rs, ll, lh, m1, alpha = [float(np.clip(res.x[i], BOUNDS_LO[i], BOUNDS_HI[i])) for i in range(6)]

    ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh, m1=m1)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    x_, v_, v_raw_ = E._sim_esn_with_params(888, T_CONFIRM, OPT_ARCH, mat, b0_delta=0.0, scale=1.0, alpha=alpha)
    confirm = all_stats_one_path(x_, v_, v_raw_)

    return dict(H_fit=H, rough_scale_fit=rs, lam_lo_fit=ll, lam_hi_fit=lh, m1_fit=m1,
                alpha_fit=alpha,
                confirm=dict(
                    Hhat=confirm["Hhat"],
                    leverage=confirm["leverage"].tolist(),
                    abs_acf=confirm["abs_acf"].tolist(),
                    sq_acf=confirm["sq_acf"].tolist(),
                    zumbach=confirm["zumbach"].tolist(),
                    kurtosis=confirm["kurtosis"],
                ),
                nfev=res.nfev, cost=float(res.cost))


def build_confidence_band(H, rs, ll, lh, m1, n_reps=N_REPLICATES_CI, T_rep=T_REP_CI,
                           seed_offset=0, alpha=1.0):
    ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh, m1=m1)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    reps = {k: [] for k in ["Hhat", "leverage", "abs_acf", "sq_acf", "zumbach", "kurtosis"]}
    for rep in range(n_reps):
        seed = 16_100_000 + seed_offset + rep
        x_, v_, v_raw_ = E._sim_esn_with_params(seed, T_rep, OPT_ARCH, mat, b0_delta=0.0, scale=1.0, alpha=alpha)
        st = all_stats_one_path(x_, v_, v_raw_)
        for k in reps:
            reps[k].append(st[k])
    out = {}
    for k, vals in reps.items():
        arr = np.array(vals)
        out[k] = dict(mean=np.nanmean(arr, axis=0).tolist() if arr.ndim > 1 else float(np.nanmean(arr)),
                      p5=np.nanpercentile(arr, 5, axis=0).tolist() if arr.ndim > 1 else float(np.nanpercentile(arr, 5)),
                      p95=np.nanpercentile(arr, 95, axis=0).tolist() if arr.ndim > 1 else float(np.nanpercentile(arr, 95)))
    return out


if __name__ == "__main__":
    real_hurst = json.load(open("bench9_real_hurst.json"))
    real_lev = json.load(open("bench9_real_leverage_wide.json"))
    real_abs = json.load(open("bench9_real_abs_acf_40.json"))
    real_sq = json.load(open("bench9_real_sq_acf_40.json"))
    real_zum = json.load(open("bench9_real_zumbach_profile.json"))
    real_kurt = json.load(open("bench9_real_kurtosis.json"))

    # Per-asset starting H, derived by inverting a direct H-to-Hhat
    # sweep (rough_scale=0.4, lam_lo=0.03, lam_hi=3.5) against each
    # asset's own real Hurst target. 2 assets need lam_hi0=5.0 instead
    # (targets below what lam_hi=3.5 reaches at H=0.01).
    H_MAP_H = [0.01, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.50]
    H_MAP_HHAT = [0.0458, 0.0473, 0.0490, 0.0516, 0.0553, 0.0615, 0.0695, 0.0827, 0.1063]

    def starting_point(asset_name):
        target = real_hurst[asset_name]
        if target < H_MAP_HHAT[0]:
            H0, lh0 = 0.01, 5.0
        else:
            H0 = float(np.interp(target, H_MAP_HHAT, H_MAP_H))
            lh0 = 3.5
        return (H0, 0.4, 0.03, lh0, 0.0)

    # Stage 1: joint 5-parameter search per asset.
    stage1 = {}
    t0 = time.time()
    for i, name in enumerate(NAMES9):
        targets = build_target_vector(
            real_hurst[name], real_lev[name], real_abs[name], real_sq[name],
            real_zum, real_kurt[name], name)
        x0 = starting_point(name)
        print(f"{name}: x0={x0}", flush=True)
        best = calibrate_joint(targets, x0, n_paths=8, max_nfev=35, seed_offset=i * 110000)
        stage1[name] = best
        print(f"  H={best['H_fit']:.3f} rs={best['rough_scale_fit']:.3f} "
              f"lam_lo={best['lam_lo_fit']:.5f} lam_hi={best['lam_hi_fit']:.4f} "
              f"m1={best['m1_fit']:+.4f} cost={best['cost']:.3f} "
              f"nfev={best['nfev']} ({time.time()-t0:.0f}s)", flush=True)
        json.dump(stage1, open("joint_calibration_stage1.json", "w"), indent=2)

    # Stage 2: kurtosis-targeted bisection for m1, holding stage-1's
    # own H/rough_scale/lam_lo/lam_hi fixed.
    final = {}
    for i, name in enumerate(NAMES9):
        c = stage1[name]
        m1_star, k_confirm = bisect_m1_for_kurtosis(
            c["H_fit"], c["rough_scale_fit"], c["lam_lo_fit"], c["lam_hi_fit"],
            real_kurt[name], n_paths=24, seed_offset=i * 1000)
        final[name] = dict(H_fit=c["H_fit"], rough_scale_fit=c["rough_scale_fit"],
                            lam_lo_fit=c["lam_lo_fit"], lam_hi_fit=c["lam_hi_fit"],
                            m1_fit=m1_star)
        print(f"{name}: m1={m1_star:+.4f} confirm_kurt={k_confirm:.3f} "
              f"target={real_kurt[name]:.3f}", flush=True)
        json.dump(final, open("joint_calibration_final.json", "w"), indent=2)

    print("\nDone. Final calibrated parameters in joint_calibration_final.json.")

