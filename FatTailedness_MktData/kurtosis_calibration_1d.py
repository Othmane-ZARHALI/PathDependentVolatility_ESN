"""
kurtosis_calibration_1d.py
======================================
Calibrates the ESN-A2 architecture's excess kurtosis of daily
log-returns against a random 9-asset bench (6 stocks + 3 indices),
fixing H, rough_scale, lam_lo, lam_hi at the SAME "good corner" that
resolved the leverage-effect, Taylor-effect, and Zumbach-effect
notes' own remaining gaps (H=0.3, rough_scale=0.7, lam_lo=1e-5,
lam_hi=1.0), and searching only m1 (the quadratic reservoir-norm
coupling in eta) per asset via a 1D bisection.

WHY 1D, NOT THE FULL 5-PARAMETER SEARCH: a first attempt jointly
calibrated all 5 parameters (H, rough_scale, lam_lo, lam_hi, m1) via
bounded Levenberg-Marquardt, starting both at m1=0 and at an informed
m1=-0.05. Both attempts got stuck: every asset converged to nearly the
SAME point (the other 4 parameters already at their own bounds from
the start, leaving m1 as the only free direction), with m1 barely
moving from its own starting value regardless of the target. This is
a real optimisation difficulty, not evidence the target is
unreachable: kurtosis is a high-variance statistic (dominated by rare,
extreme observations), so the finite-difference gradient LM estimates
for m1 is noisy relative to the true underlying sensitivity, and 5D
gradient search with only 30 evaluations does not reliably resolve it.
A manual sweep (see kurtosis_calibration.py's own docstring) already
established a clean, MONOTONIC relationship between m1 and kurtosis at
this fixed corner -- exactly the kind of well-behaved 1D problem
bisection solves reliably where gradient descent in 5D did not.

Requires: esn_base.py (core simulation module) and
bench9_real_kurtosis.json (real per-asset targets, from
kurtosis_calibration.py).
"""
import sys, json, time
sys.path.insert(0, ".")
import numpy as np
from scipy.stats import kurtosis
from scipy.optimize import brentq
import esn_base as E

OPT_ARCH = dict(
    rough_orientation=-1.0, Nr=22, Nz=29, matrix_seed=E.MATRIX_SEED,
    z_strength=0.27136029620674146, even_strength=3.1842560105665587,
    linear_strength=0.18199009132966565, gamma_norm=1.308098702657785,
    local_z_strength=0.06081431502674114, zz_scale=0.03120698941082984,
    sign_prob_neg=0.222781582556174,
)
# The fixed "good corner" (H, rough_scale, lam_lo, lam_hi), shared by every asset.
H_FIXED = 0.3
RS_FIXED = 0.7
LL_FIXED = 1e-5
LH_FIXED = 1.0
M1_BOUNDS = (-0.1, 0.1)

P_FIXED = dict(az_lo=1/400, az_hi=1/7, zr_lo=0.055, zr_hi=0.25,
               m2=0.0, b0_delta=0.0, scale=1.0,
               H=H_FIXED, rough_scale=RS_FIXED, lam_lo=LL_FIXED, lam_hi=LH_FIXED)

N_PATHS_EVAL = 20
T_EVAL = 4000
T_CONFIRM = 20000
N_REPLICATES_CI = 60
T_REP_CI = 4000
NAMES9 = ["LHX", "BK", "HUM", "TSN", "IEX", "KMB", "SP500", "RUSSELL2000", "NASDAQ100"]


def kurt_at_m1(m1, n_paths=N_PATHS_EVAL, T=T_EVAL, seed_offset=0):
    """Deterministic in m1 given seed_offset: always uses the same n_paths
    seeds, so the function bisection sees is reproducible, not stochastic."""
    ip = dict(P_FIXED, m1=m1)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    if mat["kappa0"] >= 0:
        return 1000.0
    ks = []
    for q in range(n_paths):
        seed = 15_300_000 + seed_offset + q
        x_, v_ = E._sim_esn_with_params(seed, T, OPT_ARCH, mat, b0_delta=0.0, scale=1.0)
        ks.append(kurtosis(x_, fisher=True, bias=False))
    return float(np.mean(ks))


def calibrate_m1_1d(real_target, seed_offset=0):
    def f(m1):
        return kurt_at_m1(m1, seed_offset=seed_offset) - real_target

    lo, hi = M1_BOUNDS
    f_lo, f_hi = f(lo), f(hi)
    print(f"    f({lo})={f_lo:.3f}  f({hi})={f_hi:.3f}", flush=True)

    if f_lo > 0 and f_hi > 0:
        # target below what even the most negative m1 can reach
        m1_star = lo
    elif f_lo < 0 and f_hi < 0:
        # target above what even the most positive m1 can reach
        m1_star = hi
    else:
        m1_star = brentq(f, lo, hi, xtol=0.002, maxiter=12)

    ip = dict(P_FIXED, m1=m1_star)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    # confirmation: average over several long paths too, not a single one --
    # kurtosis is too noisy for a single-path confirmation to be meaningful
    ks_confirm = []
    for q in range(8):
        x_, v_ = E._sim_esn_with_params(888 + q, T_CONFIRM, OPT_ARCH, mat, b0_delta=0.0, scale=1.0)
        ks_confirm.append(kurtosis(x_, fisher=True, bias=False))
    k_confirm = float(np.mean(ks_confirm))
    k_confirm_std = float(np.std(ks_confirm))

    return dict(H_fit=H_FIXED, rough_scale_fit=RS_FIXED, lam_lo_fit=LL_FIXED, lam_hi_fit=LH_FIXED,
                m1_fit=float(m1_star), k_confirm=k_confirm, k_confirm_std=k_confirm_std, nfev=0)


def build_confidence_band(m1, n_reps=N_REPLICATES_CI, T_rep=T_REP_CI, seed_offset=0):
    ip = dict(P_FIXED, m1=m1)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    ks = []
    for rep in range(n_reps):
        seed = 15_400_000 + seed_offset + rep
        x_, v_ = E._sim_esn_with_params(seed, T_rep, OPT_ARCH, mat, b0_delta=0.0, scale=1.0)
        ks.append(kurtosis(x_, fisher=True, bias=False))
    ks = np.array(ks)
    return dict(mean=float(np.mean(ks)), std=float(np.std(ks)),
                p5=float(np.percentile(ks, 5)), p95=float(np.percentile(ks, 95)),
                all=ks.tolist())


if __name__ == "__main__":
    real_kurt = json.load(open("bench9_real_kurtosis.json"))

    print("=== 1D (m1-only) calibration, all 9 assets ===")
    results = {}
    t0 = time.time()
    for i, name in enumerate(NAMES9):
        print(f"{name} (target={real_kurt[name]:.3f}):", flush=True)
        best = calibrate_m1_1d(real_kurt[name], seed_offset=i * 10000)
        results[name] = best
        print(f"  m1={best['m1_fit']:+.4f} nfev={best['nfev']} "
              f"confirm={best['k_confirm']:.3f} real={real_kurt[name]:.3f} "
              f"({time.time()-t0:.0f}s)", flush=True)
        json.dump(results, open("kurtosis_calibration_1d.json", "w"), indent=2)

    print("\n=== Building 90% confidence bands ===")
    ci_results = {}
    for i, name in enumerate(NAMES9):
        c = results[name]
        ci_results[name] = build_confidence_band(c["m1_fit"], seed_offset=i * 1000)
        print(f"  {name}: mean={ci_results[name]['mean']:.3f} "
              f"[{ci_results[name]['p5']:.3f}, {ci_results[name]['p95']:.3f}] "
              f"({time.time()-t0:.0f}s)", flush=True)
        json.dump(ci_results, open("kurtosis_ci_1d.json", "w"), indent=2)

    print("\nDone.")
