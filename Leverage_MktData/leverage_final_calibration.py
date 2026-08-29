"""
leverage_final_calibration.py
======================================
Final leverage-effect calibration for the ESN-A2 architecture: 5 free
parameters (H, rough_scale, lam_lo, lam_hi, m1) calibrated per asset
against the leverage lag profile, corr(x_t, v_{t+L}), lags L=1..40
(lag 0 excluded), weighted by 1/sqrt(L+1).

FULL JOURNEY (see leverage_effect_protocol.tex for the complete
write-up with equations):
1. H's bound widened from the standard [0.01,0.3] to [0.01,2.5] --
   q_i = lambda_i^(0.5-H) only reaches the reservoir's slow, genuinely
   persistent modes once H exceeds 0.5; below that q always favours
   the fastest mode, which the simulator's own shared-noise coupling
   mechanically dominates.
2. Lag 0 (contemporaneous correlation) excluded from the objective,
   and the weighting relaxed from 1/(L+1) to 1/sqrt(L+1): an earlier,
   lag-0-heavy objective risked overfitting a single, comparatively
   noisy point rather than the market's average covariance behaviour.
   Uniform weighting and a 2-summary-number (mean+decay) objective
   were also tried and rejected -- the former overshoots badly, the
   latter is under-determined with 4 free parameters and can converge
   to a degenerate spike solution.
3. Lags extended from 20 to 40 to reduce noise in both the real
   target and the simulated estimate.
4. A scalar (rather than vector) MSE objective was tried via
   Nelder-Mead and performed WORSE (converged to a shallow H~0.89) --
   minimising the sum of squared residuals and their mean are the same
   problem mathematically, but the optimiser matters in this noisy,
   non-convex landscape; bounded Levenberg-Marquardt (used throughout)
   is more reliable here than a derivative-free simplex method.
5. A persistent negative offset in the calibrated curve's own MEAN
   LEVEL remained after 1-4. Tested and REJECTED: freeing b0_delta
   (barely moves the mean -- correlation is largely insensitive to an
   additive variance shift). Tested and ADOPTED: freeing m1 (the
   quadratic reservoir-norm coupling in eta) -- a direct sweep showed
   m1<0 raises the leverage curve's mean substantially.
6. m1's FIRST bound, [-0.3,0.1] (reused from the roughness-exercise
   note's own sensitivity checks), left every asset's own optimum
   sitting AT or within 0.01 of the bound -- the bound itself, not the
   model, was limiting the fit, worse for stocks (smaller real
   targets) than indices. Widening m1's own lower bound to -1.0 and
   recalibrating closes the gap for every one of the 9 assets.

RESULT: every calibrated mean now within a few hundredths of its real
target (vs up to 3x off under the [-0.3,0.1] bound). A 90% confidence
band (60 replicates) shows the real profile falls inside it at 25-37
of 40 lags for every asset -- remaining visual gaps for LHX/BK/KMB are
consistent with sampling noise, not further systematic bias.

Requires: esn_base.py (core simulation module) and
bench9_real_leverage_wide.json (real asset targets, lags 0-40).
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
LH_BOUNDS = (1.0, 5.0)
M1_BOUNDS = (-1.0, 0.1)
BOUNDS_LO = [H_BOUNDS[0], RS_BOUNDS[0], LL_BOUNDS[0], LH_BOUNDS[0], M1_BOUNDS[0]]
BOUNDS_HI = [H_BOUNDS[1], RS_BOUNDS[1], LL_BOUNDS[1], LH_BOUNDS[1], M1_BOUNDS[1]]

CAL_LAGS = list(range(1, 41))  # lag 0 EXCLUDED from the objective
FULL_LAGS = list(range(0, 41))  # lag 0 still computed for plotting
WEIGHTS = np.array([1.0 / np.sqrt(l + 1) for l in CAL_LAGS])
X0 = (0.8, 0.45, 0.03, 1.5, -0.1)
N_PATHS_EVAL = 4
T_EVAL = 3000
MAX_NFEV = 40
T_CONFIRM = 12000
NAMES9 = ["LHX", "BK", "HUM", "TSN", "IEX", "KMB", "SP500", "RUSSELL2000", "NASDAQ100"]


def leverage_profile(x, v, lags):
    out = []
    for lag in lags:
        if lag == 0:
            out.append(np.corrcoef(x, v)[0, 1])
        else:
            out.append(np.corrcoef(x[:-lag], v[lag:])[0, 1])
    return np.array(out)


def calibrate_nolag0(real_profile_full, x0=X0, n_paths=N_PATHS_EVAL, T_eval=T_EVAL,
                      max_nfev=MAX_NFEV, seed_offset=0):
    """real_profile_full: real leverage correlation at lags 0..40 (41 values)."""
    targets = real_profile_full[1:]  # drop lag 0
    counter = [0]

    def resid(x):
        H, rs, ll, lh, m1 = x
        H = float(np.clip(H, *H_BOUNDS)); rs = float(np.clip(rs, *RS_BOUNDS))
        ll = float(np.clip(ll, *LL_BOUNDS)); lh = float(np.clip(lh, *LH_BOUNDS))
        lh = max(lh, ll * 2)
        m1 = float(np.clip(m1, *M1_BOUNDS))
        counter[0] += 1
        seed_base = 10_500_000 + seed_offset + counter[0] * 50
        ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh, m1=m1)
        mat = E.build_esn_matrices(OPT_ARCH, ip)
        if mat["kappa0"] >= 0:
            return np.full(len(CAL_LAGS), 1.0)
        profs = []
        for q in range(n_paths):
            x_, v_ = E._sim_esn_with_params(seed_base + q, T_eval, OPT_ARCH, mat,
                                             b0_delta=0.0, scale=1.0)
            profs.append(leverage_profile(x_, v_, CAL_LAGS))
        avg = np.mean(profs, axis=0)
        return (avg - targets) * WEIGHTS

    res = least_squares(resid, np.array(x0), bounds=(BOUNDS_LO, BOUNDS_HI),
                         method="trf", max_nfev=max_nfev,
                         xtol=1e-13, ftol=1e-13, gtol=1e-13, diff_step=0.15)
    H, rs, ll, lh, m1 = [float(np.clip(res.x[i], BOUNDS_LO[i], BOUNDS_HI[i])) for i in range(5)]

    ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh, m1=m1)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    x_, v_ = E._sim_esn_with_params(888, T_CONFIRM, OPT_ARCH, mat, b0_delta=0.0, scale=1.0)
    prof_confirm = leverage_profile(x_, v_, FULL_LAGS)  # lag 0 included for plotting only

    return dict(H_fit=H, rough_scale_fit=rs, lam_lo_fit=ll, lam_hi_fit=lh, m1_fit=m1,
                confirm_profile=prof_confirm.tolist(), nfev=res.nfev)


def build_confidence_band(H, rs, ll, lh, m1, n_reps=60, T_rep=3000, seed_offset=0):
    """60-replicate 90% confidence band (5th-95th percentile) at a
    calibrated point, over the full lag range 0-40."""
    ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh, m1=m1)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    profs = []
    for rep in range(n_reps):
        seed = 11_000_000 + seed_offset + rep
        x_, v_ = E._sim_esn_with_params(seed, T_rep, OPT_ARCH, mat, b0_delta=0.0, scale=1.0)
        profs.append(leverage_profile(x_, v_, FULL_LAGS))
    profs = np.array(profs)
    return dict(mean=np.nanmean(profs, axis=0).tolist(),
                p5=np.nanpercentile(profs, 5, axis=0).tolist(),
                p95=np.nanpercentile(profs, 95, axis=0).tolist())


if __name__ == "__main__":
    real_profiles = json.load(open("bench9_real_leverage_wide.json"))

    print("=== Final calibration (H, rough_scale, lam_lo, lam_hi, m1 free; "
          "lag 0 excluded), all 9 assets ===")
    results = {}
    t0 = time.time()
    for i, name in enumerate(NAMES9):
        real_prof = np.array(real_profiles[name])
        print(f"{name}:", flush=True)
        best = calibrate_nolag0(real_prof, seed_offset=i * 10000)
        results[name] = best
        confirm_prof = np.array(best["confirm_profile"])
        print(f"  H={best['H_fit']:.3f} rs={best['rough_scale_fit']:.3f} "
              f"lam_lo={best['lam_lo_fit']:.4f} lam_hi={best['lam_hi_fit']:.3f} "
              f"m1={best['m1_fit']:+.3f} nfev={best['nfev']}  "
              f"mean(lags1-40): confirm={confirm_prof[1:].mean():.4f} real={real_prof[1:].mean():.4f} "
              f"({time.time()-t0:.0f}s)", flush=True)
        json.dump(results, open("leverage_final_calibration.json", "w"), indent=2)

    print("\n=== Building 90% confidence bands (60 replicates each) ===")
    ci_results = {}
    for i, name in enumerate(NAMES9):
        c = results[name]
        ci_results[name] = build_confidence_band(
            c["H_fit"], c["rough_scale_fit"], c["lam_lo_fit"], c["lam_hi_fit"], c["m1_fit"],
            seed_offset=i * 1000)
        print(f"  {name}: done ({time.time()-t0:.0f}s)", flush=True)
        json.dump({"lags": FULL_LAGS, "by_asset": ci_results},
                  open("leverage_final_ci.json", "w"), indent=2)

    print("\nDone.")
