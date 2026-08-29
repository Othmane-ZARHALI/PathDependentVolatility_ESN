"""
zumbach_widened_lh_calibration.py
======================================
Calibrates the ESN-A2 architecture's Zumbach effect (Z(L) =
corr(pastReturn^2, futureVol) - corr(pastVol, futureReturn^2)) against
the 9-asset bench, with lambda_hi's own lower bound widened from the
standard 1.0 to 0.001.

WHY: the original shape calibration (zumbach_shape_calibration.py)
found every asset's own optimum sitting AT or against every one of the
standard bounds (H=0.30 upper, rough_scale=0.20 lower, lam_lo~0 lower,
lam_hi=1.0 lower) -- exactly the "bound was too narrow, not the model"
pattern found for leverage's H and Taylor's lambda_hi. A direct sweep
at RUSSELL2000's own calibrated point confirmed: widening H further
(beyond 0.3) makes Z(L) WORSE (values shrink toward 0, same as Taylor's
own H-widening test), but widening lambda_hi BELOW 1.0 genuinely helps
at longer lags (lam_hi=0.10 raises Z(30) from 0.173 to 0.235, against
a real target of 0.308) at some cost to the shortest lag.

This script recalibrates all 9 assets with lambda_hi's bound widened
to [0.001, 5.0], using the same multi-lag objective as
zumbach_shape_calibration.py.

Requires: esn_base.py (core simulation module) and
bench9_real_zumbach_profile.json (real per-asset Z(L) targets, lags
2..40 step 2, from zumbach_effect_evaluation.py).
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
               m1=0.0, m2=0.0, b0_delta=0.0, scale=1.0)

H_BOUNDS = (0.01, 0.3)
RS_BOUNDS = (0.2, 0.7)
LL_BOUNDS = (1e-5, 0.1)
LH_BOUNDS = (0.001, 5.0)  # WIDENED lower bound (was 1.0)
BOUNDS_LO = [H_BOUNDS[0], RS_BOUNDS[0], LL_BOUNDS[0], LH_BOUNDS[0]]
BOUNDS_HI = [H_BOUNDS[1], RS_BOUNDS[1], LL_BOUNDS[1], LH_BOUNDS[1]]

TARGET_LAG_IDX = [0, 4, 9, 14, 19]  # lags 2, 10, 20, 30, 40
X0 = (0.3, 0.2, 1e-5, 0.1)  # informed start near the hand-found good region
N_PATHS_EVAL = 4
T_EVAL = 2500
MAX_NFEV = 35
T_CONFIRM = 8000
NAMES9 = ["LHX", "BK", "HUM", "TSN", "IEX", "KMB", "SP500", "RUSSELL2000", "NASDAQ100"]


def _corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() < 1e-14 or b.std() < 1e-14:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


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
    return np.array(out)


def calibrate_widened_lh(real_Z_full, lags_all, x0=X0, n_paths=N_PATHS_EVAL, T_eval=T_EVAL,
                          max_nfev=MAX_NFEV, seed_offset=0):
    target_lags = [lags_all[i] for i in TARGET_LAG_IDX]
    targets = real_Z_full[TARGET_LAG_IDX]
    counter = [0]

    def resid(x):
        H, rs, ll, lh = x
        H = float(np.clip(H, *H_BOUNDS)); rs = float(np.clip(rs, *RS_BOUNDS))
        ll = float(np.clip(ll, *LL_BOUNDS)); lh = float(np.clip(lh, *LH_BOUNDS))
        lh = max(lh, ll * 1.5)
        counter[0] += 1
        seed_base = 14_000_000 + seed_offset + counter[0] * 50
        ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh)
        mat = E.build_esn_matrices(OPT_ARCH, ip)
        if mat["kappa0"] >= 0:
            return np.full(len(TARGET_LAG_IDX), 1.0)
        Zs = []
        for q in range(n_paths):
            x_, v_ = E._sim_esn_with_params(seed_base + q, T_eval, OPT_ARCH, mat,
                                             b0_delta=0.0, scale=1.0)
            Zs.append(zumbach_Z(x_, v_, target_lags))
        avg_Z = np.mean(Zs, axis=0)
        return avg_Z - targets

    res = least_squares(resid, np.array(x0), bounds=(BOUNDS_LO, BOUNDS_HI),
                         method="trf", max_nfev=max_nfev,
                         xtol=1e-13, ftol=1e-13, gtol=1e-13, diff_step=0.15)
    H, rs, ll, lh = [float(np.clip(res.x[i], BOUNDS_LO[i], BOUNDS_HI[i])) for i in range(4)]

    ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    x_, v_ = E._sim_esn_with_params(888, T_CONFIRM, OPT_ARCH, mat, b0_delta=0.0, scale=1.0)
    Z_confirm = zumbach_Z(x_, v_, lags_all)

    return dict(H_fit=H, rough_scale_fit=rs, lam_lo_fit=ll, lam_hi_fit=lh,
                Z_confirm=Z_confirm.tolist(), nfev=res.nfev)


def build_confidence_band(H, rs, ll, lh, lags_all, n_reps=60, T_rep=3000, seed_offset=0):
    ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    Z_reps = []
    for rep in range(n_reps):
        seed = 14_100_000 + seed_offset + rep
        x_, v_ = E._sim_esn_with_params(seed, T_rep, OPT_ARCH, mat, b0_delta=0.0, scale=1.0)
        Z_reps.append(zumbach_Z(x_, v_, lags_all))
    Z_reps = np.array(Z_reps)
    return dict(mean=np.nanmean(Z_reps, axis=0).tolist(),
                p5=np.nanpercentile(Z_reps, 5, axis=0).tolist(),
                p95=np.nanpercentile(Z_reps, 95, axis=0).tolist())


if __name__ == "__main__":
    real_prof = json.load(open("bench9_real_zumbach_profile.json"))
    lags_all = real_prof["lags"]

    results = {}
    t0 = time.time()
    for i, name in enumerate(NAMES9):
        real_Z = np.array(real_prof["pR2_fV"][name]) - np.array(real_prof["pV_fR2"][name])
        print(f"{name}:", flush=True)
        best = calibrate_widened_lh(real_Z, lags_all, seed_offset=i * 10000)
        results[name] = best
        print(f"  H={best['H_fit']:.3f} rs={best['rough_scale_fit']:.3f} "
              f"lam_lo={best['lam_lo_fit']:.5f} lam_hi={best['lam_hi_fit']:.4f} "
              f"nfev={best['nfev']} ({time.time()-t0:.0f}s)", flush=True)
        json.dump(results, open("zumbach_widened_lh_calibration.json", "w"), indent=2)

    print("\n=== Building 90% confidence bands ===")
    ci_results = {}
    for i, name in enumerate(NAMES9):
        c = results[name]
        ci_results[name] = build_confidence_band(
            c["H_fit"], c["rough_scale_fit"], c["lam_lo_fit"], c["lam_hi_fit"],
            lags_all, seed_offset=i * 1000)
        print(f"  {name}: done ({time.time()-t0:.0f}s)", flush=True)
        json.dump({"lags": lags_all, "by_asset": ci_results},
                  open("zumbach_widened_lh_ci.json", "w"), indent=2)

    print("\nDone.")
