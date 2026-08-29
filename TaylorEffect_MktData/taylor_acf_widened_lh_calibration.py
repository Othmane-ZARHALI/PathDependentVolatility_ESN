"""
taylor_acf_widened_lh_calibration.py
======================================
Extends taylor_acf_lagbylag_calibration.py: same lag-by-lag objective
(match acf(|x|) and acf(x^2) jointly at 6 lags, 1,3,5,10,15,20), but
with lambda_hi's bound widened from the standard [1.0, 5.0] to
[0.001, 5.0].

WHY: a direct sweep at SP500's own hand-checked "extreme" point
(H=0.3, rough_scale=0.7, lam_lo=1e-5) showed that lambda_hi=1.0 (the
standard lower bound used throughout this whole note family) was
ITSELF limiting how much persistence the architecture could produce --
pushing lambda_hi down to 0.05 raised acf(|x|) at lag 40 from 0.117 to
0.252, close to real SP500's own 0.228 at lag 20. This is the same
"a parameter's search bound, not the model, was the limiting factor"
pattern found for H in the leverage-effect note (there, H's bound of
0.3 was too narrow; here, lambda_hi's own LOWER bound of 1.0 is too
high).

Requires: esn_base.py (core simulation module) and
bench9_real_abs_acf.json / bench9_real_sq_acf.json (real per-asset ACF
arrays, from taylor_effect_evaluation.py).
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

FULL_LAGS = list(range(1, 41))
TARGET_LAGS = [1, 3, 5, 10, 15, 20]
X0 = (0.3, 0.7, 1e-5, 0.05)  # informed start at the hand-found good region
N_PATHS_EVAL = 4
T_EVAL = 3000
MAX_NFEV = 40
T_CONFIRM = 12000
NAMES9 = ["LHX", "BK", "HUM", "TSN", "IEX", "KMB", "SP500", "RUSSELL2000", "NASDAQ100"]


def acf(series, lags):
    return np.array([np.corrcoef(series[:-lag], series[lag:])[0, 1] for lag in lags])


def calibrate_widened_lh(name, real_abs_full, real_sq_full, x0=X0,
                          n_paths=N_PATHS_EVAL, T_eval=T_EVAL, max_nfev=MAX_NFEV,
                          seed_offset=0):
    idx = [l - 1 for l in TARGET_LAGS]
    targets = np.concatenate([real_abs_full[idx], real_sq_full[idx]])
    counter = [0]

    def resid(x):
        H, rs, ll, lh = x
        H = float(np.clip(H, *H_BOUNDS)); rs = float(np.clip(rs, *RS_BOUNDS))
        ll = float(np.clip(ll, *LL_BOUNDS)); lh = float(np.clip(lh, *LH_BOUNDS))
        lh = max(lh, ll * 1.5)
        counter[0] += 1
        seed_base = 12_500_000 + seed_offset + counter[0] * 50
        ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh)
        mat = E.build_esn_matrices(OPT_ARCH, ip)
        if mat["kappa0"] >= 0:
            return np.full(len(targets), 1.0)
        abs_list, sq_list = [], []
        for q in range(n_paths):
            x_, v_ = E._sim_esn_with_params(seed_base + q, T_eval, OPT_ARCH, mat,
                                             b0_delta=0.0, scale=1.0)
            abs_list.append(acf(np.abs(x_), TARGET_LAGS))
            sq_list.append(acf(x_ ** 2, TARGET_LAGS))
        avg = np.concatenate([np.mean(abs_list, axis=0), np.mean(sq_list, axis=0)])
        return avg - targets

    res = least_squares(resid, np.array(x0), bounds=(BOUNDS_LO, BOUNDS_HI),
                         method="trf", max_nfev=max_nfev,
                         xtol=1e-13, ftol=1e-13, gtol=1e-13, diff_step=0.15)
    H, rs, ll, lh = [float(np.clip(res.x[i], BOUNDS_LO[i], BOUNDS_HI[i])) for i in range(4)]

    ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    x_, v_ = E._sim_esn_with_params(888, T_CONFIRM, OPT_ARCH, mat, b0_delta=0.0, scale=1.0)
    abs_confirm = acf(np.abs(x_), FULL_LAGS)
    sq_confirm = acf(x_ ** 2, FULL_LAGS)

    return dict(H_fit=H, rough_scale_fit=rs, lam_lo_fit=ll, lam_hi_fit=lh,
                abs_confirm=abs_confirm.tolist(), sq_confirm=sq_confirm.tolist(),
                nfev=res.nfev)


def build_confidence_band(H, rs, ll, lh, n_reps=60, T_rep=3000, seed_offset=0):
    ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    abs_profs, sq_profs = [], []
    for rep in range(n_reps):
        seed = 12_600_000 + seed_offset + rep
        x_, v_ = E._sim_esn_with_params(seed, T_rep, OPT_ARCH, mat, b0_delta=0.0, scale=1.0)
        abs_profs.append(acf(np.abs(x_), FULL_LAGS))
        sq_profs.append(acf(x_ ** 2, FULL_LAGS))
    abs_profs = np.array(abs_profs); sq_profs = np.array(sq_profs)
    return dict(
        abs_mean=np.nanmean(abs_profs, axis=0).tolist(),
        abs_p5=np.nanpercentile(abs_profs, 5, axis=0).tolist(),
        abs_p95=np.nanpercentile(abs_profs, 95, axis=0).tolist(),
        sq_mean=np.nanmean(sq_profs, axis=0).tolist(),
        sq_p5=np.nanpercentile(sq_profs, 5, axis=0).tolist(),
        sq_p95=np.nanpercentile(sq_profs, 95, axis=0).tolist(),
    )


if __name__ == "__main__":
    real_abs = json.load(open("bench9_real_abs_acf.json"))
    real_sq = json.load(open("bench9_real_sq_acf.json"))

    results = {}
    t0 = time.time()
    for i, name in enumerate(NAMES9):
        print(f"{name}:", flush=True)
        best = calibrate_widened_lh(name, np.array(real_abs[name]), np.array(real_sq[name]),
                                     seed_offset=i * 10000)
        results[name] = best
        print(f"  H={best['H_fit']:.3f} rs={best['rough_scale_fit']:.3f} "
              f"lam_lo={best['lam_lo_fit']:.5f} lam_hi={best['lam_hi_fit']:.4f} "
              f"nfev={best['nfev']} ({time.time()-t0:.0f}s)", flush=True)
        json.dump(results, open("taylor_acf_widened_lh_calibration.json", "w"), indent=2)

    print("\n=== Building 90% confidence bands ===")
    ci_results = {}
    for i, name in enumerate(NAMES9):
        c = results[name]
        ci_results[name] = build_confidence_band(
            c["H_fit"], c["rough_scale_fit"], c["lam_lo_fit"], c["lam_hi_fit"],
            seed_offset=i * 1000)
        print(f"  {name}: done ({time.time()-t0:.0f}s)", flush=True)
        json.dump({"lags": FULL_LAGS, "by_asset": ci_results},
                  open("taylor_acf_widened_lh_ci.json", "w"), indent=2)

    print("\nDone.")
