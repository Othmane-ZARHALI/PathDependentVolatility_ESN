"""
zumbach_components_ci_compute.py
======================================
Builds a 90% confidence band (60 replicate T=3,000-day paths) for the
calibrated ESN's own 2 correlation components separately --
corr(pastReturn^2, futureVol) and corr(pastVol, futureReturn^2) -- the
2 terms whose difference defines Z(L) (Eq. 1 in
zumbach_effect_protocol.tex). Uses the same widened-lambda_hi
calibrated parameters as zumbach_widened_lh_calibration.py.

Requires: esn_base.py (core simulation module) and
zumbach_widened_lh_calibration.json (calibrated parameters).
"""
import sys, json, time
sys.path.insert(0, ".")
import numpy as np
import esn_base as E

OPT_ARCH = dict(
    rough_orientation=-1.0, Nr=22, Nz=29, matrix_seed=E.MATRIX_SEED,
    z_strength=0.27136029620674146, even_strength=3.1842560105665587,
    linear_strength=0.18199009132966565, gamma_norm=1.308098702657785,
    local_z_strength=0.06081431502674114, zz_scale=0.03120698941082984,
    sign_prob_neg=0.222781582556174,
)
P_FIXED = dict(az_lo=1/400, az_hi=1/7, zr_lo=0.055, zr_hi=0.25, m1=0.0, m2=0.0, b0_delta=0.0, scale=1.0)


def _corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() < 1e-14 or b.std() < 1e-14:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def zumbach_components(x, v, lags):
    """Returns the 2 components separately: pR2_fV, pV_fR2."""
    n = len(x)
    pR2_fV, pV_fR2 = [], []
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
            pR2_fV.append(_corr(pR ** 2, fV))
            pV_fR2.append(_corr(pV, fR ** 2))
        else:
            pR2_fV.append(np.nan); pV_fR2.append(np.nan)
    return np.array(pR2_fV), np.array(pV_fR2)


def build_components_ci(H, rs, ll, lh, lags, n_reps=60, T_rep=3000, seed_offset=0):
    ip = dict(P_FIXED, H=H, rough_scale=rs, lam_lo=ll, lam_hi=lh)
    mat = E.build_esn_matrices(OPT_ARCH, ip)
    pR2_fV_reps, pV_fR2_reps = [], []
    for rep in range(n_reps):
        seed = 14_200_000 + seed_offset + rep
        x_, v_ = E._sim_esn_with_params(seed, T_rep, OPT_ARCH, mat, b0_delta=0.0, scale=1.0)
        a, b = zumbach_components(x_, v_, lags)
        pR2_fV_reps.append(a); pV_fR2_reps.append(b)
    pR2_fV_reps = np.array(pR2_fV_reps); pV_fR2_reps = np.array(pV_fR2_reps)
    return dict(
        pR2_fV_mean=np.nanmean(pR2_fV_reps, axis=0).tolist(),
        pR2_fV_p5=np.nanpercentile(pR2_fV_reps, 5, axis=0).tolist(),
        pR2_fV_p95=np.nanpercentile(pR2_fV_reps, 95, axis=0).tolist(),
        pV_fR2_mean=np.nanmean(pV_fR2_reps, axis=0).tolist(),
        pV_fR2_p5=np.nanpercentile(pV_fR2_reps, 5, axis=0).tolist(),
        pV_fR2_p95=np.nanpercentile(pV_fR2_reps, 95, axis=0).tolist(),
    )


NAMES9 = ["LHX", "BK", "HUM", "TSN", "IEX", "KMB", "SP500", "RUSSELL2000", "NASDAQ100"]

if __name__ == "__main__":
    real_prof = json.load(open("bench9_real_zumbach_profile.json"))
    lags = real_prof["lags"]
    cal = json.load(open("zumbach_widened_lh_calibration.json"))

    results = {}
    t0 = time.time()
    for i, name in enumerate(NAMES9):
        c = cal[name]
        results[name] = build_components_ci(
            c["H_fit"], c["rough_scale_fit"], c["lam_lo_fit"], c["lam_hi_fit"],
            lags, seed_offset=i * 1000)
        print(f"{name}: done ({time.time()-t0:.0f}s)", flush=True)
        json.dump({"lags": lags, "by_asset": results},
                  open("zumbach_components_ci.json", "w"), indent=2)

    print("\nDone.")
