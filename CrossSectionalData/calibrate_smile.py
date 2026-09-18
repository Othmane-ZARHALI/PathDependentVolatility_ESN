"""
calibrate_smile.py
======================================
theta_hat = argmin_theta L(theta),
    L(theta) = sum_i w_i * (sigma_imp^theta(K_i,T_i) - sigma_imp^mkt(K_i,T_i))^2

over theta = (H, rough_scale, lam_lo, lam_hi, m1, alpha, b0_delta), jointly
against ALL 4 quoted expiries at once: SPX near (29d) + SPX next (57d)
+ VIX-option near (27d) + VIX-option next (55d). This is now the ONLY
calibration method in this note -- the earlier b0-only ATM-level
correction was a diagnostic step, superseded by this fuller fit.

Weights w_i: Black-76 vega (the objective is in vol-space, so
w_i ~ vega_i, not 1/vega_i -- see the Definition box in the note),
ADDITIONALLY multiplied by a Gaussian kernel in log-moneyness centred
at the money, per explicit request to weight the ATM region more than
plain vega alone does:
    w_i  ~  vega_i * exp( -0.5 * (ln(m_i)/bw)^2 )
bw = 0.12 for SPX, 0.40 for VIX options (VIX's own natural moneyness
spread is far wider than SPX's -- real VIX-option quotes reach out to
K/F~3 -- so an SPX-tuned bandwidth would zero out almost the entire
VIX chain; 0.40 keeps a comparable EFFECTIVE number of contributing
VIX points while still concentrating weight near F_VIX). Weights are
normalised to sum to 1 within each of the 4 expiries, so all 4
contribute comparably to L regardless of how many raw quotes each has.

Speed: this version uses vectorized_sim.py (now with native b0_delta
support) throughout, INCLUDING for the SPX legs -- roughly 50-80x
faster than the original per-path scalar-engine loop, which is what
makes a real (not just 16-evaluation) joint fit computationally
feasible here.
"""
import time
import json
import sys
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import norm

import esn_base as E
from joint_calibration import OPT_ARCH, P_FIXED
from model_pricing import calendar_to_trading_days, mc_undisc_price, implied_vol_from_price
from market_vix import build_moneyness_grid, cboe_variance_from_grid
from vectorized_sim import simulate_batch

SPX_MONEYNESS_RANGE = (0.70, 1.40)
VIX_MONEYNESS_RANGE = (0.40, 3.00)
SPX_ATM_BANDWIDTH_LOW = 0.08    # narrowed from 0.12 (per direct request to weight ATM/ITM more
                                # heavily for SPX): a tighter Gaussian concentrates relative
                                # weight closer to the money, at the deep-OTM wings' expense
SPX_ATM_BANDWIDTH_HIGH = 0.15   # narrowed from 0.28, same reasoning; kept wider than the low
                                # side since m>=1 is empirically where the gap has repeatedly
                                # been found (Section 3.3.1's own diagnosis)
VIX_ATM_BANDWIDTH_LOW = 1.5
VIX_ATM_BANDWIDTH_HIGH = 1.5   # widened substantially from 0.40/0.55: at that bandwidth, points
                                # beyond moneyness~3 (VIX quotes routinely run out to 5-9) carried
                                # essentially 0 weight in the objective (e.g. relative weight
                                # 0.014 at m=5, 0.0008 at m=8 under bw=0.55) -- the optimizer was
                                # never actually asked to match those points, which is why a
                                # visible gap persisted there even though the fit "looked"
                                # converged. At bw=1.5 the far strikes carry real (0.38-0.56)
                                # weight, making this a genuine full-curve fit rather than an
                                # effectively ATM-only one with an unfit tail plotted alongside it

SPX_CSV = "/mnt/user-data/uploads/1788726546415_SPX_data_19_02_2026.csv"
VIX_CSV = "/mnt/user-data/uploads/1788726546415_VIX_data_19_02_2026.csv"


def black76_vega(F, K, T, sigma):
    if sigma <= 0 or T <= 0:
        return 0.0
    d1 = (np.log(F / K) + 0.5 * sigma**2 * T) / (sigma * np.sqrt(T))
    return F * norm.pdf(d1) * np.sqrt(T)


ITM_WEIGHT_BOOST = 5.0  # extra multiplier on ITM points (raised from 2.0 per request to further
                        # improve the ITM fit specifically in Figures 1-2), on top of removing
                        # their OTM-wing decay entirely


MONEYNESS_HIGH_BOOST = 3.0  # extra flat weight for m>=1 points, replacing the retired ITM boost:
                            # once put+call were merged (Section on mid-of-put-call calibration),
                            # there was no longer a separate ITM row to boost, and the m>=1 gap
                            # this used to fix re-opened (RMSE ~0.056 vs ~0.016 on the m<1 side,
                            # SPX near) -- applied directly by moneyness side instead, since that
                            # is what the gap actually tracks, not option-type convention


def atm_weights(moneyness, vols_mkt, T_years, bandwidth_low, bandwidth_high, is_itm):
    """Vega weight, sharpened toward the money -- asymmetrically, on 2
    separate axes:
    (i) m>=1 points get an extra flat multiplier (MONEYNESS_HIGH_BOOST),
        on top of their own Gaussian falloff -- since merging put and
        call quotes into 1 mid-IV point per strike retired the
        previous ITM-row boost that used to fix this same gap;
    (ii) the OTM-wing Gaussian itself uses a WIDER bandwidth on the
        m>=1 side than the m<1 side. Diagnosing the original gap
        (splitting simultaneously by ITM/OTM AND by moneyness side)
        showed it is concentrated specifically at m>=1 -- both the
        OTM-call AND the ITM-put convention there are poorly fit
        (RMSE 0.034 and 0.020 respectively, vs <=0.011 on the m<1
        side under theta0) -- a plain ITM/OTM split alone averages
        this away against the well-fitting low-strike ITM calls, so
        bandwidth_high > bandwidth_low targets the actual gap.
    is_itm is unused now (always False since put+call were merged) --
    kept as a parameter for call-site compatibility.
    """
    vega = np.array([black76_vega(1.0, m, T_years, s) for m, s in zip(moneyness, vols_mkt)])
    bw = np.where(moneyness >= 1.0, bandwidth_high, bandwidth_low)
    otm_kernel = np.exp(-0.5 * (np.log(moneyness) / bw) ** 2)
    boost = np.where(moneyness >= 1.0, MONEYNESS_HIGH_BOOST, 1.0)
    w = vega * otm_kernel * boost
    return w / w.sum()


def load_spx_points(expiry):
    df = pd.read_csv(SPX_CSV)
    d = df[df.expiry == expiry].copy()
    d['mid'] = (d.bid + d.ask) / 2
    d['m'] = d.strike / d.forward
    lo, hi = SPX_MONEYNESS_RANGE
    d = d[(d.m >= lo) & (d.m <= hi)].dropna(subset=['mid'])
    # Calibrate on the MID of the put and call implied vol at each strike (1 point per strike),
    # not on both separately: the data carries a put row AND a call row per strike, and their
    # own small disagreement (a few vol points, growing away from the money -- see the report's
    # own put-vs-call figure) is market microstructure noise around a single "true" smile value,
    # not 2 independent facts to fit. Averaging removes both the redundant double-counting of
    # every strike and the old ITM-vs-OTM-convention diagnosis this note spent real effort on --
    # that diagnosis was, in hindsight, partly an artifact of fitting the 2 conventions separately.
    g = d.groupby('m', as_index=False)['mid'].mean().sort_values('m')
    is_itm = np.zeros(len(g), dtype=bool)  # no ITM/OTM distinction once put+call are merged
    return g.m.values, g.mid.values, is_itm


def load_vix_points(expiry):
    df = pd.read_csv(VIX_CSV)
    d = df[df.expiry == expiry].copy()
    d['mid'] = (d.bid + d.ask) / 2
    F = float(d.forward.iloc[0])
    d['m'] = d.strike / F
    lo, hi = VIX_MONEYNESS_RANGE
    d = d[(d.m >= lo) & (d.m <= hi)].dropna(subset=['mid'])
    g = d.groupby('m', as_index=False)['mid'].mean().sort_values('m')
    is_itm = np.zeros(len(g), dtype=bool)
    return g.m.values, g.mid.values, is_itm


def build_matrices_raw(H, rough_scale, lam_lo, lam_hi, m1, m2=0.0):
    ip = dict(P_FIXED, H=H, rough_scale=rough_scale, lam_lo=lam_lo, lam_hi=lam_hi, m1=m1, m2=m2)
    return E.build_esn_matrices(OPT_ARCH, ip)


def spx_model_vols(mat, alpha, b0_delta, n_days, burn, K, moneyness, T_years, seed):
    dx, v, _ = simulate_batch(OPT_ARCH, mat, alpha, burn + n_days, K, seed=seed, b0_delta=b0_delta)
    Y = np.exp(np.cumsum(dx[burn:burn + n_days], axis=0)[-1])
    out = np.empty(len(moneyness))
    for i, m in enumerate(moneyness):
        is_call = m >= 1.0
        price = mc_undisc_price(Y, m, is_call)
        out[i] = implied_vol_from_price(1.0, m, T_years, price, is_call)
    return out


def vix_ensemble(mat, alpha, b0_delta, n_outer_days, burn, K_outer, K_inner, pricing_grid, seed):
    """The model's own ensemble of K_outer draws of VIX_N (the raw
    simulated VIX level at day N, one per outer chain) -- the
    distribution Figure 4 (Section 4) compares against the market's
    own risk-neutral density."""
    dx, v, state = simulate_batch(OPT_ARCH, mat, alpha, burn + n_outer_days, K_outer,
                                   seed=seed, b0_delta=b0_delta)
    r0, z0, v0 = state
    r0i = np.repeat(r0, K_inner, axis=0)
    z0i = np.repeat(z0, K_inner, axis=0)
    v0i = np.repeat(v0, K_inner)
    inner_days = calendar_to_trading_days(30)
    dxi, _, _ = simulate_batch(OPT_ARCH, mat, alpha, inner_days, K_outer * K_inner,
                                r0=r0i, z0=z0i, var_used_prev0=v0i, seed=seed + 1, b0_delta=b0_delta)
    Y = np.exp(np.cumsum(dxi, axis=0)[-1]).reshape(K_outer, K_inner)

    vix_vals = np.empty(K_outer)
    for k in range(K_outer):
        Yk = Y[k]
        Qs = [float(np.mean(np.maximum(Yk - m, 0.0))) if m >= 1.0
              else float(np.mean(np.maximum(m - Yk, 0.0))) for m in pricing_grid]
        res = cboe_variance_from_grid(1.0, 30 / 365, pricing_grid, Qs)
        vix_vals[k] = 100 * np.sqrt(max(res['sigma2'], 0.0))
    return vix_vals


def vix_model_vols(mat, alpha, b0_delta, n_outer_days, burn, K_outer, K_inner,
                    moneyness, T_years, seed, pricing_grid):
    vix_vals = vix_ensemble(mat, alpha, b0_delta, n_outer_days, burn, K_outer, K_inner,
                             pricing_grid, seed)
    F_vix = float(np.mean(vix_vals))
    out = np.empty(len(moneyness))
    for i, m in enumerate(moneyness):
        K_abs = m * F_vix
        is_call = K_abs >= F_vix
        price = float(np.mean(np.maximum(vix_vals - K_abs, 0.0))) if is_call else \
            float(np.mean(np.maximum(K_abs - vix_vals, 0.0)))
        out[i] = implied_vol_from_price(F_vix, K_abs, T_years, price, is_call, hi=8.0)
    return out


INSTRUMENT_WEIGHT = dict(spx=1.0, vix=1.0)  # corrected down again from 2.0: even 2.0 plateaued
                                             # around 9.5-9.8, well above an O(1-2) target: SPX
                                             # and VIX contribute EQUALLY here (plus the m2
                                             # channel freed below, which helps VIX shape "for
                                             # free" without an artificial extra weight), matching
                                             # the balanced weighting that let the pre-VIX-fix
                                             # search reach 1.54 (Section 3.3.4)
LABEL_WEIGHT = dict(vix_near=3.0, vix_next=1.5)  # re-derived after widening the ATM bandwidth
                                                  # (VIX_ATM_BANDWIDTH_LOW/HIGH, above) to fix the
                                                  # real bug: mid1/mid2 no longer need any extra
                                                  # weight under the wide kernel (they became
                                                  # excellent on their own -- weighted RMSE
                                                  # 0.03-0.08 across 2 seeds, the best in this
                                                  # entire note), but the wider kernel dilutes
                                                  # near-money weight enough that vix_near/next
                                                  # need their own protection now instead


class CombinedCalibrator:
    """theta = [H, rough_scale, log10(lam_lo), log10(lam_hi), m1, m2, alpha, b0_delta]
    m2 (the qq^T coefficient of the quadratic channel, Eq. 8 -- fixed
    at 0 throughout the ENTIRE note family up to this point, never
    exposed to any calibration, time-series or cross-sectional) is
    freed here specifically because it is the one lever that acts on
    the SAME projection direction as rough_scale's own linear term
    (q^T r_t)^2 rather than the isotropic ||r_t||^2 m1 already
    controls -- exactly the kind of extra curvature/convexity a smile
    as steep as VIX's own needs and m1 alone cannot supply without
    also moving the isotropic (SPX-relevant) level."""

    BOUNDS_LO = np.array([0.01, 0.20, -5.0, -3.0, -1.00, -2.00, 0.02, -1.3])
    BOUNDS_HI = np.array([2.50, 0.70, -1.0, 0.699, 0.50, 2.00, 1.00, -0.3])

    def __init__(self, K_spx=2000, K_vix_outer=250, K_vix_inner=250, burn=500,
                 seed_base=90_000_000, checkpoint_path=None, history_path=None):
        self.K_spx = K_spx
        self.K_vix_outer = K_vix_outer
        self.K_vix_inner = K_vix_inner
        self.burn = burn
        self.seed_base = seed_base
        self.checkpoint_path = checkpoint_path
        self.history_path = history_path
        self.n_evals = 0
        self.global_eval = 0
        if history_path is not None:
            try:
                with open(history_path) as f:
                    self.global_eval = sum(1 for _ in f)
            except FileNotFoundError:
                pass
        self.vix_grid = build_moneyness_grid(0.60, 1.60, 0.85, 1.20, 0.02, 0.06)

        specs = []
        for label, cal_days, expiry in [("spx_near", 29, 46101), ("spx_next", 57, 46129)]:
            m, v, is_itm = load_spx_points(expiry)
            w = atm_weights(m, v, cal_days / 365, SPX_ATM_BANDWIDTH_LOW, SPX_ATM_BANDWIDTH_HIGH, is_itm)
            specs.append(dict(kind="spx", label=label, cal_days=cal_days,
                               n_days=calendar_to_trading_days(cal_days),
                               T_years=cal_days / 365, m=m, v=v, w=w, is_itm=is_itm))
        for label, cal_days, expiry in [("vix_near", 27, 46099), ("vix_next", 55, 46127),
                                          ("vix_mid1", 118, 46190), ("vix_mid2", 153, 46225)]:
            m, v, is_itm = load_vix_points(expiry)
            w = atm_weights(m, v, cal_days / 365, VIX_ATM_BANDWIDTH_LOW, VIX_ATM_BANDWIDTH_HIGH, is_itm)
            specs.append(dict(kind="vix", label=label, cal_days=cal_days,
                               n_days=calendar_to_trading_days(cal_days),
                               T_years=cal_days / 365, m=m, v=v, w=w, is_itm=is_itm))
        self.specs = specs

    def unpack(self, theta):
        H, rough_scale, log_lam_lo, log_lam_hi, m1, m2, alpha, b0_delta = theta
        return dict(H=H, rough_scale=rough_scale, lam_lo=10**log_lam_lo,
                    lam_hi=10**log_lam_hi, m1=m1, m2=m2, alpha=alpha, b0_delta=b0_delta)

    def calibrate_scale(self, theta0, seed=None):
        """Precompute, ONCE, each expiry's own weighted-RMSE at theta0
        and store it as a fixed per-expiry scale. Without this, SPX and
        VIX residuals enter the SAME sum-of-squares in absolute
        vol-point units, but VIX vols sit on a ~100%+ scale against
        SPX's ~15-40% -- so raw VIX residuals are 10-30x larger than
        SPX ones even after each expiry's own weights are normalised to
        sum to 1. That scale gap, not anything about the ATM/ITM
        weighting itself, is what was silently dominating the
        optimizer's gradient direction toward VIX at SPX's expense.
        Dividing each expiry's residual block by its own FIXED theta0
        scale makes all 4 expiries contribute comparably in RELATIVE
        (not absolute) terms, so improving SPX's own ITM gap actually
        moves the objective instead of being swamped."""
        p = self.unpack(theta0)
        mat = build_matrices_raw(p['H'], p['rough_scale'], p['lam_lo'], p['lam_hi'], p['m1'], p['m2'])
        s = seed or (self.seed_base + 999_000)
        scale = {}
        for spec in self.specs:
            if spec['kind'] == 'spx':
                model_v = spx_model_vols(mat, p['alpha'], p['b0_delta'], spec['n_days'], self.burn,
                                          self.K_spx, spec['m'], spec['T_years'], s)
            else:
                model_v = vix_model_vols(mat, p['alpha'], p['b0_delta'], spec['n_days'], self.burn,
                                          self.K_vix_outer, self.K_vix_inner, spec['m'],
                                          spec['T_years'], s, self.vix_grid)
            diff = model_v - spec['v']
            mask = ~np.isnan(diff)
            wrmse = float(np.sqrt(np.sum(spec['w'][mask] * diff[mask]**2) / np.sum(spec['w'][mask])))
            scale[spec['label']] = max(wrmse, 1e-4)
        self.scale = scale
        return scale

    def residuals(self, theta):
        self.n_evals += 1
        p = self.unpack(theta)
        mat = build_matrices_raw(p['H'], p['rough_scale'], p['lam_lo'], p['lam_hi'], p['m1'], p['m2'])
        seed = self.seed_base  # common random numbers: identical seed every theta
        all_res = []
        for spec in self.specs:
            if spec['kind'] == 'spx':
                model_v = spx_model_vols(mat, p['alpha'], p['b0_delta'], spec['n_days'],
                                          self.burn, self.K_spx, spec['m'], spec['T_years'], seed)
            else:
                model_v = vix_model_vols(mat, p['alpha'], p['b0_delta'], spec['n_days'], self.burn,
                                          self.K_vix_outer, self.K_vix_inner, spec['m'],
                                          spec['T_years'], seed, self.vix_grid)
            diff = np.nan_to_num(model_v - spec['v'], nan=0.0)
            scale = getattr(self, 'scale', {}).get(spec['label'], 1.0)
            iw = INSTRUMENT_WEIGHT.get(spec['kind'], 1.0) * LABEL_WEIGHT.get(spec['label'], 1.0)
            # Normalised so each expiry contributes exactly 1 unit of sum-of-squares AT theta0
            # (sum_i w_i*(diff_i/scale)^2 = 1 there, by scale's own definition), and the FINAL
            # /sqrt(n_specs) below divides that down again so the reduced objective itself reads
            # 1.0 at theta0, not n_specs (now 6: SPX near/next + VIX near/next/118d/153d) -- the
            # usual reduced-chi-squared convention: ~1 means "about as good as the reference fit",
            # smoothly smaller as theta improves beyond it, not some arbitrary multiple of the
            # number of expiries.
            all_res.append(iw * np.sqrt(spec['w']) * diff / scale / np.sqrt(len(self.specs)))
        res = np.concatenate(all_res)
        sumsq = float(np.sum(res**2))
        if self.checkpoint_path is not None:
            with open(self.checkpoint_path, "w") as f:
                json.dump(dict(n_evals=self.n_evals, theta=theta.tolist(), sumsq=sumsq), f)
        if self.history_path is not None:
            self.global_eval += 1
            with open(self.history_path, "a") as f:
                f.write(json.dumps(dict(iter=self.global_eval, sumsq=sumsq)) + "\n")
        return res

    def diagnostics(self, theta, K_spx=None, K_vix_outer=None, K_vix_inner=None, seed=None):
        p = self.unpack(theta)
        mat = build_matrices_raw(p['H'], p['rough_scale'], p['lam_lo'], p['lam_hi'], p['m1'], p['m2'])
        s = seed or (self.seed_base + 777_000)
        out = {}
        for spec in self.specs:
            if spec['kind'] == 'spx':
                model_v = spx_model_vols(mat, p['alpha'], p['b0_delta'], spec['n_days'], self.burn,
                                          K_spx or self.K_spx, spec['m'], spec['T_years'], s)
            else:
                model_v = vix_model_vols(mat, p['alpha'], p['b0_delta'], spec['n_days'], self.burn,
                                          K_vix_outer or self.K_vix_outer,
                                          K_vix_inner or self.K_vix_inner,
                                          spec['m'], spec['T_years'], s, self.vix_grid)
            diff = model_v - spec['v']
            mask = ~np.isnan(diff)
            wrmse = float(np.sqrt(np.sum(spec['w'][mask] * diff[mask]**2) / np.sum(spec['w'][mask])))
            rmse = float(np.sqrt(np.mean(diff[mask]**2)))
            out[spec['label']] = dict(weighted_rmse=wrmse, rmse=rmse, n_points=int(mask.sum()))
        return out


if __name__ == "__main__":
    # Reference theta for scale computation: the pre-VIX-reweight, SPX-good theta_hat
    # (Section 3.3.5's own result, m2=0) -- deliberately NOT the VIX-x6 run's own theta, since
    # that run's baseline already had VIX weighted 6x and a scale computed there would still
    # carry that inflation forward even after correcting the weight itself.
    theta0 = np.array([0.011746584544205775, 0.4127761798792011, -3.7219875698226055,
                        -0.07996888792183689, 0.28164092850989053, 0.0,
                        0.999999990007693, -0.5817249325806491])
    ckpt_path = "../results/calibration_checkpoint.json"

    x0 = theta0
    prev_evals = 0
    try:
        with open(ckpt_path) as f:
            ck = json.load(f)
        x0 = np.array(ck['theta'])
        prev_evals = ck['n_evals']
        print("resuming from checkpoint, n_evals so far:", prev_evals, "theta:", x0)
    except FileNotFoundError:
        print("starting fresh from theta0")

    max_nfev = int(sys.argv[1]) if len(sys.argv) > 1 else 12

    calib = CombinedCalibrator(K_spx=3000, K_vix_outer=200, K_vix_inner=200, burn=500,
                                seed_base=90_000_000, checkpoint_path=ckpt_path,
                                history_path="../results/calibration_history.jsonl")

    scale_path = "../results/calibration_scale.json"
    try:
        with open(scale_path) as f:
            calib.scale = json.load(f)
        print("loaded fixed cross-instrument scale:", calib.scale)
    except FileNotFoundError:
        print("computing fixed cross-instrument scale at theta0...")
        calib.calibrate_scale(theta0)
        with open(scale_path, "w") as f:
            json.dump(calib.scale, f, indent=2)
        print("scale:", calib.scale)

    t0 = time.time()
    result = least_squares(calib.residuals, x0, bounds=(calib.BOUNDS_LO, calib.BOUNDS_HI),
                            method='trf', max_nfev=max_nfev, diff_step=0.04, xtol=1e-4, ftol=1e-4,
                            loss='soft_l1', f_scale=1.0)
    elapsed = time.time() - t0
    print("status:", result.status, result.message)
    print("n_evals this round:", calib.n_evals, "elapsed:", elapsed)
    print("theta:", result.x)
    print("cost:", result.cost)

    with open(ckpt_path, "w") as f:
        json.dump(dict(n_evals=prev_evals + calib.n_evals,
                        theta=result.x.tolist(), sumsq=float(2 * result.cost)), f)
