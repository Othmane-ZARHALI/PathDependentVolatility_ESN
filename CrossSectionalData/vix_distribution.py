"""
vix_distribution.py
======================================
Builds the data for "VIX empirical distribution: market vs. model"
(replaces the Test 1 term-structure bar chart), all 4 quoted VIX
maturities (27d/55d calibrated, 118d/153d held out).

Market side: the risk-neutral density of VIX_T implied by the VIX
option chain, via Breeden-Litzenberger,
    f_Y(K) = d^2 C(K) / dK^2
(exact identity for undiscounted/forward-measure prices, the same
convention used throughout this note -- no discount factor needed).

2 smoothing steps are needed, not 1: (i) a low-degree polynomial fit
to the market's own mid-IV quotes in log-moneyness (raw
finite-differencing the scattered quotes directly would be far too
noisy for a second derivative); (ii) pricing every strike from a
SINGLE smooth branch (calls throughout the whole strike range, not
the usual put-below/call-above OTM splice) -- put-call parity makes
the 2 choices equal in theory (C-P=F-K is linear, so d^2C/dK^2 =
d^2P/dK^2 exactly), but the OTM splice's first derivative genuinely
jumps by 1 at the crossover strike (differentiating the parity
identity), which is a real kink, not a numerical artifact -- and its
second derivative is correspondingly a spike. Using 1 consistent
branch removes it. The resulting smooth price curve's exact 2nd
derivative is then read off a cubic-spline interpolant
(`CubicSpline(...).derivative(2)`, an analytic derivative of a smooth
curve) rather than a 2nd-order finite difference, avoiding the
separate, smaller noise a double `np.gradient` would add on top.

Model side: the K_outer-draw ensemble of simulated VIX_N values
(vix_ensemble(), already used internally by the smile-fitting code) --
a direct Monte-Carlo sample of the model's own distribution, no
smoothing needed, just a histogram/KDE.
"""
import json
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

from calibrate_smile import build_matrices_raw, vix_ensemble
from model_pricing import calendar_to_trading_days
from market_vix import black76_undisc, build_moneyness_grid

VIX_CSV = "/mnt/user-data/uploads/1788726546415_VIX_data_19_02_2026.csv"


def market_vix_density(expiry, T_years, poly_degree=5, K_grid_step=0.25, pad=1.0):
    df = pd.read_csv(VIX_CSV)
    d = df[df.expiry == expiry].copy()
    d['mid'] = (d.bid + d.ask) / 2
    d = d.dropna(subset=['mid']).sort_values('strike')
    F = float(d.forward.iloc[0])

    # (i) smooth IV(log-moneyness) fit -- a single global polynomial,
    # not a knot spline, so there is no risk of per-knot wiggle.
    x = np.log(d.strike.values / F)
    y = d.mid.values
    coefs = np.polyfit(x, y, poly_degree)

    K_lo, K_hi = d.strike.min() + pad, d.strike.max() - pad
    Ks = np.arange(K_lo, K_hi, K_grid_step)
    ivs = np.clip(np.polyval(coefs, np.log(Ks / F)), 0.05, 5.0)

    # (ii) ONE smooth branch (calls) across the whole grid -- avoids the
    # real put/call-splice kink at F; then an analytic 2nd derivative
    # of a spline through the (already smooth) price curve.
    prices = np.array([black76_undisc(F, K, T_years, iv, True) for K, iv in zip(Ks, ivs)])
    cs = CubicSpline(Ks, prices)
    density = np.clip(cs(Ks, 2), 0.0, None)
    area = np.trapezoid(density, Ks)
    if area > 0:
        density = density / area
    return dict(K=Ks.tolist(), density=density.tolist(), F=F, T_years=T_years)


if __name__ == "__main__":
    theta_hat = json.load(open("../results/theta_hat_final.json"))['unpacked']
    mat = build_matrices_raw(theta_hat['H'], theta_hat['rough_scale'], theta_hat['lam_lo'],
                              theta_hat['lam_hi'], theta_hat['m1'], theta_hat.get('m2', 0.0))
    pricing_grid = build_moneyness_grid(0.60, 1.60, 0.85, 1.20, 0.02, 0.06)

    out = {}
    for label, cal_days, expiry, seed in [("near", 27, 46099, 130_000_000),
                                            ("next", 55, 46127, 131_000_000),
                                            ("mid1", 118, 46190, 132_000_000),
                                            ("mid2", 153, 46225, 133_000_000)]:
        n = calendar_to_trading_days(cal_days)
        mkt = market_vix_density(expiry, cal_days / 365)
        ens = vix_ensemble(mat, theta_hat['alpha'], theta_hat['b0_delta'], n, 500,
                            1000, 1000, pricing_grid, seed)
        out[label] = dict(market=mkt, model_ensemble=ens.tolist(),
                           model_mean=float(ens.mean()), model_median=float(np.median(ens)))
        print(label, "model mean/median:", out[label]['model_mean'], out[label]['model_median'],
              "market F:", mkt['F'])
        with open("../results/vix_distribution.json", "w") as f:
            json.dump(out, f, indent=2)

