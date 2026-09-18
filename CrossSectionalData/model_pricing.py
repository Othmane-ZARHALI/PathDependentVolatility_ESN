"""
model_pricing.py
======================================
Turns the SP500-calibrated ESN into an option-pricing engine so its
own cross-sectional output (implied-vol smile, 30-day VIX) can be
compared directly against the real SPX/VIX chains -- the only 2 tests
kept in this note.

Why this is a legitimate pricing measure without further adjustment
-----------------------------------------------------------------
At every simulated day, conditional on the state r_t (hence on
sigma_t^used), the one-day log-return is
    dx_{t+1} = -0.5*sigma_t^2 dt + sigma_t sqrt(dt) * eps_t,   eps_t ~ N(0,1)
so E[e^{dx_{t+1}} | F_t] = 1 exactly: the simulator is already a
discrete-time exponential martingale under its own measure. Chaining
this over N days, S_t = S_0 * exp(cumsum dx) is itself a martingale,
so it can be used directly as a "risk-neutral" path simulator for
Monte-Carlo option pricing with forward F_model = 1 (in moneyness
units) -- no risk-premium adjustment, no re-weighting, is needed or
assumed.

Method
------
1.  Burn in each of K independent reservoir chains for BURN days.
2.  Continue simulating; record cumulative log-return at 2 checkpoints
    N1, N2 (trading-day equivalents of the market's near/next expiries).
3.  Y_T = exp(cumulative log-return) is then an unconditional (i.e.
    integrated over the reservoir's own stationary state distribution)
    Monte-Carlo sample of the terminal price ratio S_T/S_0.
4.  Undiscounted OTM option prices on a moneyness grid are plain MC
    averages: C(m) = mean(max(Y_T-m,0)), P(m) = mean(max(m-Y_T,0)).
5.  These feed EXACTLY the same cboe_variance_from_grid() used on the
    market's own quotes (market_vix.py) -- for the VIX number -- and,
    strike-by-strike, a Black-76 inversion -- for the smile plot.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

import esn_base as E
from joint_calibration import OPT_ARCH, P_FIXED
from market_vix import black76_undisc, cboe_variance_from_grid, build_moneyness_grid

TRADING_DAYS_PER_YEAR = 252.0
CALENDAR_DAYS_PER_YEAR = 365.0


def calendar_to_trading_days(n_calendar):
    return max(1, round(n_calendar * TRADING_DAYS_PER_YEAR / CALENDAR_DAYS_PER_YEAR))


def build_matrices(params):
    ip = dict(P_FIXED, H=params['H'], rough_scale=params['rough_scale'],
              lam_lo=params['lam_lo'], lam_hi=params['lam_hi'], m1=params['m1'])
    return E.build_esn_matrices(OPT_ARCH, ip)


def simulate_terminal_ratios(mat, alpha, checkpoints, K, burn, seed0):
    """K independent chains, each burned in for `burn` days then run
    to max(checkpoints) more days; returns {n: Y_T array of length K}
    for every n in checkpoints, Y_T = exp(cumulative log-return over
    the first n post-burn-in days)."""
    n_max = max(checkpoints)
    out = {n: np.empty(K) for n in checkpoints}
    for k in range(K):
        x, _, _ = E._sim_esn_with_params(seed0 + k, burn + n_max, OPT_ARCH, mat,
                                          b0_delta=0.0, scale=1.0, alpha=alpha)
        cum = np.cumsum(x[burn:burn + n_max])
        for n in checkpoints:
            out[n][k] = np.exp(cum[n - 1])
    return out


def mc_undisc_price(Y, m, is_call):
    if is_call:
        return float(np.mean(np.maximum(Y - m, 0.0)))
    return float(np.mean(np.maximum(m - Y, 0.0)))


def implied_vol_from_price(F, K, T, price, is_call, lo=1e-4, hi=6.0):
    intrinsic = max(0.0, (F - K) if is_call else (K - F))
    if price <= intrinsic + 1e-12:
        return np.nan
    f = lambda s: black76_undisc(F, K, T, s, is_call) - price
    try:
        if f(lo) > 0 or f(hi) < 0:
            return np.nan
        return brentq(f, lo, hi, xtol=1e-6)
    except ValueError:
        return np.nan


def model_smile(Y, T_years, moneyness_grid):
    """Model-implied vol smile at one horizon: OTM MC price -> Black-76
    inversion, per strike, using F=1 (moneyness units)."""
    vols = []
    for m in moneyness_grid:
        is_call = m >= 1.0
        price = mc_undisc_price(Y, m, is_call)
        vols.append(implied_vol_from_price(1.0, m, T_years, price, is_call))
    return np.array(vols)


def model_cboe_sigma2(Y, T_years, moneyness_grid):
    """Model's own sigma^2 at one horizon, via the IDENTICAL CBOE
    discretized-strike replication used on the market's own quotes,
    fed directly with MC undiscounted OTM prices (no IV inversion
    needed for this number -- only for the smile plot)."""
    F = 1.0
    Qs = []
    for m in moneyness_grid:
        is_call = m >= F
        Qs.append(mc_undisc_price(Y, m, is_call))
    return cboe_variance_from_grid(F, T_years, moneyness_grid, Qs)


def run(params, near_calendar_days, next_calendar_days, K=4000, burn=500,
        seed0=50_000_000, moneyness_grid=None):
    if moneyness_grid is None:
        moneyness_grid = build_moneyness_grid(lo=0.60, hi=1.60, atm_lo=0.85, atm_hi=1.20,
                                               atm_step=0.005, wing_step=0.02)

    mat = build_matrices(params)
    alpha = params['alpha']
    n1 = calendar_to_trading_days(near_calendar_days)
    n2 = calendar_to_trading_days(next_calendar_days)
    T1_years = near_calendar_days / CALENDAR_DAYS_PER_YEAR
    T2_years = next_calendar_days / CALENDAR_DAYS_PER_YEAR

    draws = simulate_terminal_ratios(mat, alpha, [n1, n2], K, burn, seed0)

    smile1 = model_smile(draws[n1], T1_years, moneyness_grid)
    smile2 = model_smile(draws[n2], T2_years, moneyness_grid)

    sig1 = model_cboe_sigma2(draws[n1], T1_years, moneyness_grid)
    sig2 = model_cboe_sigma2(draws[n2], T2_years, moneyness_grid)

    N_T1, N_T2, N_30, N_365 = near_calendar_days, next_calendar_days, 30, 365
    vix2 = (T1_years * sig1['sigma2'] * (N_T2 - N_30) / (N_T2 - N_T1)
            + T2_years * sig2['sigma2'] * (N_30 - N_T1) / (N_T2 - N_T1)) * (N_365 / N_30)
    vix_model = 100 * np.sqrt(vix2)

    return dict(
        params=params, kappa0=float(mat['kappa0']), stability_ok=bool(mat['kappa0'] < 0),
        K=K, burn=burn, seed0=seed0,
        near=dict(calendar_days=near_calendar_days, trading_days=n1, T_years=T1_years,
                  sigma2=sig1['sigma2'], implied_vol=100 * np.sqrt(sig1['sigma2'])),
        next=dict(calendar_days=next_calendar_days, trading_days=n2, T_years=T2_years,
                  sigma2=sig2['sigma2'], implied_vol=100 * np.sqrt(sig2['sigma2'])),
        vix_30d=float(vix_model),
        moneyness_grid=moneyness_grid.tolist(),
        smile_near=smile1.tolist(), smile_next=smile2.tolist(),
    )


if __name__ == "__main__":
    import json
    SP500_PARAMS = dict(H=0.386, rough_scale=0.672, lam_lo=0.00001,
                         lam_hi=4.8809, m1=0.1745, alpha=0.220)
    out = run(SP500_PARAMS, near_calendar_days=29, next_calendar_days=57)
    print("model VIX 30d:", round(out['vix_30d'], 3))
    print("near-term implied vol:", round(out['near']['implied_vol'], 3))
    print("next-term implied vol:", round(out['next']['implied_vol'], 3))
    with open("../results/model_pricing_results.json", "w") as f:
        json.dump(out, f, indent=2)
