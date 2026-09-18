"""
vix_smile.py
======================================
Model-implied VIX-option smile via NESTED Monte Carlo.

Real VIX options are options on the 30-day model-free variance swap
rate itself (Section 3 of the companion note), observed at a future
date. Reproducing their smile from the ESN therefore needs a
genuinely different, 2-level simulation from the SPX smile
(Section 3): an OUTER ensemble of possible reservoir states at the
option's own expiry date N, and, FROM EACH outer draw, an INNER
ensemble of 30-day-forward paths used to compute THAT day's own
model-implied VIX (by applying the identical CBOE replication used
throughout this note to the inner ensemble's own terminal
distribution). The outer ensemble of resulting VIX_N values is then
itself priced (calls/puts, Black-76 inversion) exactly as the SPX
smile was.

This requires conditioning a simulation on an arbitrary MID-PATH
state (the outer draw's own reservoir state at day N) -- something
esn_base._sim_esn_with_params cannot do (it always starts at
r=z=var_used=0) -- hence vectorized_sim.py.
"""
import numpy as np
import esn_base as E
from joint_calibration import OPT_ARCH, P_FIXED
from market_vix import black76_undisc, cboe_variance_from_grid, build_moneyness_grid
from vectorized_sim import simulate_batch
from model_pricing import calendar_to_trading_days, build_matrices, implied_vol_from_price

INNER_WINDOW_CALENDAR_DAYS = 30
INNER_WINDOW_TRADING_DAYS = calendar_to_trading_days(INNER_WINDOW_CALENDAR_DAYS)  # ~21


def outer_states(mat, alpha, n_outer_days, K_outer, burn, seed):
    """Burn in K_outer independent chains, continue n_outer_days
    further, return the terminal state of each -- the ensemble of
    'possible future reservoir states on the option's own expiry
    date'."""
    dx, v, state = simulate_batch(OPT_ARCH, mat, alpha, burn + n_outer_days, K_outer, seed=seed)
    return state   # (r [K_outer,Nr], z [K_outer,Nz], var_used_prev [K_outer])


def conditional_vix(mat, alpha, outer_state, K_inner, moneyness_grid, seed):
    """From EACH outer state, run K_inner inner paths INNER_WINDOW
    trading days forward; return one model-implied VIX value per
    outer draw (K_outer,)."""
    r0, z0, v0 = outer_state
    K_outer = r0.shape[0]
    r0_rep = np.repeat(r0, K_inner, axis=0)
    z0_rep = np.repeat(z0, K_inner, axis=0)
    v0_rep = np.repeat(v0, K_inner, axis=0)

    dx_inner, _, _ = simulate_batch(OPT_ARCH, mat, alpha, INNER_WINDOW_TRADING_DAYS,
                                     K_outer * K_inner, r0=r0_rep, z0=z0_rep,
                                     var_used_prev0=v0_rep, seed=seed)
    Y = np.exp(np.cumsum(dx_inner, axis=0)[-1])           # (K_outer*K_inner,)
    Y = Y.reshape(K_outer, K_inner)

    T_years = INNER_WINDOW_CALENDAR_DAYS / 365.0
    vix_vals = np.empty(K_outer)
    for k in range(K_outer):
        Yk = Y[k]
        Qs = [float(np.mean(np.maximum(Yk - m, 0.0))) if m >= 1.0
              else float(np.mean(np.maximum(m - Yk, 0.0))) for m in moneyness_grid]
        res = cboe_variance_from_grid(1.0, T_years, moneyness_grid, Qs)
        vix_vals[k] = 100 * np.sqrt(max(res['sigma2'], 0.0))
    return vix_vals


def price_vix_options(vix_vals, F, T_years, strike_grid):
    """MC price + Black-76 invert a VIX-value ensemble against a real
    strike grid, using the ensemble's own forward F = mean(vix_vals)."""
    vols = []
    for K in strike_grid:
        is_call = K >= F
        price = float(np.mean(np.maximum(vix_vals - K, 0.0))) if is_call else \
                float(np.mean(np.maximum(K - vix_vals, 0.0)))
        vols.append(implied_vol_from_price(F, K, T_years, price, is_call, hi=8.0))
    return np.array(vols)


def run(params, near_calendar_days, next_calendar_days, strike_grid_near, strike_grid_next,
        K_outer=800, K_inner=800, burn=500, seed0=70_000_000):
    mat = build_matrices(params)
    alpha = params['alpha']
    moneyness_grid = build_moneyness_grid(lo=0.60, hi=1.60, atm_lo=0.85, atm_hi=1.20,
                                           atm_step=0.01, wing_step=0.04)  # coarser than SPX (nested MC is costlier)

    out = {}
    for label, cal_days, strikes in [("near", near_calendar_days, strike_grid_near),
                                      ("next", next_calendar_days, strike_grid_next)]:
        n_outer = calendar_to_trading_days(cal_days)
        ostate = outer_states(mat, alpha, n_outer, K_outer, burn, seed0 + (0 if label == "near" else 1))
        vix_vals = conditional_vix(mat, alpha, ostate, K_inner, moneyness_grid,
                                    seed0 + 1000 + (0 if label == "near" else 1))
        F = float(np.mean(vix_vals))
        T_years = cal_days / 365.0
        vols = price_vix_options(vix_vals, F, T_years, strikes)
        out[label] = dict(calendar_days=cal_days, trading_days=n_outer, forward=F,
                           vix_ensemble_mean=F, vix_ensemble_std=float(np.std(vix_vals)),
                           vix_ensemble_p5=float(np.percentile(vix_vals, 5)),
                           vix_ensemble_p95=float(np.percentile(vix_vals, 95)),
                           strikes=list(strikes), implied_vols=vols.tolist())
    out['K_outer'] = K_outer
    out['K_inner'] = K_inner
    out['inner_window_trading_days'] = INNER_WINDOW_TRADING_DAYS
    return out


if __name__ == "__main__":
    import json
    import pandas as pd
    SP500_PARAMS = dict(H=0.386, rough_scale=0.672, lam_lo=0.00001,
                         lam_hi=4.8809, m1=0.1745, alpha=0.220)

    vix_df = pd.read_csv("/mnt/user-data/uploads/1788726546415_VIX_data_19_02_2026.csv")
    near_strikes = sorted(vix_df[vix_df.expiry == 46099].strike.unique())
    next_strikes = sorted(vix_df[vix_df.expiry == 46127].strike.unique())

    out = run(SP500_PARAMS, 27, 55, near_strikes, next_strikes, K_outer=800, K_inner=800)
    print("near forward VIX:", round(out['near']['forward'], 2),
          " [p5,p95]=", round(out['near']['vix_ensemble_p5'],1), round(out['near']['vix_ensemble_p95'],1))
    print("next forward VIX:", round(out['next']['forward'], 2),
          " [p5,p95]=", round(out['next']['vix_ensemble_p5'],1), round(out['next']['vix_ensemble_p95'],1))
    with open("../results/vix_smile_results.json", "w") as f:
        json.dump(out, f, indent=2)
