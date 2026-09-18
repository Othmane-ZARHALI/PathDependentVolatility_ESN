"""
make_figures.py
======================================
Regenerates every figure used in the final report:
  1. calibration_convergence.png       -- L(theta) vs. evaluation count, best
                                           (final) configuration only (Section 3.3)
  2. spx_smile_final.png               -- market vs. theta_hat, 4 maturities, curve
                                           clipped to each panel's own market range
  3. vix_smile_final.png               -- market vs. theta_hat, 4 maturities, curve
                                           clipped to each panel's own market range
  4. vix_distribution_comparison.png   -- VIX empirical distribution, market
                                           (Breeden-Litzenberger) vs. model ensemble (Test 1)

Run model_pricing.py, vix_smile.py, calibrate_smile.py, and
vix_distribution.py first to (re)produce the JSON results this script
reads.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 10, "figure.dpi": 150})
FIGDIR = "../figures/"
SPX_CSV = "/mnt/user-data/uploads/1788726546415_SPX_data_19_02_2026.csv"
VIX_CSV = "/mnt/user-data/uploads/1788726546415_VIX_data_19_02_2026.csv"

spx = pd.read_csv(SPX_CSV)
vix = pd.read_csv(VIX_CSV)

with open("../results/model_pricing_results.json") as f:
    model = json.load(f)
with open("../results/market_results.json") as f:
    market = json.load(f)

# ---------- Fig 1: L(theta) vs. evaluation count, multiple random inits ----------
# Plotted as L(theta)/N_points -- the mean weighted SQUARED residual per point (N_points=858,
# the total across all 6 specs) -- rather than the raw sum, so the true optimum's own floor
# reads as a genuinely small number close to 0, not an arbitrary-looking ~4.6. This rescaling
# changes nothing about the shape of the curve or the relative comparison between points, only
# the y-axis's own units.
#
# This single evaluation sequence actually concatenates 10 separate optimizer runs, back to
# back, with NO visual marker on the plot itself distinguishing where one ends and the next
# begins (by request) -- documented here instead:
#   evals   1-64   : a deliberately extreme, near-worst-case corner of bounds Theta (H, lam_hi,
#                     m1, m2 all at their own upper bound; alpha and Delta b0 at their lower
#                     bound) -- not a random draw, the single highest-loss point found in this
#                     note (L/N approx 1.42, an order of magnitude above every random draw
#                     tried elsewhere here) -- converges on its own, quickly, down to ~0.0089
#   evals  65-160  : random theta (uniform draw within bounds Theta, seed 42)
#                     -> gets stuck; plateaus ~26-33, alpha and m2 both pin near their own bounds
#   evals 161-190  : restarted from the note's own informed theta0 (Table 3)
#                     -> converges immediately to ~4.6-4.8, the true optimum
#   evals 191-244  : random theta (seed 123)
#                     -> finds its own way to the true optimum (~5.0) unaided, just less directly
#   evals 245-335  : random theta (seed 7)
#                     -> stuck again, at a 3rd, different plateau (~7.5-8.9)
#   evals 336-353  : restarted again from theta0 -> converges to ~4.6 again
#   evals 354-459  : random theta (seed 2024)
#                     -> stuck ~25-31 for a long stretch (>60 evaluations), then breaks free and
#                        also finds its own way to the true optimum (~5.4) unaided
#   evals 460-477  : restarted a 3rd time from theta0 -> converges to ~4.6 again
#   evals 478-533  : random theta (seed 999)
#                     -> stuck ~31-52 throughout, never finds the true optimum in this budget
#   evals 534-590  : restarted a 4th time from theta0 -> converges to ~4.6 again, and further
#                     evaluations spent there (57 total) do not push it any lower -- this is a
#                     genuine, repeatable floor for this objective, not a slow approach to 0
#
# The honest takeaway, precisely because nothing is marked: random initialization is
# unreliable across seeds (2 of 5 tried here eventually reached the true optimum unaided, 2
# got badly stuck throughout, 1 got moderately stuck throughout), while every one of the 4
# informed-start attempts converges to the same point immediately, in well under 30
# evaluations every time, and no amount of further search there improves on it. That is the
# actual argument for warm-starting from theta0 throughout this note -- and the reason this
# floor, not 0, is what a correctly-posed version of this objective actually looks like.
N_POINTS = 858
hist = [json.loads(l) for l in open('../results/calibration_history.jsonl')]
iters = np.array([h['iter'] for h in hist])
sumsq = np.array([h['sumsq'] for h in hist]) / N_POINTS
running_min = np.minimum.accumulate(sumsq)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axes[0]
ax.plot(iters, sumsq, '.', ms=4, color='#888888', alpha=0.4, label=r'$\mathcal{L}(\theta)/N$ per evaluation')
ax.step(iters, running_min, where='post', color='#c44e52', lw=2, label='running minimum')
ax.set_xlabel('residual evaluation')
ax.set_ylabel(r'$\mathcal{L}(\theta)/N$ (mean weighted squared residual, $N=858$)')
ax.set_title('Full range')
ax.set_ylim(0, None)
ax.legend(fontsize=8, loc='upper right')

ax = axes[1]
ax.plot(iters, sumsq, '.', ms=4, color='#888888', alpha=0.4, label=r'$\mathcal{L}(\theta)/N$ per evaluation')
ax.step(iters, running_min, where='post', color='#c44e52', lw=2, label='running minimum')
ax.set_xlabel('residual evaluation')
ax.set_ylabel(r'$\mathcal{L}(\theta)/N$')
ax.set_title('Same data, zoomed to the 0-0.06 range')
# A handful of trf's own finite-difference Jacobian probes reach L/N > 0.06, off the top of
# this panel's own range, and are not shown -- see calibration_history.jsonl for every point.
ax.set_ylim(0, 0.06)
ax.legend(fontsize=8, loc='upper right')

fig.suptitle('Joint calibration objective vs. evaluation count')
ax.legend(fontsize=9, loc='upper right')
fig.tight_layout()
fig.savefig(FIGDIR + "calibration_convergence.png")
plt.close(fig)

# ---------- Fig 2: SPX final, market vs. mean model IV +/- 90% CI, 4 maturities ----------
# Plots the MEAN curve across N_REPEATS independent MC seeds, with a shaded 5th-95th percentile
# band, rather than 1 single noisy realization -- single realizations were shown (Section on
# 153d) to be seed-dependent enough to materially mislead, so this is the honest default now.
# Market quotes carry BOTH a put and a call row per strike (columns callput=-1/+1) -- both are
# shown, distinguished by marker, since together they are what the calibration itself fits
# against (Section 3.3.1's own ITM/OTM split), and their own small, real disagreement (put IV
# vs call IV at the identical strike -- not a bid/ask spread, both are already bid/ask mid
# prices) is itself informative, not a plotting duplicate.
d = json.load(open('../results/spx_smile_final_curves.json'))
grid = np.array(d['grid'])
n_rep_spx = d.get('n_repeats', 1)
panels = [('near', 46101, 29, 'calibrated'), ('next', 46129, 57, 'calibrated'),
          ('mid1', 46191, 119, 'out-of-sample'), ('mid2', 46283, 211, 'out-of-sample')]
fig, axes = plt.subplots(2, 2, figsize=(11, 8.4), sharey=True)
for ax, (label, expiry, days, tag) in zip(axes.flat, panels):
    dd = spx[spx.expiry == expiry].copy()
    dd['mid'] = (dd.bid + dd.ask) / 2
    dd['m'] = dd.strike / dd.forward
    dd = dd[(dd.m >= 0.6) & (dd.m <= 1.6)].dropna(subset=['mid'])
    dd_put = dd[dd.callput == -1]
    dd_call = dd[dd.callput == 1]
    c = d['curves']['theta_hat'][label]
    mean, lo, hi = np.array(c['mean']), np.array(c['ci_lo']), np.array(c['ci_hi'])
    m_lo, m_hi = dd.m.min(), dd.m.max()
    clip = (grid >= m_lo) & (grid <= m_hi)
    pad = 0.02 * (m_hi - m_lo)
    ax.plot(dd_put.m, dd_put.mid, 'o', ms=3, mfc='none', mec='#888888', alpha=0.7, label='market put (mid IV)')
    ax.plot(dd_call.m, dd_call.mid, 'x', ms=4, color='#555555', alpha=0.7, label='market call (mid IV)')
    ax.plot(grid[clip], mean[clip], '-', lw=2.2, color='#1f4fd6', label=f'model IV (mean of {n_rep_spx})')
    ax.fill_between(grid[clip], lo[clip], hi[clip], color='#1f4fd6', alpha=0.2, label='90% CI (5th-95th pctile)')
    if m_lo <= 1.0 <= m_hi:
        ax.axvline(1.0, color='gray', ls='--', lw=0.6)
    ax.set_xlim(m_lo - pad, m_hi + pad)
    ax.set_ylim(0.0, 0.65)
    ax.set_xlabel("moneyness $K/F$")
    ax.set_title(f"{days}d ({tag})")
    ax.legend(fontsize=7.5)
for ax in axes[:, 0]:
    ax.set_ylabel("implied volatility")
fig.suptitle(f"SPX smile: market put (circles) and call (crosses) mid implied vol vs.\n"
              f"mean model IV +/- 90% CI (blue, {n_rep_spx} independent MC seeds), both axis "
              "range and curve limited to each maturity's own market-quoted moneyness")
fig.tight_layout()
fig.savefig(FIGDIR + "spx_smile_final.png")
plt.close(fig)

# ---------- Fig 3: VIX final, market vs. mean model IV +/- 90% CI, all 4 available maturities ----------
# y-axis is dynamic per panel (min/max of BOTH market and the model's own CI band, with padding)
# -- no fixed cutoff. Put/call quotes distinguished as in Figure 2.
dh = json.load(open('../results/vix_smile_final_curves_thetahat.json'))
grid_vix = np.array(dh['grid'])
n_rep_vix = dh.get('n_repeats', 1)
panels = [('near', 46099, 27, 'calibrated'), ('next', 46127, 55, 'calibrated'),
          ('mid1', 46190, 118, 'calibrated'), ('mid2', 46225, 153, 'calibrated')]

fig, axes = plt.subplots(2, 2, figsize=(11, 8.4), sharey=False)
for ax, (label, expiry, days, tag) in zip(axes.flat, panels):
    mkt = vix[vix.expiry == expiry].copy()
    mkt['mid'] = (mkt.bid + mkt.ask) / 2
    F_mkt = mkt.forward.iloc[0]
    mkt['m'] = mkt.strike / F_mkt
    mkt = mkt.dropna(subset=['mid'])
    mkt_put = mkt[mkt.callput == -1]
    mkt_call = mkt[mkt.callput == 1]
    c = dh['curves'][label]
    mean, lo, hi = np.array(c['mean'], dtype=float), np.array(c['ci_lo'], dtype=float), np.array(c['ci_hi'], dtype=float)
    m_lo, m_hi = mkt.m.min(), mkt.m.max()
    clip = (grid_vix >= m_lo) & (grid_vix <= m_hi)
    pad = 0.02 * (m_hi - m_lo)
    mean_c, lo_c, hi_c = mean[clip], lo[clip], hi[clip]
    valid = ~np.isnan(mean_c)
    y_lo = min(mkt.mid.min(), np.nanmin(lo_c)) if valid.any() else mkt.mid.min()
    y_hi = max(mkt.mid.max(), np.nanmax(hi_c)) if valid.any() else mkt.mid.max()
    y_pad = 0.05 * (y_hi - y_lo)
    ax.plot(mkt_put.m, mkt_put.mid, 'o', ms=4, mfc='none', mec='#888888', alpha=0.7, label='market put (mid IV)')
    ax.plot(mkt_call.m, mkt_call.mid, 'x', ms=5, color='#555555', alpha=0.7, label='market call (mid IV)')
    ax.plot(grid_vix[clip], mean_c, '-', lw=2.2, color='#1f4fd6', label=f'model IV (mean of {n_rep_vix})')
    ax.fill_between(grid_vix[clip], lo_c, hi_c, color='#1f4fd6', alpha=0.2, label='90% CI (5th-95th pctile)')
    if m_lo <= 1.0 <= m_hi:
        ax.axvline(1.0, color='gray', ls='--', lw=0.6)
    ax.set_xlim(m_lo - pad, m_hi + pad)
    ax.set_ylim(y_lo - y_pad, y_hi + y_pad)
    ax.set_xlabel("moneyness $K/F$")
    ax.set_title(f"{days}d ({tag})")
    ax.legend(fontsize=7.5)
for ax in axes[:, 0]:
    ax.set_ylabel("implied volatility (vol-of-VIX)")
fig.suptitle(f"VIX-option smile: market put (circles) and call (crosses) mid implied vol vs.\n"
              f"mean model IV +/- 90% CI (blue, {n_rep_vix} independent MC seeds) -- axis range "
              "shows the complete curves")
fig.tight_layout()
fig.savefig(FIGDIR + "vix_smile_final.png")
plt.close(fig)

# ---------- Fig 4: VIX empirical distribution, market (Breeden-Litzenberger) vs. model ensemble, 4 maturities ----------
dist = json.load(open('../results/vix_distribution.json'))
panels = [('near', 27, 'calibrated'), ('next', 55, 'calibrated'),
          ('mid1', 118, 'calibrated'), ('mid2', 153, 'calibrated')]
fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))
for ax, (label, days, tag) in zip(axes.flat, panels):
    mkt = dist[label]['market']
    K = np.array(mkt['K'])
    dens = np.array(mkt['density'])
    ens = np.array(dist[label]['model_ensemble'])

    ax.plot(K, dens, '-', lw=2, color='#888888',
            label='market (risk-neutral density,\nBreeden-Litzenberger)')
    ax.hist(ens, bins=60, range=(0, 100), density=True, color='#1f4fd6', alpha=0.55,
            label=r'model ($\hat\theta$, ensemble of' + f'\n{len(ens)} simulated draws)')
    ax.axvline(mkt['F'], color='#888888', ls='--', lw=1, label=f"market forward ({mkt['F']:.1f})")
    ax.axvline(dist[label]['model_mean'], color='#1f4fd6', ls='--', lw=1,
               label=f"model mean ({dist[label]['model_mean']:.1f})")
    ax.set_xlim(0, 100)
    ax.set_xlabel("VIX level")
    ax.set_title(f"{days}d ({tag})")
    ax.legend(fontsize=7.5)
for ax in axes[:, 0]:
    ax.set_ylabel("probability density")
fig.suptitle("VIX empirical distribution: market (grey curve, risk-neutral density) vs.\n"
              "model (blue histogram, $\\hat\\theta$ simulated ensemble), all 4 quoted maturities")
fig.tight_layout()
fig.savefig(FIGDIR + "vix_distribution_comparison.png")
plt.close(fig)

print("all figures written to", FIGDIR)
