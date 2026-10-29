# Warm start and continuous controls — 24 October 2026

The user asks to start from the previous successful parameter vector and remove
strike jitter. This is a pricing-method follow-up, not a new calibration.

## Optimizer and starting point

The implementation uses scipy.optimize.least_squares with its default
trust-region reflective (TRF) method, box bounds and analytic derivatives.
This is a local least-squares optimizer, not a global minimization guarantee.
The eleven-parameter run started from the preceding all-maturity vector except
beta_J was set to 0.15 (second start also moved other parameters). A cleaner
new first start preserves all ten other entries exactly and changes only
beta_J from zero to its declared positive lower bound 0.01. The helper
`continuous_strike_controls.warm_start` constructs that vector and rejects
unrelated bound violations rather than silently moving other coordinates.
The exact proposed vector is retained in the diagnostic run's warm_start.json.
Future search must retain that feasible starting vector as an incumbent, use
the SAME price oracle for its baseline and every subsequent candidate, and
report searching and separate numerical checks distinctly.

## Continuity construction

Let Z_m be an existing frozen, known-zero control correction at anchor strike
K_m, including its independently fitted frozen coefficient and importance
weight. Define deterministic piecewise-linear weights a_m(K) on anchor strikes,
with constant endpoint extrapolation, and form the correction PATH BY PATH:

    Z(K) = sum_m a_m(K) Z_m
    price_hat(K) = mean(weight * payoff(K) - Z(K)).

The weights depend only on strikes; E[Z(K)] = sum_m a_m(K) E[Z_m] = 0 under
the original target-law/control assumptions. Z(K) and the empirical price are
continuous in K. This interpolates zero-mean random controls before averaging;
it does not interpolate or smooth model or market prices after pricing.
Anchors are the nearest quotes to scaled moneyness 0.5/1.5/2.5 on each option
side, selected from market-only information before inspecting the results.
All seven agreed SPX dates, 1294 cutoff-3 OTM quotes and both existing 512-path
search/check ensembles enter the saved-draw experiment. Parameters are the
preceding all-maturity selected zero-shift reference throughout this diagnostic.

## Results and limitations

| Checking ensemble method | Strike-order violations | Convexity violations | Negative prices |
|---|---:|---:|---:|
| Historical nearest-donor controls | 11 | 381 | 0 |
| Continuous piecewise-linear controls | 0 | 15 | 0 |
| Global affine control comparison | 4 | 0 | 5 |
| Raw positive weighted payoffs | 0 | 0 | 0 |

Searching ensemble: nearest-donor controls have 10 strike-order violations;
continuous controls have zero, with 13 remaining convexity violations and no
negative prices. The affine comparison has zero convexity violations by
construction but 13 strike-order violations and one negative price.
Subtracting a global affine function preserves convexity, but does not ensure
positivity or monotonicity. These failures are retained, not clipped.

The continuous construction keeps useful measured precision. At 119/211/666
days, median price SE is respectively 4.80/5.62/10.71, versus 4.81/5.72/12.39
with historical controls and 19.10/17.31/29.06 for raw payoffs on THESE SAME
samples. This is not a comparison with the newer low-budget raw run.

Thus the construction removes discontinuous donor switching and all observed
strike-order failures in this bounded experiment, but piecewise-linear anchor
changes can still introduce convexity failures. It is a promising control
repair, not a fully qualified arbitrage-consistent pricing surface. The affine
alternative is not adopted either. Do not silently substitute either into a
new calibration or smooth plotted IVs to hide those failures. Smoothness,
shape validity, precision and population accuracy are separate requirements.

Six focused tests cover known-zero expectation, anchor continuity, interpolation
weights, affine convexity and exact preservation of the warm-start coordinates.
No new simulation, normalizer, market fit, model price or theorem was produced.
Only alternative statistical estimates from retained payoff/control draws were
computed. Original runs and pictures are unchanged.

## Retained evidence

Final diagnostic: runs/m2_continuous_strike_controls_003/run.json and
shape_and_precision.csv, control_anchors.csv, warm_start.json and
saved_estimator_draws.npz (all raw/old/continuous/affine pathwise estimates).
Each diagnostic captures its source/input bytes and hashes. Run 001 is the
first piecewise-linear experiment; run 002 adds the affine comparison; run 003
retains all comparison draw arrays. These are not three independent studies.
The initial module-import startup failure is retained in startup_failure.log;
fixing import order did not alter numerical results. No calibration was launched.
