# Dated pricing implementation status — 21 October 2026

This is an implementation update to [CALIBRATION_PLAN.md](CALIBRATION_PLAN.md),
not a change to the model or its historical research claims.

## Implemented

- `src/esn_eq/market_targets.py` loads the existing processed price intervals,
  frozen scales, split labels, VIX forward targets, curve candidates and audit
  flags. Training selection is explicit. It validates contract identities,
  forwards, moneyness and units. Failed/nonfinite model prices raise an error
  instead of vanishing from the residual.
- `src/esn_eq/pricing.py` supports Gaussian M2 (zero echo and orthogonal
  loadings), fixed full dated state and the minimum-prior-norm full-state map
  for supplied projections. Its API orders those projections as **(C,Y,F)**.
  No filtered state or inferred variance-risk premium is supplied by this API.
- Exact joint OU/equity transitions reuse the existing covariance functions.
  Q drift follows `states.q_state_paths()`'s negative risk-price shift. Rates
  and time are annualized; quoted maturities are not rounded to trading days.
- A reusable state/projection/equity cache separates expensive simulation from
  readout changes. The optional checkpoint states support later conditional
  valuation. Incompatible architecture/risk-price caches are rejected.
- The date-zero normalizer uses the three-dimensional conditional Gaussian
  law of (C,Y,F) and tensor Gauss-Hermite quadrature. It is deterministic given
  X0 and parameters; no stock sample is used to normalize itself. Quadrature
  remains approximate and needs independent refinement.
- SPX terminal ratios use adapted exponential Euler. Absolute undiscounted
  prices are the quoted SPX forward times the normalized payoff expectation.
  This uses deterministic quoted carry; cash premiums/settlement precision
  require further inputs.
- Curve jumps and expiries are exact grid events. A coarsening helper sums
  the same Brownian increments and selects the same exact OU states, enabling
  paired stock-integration refinement rather than independent noisy runs.
- Price MC standard errors, simulated variance means, normalizer changes,
  cap use and paired timestep changes are retained explicitly.
- `src/esn_eq/vix_pricing.py` now prices the conditional continuous 30/365-year
  VIX functional by Gaussian quadrature or nested full-state OU simulation.
  A read-only window stores the date-zero normalizer and rejects stale model
  parameters. Midpoint integration splits at curve jumps. Independent inner
  streams per outer state preserve prefixes across inner-count refinements.
  Absolute VIX strikes are priced; quoted forwards remain separate targets.
- `src/esn_eq/variance_reduction.py` now supplies SPX controls using the stock
  martingale and a deterministic-volatility reference option with known Black
  expectation. Regression coefficients are trained on an independent pilot,
  then frozen for pricing and paired time-step diagnostics. Negative adjusted
  estimates are retained, never clipped.
- `src/esn_eq/calibration_loss.py` implements the frozen common family/expiry
  interval objective and separate VIX forward residuals, with explicit failure
  handling and held-out boundaries. Synthetic parameter recovery is tested.
- Conditional VIX supports Gaussian time rules on geometrically growing lag
  panels split at curve knots. Order 8/16 checks against the fine midpoint
  reference and sampled states support fewer time nodes at the frozen anchor.
  Stock-path resolution and outer sampling requirements are unchanged.
- `src/esn_eq/outer_sampling.py` supplies full-state Gaussian iid/Sobol
  sampling with projection-ordered PCA and uncertainty across independent
  complete replications. A matched-cost 55-day pilot supports a substantial
  sampling improvement, but tail and quadrature accuracy remain unresolved.
- `src/esn_eq/market_calibration.py` constructs the declared three-coefficient
  nuisance curve, penalized `(C,Y,F)` initial state and five primitive risk-price
  scenarios. These are calibration inputs, not fitted results.
- `src/esn_eq/vix_curve_profile.py` caches conditional curve-segment exposures.
  At fixed model/state, integrated variance is affine in segment levels;
  candidate VIX prices retain the square root and nonzero variance floor.
  Direct-repricing, affine, constant-variance and stale-cache checks pass.
- `src/esn_eq/hedging_controls.py` adds two predictable Black-inspired trading
  gains as zero-mean SPX controls. Independent pilot coefficients combine them
  with the existing static controls. Matched-path tests and a three-grid,
  two-trial diagnostic show variance reductions, with tails still unresolved.

## Evidence and limits

Analytic-oracle tests cover Black prices with deterministic integrated
variance, shared-path put-call parity, nonconstant-volatility stock-forward
checks, exact state moments, initial-state projection recovery, the Q drift,
and independent Monte Carlo checks of Gaussian normalization. A coupled
constant-volatility grid test preserves stock paths up to floating-point
error. These are numerical/implementation checks, not formal proofs.

The market runner lives at
[`../../model_comparison/option_pricing/m2_spx.py`](../../model_comparison/option_pricing/m2_spx.py).
It uses the frozen Stage-3 candidate 18/scenario 0, disables the unused echo
calculation, sets orthogonal/echo readout loadings to zero, and records the
complete resulting configuration. It reads training quotes only at 29/57 days.
It does not change the historical anchor or optimize any parameter.

The first experiments use zero projected X0 and the diagnostic training-mid
curve. The preprocessing audit flags that curve's incompatibility with later
VIX quotes; this flag is copied into every run. It may be used to check the
engine, but not as a certified feasible joint-market input. Plain MC prices
at 4096 paths fail the planned quote-relative accuracy gate. Neither weak
initial interval coverage nor apparent smile differences establish a model
failure before fitting and numerical resolution.

See [the market diagnostic record](../../model_comparison/option_pricing/M2_SPX_RESULTS.md)
for run configurations and outcomes. Original physical experiments and the
legacy ESN sources remain preserved.

The next conditional-VIX and SPX-control diagnostic record is
[M2_VIX_RESULTS.md](../../model_comparison/option_pricing/M2_VIX_RESULTS.md).
That checkpoint recorded 60 passing package tests and six legacy checks. Its new
VIX tests include a constant-volatility oracle, state-dependent conditional
values, date-zero normalization, nested/quadrature correspondence, the outer
second moment and stale-normalizer rejection. These tests do not supply
quote-relative convergence guarantees. Initial nominal SE screens with zero
sampled payoffs are corrected by additive `tail_diagnostics.json` records;
those original counters must not be read as successful tail-error gates.

The requested coupled 64/128 SPX and fixed-state VIX checks have now run; see
[FINE_GRID_RESULTS.md](../../model_comparison/option_pricing/FINE_GRID_RESULTS.md).
SPX has not met quote-relative convergence tolerances. VIX quadrature window
integration is stable across these grids at the three tested states, while
nested inner sampling remains noisy. Standalone runner defaults are now 64
steps/day, and a batched driver supports the larger SPX path caches.
The 128/256 SPX follow-up has also completed with the same fixed configuration
and path budget. Its largest controlled paired changes are 0.576/0.905 SPX
points (paired SEs 0.308/0.512); quote-relative numerical gates remain open.
The following work targeted sampling efficiency and the common calibration objective,
with the nuisance-curve/state profiles retained from the existing plan.

The common objective and faster time-integration checkpoint is recorded in
[CALIBRATION_PROGRESS.md](../../model_comparison/option_pricing/CALIBRATION_PROGRESS.md).
That checkpoint recorded 68 passing package tests and eight legacy/driver checks.
A restricted joint anchor score combines archived 256/day SPX prices at 29/57
days with fresh VIX samples at all four maturities. It remains unfitted,
numerically unresolved and subject to the original input-curve conflict.

The subsequent sampling and constrained curve-initializer checkpoint is in
[SAMPLING_PROGRESS.md](../../model_comparison/option_pricing/SAMPLING_PROGRESS.md).
That checkpoint recorded 77 passing package tests and eight legacy/driver checks.
The initializer satisfies necessary moment inequalities at quoted VIX means
but substantially changes short-end stripped total variance. It is a starting
point, not a joint option fit. The next checkpoint below implements curve
exposures and the preliminary VIX-only market optimization.

The exposure cache and first exploratory two-coefficient VIX curve fit are
now implemented; see [CURVE_PROFILE_RESULTS.md](../../model_comparison/option_pricing/CURVE_PROFILE_RESULTS.md).
All eight starts solve the fixed-sample problem, but the apparent improvement
reverses on fresh checking scrambles. No reliable fit improvement or joint
feasibility is established. Numerical comparison stability remains a priority.
That checkpoint recorded 81 passing package tests and eight legacy/driver checks.
The resulting curve was also checked against 29/57-day SPX training prices
with coupled 128/256 grids and independent-pilot controls. Its descriptive
SPX loss is 384.11 versus 963.93 for the earlier stripped curve; much of the
change comes from the nuisance family's higher near variance. Numerical gates
remain unresolved, and the VIX sampling reversal prevents a reliable joint
improvement claim.

The subsequent diagnosis is in
[NUMERICAL_STABILITY_RESULTS.md](../../model_comparison/option_pricing/NUMERICAL_STABILITY_RESULTS.md).
At the reversal's 153-day VIX block, increasing conditional GH 32/48 and
normalization 48/64 changes the fixed-candidate comparison by at most about
0.0031 loss units, whereas changing the outer sample moves it about 10.8.
Rare-event outer sampling is the main supported explanation at that block.
Predictable SPX hedging controls reduce estimated variance by a median factor
5.20 over twelve representative contracts in two trials, but none meets the
sampling tolerance. The new 64/128/256 comparisons share one finest path set
per trial; time bias and normalizer accuracy remain unresolved.
There are now 84 passing package tests and nine passing legacy/driver checks.

## Remaining work

1. Carry forward the pre-existing time-step convergence decision: 64 steps per
   trading day is the working baseline for the frozen fast-rate architecture;
   compare coupled 64/128 grids and extend to 256 as needed for option-price
   tolerances. The initial 2/8-step pricing runs are coarse diagnostics only.
   Reduce SPX Monte Carlo variance and assess stock-integration bias against
   quote scales using independent final seeds. Complete normalizer quadrature
   refinement in cap/tail regimes. Keep model/state fixed across resolutions.
   Extend the tested hedging controls with tail-aware pilot/evaluation sampling
   and measure multilevel correction variance/cost before scaling path counts.
2. Refine the implemented conditional VIX engines, especially outer/tail
   sampling. Conditional quadrature and nested integration agree in controlled
   tests, but market quote-relative sampling/bias gates are still unresolved.
3. Profile the documented nuisance-curve family to address input conflicts,
   keeping the original stripped curves as diagnostics. Do not label this as
   repairing raw quotes or infer a unique curve from this one snapshot.
4. Build the common-protocol ESN evaluator, state/risk-price scenario layer
   and joint fitter, then physical compatibility and held-out evaluation as
   specified in the calibration plan. No held-out market performance has been
   examined in the implementation diagnostics.
