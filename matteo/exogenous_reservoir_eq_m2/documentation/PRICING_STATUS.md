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
There are now 60 passing package tests and six passing legacy checks. The new
VIX tests include a constant-volatility oracle, state-dependent conditional
values, date-zero normalization, nested/quadrature correspondence, the outer
second moment and stale-normalizer rejection. These tests do not supply
quote-relative convergence guarantees. Initial nominal SE screens with zero
sampled payoffs are corrected by additive `tail_diagnostics.json` records;
those original counters must not be read as successful tail-error gates.

## Remaining work

1. Reduce SPX Monte Carlo variance and assess stock-integration bias against
   quote scales, using coupled grids and independent final seeds. Complete
   normalizer quadrature refinement in cap/tail regimes.
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
