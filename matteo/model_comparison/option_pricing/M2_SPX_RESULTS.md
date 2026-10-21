# Dated M2 SPX diagnostic — 21 October 2026

## What is now available

The common price-target interface and dated Gaussian-M2 SPX prototype are
implemented. The pricer uses an explicit full X0, exact correlated OU/equity
transitions, adapted stock integration, exact expiry/curve events, and an
independent deterministic Gaussian-quadrature normalizer. It prices absolute
SPX strikes in undiscounted forward-value units. All model prices enter the
interval diagnostic; nonfinite prices cause an error rather than being masked.

Validation: **49 M2 package tests and six legacy benchmark tests pass**.
The new tests include analytic Black prices, common-path parity, state moments,
the Q shift, independent normalization checks, projected-state recovery,
split preservation and a coupled-grid constant-volatility oracle.

This is an unfitted implementation experiment, **not an achieved option fit**,
a historical-versus-pricing compatibility result, or a model-family ranking.

## Frozen experimental choices

- Existing Stage-3 M2 candidate 18, structural scenario 0, read from the frozen
  manifest/table. Its slow spike grid is retained; default architecture values
  are not substituted. The unused echo calculation is disabled and both echo
  and orthogonal readout loadings are zero.
- X0 is the minimum-stationary-prior-norm state for **(C0,Y0,F0)=(0,0,0)**.
  It is a declared reference state, not a filtered February 2026 state.
- Risk-price entries are inherited from that frozen row: equity 0.35, both
  unspanned feedback/spike prices zero. The effective Q drifts retain the
  nonzero equity-correlation contributions.
- The training-mid piecewise forward-variance curve is used only for engine
  diagnostics. The existing preprocessing audit flags a conflict with later
  VIX quotes and requires nuisance-curve profiling; every run preserves that
  status. No quote, curve or audit flag is repaired here.
- Training price intervals only: 328 SPX quotes at 29 days and 320 at 57 days.
  No held-out quote performance is evaluated. The legacy benchmark instead
  uses pooled mid-IV targets, so its RMSE columns are not directly comparable.
- 4096 independent stock paths per run, variance floor 0.0001, ACT/365 quoted
  maturities and inverse-year OU rates. Quoted deterministic SPX forwards scale
  normalized payoffs; cash discount factors are not inferred.

## Run ledger

All directories below are under `runs/`. Each contains `run.json`,
`prices.csv` and `normalizer.csv`, with source/target hashes, complete model
configuration, numerical controls, quote-level residuals and uncertainty.

| Run | Seed | Time grid | Normalizer orders |
| --- | ---: | --- | --- |
| `m2_spx_reference_001` | 20261021 | max 1/(252×2) years | 20, 32 |
| `m2_spx_fresh_seed_001` | 20261022 | max 1/(252×2) years | 20, 32 |
| `m2_spx_timestep_001` | 20261021 | max 1/(252×8) years; separate random allocation | 20, 32 |
| `m2_spx_coupled_grid_001` | 20261021 | exact 4-way refinement of the 2-step/day event grid | 20, 32 |
| `m2_spx_quadrature_001` | 20261021 | same coupled coarse/fine grid | 32, 48 |

The separate 2- and 8-step/day runs are not paired despite their equal seed:
changing step counts changes normal allocation. In the coupled runs, the
coarse stock calculation instead uses exact fine-grid OU states at coarse
times and sums the **same** fine equity increments. Paired payoff differences
therefore isolate stock-integration changes conditional on the numerical
normalizer; their standard errors use independent path-level differences.

## Findings and accuracy status

1. **Gaussian normalization can be assessed independently of stock draws.**
   The largest relative order-20-to-32 normalizer change is about 0.374%;
   order 32 to 48 reduces this to about 0.112%. At the fine grid, the largest
   corresponding 32-to-48 price changes are 0.00396 and 0.00547 SPX index
   points for 29/57 days, at most about 0.011 frozen price scales. This is a
   useful refinement diagnostic, not a certified error bound.

2. **Stock sampling error dominates quote-relative accuracy at this budget.**
   In the final diagnostic, the maximum option-price MC SE is 4.52 points
   (29d) and 6.46 points (57d). Zero quotes pass the plan's sampling gate
   `2*SE <= price_scale/4`. These are standard errors conditional on the chosen
   quadrature and stock discretization; they do not include bias.

3. **The stock-forward checks are consistent with one at this budget.**
   Fine-grid terminal-ratio means are 1.000057 (SE 0.000832) and 0.999594
   (SE 0.001174). This supports the adapted implementation check; a sample
   expectation check is not a proof of the continuous-time model's martingality.

4. **Time integration remains unresolved.** At the strikes with the largest
   absolute paired coarse/fine price differences, changes are +4.50 points
   (29d, put K=6630; paired SE 2.28) and -7.87 points (57d, call K=7180;
   paired SE 3.51). These are large relative to quote scales, with substantial
   sampling uncertainty. They do not establish a converged price or an exact
   bias. The frozen feedback grid reaches 5000 inverse years, so exact OU
   endpoint simulation alone does not remove stock-integral discretization.

5. **The unfitted anchor does not reproduce the quote intervals.** In the
   fine-grid/order-48 run, midpoint-price RMSE is 9.07/14.06 points and interval
   coverage is 1.22%/0.94% at 29/57 days. These low coverages describe one
   synthetic anchor, fixed reference state and diagnostic curve with unresolved
   numerical error; they do not reject M2's fitting capacity.

No numerical gates are relaxed to improve these diagnostics. Optimization and
comparative-fit claims remain premature at the current accuracy.

## Next work

Implement SPX variance reduction and further coupled refinement before expensive
optimization. In parallel in the implementation sequence, add conditional VIX
valuation using the **date-zero** normalizer and full conditional OU state.
Then profile the declared nuisance curves and build the common-protocol ESN
price evaluator. The physical-compatibility and held-out rules in the existing
calibration plan remain unchanged.

Reproduce the most refined diagnostic from the repository root with fresh output:

```sh
python -m matteo.model_comparison.option_pricing.m2_spx \
  --refinement-factor 4 --quadrature-order 32 --check-order 48 \
  --output /tmp/m2_spx_quadrature_new
```

Run the package checks with its source directory on the Python path:

```sh
PYTHONPATH=matteo/exogenous_reservoir_eq_m2/src \
  python -m unittest discover -s matteo/exogenous_reservoir_eq_m2/tests -p 'test_*.py'
python -m unittest matteo.model_comparison.option_pricing.test_benchmark
```
