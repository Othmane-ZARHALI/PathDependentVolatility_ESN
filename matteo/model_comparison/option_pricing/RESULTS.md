# Initial option-engine diagnostics — 21 October 2026

## Scope and conclusion

The legacy ESN pricing engine can be exercised on the local snapshot without
editing its sources. The new same-seed scalar/batch correspondence tests pass.
The benchmark retains 858 pooled-IV targets across six expiries, matching the
report's total count, but does not include VIX targets beyond moneyness 3.

**No final-fit reproduction or M2 option fit is established.** These runs use
the report's rounded parameter vector, not an archived exact optimizer result.
The supplied archives contain code and figures, but no final parameter/scales
ledger. See [AUDIT.md](AUDIT.md) for the code/report discrepancies.

## Executed runs

All use 500 burn-in days and the executable legacy quote selection, weights,
daily simulation, finite strike grid, and model-forward moneyness convention.
No optimization is performed.

| Run directory under `runs/` | SPX paths | VIX outer | VIX inner | Seed |
| --- | ---: | ---: | ---: | ---: |
| `report_rounded_smoke_001` | 512 | 32 | 64 | 90000000 |
| `report_rounded_refinement_001` | 2000 | 128 | 256 | 90000000 |
| `report_rounded_fresh_seed_001` | 2000 | 128 | 256 | 91000000 |

Each contains `run.json` (legacy source hashes, archive correspondence, input
hashes, versions, budgets, parameters and diagnostics) and `points.csv`
(all targets and failures). The two larger runs use the same numerical budget
but independent seeds. Changing path counts changes random-number allocation;
the first two runs are not a paired nested-sample bias estimator. Both inner
and outer counts change together, so this is not an inner-bias convergence study.

## Larger-budget outcomes

RMSE below is weighted IV error on **finite model-IV points only**, in decimal
volatility units. Thus 0.03 means three percentage points of annualised IV.
Weights are renormalized over surviving points as in the legacy diagnostic;
this is not the scaled/robust optimization objective. Columns with failures
must not be interpreted as complete-smile scores.

| Expiry | Targets | RMSE, seed 90000000 | Failed IVs | RMSE, seed 91000000 | Failed IVs |
| --- | ---: | ---: | ---: | ---: | ---: |
| SPX 29d | 362 | 0.03445 | 13 | 0.02224 | 0 |
| SPX 57d | 337 | 0.03248 | 0 | 0.02089 | 0 |
| VIX 27d | 41 | 0.17798 | 0 | 0.19411 | 2 |
| VIX 55d | 43 | 0.25777 | 0 | 0.30900 | 2 |
| VIX 118d | 37 | 0.36299 | 9 | 0.16410 | 2 |
| VIX 153d | 38 | 0.29921 | 5 | 0.11800 | 2 |

For VIX 118d, failed points carry 12.75% of the original weight in the first
larger run versus 0.33% in the fresh-seed run. At the smallest smoke budget,
49.46% disappears. This demonstrates why surviving-point RMSE needs an explicit
coverage diagnostic; it does not establish whether a failure disappears in
the infinite-sample limit.

Larger-run model VIX forwards range from 14.71 to 15.10 for seed 90000000 and
15.11 to 16.21 for seed 91000000, versus the snapshot's 20.10–21.72. These
are preliminary estimates with unresolved inner-MC and strike-grid effects,
not converged forward errors. They nevertheless show that a smile-only
comparison in each model's own moneyness is not a forward-level fit.

## Implementation checks

Six focused unit tests pass:

- Scalar versus vectorized same-seed paths for two smoothing factors and two
  bias offsets, with nonzero directional quadratic loading: maximum return
  discrepancy approximately 2.1e-17 and variance discrepancy 6.5e-19.
- Restart from a nonzero full state agrees with explicit continuation steps;
  caller-owned state arrays are preserved.
- Black price/IV round trips recover a known constant volatility.
- Constant-volatility finite-grid variance replication approaches its oracle
  under refinement (details in the audit).
- Partial and complete IV failures remain visible in the benchmark metrics.

These are implementation correspondence and numerical checks. They are not
formal proofs, distribution-wide equivalence guarantees or pricing convergence
certificates. Original `CrossSectionalData/` files remain unchanged.

## Next implementation step

Freeze a common target interface using our existing audited bid/ask price
pipeline and actual VIX forward targets. Add the dated M2 SPX pricer and its
independent normalizer checks, then conditional Gaussian/nested-MC VIX pricing.
The original ESN benchmark remains available as a separately labeled diagnostic.
Do not use this six-expiry in-sample comparison as held-out validation; retain
the original target/split boundaries from the M2 calibration plan for the
subsequent common-protocol study.

Reproduce the larger runs from the repository root, choosing fresh output paths:

```sh
python -m matteo.model_comparison.option_pricing.benchmark \
  --spx-paths 2000 --outer 128 --inner 256 --seed 90000000 \
  --output /tmp/esn_report_rounded_refinement_new
python -m matteo.model_comparison.option_pricing.benchmark \
  --spx-paths 2000 --outer 128 --inner 256 --seed 91000000 \
  --output /tmp/esn_report_rounded_fresh_seed_new
```
