# ESN/M2 option-pricing benchmark

First milestone: review and diagnose the preserved ESN engine before building
M2 option pricing. See [AUDIT.md](AUDIT.md) for code/report inconsistencies.
See [RESULTS.md](RESULTS.md) for the three initial diagnostic runs and their limits.
No original ESN sources or historical experiment records are changed.

The dated M2 SPX prototype, common price targets, conditional VIX engines and
SPX independent-pilot controls are now implemented. See
[M2_SPX_RESULTS.md](M2_SPX_RESULTS.md) and [M2_VIX_RESULTS.md](M2_VIX_RESULTS.md)
for unfitted diagnostics. Numerical accuracy gates and joint fitting remain open.

```sh
python -m matteo.model_comparison.option_pricing.m2_spx \
  --refinement-factor 4 --output /tmp/m2_spx_coupled_new
```

This reads training quotes at 29/57 days only. It does not calibrate or access
held-out prices for scoring. The copied preprocessing flags and numerical
diagnostics must accompany any interpretation of its output.

```sh
python -m matteo.model_comparison.option_pricing.m2_vix \
  --output /tmp/m2_vix_new
python -m matteo.model_comparison.option_pricing.m2_spx \
  --pilot-paths 2048 --refinement-factor 4 \
  --quadrature-order 32 --check-order 48 --output /tmp/m2_spx_controls_new
```

VIX uses actual training strikes and separate market forward targets. For old
VIX runs, read `tail_diagnostics.json` with the original nominal SE counters;
zero sampled payoffs cannot establish tail-price accuracy.
New M2 runs capture executable sources before calculation, reject mid-run
source changes, and save the snapshot under each run's `source/` directory.

The 64/128-step-per-day checks are recorded in
[FINE_GRID_RESULTS.md](FINE_GRID_RESULTS.md). Use the bounded-memory driver:

```sh
python -m matteo.model_comparison.option_pricing.fine_grid --mode spx \
  --output /tmp/m2_spx_64_128_new
```

Standalone SPX and VIX runners now default to 64 steps/day. To reproduce old
coarse runs, specify their original step counts explicitly or use their saved
source snapshots. Finer grids do not remove Monte Carlo sampling error.

From the repository root, with NumPy, SciPy and pandas installed:

```sh
python -m unittest matteo.model_comparison.option_pricing.test_benchmark -v
python -m matteo.model_comparison.option_pricing.benchmark \
  --output matteo/model_comparison/option_pricing/runs/report_rounded_smoke_001
```

The default budget (512 SPX paths, 32 VIX outer draws, 64 inner continuations,
500 burn-in days) is a smoke diagnostic, not calibration-quality evidence.
Use `--parameters reference` for the exact entry-point starting vector;
`report-rounded` uses the approximate final vector printed in the report.
Both preserve the actual executable target selection and weights. Explicit
CSV paths are available through `--spx` and `--vix`.

Each new run directory contains `run.json` and `points.csv`. Existing output
directories are refused. No optimization is launched and no claims about
matching the published optimum follow from this command. Larger-budget runs
must refine outer/inner counts and strike integration separately, use fresh
seeds, and report failed prices as well as surviving-point errors.

The implementation sequence remains the one in the
[M2 calibration plan](../../exogenous_reservoir_eq_m2/documentation/CALIBRATION_PLAN.md):
reproduced ESN benchmark; common audited targets; dated M2 SPX pricer;
conditional M2 VIX pricer; joint fitting; physical-fit and held-out evaluation.

The shared price loss, synthetic optimizer check, Gaussian VIX time rule and
first restricted joint anchor score are recorded in
[CALIBRATION_PROGRESS.md](CALIBRATION_PROGRESS.md). Reproduce with fresh output
directories using `vix_time_rule` and `joint_diagnostic` as described there.
Neither diagnostic performs market calibration or establishes price accuracy.
