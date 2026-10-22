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

The next sampling/profile checkpoint is [SAMPLING_PROGRESS.md](SAMPLING_PROGRESS.md).
It records the matched-cost iid/Sobol experiment and a constrained nuisance-curve
starting point, with commands for `vix_outer_sampling` and `curve_initializer`.

Cached exposures, the exploratory VIX curve fit and fresh-scramble checking are
recorded in [CURVE_PROFILE_RESULTS.md](CURVE_PROFILE_RESULTS.md). Commands are
`vix_curve_fit` and `curve_spx_check`. The fixed-sample fit improvement reverses
on checking scrambles; it must not be promoted as a reliable market fit.

[NUMERICAL_STABILITY_RESULTS.md](NUMERICAL_STABILITY_RESULTS.md) investigates
that reversal and the SPX sampling bottleneck. It includes the paired audit,
matched-state quadrature check, independent pilot coverage audit and a
three-grid SPX hedging-control experiment. No full-surface convergence is claimed.

[VIX_SCOPE_DECISION.md](VIX_SCOPE_DECISION.md) records the scope history.
**The active analysis scope is again all four VIX maturities: 27/55/118/153
days, 175 training options and four forwards**, following the user's latest
instruction. SPX scope is unchanged. `active_scope.json` is authoritative.
The earlier near-only refit in [VIX_NEAR_RESULTS.md](VIX_NEAR_RESULTS.md) is
historical; its gain was smaller than numerical uncertainty. Reproducing it
requires `vix_near_fit --scope .../near_scope_historical.json` explicitly.
New full-scope sampling diagnostics use `vix_importance_sampling`.

[VIX_PILOT_SAMPLING_RESULTS.md](VIX_PILOT_SAMPLING_RESULTS.md) follows up with
pilot-learned defensive proposals, frozen before larger fresh all-four-maturity
checks (`vix_pilot_sampling`). `vix_sampling_rank_audit` then evaluates the two
original frozen curves on the same weighted observations, without refitting.

[VIX_STRATIFIED_SAMPLING_RESULTS.md](VIX_STRATIFIED_SAMPLING_RESULTS.md) compares
fresh random and fixed component allocations under the identical frozen
proposal (`vix_stratified_sampling`). The exact allocation correction and
component-wise prefix handling are recorded and checked by
`verify_vix_stratified_sampling`; all four maturities remain included.

[VIX_MOMENT_DIAGNOSIS.md](VIX_MOMENT_DIAGNOSIS.md) attributes the new consistency
warning using direct Gaussian expectations, exact old-state refinements and
multiplier-tail diagnostics. It favors outer sampling sensitivity at the frozen
model; it does not patch normalization or certify prices. All four maturities
remain included. `verify_vix_moment_diagnostics` replays artifact checks.

[VIX_SCORE_SAMPLING_RESULTS.md](VIX_SCORE_SAMPLING_RESULTS.md) tests a deterministic
model-score subspace proposal against the frozen pilot-learned mean mixture on
fresh matched-budget R16 x 64 draws, with all four maturities retained. The
proposal broadens four directions while retaining the full 36-dimensional state
and exact weights. `verify_vix_score_sampling` replays full states, weights and
statistics; `audit_vix_score_sampling` produces descriptive method comparisons.
This experiment does not refit the model or impose the moment identity.

[VIX_TARGETED_PILOT_RESULTS.md](VIX_TARGETED_PILOT_RESULTS.md) compares four
independent payoff-learning pilots, explicitly repurposing old checks as design
data, then checks one frozen mean-shift proposal on fresh equal-budget draws.
The new proposal helps some option families but worsens every forward SE and
has a new 118-day moment warning. It remains experimental. Reproduction uses
`vix_targeted_pilots`, `vix_targeted_check` and their corresponding verifier
modules; all four maturities and forwards are retained.

[VIX_ALLOCATION_REFINEMENT_RESULTS.md](VIX_ALLOCATION_REFINEMENT_RESULTS.md)
replays the exact 118-day warning states at higher conditional and normalizer
orders, then compares random and fixed allocations under the identical frozen
17-component density on fresh R12 x 128 draws. The density weights and allocation
corrections are recorded separately; all four maturities remain active.
Reproduce with `vix_targeted_refinement`, `vix_targeted_allocation` and their
corresponding verifiers; the same-density comparator is random targeted sampling.

[VIX_GAUSSIAN_CONTROL_RESULTS.md](VIX_GAUSSIAN_CONTROL_RESULTS.md) trains nine
Gaussian controls on independent complete pilot replications, then compares
raw/mass/pooled/half-pilot estimates on identical fresh R16 x 128 paths at all
four maturities. Primary forward SE decreases everywhere; 27/55-day forward
SE screens pass, but longer maturities, pilot stability and most option quotes
remain unresolved. Commands: `vix_gaussian_controls`,
`verify_vix_gaussian_controls`, `audit_vix_gaussian_controls`.

[VIX_GAUSSIAN_STABILITY_RESULTS.md](VIX_GAUSSIAN_STABILITY_RESULTS.md) increases
compatible independent pilot replications from 12 to 28, retains the old fit,
and checks pooled/balanced pilot fits on fresh nested 128/256-point nets at
all four maturities. Reproduce with `vix_gaussian_stability`,
`verify_vix_gaussian_stability` and `audit_vix_gaussian_stability`. Historical
observations used for the larger pilot are explicitly design data.

[NEXT_STEPS.md](NEXT_STEPS.md) records the user's direction to leave the VIX
tail-control extension pending, retain all four maturities/forwards, and move
forward. It records SPX numerical diagnostics as the working next focus, with
all existing scope and accuracy limitations preserved.

[SPX_PATH_IMPORTANCE_RESULTS.md](SPX_PATH_IMPORTANCE_RESULTS.md) records the
SPX full-path Gaussian shift pilot, exact likelihood weights and independent
hedging-control pilots, with coupled 64/128/256 grids and separate GH48/64
normalizer changes. Reproduce with `spx_path_importance`; replay checking draws
and source integrity with `verify_spx_path_importance`. The diagnostic preserves
all maturity/split scopes and does not authorize market refitting.
