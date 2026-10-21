# ESN option pricing: code and report correspondence

Audit started 21 October 2026. Scope: `CrossSectionalData/`, relative to the
repository root. The original files are preserved. This is a source review
and an initial numerical diagnostic, not a reproduced optimization or an M2 fit.
The historical calibration audit lives separately in
[`../othmane_replication/OTHMANE_APPROACH_ISSUES.md`](../othmane_replication/OTHMANE_APPROACH_ISSUES.md).

## Recorded discrepancies

1. **VIX fitting range differs from the report's wing discussion.**
   `calibrate_smile.py`, `VIX_MONEYNESS_RANGE` and `load_vix_points()`, retain
   only `0.40 <= K/F <= 3.00`. The report's final VIX-smile discussion describes
   calibration and weights at moneyness 5–9. Increasing the weighting bandwidth
   cannot include quotes already removed by the loader. The runner records
   actual retained targets and raw rows outside the code's range. The figures
   may come from another code revision; this has not been established.

2. **Failed implied-volatility inversions disappear from the loss.**
   `calibrate_smile.py`, `residuals()`, sets NaN differences to zero;
   `diagnostics()` removes NaNs before computing errors. This is acknowledged
   in the report's final VIX figure caption, but limits any claim to fit all
   quoted strikes. The runner retains failures and their original weight,
   labels surviving-point RMSE explicitly, and does not optimize this loss.

3. **The plotted/labeled objective is not the executed optimizer objective.**
   The report writes a weighted IV sum of squares and describes a convergence
   plot as that sum divided by 858 points. Executed `residuals()` additionally
   divides by fixed reference-fit scales and `sqrt(number of expiry specs)`,
   applies expiry multipliers (`vix_near=3`, `vix_next=1.5`), and the entry point
   uses `least_squares(loss='soft_l1', f_scale=1)`. Multipliers in the residual
   are squared in the ordinary sum of squares. The history logger writes the
   raw residual sum of squares; the checkpoint after fitting stores twice the
   robust optimizer cost under the same `sumsq` key. These quantities must not
   be treated as interchangeable. Exact original scales/history are not in
   either supplied ZIP archive.

4. **Several headers describe an older calibration.**
   `calibrate_smile.py`'s opening documentation says seven parameters and four
   expiries; `unpack()` and `CombinedCalibrator.__init__()` implement eight
   parameters and six expiries. Old bandwidth descriptions also differ from
   current constants. The report's eight-parameter/six-expiry description
   agrees with the executable class on those two points.

5. **Printed final parameters cannot reproduce the exact optimizer output.**
   The report prints approximate final parameters. The entry point contains
   an exact reference `theta0`, not the printed final solution. Neither ZIP
   archive contains final parameter, scale, or history JSONs. The runner
   labels its `report-rounded` and `reference` presets distinctly. Agreement
   or disagreement at rounded parameters is not an exact reproduction claim.

6. **Standalone execution depends on missing original filesystem paths.**
   Several scripts hard-code `/mnt/user-data/uploads/...`; output paths expect
   sibling `results/` directories. The runner substitutes explicit local CSV
   paths in memory and writes only to a fresh output directory. It does not
   change numerical equations or target-processing rules.

7. **Archive and checkout are not wholly byte-identical.**
   All ten top-level Python files were compared with `esn_market_test_code.zip`.
   Nine match; `joint_calibration.py` differs in historical target paths and a
   leverage-target index correction. The option engine only imports its fixed
   architecture and parameter constants, which match. These differences do not
   explain the option code/report discrepancies above.

## Pricing conventions and limitations to keep separate from bugs

- The daily Gaussian log-return correction makes the discrete stock ratio a
  martingale under the simulator's chosen law. It does not identify that law
  with an empirically estimated change from the historical measure.
- Burn-in averages over an ensemble of possible current states. It does not
  condition prices on a state filtered for 19 February 2026.
- VIX strikes are evaluated at market moneyness times the **model's** VIX
  forward. A shape fit therefore does not enforce market forward levels or
  prices at the actual absolute strikes. The report's model means around
  15–16 versus market forwards around 20–22 remain material level differences.
- The conditional VIX uses a finite SPX strike range (0.60–1.60) and finite
  inner simulation. Wing truncation, integration error, inner bias and outer
  sampling error must be measured separately before comparing model families.
  A deterministic 20% constant-volatility oracle at 30/365 years gives variance
  0.0383673 on the calibration's 29-point grid versus the exact 0.04; refining
  to 284 points gives 0.0399595. This isolates a real grid effect without any
  Monte Carlo noise, but is not a bias estimate for the calibrated ESN.
- Maturities are rounded to trading days for diffusion, but exact calendar-day
  fractions enter Black inversion. The runner records both clocks.
- The original validation script compares independent samples. The new check
  also compares scalar and batched single-path outputs with the same seed,
  across smoothing factors and bias offsets, and checks conditional restart.
- Put/call mid-IV pooling is reproduced for the legacy benchmark only. Our
  common comparison will use the existing audited bid/ask price pipeline.

## Reproduction status

`benchmark.py` records source and input hashes, archive byte correspondence,
versions, parameters, seed, budgets, target counts, every pricing point,
failed-IV weight, and VIX forward estimates. Its mean standard errors do not
include finite-inner or strike-truncation bias. It does not reconstruct the
legacy robust/scaled optimization objective or establish global optimality.

Next: assess independent-seed and sampling/grid refinement; recover the exact
final parameters and scales if available; then freeze the common ESN/M2 target
contract before optimization. Missing provenance does not block a new,
explicitly versioned benchmark fit.

## Reporting correction in our own VIX diagnostics

The first M2 VIX runner counted `2*sample_SE <= price_scale/4` as an MC gate
without checking whether any positive payoff had been observed. A deep OTM
option with zero payoffs in the outer sample has zero empirical SE; that does
not bound its rare-event expectation. This is a diagnostic reporting problem
in **our runner**, not a newly identified defect in Othmane's source.

The revised runner reports payoff-hit counts, separates the nominal SE screen,
and excludes fewer-than-two-hit cases from even the minimal screen. Two hits
still do not certify tail accuracy. Earlier `run.json` files are preserved;
`audit_tail_samples.py` adds `tail_diagnostics.json` with input/source hashes.
In the 256-outer-state run, all 21 nominal screen passes are zero-hit cases
(1/6/7/7 across the four expiries). No quote passes the corrected minimal
screen. Model prices, residuals and original records are not rewritten.

The first development runs recorded source hashes at completion. That is not
an immutable record of loaded runtime sources if files change during a long
calculation. The new runners capture sources before calculation, reject
mid-run changes, and save the complete captured files under `source/` in each
run. The `m2_vix_frozen_*` and `m2_spx_frozen_controls_001` runs replay the main
diagnostics with this stronger provenance. Initial runs remain preliminary
records; no old hash or numerical result is silently updated.

## Existing time-step decision not carried into initial pricing grids

On 21 October the user pointed back to the model's high mean-reversion rates.
The pre-existing `experiments/timestep_convergence_001/RESULTS.md` in the M2
package explicitly rejects 2/4/8 steps per day for final evidence, sets 64 as
the working baseline and requests a targeted 64-versus-128 confirmation. The
model source's later fine-grid appendix also records 64-step experiments.
The new SPX diagnostics nevertheless began at 2/8 steps per day. They were
labeled numerically unresolved, but the initial implementation work did not
carry forward that existing resolution decision. Those runs are coarse engine
and variance-reduction diagnostics only; they cannot establish accurate
continuous-time option prices.

At the frozen maximum rate 5000/year, lambda_max*dt is approximately
9.92/2.48/0.310/0.155/0.0775 at 2/8/64/128/256 steps per trading day. Exact OU
endpoints do not make the nonlinear stock integral exact between endpoints.
The next option-pricing convergence study must keep the model/state fixed and
compare coupled 64/128 grids, extending to 256 if price-relative tolerances
remain unresolved. The old physical-statistic convergence study motivates this
starting resolution; it does not itself establish option-price convergence.

The calibration plan's numerical-budget paragraph also still prescribed a
2/8/32 starting ladder after that decision. It has now been updated to the
64-step working stock baseline and coupled 128/256 refinements. This is a
correction to the implementation plan, not a retrospective change to the
recorded budgets or conclusions of earlier runs.

## Faster VIX time integration and common loss boundaries

The new Gaussian **time** rule is distinct from the three-dimensional
Gauss-Hermite conditional-expectation rule. It resolves lag panels starting at
1/lambda_max, growing geometrically, with every curve jump included. Its
smaller node count must not be described as an equivalent coarse SPX stock
timestep. SPX integration is unchanged. Fixed-state and sampled-state order
checks are recorded separately in `CALIBRATION_PROGRESS.md`.

`JointPriceObjective` now implements the planned family/expiry weighting and
separate VIX forward family. Its first anchor score uses only the two already
priced SPX maturities, plus all four VIX maturities. It is therefore a
restricted-maturity diagnostic rather than the complete planned training loss.
Its inherited input conflict and sampling errors remain open. No optimization,
P-feasibility or ESN/M2 ranking follows from a small scalar score.

## Dated-state coordinate order

The calibration plan's state paragraph listed `(Y,F,C)`, whereas the dated
pricing API and its projection matrix use `(C,Y,F)`. Either ordering is
mathematically valid with consistently permuted covariance and H, but passing
the former directly into this API would condition on the wrong state. The plan
now explicitly uses `(C,Y,F)`; the restricted-state constructor and tests name
that order. Existing runs used zero projections, so this documentation mismatch
did not change their inputs or results.

## First curve fit: fixed-sample improvement reverses on fresh scrambles

The exploratory VIX-only nuisance-curve run uses 256 outer points per maturity
for fitting and another independent 256 for checking, below the calibration
plan's 4096-point screening budget. This is explicitly a preliminary machinery
and sensitivity study, not final screening or price qualification. The original
plan's budgets and accuracy gates have not been relaxed.

Its selected curve reduces fitting loss from 72.6724 to 62.3668, but on fresh
scrambles scores 72.2576 against the initializer's 66.8887. All eight optimizer
starts converge to essentially the same fixed-sample solution. Optimizer
convergence must not be reported as reliable market-fit improvement. Checking
scrambles reuse training contracts and do not constitute held-out market data.
The reversal and remaining tail/quadrature errors are preserved in
`CURVE_PROFILE_RESULTS.md` and `runs/m2_vix_curve_fit_001`; no candidate was
reselected using checking loss.

The subsequent paired audit reproduces those prices and isolates the main
153-day block. Its comparison reverses by about 10.8 loss units between outer
samples, while conditional/normalizer quadrature refinement shifts it by only
about 0.0031. This supports outer-sample sensitivity at that block; it does not
justify declaring quadrature accurate everywhere or ranking the two curves.
See `NUMERICAL_STABILITY_RESULTS.md` for loss-scale, tail-concentration and
leave-one-scramble-out diagnostics.

The SPX control experiment also retains a zero-raw-payoff case whose empirical
plain SE is zero. Its augmented estimate has nonzero uncertainty; neither is
an accurate tail-price certificate. Entirely zero pilot targets can produce
zero regression coefficients, so tail coverage must be monitored during
pilot training as well as evaluation. The new three-level stock study couples
all 64/128/256 grids within a single finest sample, avoiding the earlier
cross-pair limitation while leaving continuous-time accuracy unresolved.

## Explicit two-maturity VIX scope and identifiable curve dimension

The user subsequently authorized deferring both 118/153-day VIX options and
their forwards. `active_scope.json` now governs the exploratory 27/55-day
objective; SPX scope is unchanged. Earlier four-maturity scope decisions are
historical and superseded for active fitting, while their data/results remain
preserved. The deferred training maturities are not relabeled as pristine
held-out observations.

Both retained VIX windows lie before the first nuisance knot at 0.25 years.
The active fit therefore identifies only one variance coefficient. Inactive
segments are equal-level computational placeholders, not fitted longer-term
values; their smoothness penalty is zero. This prevents implying that a
two-parameter curve has been estimated from the reduced data. Current model
parameters and scope cannot be ranked as full-maturity fits.

The new near-only fit's fresh-sample gain is about 0.0352 loss units, versus
descriptive paired jackknife uncertainty about 0.2760. The large prior reversal
is absent in these point comparisons, but reliable improvement and convergence
are still unestablished. See `VIX_NEAR_RESULTS.md` for full numerical budgets,
scope limitations and separate sampling/integration checks.

## Scope restoration and importance-sampling boundary — 21 October 2026

The latest user instruction restores both previously deferred VIX maturities
(118/153 days) and their forwards for new analyses. The current scope is
27/55/118/153 days: 175 training options and four forwards. This supersedes
the earlier same-day near-only decision; the prior exclusion and its results
remain historical records. `active_scope.json` is authoritative and the old
policy is saved as `near_scope_historical.json`. The one-level near-only
driver rejects the restored scope and requires an explicit historical policy
for reproduction, preventing a silent omission of longer maturities.

Restoration does not clear input conflicts, approximation errors or feasibility
gates. The new Gaussian-mixture diagnostic changes only the sampling law,
using exact likelihood weights to preserve the original model expectation.
No market-forward substitution, self-normalized weighting or parameter
refitting is used. Its outer standard errors exclude conditional GH, time
integration and normalizer biases. See `VIX_IMPORTANCE_SAMPLING_RESULTS.md`.

## Pilot reuse and weighted fixed-candidate comparison — 21 October 2026

The previous defensive-mixture observations are now explicitly pilot data for
`vix_pilot_sampling`; they cannot also serve as fresh validation of the learned
proposal. Prior ordinary-Sobol outcomes remain historical evidence. Fresh
ordinary-Sobol and pilot-mixture seed streams are generated after all four
proposals are frozen and hashed. One fitted pilot cannot establish robustness
across repeated independent pilot fits. The far-call pilot contribution ESS
is only about 5–12, so proposal overfitting remains a recorded concern.

Pilot target normalization learns proposal means/covariances; evaluation still
uses the exact original Gaussian / proposal density ratio, without dividing
by its sample mean. First-eight-coordinate selection and the norm-3 cap are
heuristics declared before evaluation. All four maturities and forwards remain
in scope. The new run holds conditional, time and normalizer orders fixed;
selected-state GH refinements cannot certify complete-price convergence.

The paired ranking follow-up changes only the two previously frozen variance
curves. Its algebraic VIX transformation requires an entire window within one
curve segment and a curve-independent date-zero normalizer. The implementation
rejects curve-knot crossings and checks against the archived conditional-
exposure cache at each expiry. No nonlinear-loss ranking is certified by an
eight-scramble descriptive jackknife. See `VIX_PILOT_SAMPLING_RESULTS.md`.

## Fixed component allocation and density-weight distinction — 21 October 2026

The all-four-maturity follow-up freezes the preceding learned Gaussian proposal
and compares fresh random component selection with fixed per-component counts.
The integers (32,6,6,5,5,5,5) at N=64 are only approximately proportional to the
probabilities (1/2,1/12,...,1/12). A plain average using only f/g under these
fixed counts would generally estimate the wrong expectation. The implemented
estimator is sum_k pi_k mean_k[(f/g)H], equivalently a plain average with factor
(N*pi_k/n_k)*(f/g). Density weights and deterministic allocation factors are
saved separately. This is exact allocation correction, not self-normalization.

Coarse prefixes must be taken within each component, not from global row order.
Complete sets of independent component scrambles define a replication. Tests
cover the allocation correction under deliberately unequal component counts and
the prefix indexing distinction. Non-power-of-two component counts use partial
Sobol nets; no full-net balance claim is made. No control from the previous
negative-price experiment is applied. All original scope and numerical-error
limitations remain; see `VIX_STRATIFIED_SAMPLING_RESULTS.md`.

The QMC layouts also differ: one global (d+1)-dimensional net with a component
selection coordinate versus seven independent d-dimensional component nets.
Any observed improvement belongs to the complete stratified sampling design;
it cannot be attributed solely to reduced random component counts. Equal
conditional-pricing counts do not imply identical random-number generation cost.

## New numerical identity discrepancy, not a cleared normalizer — 21 October 2026

The fixed-allocation study retains all four VIX maturities and forwards but
does not establish a uniform precision gain. The 27-day E[VIX^2] target under
the exact normalized model is 542.49611; fresh random and stratified estimates
are 556.16551 ± 3.23481 and 554.03172 ± 3.10259 (reported scramble SE). The
exact-mass-centered residual w*(VIX^2-target) remains +15.66458 ± 2.75395 and
+12.60904 ± 3.46012. This discrepancy is preserved in
`runs/m2_vix_stratified_moment_audit_001`. It is not corrected by imposing the
moment identity or renormalizing prices. Its interpretation remains open: outer
tails/SE estimation, conditional GH composition, global M_Q quadrature and time
integration must be separated. The exploratory audit was prompted by the
observed second moments; no predeclared significance test is claimed.

Earlier noisy moment agreement is a historical observation, not normalization
clearance. Likewise the preceding pilot's small selected-state GH changes do
not cover every new tail state: the 55-day GH32/48 tail-state difference is now
0.022312 VIX points. Known-Gaussian likelihood-mass checks concern sampling
weights and do not establish accuracy of the model's separate M_Q normalizer.

Fresh fixed-curve loss changes also differ in sign between random mixture
(-0.21042) and stratification (+0.93316), with descriptive jackknife SEs about
1.45 and 1.12. Preserve the previous positive changes around +2.7/+2.9 as
historical results; no final ranking follows from any of these small studies.
See `VIX_STRATIFIED_SAMPLING_RESULTS.md` for full records and next diagnostics.
