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

## Moment-warning attribution and rare multiplier concentration — 21 October 2026

The original warning remains a valid observed discrepancy, but new diagnostics
favor outer-sample/tail sensitivity over a large normalization error at this
frozen model. The independent direct Gaussian expectation gives 542.49136 ±
0.09509 at 27 days against target 542.49611 and agrees at all other maturities.
All original 1024 random-mixture states retain a residual +13.65296 after
conditional32→48 and normalizer48→64; the paired rule change is only -0.01643
± 0.01193. Limited time/higher-order checks and full Gaussian-law composition
checks support this attribution without globally certifying quadrature errors.

The separate population-composition calculation remains noisy. Its outputs,
including the large paired GH48 comparison, are retained as inconclusive, not
used to claim conditional quadrature convergence. The initial composed-run
startup failure was a relative-path source-snapshot issue, corrected before
numerical work in run 002; failed directory 001 is preserved with its status.

At representative middle dates, future G states with probability about
0.056%–0.072% supply about 49%–53% of E[G]. This concerns the Gaussian variance
multiplier, not directly the expiry VIX payoff distribution. The smooth cap
bounds G at about 188619 but does not prevent severe sampling concentration.
No normalization, parameters, prices or market targets are changed or rescaled.
All four maturities and forwards remain in scope; detailed evidence, limitations,
113 passing tests and replay checks are in `VIX_MOMENT_DIAGNOSIS.md`. The
lean-numerics-bridge workflow is used without claiming Lean certification.

## Score sensitivity is not payoff variance capture — 22 October 2026

The new all-four-maturity frozen VIX experiment chooses four terminal Gaussian
score directions capturing 99.74%–99.77% of its deterministic sensitivity-matrix
trace. On fresh matched-budget checking draws, covariance inflation in those
directions worsens median option SE at every maturity (ratios 2.242/1.858/
1.358/1.806 versus the earlier learned mixture), and worsens every forward SE.
It must not be described as retaining 99.7% of pricing variance or as a
successful sampling improvement. The design and negative result are retained in
[VIX_SCORE_SAMPLING_RESULTS.md](VIX_SCORE_SAMPLING_RESULTS.md); it is not adopted.

The fresh learned-mixture second moments are 543.724/540.539/600.978/593.863,
against targets 542.496/542.496/600.319/600.319, with SEs 5.75/6.90/6.89/6.18.
The old 27-day +13.65 discrepancy does not persist in this new sample; this
supports the prior sampling-sensitivity attribution, not an identity-enforcing
normalizer correction. Only 28/175 learned and 11/175 inflated option estimates
pass the observed SE screen; none of these counts certifies accuracy. No
maturity, forward, model parameter or price normalization has been changed.

The new artifact verifier initially compared Python risk-premium tuples to
JSON lists and stopped at metadata correspondence. Canonical JSON comparison
corrected that verifier-only mismatch; the completed pricing artifacts were
unchanged. Replay now passes for all 8192 full states/weights and 732 summary
rows. Negative results and incomplete numerical gates remain explicit.

## Targeted pilot proxy versus fresh precision — 22 October 2026

Four disjoint pilots favor a payoff-learned score-subspace mean-shift proposal
on a historical independent screening set. The screen estimates an IID
importance-estimator second moment, not variance or randomized-Sobol SE.
All historical observations used to choose it are explicitly promoted to
design data. Fresh R8 x 64 pricing reduces median option SE at 55/153 days
(ratios 0.592/0.500) but worsens it at 27/118 days (1.165/1.066), and worsens
all four forward SEs. Only 24/175 targeted and 19/175 baseline quotes pass the
observed SE screen. This is a recorded limit of the design proxy, not a full
sampling success or calibration result.

**New 118-day moment warning:** targeted E[VIX^2]=588.51408 ± 4.30771 versus
600.31899, residual -11.80491 (about -2.74 reported SE). The baseline 55-day
residual is also +15.72785 ± 7.94602. Previous fresh moment agreement did not
establish convergence. Preserve this warning alongside limited selected-state
GH checks; it does not prove a normalization bug or negligible quadrature
bias. No moment rescaling or maturity exclusion is applied. See
[VIX_TARGETED_PILOT_RESULTS.md](VIX_TARGETED_PILOT_RESULTS.md) and the saved
`m2_vix_targeted_check_001` artifacts. All 119 tests and correspondence replays
pass, without certifying statistical accuracy.

## Exact-state 118-day warning attribution — 22 October 2026

On all 512 original targeted 118-day states and their unchanged likelihood
weights, conditional GH32/M48 to GH48/M64 changes E[VIX squared] by
-0.00603 ± 0.02206, leaving residual -11.81094 versus the 600.31899 model target.
The corresponding forward change is -0.000205 ± 0.000233; the largest 118-day
option price change divided by its quote scale is 0.00464. The observed warning
is therefore not explained by these two refinements on this sample. Small
higher-order/time subset checks remain limited, not global price-bias bounds.
No normalizer patch, price rescaling, maturity removal or new fit is justified.
See [VIX_ALLOCATION_REFINEMENT_RESULTS.md](VIX_ALLOCATION_REFINEMENT_RESULTS.md).

The fresh R12 x 128 same-density allocation comparison is complete. Fixed
17-component allocation lowers median option SE at 27/118/153 days but worsens
55 days (ratios 0.797/1.124/0.931/0.883); observed SE screens total 35/175 fixed
versus 36/175 random. No forward passes its sampling SE tolerance. Fresh moments
are within about 1.7 reported SE of all model targets under both methods,
further supporting sample sensitivity of the old warning, not a convergence
clearance. Half-replication SEs remain variable. The fixed design's effective
weights include allocation factors and can exceed 2 (bound 8/3), while likelihood
weights remain bounded by 2; the larger effective weights are not a density
ratio bug. All 121 tests and both correspondence replays pass. Fixed allocation
is retained as experimental, not adopted as a uniform remedy. All four VIX
maturities/forwards and prior warning records remain included.

## Gaussian controls: sampling gain versus unresolved accuracy — 22 October 2026

Frozen Gaussian controls trained at the complete-replication level reduce fresh
forward SE at all four maturities; primary 27/55-day SEs 0.05735/0.04649 pass
the observed 0.0625 threshold. This must not be reported as full price accuracy:
118/153-day forward SEs remain 0.08702/0.11230, only 39/175 primary option SE
screens pass (37/175 raw), and count/quadrature/pilot gates remain open.
No model refit, maturity exclusion, forced VIX moment or price rescaling occurs.

The two independent pilot-half fits are not uniformly useful: half-A worsens
118-day forward SE and raises its median option SE by about 39%. Its 27-day
forward also fails the screen. **Half-B 153-day moment warning:** 605.93740 ±
2.69630 versus target 600.31899, residual +5.61841 (about 2.08 reported SE).
The pooled primary moment is 599.78679 ± 3.75925. All fits and this diagnostic
warning are retained in
[VIX_GAUSSIAN_CONTROL_RESULTS.md](VIX_GAUSSIAN_CONTROL_RESULTS.md), without
choosing a favorable half after checking or inferring a normalizer bug. Every
reported option estimate is positive; no adjusted value was clipped. All
124 tests and pilot/state/control/statistics correspondence checks pass.

## Larger Gaussian pilots: forward screens versus surface accuracy — 22 October 2026

Twenty-eight compatible historical 128-point pilot replications, explicitly
promoted to design data, train a pooled fit and two disjoint balanced 14-pilot
fits. On fresh R12 x 256 paths, pooled forward SEs
0.04230/0.03750/0.03961/0.06157 pass the observed 0.0625 sampling screen at all
four maturities. The 153-day pass is only about 1.5% below the threshold and
fails for one balanced fit and the unchanged legacy12 fit. It must not be
reported as robust full-price convergence or successful market calibration.

Option screens remain 45/175 pooled versus 47/175 raw; the prior 27-day median
option variance gain does not survive this larger-net check. **Retained 118-day
moment warnings:** pooled residual -4.55861 ± 2.24254 (about -2.03 reported SE),
balanced A -4.21342 ± 1.98663 (about -2.12). All variants, dependent diagnostic
flags and paired count changes are preserved in
[VIX_GAUSSIAN_STABILITY_RESULTS.md](VIX_GAUSSIAN_STABILITY_RESULTS.md), with
no price rescaling, clipping, tolerance change or maturity exclusion. Larger
pilot size and lower forward SE do not settle option tails or population
conditional-quadrature bias. All 126 tests and pilot/state/control/statistics
correspondence checks pass, without certifying pricing accuracy.

## SPX full-path importance pilot — 22 October 2026

The user explicitly confirmed focusing on SPX sampling, stock integration and
normalization. A defensive Brownian path mixture shifts the correlated OU
innovations consistently, uses exact nonnormalized full-path weights and fits
weighted static/adapted/mass controls on independent pilots. Gaussian joint
moments and Black rare-payoff oracles pass; all 129 package/driver tests pass.
A fresh two-trial 512-pilot/1024-checking-path diagnostic preserves all scopes.

**Unsuccessful cases retained:** controlled variance improves in only 13/24
contract/trial cases (median target/importance ratio 1.2944), and neither
method passes any of the 24 positive-hit quote-relative sampling screens.
Target sampling has two raw zero-hit cases and two negative controlled price
estimates; no clipping or zero-hit accuracy claim occurs. Importance sampling
has payoff hits everywhere but remains inaccurate relative to quote scales.

128-to-256-step changes remain large relative to tolerance and mostly noisy;
no continuous-time convergence or bias sign is inferred. Same-path GH48-to-64
changes are below 0.0013 quote scales at this sample, a sensitivity result
rather than a population bound. Full weights include correlated OU shifts and
must multiply target-zero controls too. A raw mean-one stock/control identity
under the proposal would be incorrect. There is no self-normalization, model
refit, maturity pruning or completed MLMC estimator. All 384 estimates, 288
paired changes and 20 monitors replay from checking draws; source hashes and
seed separation pass. See [SPX_PATH_IMPORTANCE_RESULTS.md](SPX_PATH_IMPORTANCE_RESULTS.md).

## SPX volatility-path proposals and independent checking — 22 October 2026

The new full-path Gaussian basis targets equity plus conditional OU projection
excursions. Its exact f/q likelihood retains all path dimensions and an original
Gaussian component of probability 1/2. Independent training and selection choose
the 7-day proposal before checking. It improves observed controlled variance in
20/24 fresh contract/trial cases (10/12 in each trial; median gain 1.4192), and
8/8 near-ATM cases improve. **None of the 24 sampling screens passes.** All
candidate/half-fit results remain; all four VIX maturities and forwards remain
active and no market refit occurs.

**Preserved failures/discrepancies:** launch 001 fails before any market paths
because seven list initializers were assigned to six variables. Its source and
failure status are kept; corrected launch 002 adds an end-to-end collection
test. Equity-comparator trial 1 weight mass is 1.04935 ± 0.01945 (+2.54 SE),
with correlated stock/reference warnings; no mass normalization is imposed.
Vol21 half-B gives -4.77501 ± 7.84805 for the 29-day 3400 put, a 13.60-fold SE
increase relative to its pooled fit. No negative estimator is clipped.

**Time-integration warning persists despite sampling improvement:** selected
vol7 128->256 changes exceed scale/8 in 20/24 cases, five above two paired SEs.
The 57-day 5550 put change is -1.62367 ± 0.67967 against scale 0.20024. Neither
256-step convergence nor an MLMC allocation is established. Correction variance
and warmed costs are saved, with finest-cache/standalone-cost limitations.

An additional 13-control Gaussian-score/ridge dictionary is frozen on design
data while checking paths are generated, before any completed checking archive
or numerical results exist. On identical fresh paths, its pooled vol7 median
variance gain is only 1.0178; no sampling screen passes. It does not reselect the
proposal or a favorable half. All 864 original estimates, 576 paired corrections,
30 monitors and 648 additional score-control estimates replay; all 134 tests
pass. See [SPX_VOLATILITY_SAMPLING_RESULTS.md](SPX_VOLATILITY_SAMPLING_RESULTS.md).

## SPX MLMC: negative efficiency result and 512-step warnings — 22 October 2026

A four-level 64/128/256/512 prototype uses standalone fine-grid state draws,
exact level-specific vol7 likelihoods and independent correction-specific pilot
controls. Independent allocation paths freeze [2752,640,384,192] counts before
checking. The estimator targets the discretized 512-step payoff by an expectation
telescope; its SE adds independent-level variance/N terms. No continuous-model
convergence or globally optimal integer allocation is claimed.

**Negative result retained:** on the fixed four-central-contract objective,
actual variance-cost efficiency relative to direct 512-step simulation is
0.6069/0.5491 in two fresh trials. Only 5/12 and 4/12 contract variances improve;
none of 24 sampling screens passes. Pilot-predicted cost matched the budget,
but realized MLMC cost is 2.5%/20.2% higher than direct. The separate cross-strike
median SEs must not be presented as an MLMC efficiency gain. Correction controls
help relative to raw corrections, but remaining variances decay too slowly
relative to measured costs for this hierarchy.

**Retained warnings:** trial-0 final MLMC/direct differences for the 29-day 3400
put and 7030 call are +2.56/-2.50 combined SEs. Top-level importance mass is
1.09478 ± 0.03365 and 1.07058 ± 0.03476 across trials, with correlated reference
warnings; no weights are renormalized. Six trial-0 level/contract rows have no
fine/coarse payoff hits, so their small raw/controlled uncertainty is not a tail
accuracy clearance. Half-pilot SE inflation reaches 4.00 direct and 2.53 MLMC.

**512 is still unresolved:** 256->512 changes exceed scale/8 in 8/12 and 9/12
cases. The 57-day 6775 put changes -2.57035 ± 0.87963 (trial 0, scale 0.65303);
the 5550 put changes +1.69703 ± 0.69914 (trial 1, scale 0.20024). These and the
sparse tail corrections preclude convergence clearance or model refitting.
All 138 tests pass; allocation, likelihoods, fits, 192 estimates, 384 level rows
and 50 monitors replay. The verified prototype remains a comparator; direct
vol7 is not replaced. All VIX maturities and original scopes remain active.
See [SPX_MLMC_RESULTS.md](SPX_MLMC_RESULTS.md) for sources, costs and the pending
better-coupling design question.

## Endpoint reflection: valid target, negative cost result — 22 October 2026

An SVD-based orthogonal Gaussian reflection preserves coarse OU endpoints
(max observed mismatch 6.31e-14) and equity increments (5.38e-17). The original
f/q weight is used for both fine payoffs, corresponding to the transported
reflected proposal; evaluating the old density at reflected scores would be
incorrect. Each pair is one observation. Full-path Gaussian/Black and adapted
stock-forward oracles pass, and the ordinary sampler is rerun for actual cost
and numerical path correspondence.

**Negative adoption result:** pooled central variance gains 1.01–1.16 do not
cover runtime multipliers around 1.50–1.52. Variance-cost efficiencies are
0.7667/0.6970 at 128->256 and 0.6689/0.7324 at 256->512 across two fresh trials.
All comparisons and half-pilot fits remain. A same-frozen-coefficient replay
finds the averaged residual retains 86–96% of original controlled variance;
its nonzero covariance with the removed difference component is preserved.
No checking-based refit, reweighting or coupling selection occurs.

**Retained warnings:** three paired pooled mean differences exceed two SEs:
29-day 6755 put (-0.43498 ± 0.19928), 57-day 5550 put (+0.42483 ± 0.21037) at
128->256, and 29-day 7950 call (-0.00853 ± 0.00413) at 256->512, all trial 1.
They are finite diagnostics, not proof of different population targets. All
mass/stock/reference checks are within two SEs here; hit counts can still be
only two. Half-pilot SE inflation reaches 2.55 bridge and 3.59 ordinary.

Reflected 256->512 changes exceed scale/8 in 20/24 cases, despite none exceeding
two observed SEs in this check. Earlier 512-step warnings remain active;
absence of a new significance flag is not a convergence bound. This is a
correction study, not a new complete pricing result. All 141 tests pass;
384 estimates, 192 paired differences and 28 expectation monitors replay.
Direct vol7 and ordinary corrections remain the references. See
[SPX_BRIDGE_COUPLING_RESULTS.md](SPX_BRIDGE_COUPLING_RESULTS.md), including the
pending projected-OU innovation control hypothesis and its full-state
conditional-mean requirement. No model refit or maturity exclusion occurs.

## Exact OU innovation controls: direct variance signal, tail warnings — 22 October 2026

Full-state observers record eps=P*(X_next-D*X-b_Q) on nested grids without
changing the stock simulator. Nine predictable linear, square and equity-cross
martingale gains use full target covariances (not conditional-on-equity
covariances) and full-path likelihood weights. Read-only observer, Gaussian,
weighted control/Black and exact existing-stock/hedge correspondence tests pass.
Baseline and augmented ridge use the same penalty; legacy OLS is separate.

On two fresh 512-path trials, central direct-512 variance gains are
**1.8652/1.8421**, but joint-grid runtime multipliers **1.7882/1.7707** leave only
**1.0431/1.0404** measured cost efficiencies. Correction efficiencies remain
**0.62–0.68**. All price estimators/grids pass **0/24** positive-hit sampling
screens. These do not clear calibration or establish a robust speedup.

**Repeated wing warnings remain:** at 512 steps, augmented-minus-base 29-day
7900-call shifts are +0.14339 ± 0.03930 and +0.16172 ± 0.03926 (3.65/4.12 paired
SEs); 7950-call shifts are +0.11384 ± 0.03299 and +0.11893 ± 0.03386
(3.45/3.51). Do not dismiss same-sign replication by treating related strikes
and grids as independent tests. All 64 bare innovation/forward checks pass a
two-SE screen, but **126/1680 correlated contract-gain mean flags** remain,
including reference-payoff, stock-delta and OU-linear gains. Fine-root covariance
composition differs from target by only 2.78e-17 (relative diagonal 3.84e-16);
no compensator mismatch is identified, but the warning cause is not settled.

Negative adjusted estimators are preserved, including baseline 29-day 3400 put
-3.69264 ± 6.62777 at 128 steps and augmented 57-day 3450 put
-2.12469 ± 2.59609 at 256 steps (trial 0). Half-pilot SE inflation reaches
9.42 for a correction and 6.41 for a direct price. Augmented 256->512 changes
exceed scale/8 in 16/24 cases; the 57-day 6775 put changes
-1.17242 ± 0.45488, and a wing warning has just one nonzero correction.
No clipping, recentering or variance/bias conflation occurs.

All **146 tests** pass; 960 estimates, 120 paired differences, 64 bare monitors
and 1680 contract-control means replay. A post-hoc two-fold analysis of OLD
training paths finds three linear controls retain most nominal direct-price
benefit, but it is not fresh confirmation or a measured runtime saving.
See [SPX_INNOVATION_CONTROL_RESULTS.md](SPX_INNOVATION_CONTROL_RESULTS.md).
The full-nine dictionary stays experimental; the next candidate is a cheaper
linear-only kernel with fresh tail/control-mean and standalone-cost checks.
All VIX maturities and original scopes remain active, with no market refit.


## Linear-only direct512 check (22 October 2026)

[SPX_LINEAR_INNOVATION_RESULTS.md](SPX_LINEAR_INNOVATION_RESULTS.md) records the
implemented three-control branch and a fresh 512-path pilot plus two 1024-path
checks. Central variance gains 1.7953/1.7105 outweigh measured standalone pipeline
cost multipliers 1.5718/1.5030 by only 1.1422/1.1381. It remains experimental:
0/24 primary quote-scale screens pass, a far put has a negative adjusted estimate
and one path contributes 93.99% of squared paired differences, half-pilot SE
inflation reaches 30.52, and 19/192 contract-control mean flags remain. The old
full-nine far-call warning is not cleared by this different pilot/dictionary.
All 147 tests pass; saved coefficients, weights, statistics and source hashes
replay. Baseline remains reference; all four VIX maturities stay active.

## Moneyness and positivity follow-up (22 October 2026)

[SPX_MONEYNESS_POSITIVITY_AUDIT.md](SPX_MONEYNESS_POSITIVITY_AUDIT.md) checks saved
linear-control results. All 12 representatives are technically OTM; the four
about 2% from the forward are positive, with only 0.993–1.026x half-pilot SE
ratios. All eight negative estimates across pooled/half-pilot methods are puts
at K/F about 0.493–0.499. Near-ATM SE/price is 0.95–1.43%, but quote-scale screens
still fail. These facts support a central research focus, not convergence or
ITM/full-surface claims. Negative controlled estimates are unusable as quotes
and retained as signed diagnostics; no clipping or checking-selected fallback.
Raw positive-weight payoffs are nonnegative by construction but can be noisy.
All maturity/quote scopes are retained. Reproduction inputs and hashes are saved
in `runs/m2_spx_moneyness_audit_001`.

## User-confirmed near-ATM priority (22 October 2026)

Prioritize the four existing near-ATM SPX training representatives for numerical
stability/accuracy. Record exact contracts and boundaries in `active_scope.json`
and the current section of `NEXT_STEPS.md`. Wings stay in diagnostic records and
all four VIX maturities remain active. This is a work-priority decision, not a
quote-scope exclusion, sampling/bias clearance or model refit. No simulation or
estimator changes were made in this update.

## Maturity-scaled SPX central band

[SPX_SCALED_MONEYNESS_COVERAGE.md](SPX_SCALED_MONEYNESS_COVERAGE.md) implements
training-only market-IV-scaled log-moneyness. The z<=2 priority band includes
566/980 SPX training quotes (57.8%), over half at each of 29/57/211/302 days.
The 648 count used for short-maturity studies refers only to 29/57 days, not
the full training surface; this corrects earlier denominator shorthand.
Coverage does not establish pricing accuracy. Four old central contracts remain
representatives, all wing diagnostics and all VIX maturities remain active.

## User-selected provisional SPX cutoff 3

User selected z<=3 for the maturity-scaled central working band, superseding
cutoff 2. This includes **707/980 training quotes (72.1%)**. Current membership
is `runs/spx_scaled_moneyness_003/training_membership.csv`, `priority_selected`;
old cutoff-2 flags remain explicitly historical. All quote/maturity scopes and
wing warnings remain recorded; no new pricing accuracy claim follows.

## Cutoff-3 representative sampling panel completed

[SPX_CENTRAL_PANEL_RESULTS.md](SPX_CENTRAL_PANEL_RESULTS.md) freezes 24 contracts
(put/call at z near .5/1.5/2.5 at all 29/57/211/302-day maturities), then runs
128-path independent pilots and two fresh 256-path checks per expiry at 512 steps/day.
All 384 estimates are nonnegative; 46/48 primary variance comparisons improve,
but all eight full-panel variance-cost comparisons lose (0.823–0.989 gain).
All primary methods still pass 0/48 quote-relative sampling screens. Long outer
calls have only one pilot payoff and 1–4 checking payoffs; residual concentration
can exceed 99%. Half-pilot SE inflation reaches 4.20, with 3.09 also in the inner
302-day band. This limits the earlier near-ATM stability observation to its
original contracts/seeds/larger pilot; it is not a contradiction or universal
stability clearance. Seven contract-control mean flags, three bare-monitor flags
and two paired mean shifts are retained. The 211-day stock mean flags below one.
The baseline remains reference and the linear dictionary stays experimental.
No calibration, new scope pruning or positivity clipping occurs. All four VIX
maturities remain active; all 707 cutoff-3 quotes remain the priority population.

## Upward-equity sampling and pilot-size follow-up (23 October 2026)

[SPX_CALL_TILT_RESULTS.md](SPX_CALL_TILT_RESULTS.md) records a pilot-only component
analysis, then fresh original/+3 equity-tilt comparisons on the same 24 contracts.
Each proposal gets 256-path pilots and two 256-path checks at 512 steps/day. All 480 estimates
are nonnegative and all eight outer-call checking hit counts increase. However,
only 10/48 primary variances improve; aggregate efficiency wins only 4/8 trials.
Identical old128 coefficients isolate proposal effects and also show no robust
broad-panel win. Larger-pilot learning is mixed (one new fit has 49.85x old-fit
variance). All primary sampling screens fail. Preserve 30 correlated contract-mean
flags, 18/176 bare flags and paired/learning warnings; no weight self-normalization.
The original sampler/ridge5 stays reference; call3 remains experimental.
All 707 cutoff-3 quotes and all four VIX maturities remain in scope. Correctly
validated cached GH48 setup takes 1 second, avoiding recomputation without convergence
clearance. Next candidate is a source-faithful conditional-SPX Gaussian-noise
integration study, starting with mathematical/rank/covariance checks. It is
recorded as pending, not an implemented positivity or pricing fix.
