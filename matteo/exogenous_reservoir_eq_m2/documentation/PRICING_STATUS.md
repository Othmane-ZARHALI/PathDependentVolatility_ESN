# Dated pricing implementation status — 23 October 2026

This is an implementation update to [CALIBRATION_PLAN.md](CALIBRATION_PLAN.md),
not a change to the model or its historical research claims.

Current user-authorized expansion is COMPLETE:
[ALL_MATURITY_FIT_RESULTS.md](../../model_comparison/option_pricing/ALL_MATURITY_FIT_RESULTS.md).
All seven retained SPX maturities (29/57/119/211/302/666/1030 days) enter actual
recalibration, cutoff 3; exclude only 1/1401/1765/2129. The entire retained
market population is now in sample: 1294 SPX options, 217 VIX options and all
four VIX forwards. Original split labels remain metadata; no market out-of-sample
evaluation is requested or reported. Refit the three variance levels at fixed
selected state/readout/Q and show model curves at all dates. Selected levels
are [0.05762055,0.07134085,0.06321057]; separate numerical loss on the same
market observations improves 351.36→308.64. Long SPX improves substantially,
but only 103/1511 prices lie inside bid/ask and short-SPX/VIX errors persist.
This expanded curve is the current working reference for the next Q-scenario
fit on the full in-sample scope. Source/membership/score/sample-mean/SE audits
pass; an audit-reader raw-VIX shape mistake was corrected without changing
prices or parameters, retaining failed evidence. No jobs remain running.
All 182 Python tests and the thirteen scoped
Lean identities (audit 006) pass; numerical accuracy remains unqualified.

Latest completed checkpoint:
[JOINT_STATE_REFINEMENT_RESULTS.md](../../model_comparison/option_pricing/JOINT_STATE_REFINEMENT_RESULTS.md)
records seven restricted-state profiles and a nested five-piece curve test.
State plus three-piece refitting selects C0=1/Y0=F0=0 and improves pooled loss
355.46→316.93 on three new paired checks, with only 41/882 interval hits.
Five pieces add a small inconsistent 0.64% gain, worsening forwards and interval
coverage; keep the three-piece state-adjusted working reference. Failed runs
and frozen selections are retained; resumed runs 003 do not repeat optimization.
All source/score/membership audits complete, 179 Python tests pass and thirteen
scoped Lean identities pass audit 005. Sampling, stock-grid, normalizer/GH/time
errors remain separate and no qualified market/P fit or ranking is established.
Next bounded fitting hypothesis: declared Q-risk-price scenarios one at a time,
with updated normalization/state mapping/prices and P diagnostics; not launched.
The model-own same-expiry forward improves IV display consistency but changes
no price or forward error. AUDIT.md records and corrects a broader raw-parity
diagnostic's exposure of 35 reserved VIX strike keys; no held-out model scoring
or parameter selection used them. Future evaluation on those keys is quote-exposed.

Preceding curve-only fitting addition:
[JOINT_CURVE_FIT_RESULTS.md](../../model_comparison/option_pricing/JOINT_CURVE_FIT_RESULTS.md)
records actual three-level curve optimization on all 707 cutoff-3 SPX options,
all 175 VIX options and four forwards. All three fresh numerical comparisons
improve; pooled loss 597.90→347.10, but VIX fit worsens in aggregate and only
50/882 objective options lie inside bid/ask. No qualified fit or model ranking.
Fast fixed-shock stock repricing and analytical derivatives are checked against
the existing scheme. An explicit signed exploratory loss retains negatives;
production validation stays strict. There are 170 passing Python tests and
ten checked Lean identities with clean compiler/axiom audits. See the
[formal scope](../../model_comparison/option_pricing/formalization/CORRESPONDENCE.md).
Next test restricted state while updating dated normalization; current curve-only
exposures must not be reused unchanged. Numerical accuracy remains unresolved.

Preceding bounded addition: `conditional_controls.py` analytically conditions the
existing stock, deterministic-reference payoff and two adapted stock-increment
controls, with valid importance sampling on retained Gaussian coordinates.
Seven new Gaussian/Black/forward/control tests pass (124 package tests total;
42 study-driver tests also pass). The fresh market comparison loses whole-panel
variance-cost efficiency in both trials, 0.7483/0.8541. Retain original
vol7/ridge5. Tiny conditional positive values are not raw exercise hits, and
the two nominal SE passes remain tail-unresolved. No pricing or fit clearance.
See [conditional-control results](../../model_comparison/option_pricing/SPX_CONDITIONAL_CONTROL_RESULTS.md)
and [fitting reentry decision](../../model_comparison/option_pricing/FITTING_REENTRY_DECISION.md)
for the derivation, retained failures and separate accuracy gates. The previously
specified local profile was superseded by the actual joint fit above. Older
numerical test counts below are historical.

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
- `src/esn_eq/path_importance.py` adds an iid defensive mixture of full equity
  Brownian path shifts, including the correlated OU innovations and exact
  bounded target/proposal likelihood weights. Weighted controls use independent
  pilots. Gaussian joint-moment and Black tail oracles pass. See
  [SPX_PATH_IMPORTANCE_RESULTS.md](../../model_comparison/option_pricing/SPX_PATH_IMPORTANCE_RESULTS.md)
  for the bounded market diagnostic and separate grid/normalizer checks;
  these do not qualify the SPX surface or change VIX maturity scope.
- `src/esn_eq/volatility_path_importance.py` builds orthonormal full-path scores
  for equity and conditional OU projection excursions. Exact joint likelihoods
  and independent Gaussian/Black oracles cover asymmetric volatility shifts.
  The [volatility-sampling diagnostic](../../model_comparison/option_pricing/SPX_VOLATILITY_SAMPLING_RESULTS.md)
  separates proposal selection from fresh checking and retains pilot-half
  instability, negative estimates and all candidate results.
  Its selected 7-day proposal improves fresh observed variance in 20/24 cases
  but clears no quote-relative sampling screen; coupled time errors remain.
  Additional weighted Gaussian-score/ridge controls are frozen before checking
  results and give only a modest further gain. No market-fit clearance follows.
- `src/esn_eq/mlmc.py` supplies independent-level mean/SE aggregation and a
  minimum-constrained, block-rounded variance/cost allocation. The
  [SPX MLMC prototype](../../model_comparison/option_pricing/SPX_MLMC_RESULTS.md)
  targets the 512-step discretized stock price using independent base and
  coupled-correction samples. Oracle tests cover exact constant-volatility
  corrections and Black expectations; cost/accuracy claims require fresh checks.
  The first fresh study is worse than direct vol7 on the central variance-cost
  objective in both trials, retains sparse-tail/expectation warnings and does
  not clear 256->512 timestep changes. The implementation remains exploratory.
- `src/esn_eq/bridge_coupling.py` implements an endpoint-preserving orthogonal
  Gaussian reflection and streaming original/reflected exact-OU paths.
  Target-law, transported-likelihood, Black and adapted stock-forward checks
  pass. The [bridge study](../../model_comparison/option_pricing/SPX_BRIDGE_COUPLING_RESULTS.md)
  finds modest correction variance reductions that do not pay for measured
  runtime. It remains a comparator; no price/timestep clearance follows.
- `src/esn_eq/innovation_controls.py` records exact target-Q projected
  innovations from full-state conditional means on nested grids. Nine new
  predictable linear/square/equity-cross gains complement the stock hedges.
  Gaussian covariance, weighted martingale, Black and backward-correspondence
  tests pass. The [innovation study](../../model_comparison/option_pricing/SPX_INNOVATION_CONTROL_RESULTS.md)
  separates dictionary/regularization effects and charges measured recording
  and control runtime; variance reduction does not establish a bias clearance.
  Fresh central direct variance improves around 1.8–1.9, with nearly matching
  runtime overhead, while correction efficiencies worsen. Repeated far-wing
  paired mean/control-gain warnings and pilot sensitivity remain recorded;
  the full dictionary is experimental and clears no quote-relative gate.

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

The proposed last-VIX-maturity deferral is examined in
[VIX_SCOPE_DECISION.md](../../model_comparison/option_pricing/VIX_SCOPE_DECISION.md).
Dropping the 153-day options and forward from the fixed-candidate comparison
does not eliminate its fitting/checking reversal; the 118-day block also
matters. That checkpoint retained four expiries; the user's subsequent explicit
instruction initially deferred both 118/153-day options and forwards. That historical
exploratory VIX scope was 27/55 days, while SPX scope stays unchanged.
`active_scope.json` now records the subsequent restoration of all four maturities and forwards.

The historical near-only refit is in
[VIX_NEAR_RESULTS.md](../../model_comparison/option_pricing/VIX_NEAR_RESULTS.md).
Only one variance level is identifiable from the retained windows in the
declared curve family. Fresh 8x128-state checks show consistent point rankings
but only a 0.0352 loss-unit gain, smaller than its descriptive uncertainty.
Sampling and quadrature gates remain open. There are 84 package tests and
11 legacy/driver checks passing. No deferred-maturity or held-out performance
claim is made.

The user's latest same-day instruction restores **all four VIX maturities and
forwards** for new analyses (175 options); no expiry is currently deferred.
Historical near-only evidence is preserved. Numerical gates remain open.

The first full-scope defensive-mixture importance-sampling pilot is in
[VIX_IMPORTANCE_SAMPLING_RESULTS.md](../../model_comparison/option_pricing/VIX_IMPORTANCE_SAMPLING_RESULTS.md).
At 8x64 states per expiry it removes 17 zero-hit option estimates, but precision
changes are mixed: many 55/118-day estimates improve, while 27-day and 153-day
forward uncertainty worsens. Single-observation tail concentration persists;
no accuracy gate is cleared and no refit is performed. Current checks pass:
87 package tests and 11 legacy/driver tests (98 total).

The pilot-designed full-scope follow-up is in
[VIX_PILOT_SAMPLING_RESULTS.md](../../model_comparison/option_pricing/VIX_PILOT_SAMPLING_RESULTS.md).
Fresh 8x128-point checks lower median option SE at all four maturities and remove
15 ordinary-Sobol zero-hit option estimates. Raw forward SE increases at every
expiry; simple independent-pilot mass centering gives mixed results and ten
negative option means, retained as failures without clipping. Only 27/175 raw
mixture estimates meet the observed sampling-SE threshold; no full gate passes.
The original selected curve again has higher fresh loss than the starting curve:
+2.655 under ordinary Sobol and +2.884 under the learned mixture, with unresolved
uncertainty. No refit, final ranking, or held-out claim follows. All 101 tests
pass; source/weight/summary correspondence checks pass separately.

The frozen-proposal allocation follow-up is in
[VIX_STRATIFIED_SAMPLING_RESULTS.md](../../model_comparison/option_pricing/VIX_STRATIFIED_SAMPLING_RESULTS.md).
All four maturities remain included. Fixed allocation lowers forward SE at
55/153 days but worsens it at 118 days; option precision is mixed. Both methods
have no zero-hit option estimates, but no full accuracy gate passes. The fresh
fixed-curve loss changes are -0.210 (random mixture) and +0.933 (stratified),
with descriptive uncertainty around 1.1–1.4: the original ranking remains
unresolved. A 55-day conditional tail-state GH change is 0.0223 VIX points,
larger than the preceding pilot-run checks. A new 27-day second-moment identity
discrepancy persists after removing random mass fluctuation and warrants a
separate outer/conditional/normalizer/time audit before refitting. It is not
proof of a normalization bug. Current validation: 106 tests and independent
artifact correspondence checks pass; statistical/numerical gates remain open.

The follow-up identity attribution is in
[VIX_MOMENT_DIAGNOSIS.md](../../model_comparison/option_pricing/VIX_MOMENT_DIAGNOSIS.md).
The independent, variance-inflated direct Gaussian calculation agrees with all
four moment targets (27d: 542.49136 ± 0.09509 versus 542.49611). Gaussian-law
composition agrees to floating-point precision. On all 1024 original 27d
states, conditional32→48 plus normalizer48→64 changes the moment by only
-0.01643 ± 0.01193, leaving the +13.65 residual. Limited higher-order/time
checks are small; the independent population-composition estimator remains
noisy and inconclusive. Outer-sample/tail sensitivity is the leading explanation,
not a demonstrated normalization bug. At one representative future date per
maturity, roughly 0.056%–0.072% of Gaussian states contribute about 49%–53%
of mean variance multiplier, showing why sampling is difficult. All four
maturities/forwards remain in scope; no prices are rescaled or refitted. Current
validation: 113 tests plus source/statistics replay checks pass. Full price,
quadrature and calibration gates remain open.

The deterministic score-subspace sampling follow-up is in
[VIX_SCORE_SAMPLING_RESULTS.md](../../model_comparison/option_pricing/VIX_SCORE_SAMPLING_RESULTS.md).
The new `score_sampling.py` retains the full terminal state and uses exact
bounded defensive-mixture weights. On fresh R16 x 64 equal-budget checking
samples at all four maturities, the proposed rank-four covariance inflation
worsens median option SE by factors 2.24/1.86/1.36/1.81 and worsens every forward
SE versus the frozen learned mixture. It is not adopted. Fresh learned-mixture
moments lie within 1.1 reported SE of all targets, supporting sampling sensitivity
of the preceding warning, without establishing option convergence. Observed
option SE screens pass only 28/175 (learned) and 11/175 (inflated); count and
quadrature gates remain open. All four maturities and forwards are retained;
no refit or rescaling is performed. Validation: 116 tests and replay of 8192
full states/weights plus 732 summary rows pass. Next use independent design
pilots and fresh count refinements; no joint-fit ranking is accepted.

The payoff-targeted independent-pilot follow-up is in
[VIX_TARGETED_PILOT_RESULTS.md](../../model_comparison/option_pricing/VIX_TARGETED_PILOT_RESULTS.md).
Four disjoint pilot fits favor mean shifts in score-sensitive directions on
historical design screening. Fresh R8 x 64 checks improve median option SE at
55/153 days but worsen 27/118 days and all four forward SEs. Only 24/175
targeted versus 19/175 baseline observed SE screens pass. A new targeted
118-day moment residual of -11.80491 ± 4.30771 is explicitly retained; the
preceding fresh moment agreement did not clear convergence. This design is
experimental, not adopted for calibration. All four maturities/forwards remain
included; no model refit or price rescaling. Validation: 119 tests plus pilot
and full-state/weight/statistics replay checks pass. Next isolate component
allocation effects and refine the offending moment on identical states.

The saved-state attribution and fresh allocation checkpoint is in
[VIX_ALLOCATION_REFINEMENT_RESULTS.md](../../model_comparison/option_pricing/VIX_ALLOCATION_REFINEMENT_RESULTS.md).
All 512 offending 118-day states retain their original weights. Raising
conditional GH32/M48 to GH48/M64 changes their moment by only -0.00603 ± 0.02206,
leaving -11.81094; tested refinements do not explain that sample's warning.
Fresh R12 x 128 random/fixed sampling of the identical 17-component density
gives all moment residuals within about 1.7 reported SE, but only 36/175 versus
35/175 option SE screens and no forward SE gate pass. Fixed allocation is mixed,
not adopted universally; count/price-bias gates remain open. All four maturities
and forwards remain active. Validation: 121 tests, exact-state node/statistics
replay, and 12288-state/weight/allocation replay pass. Next test known-mean
Gaussian controls with independent pilot coefficients, preserving all bias gates.

The Gaussian-control checkpoint is in
[VIX_GAUSSIAN_CONTROL_RESULTS.md](../../model_comparison/option_pricing/VIX_GAUSSIAN_CONTROL_RESULTS.md).
`gaussian_controls.py` supplies exact-mean Gaussian linear/quadratic controls
and replication-level ridge fitting. Independent pilot coefficients are frozen
before fresh R16 x 128 checks; raw and controlled estimates share paths. Primary
forward SE decreases at all four maturities and 27/55-day SE screens pass, but
118/153 days still fail; only 39/175 option SE screens pass (37/175 raw).
Pilot halves show sensitivity, including a retained half-B 153-day moment
warning +5.61841 ± 2.69630. No identity is imposed or negative estimate clipped.
All four maturities/forwards remain included; no calibration or bias clearance.
Validation: 124 tests and replay of pilot fits, 8192 states/control rows and
1830 summary rows pass. Next assess independent-pilot and count stability,
prioritizing long forwards and tails while keeping the full scope.

The larger-pilot and 128/256 outer-count checkpoint is in
[VIX_GAUSSIAN_STABILITY_RESULTS.md](../../model_comparison/option_pricing/VIX_GAUSSIAN_STABILITY_RESULTS.md).
The same nine controls/ridge use 28 compatible pilot replications, with two
disjoint balanced 14-pilot fits and the original 12-pilot matrix retained.
Fresh R12 x 256 pooled forward SEs pass the observed screen at all four
maturities, but the 153-day pass is borderline and not supported by every
pilot fit. Only 45/175 pooled option screens pass (47/175 raw); count and
bias gates remain unresolved. New pooled/balanced-A 118-day moment warnings
are retained (-4.55861 ± 2.24254 and -4.21342 ± 1.98663). No refit, exclusion,
rescaling or clipping occurs. Validation: 126 tests plus replay of 12288 states,
controls, all pilot fits and 2196 summary rows pass. Next screen payoff-shaped
Gaussian tail controls on independent design pilots, preserving the full scope.

The user has asked to record the VIX tail-control extension and move forward
for now, retaining all maturities. The pending extension and next SPX diagnostic
work are recorded in [NEXT_STEPS.md](../../model_comparison/option_pricing/NEXT_STEPS.md).
All four VIX maturities, 175 training options and four forwards remain active;
`active_scope.json` is unchanged. The tail-control pilot is pending, not launched,
and the working interpretation is to resume SPX numerical diagnostics. No VIX
convergence, scope reduction or new market fit is implied.

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


The 22 October linear-only SPX follow-up reports about 1.14x central measured
variance-cost efficiency in two fresh direct512 trials. It remains experimental:
0/24 quote-relative sampling screens pass, a negative far-put estimate is driven
by a single high-leverage path, and pilot sensitivity remains severe (30.52x SE).
See `matteo/model_comparison/option_pricing/SPX_LINEAR_INNOVATION_RESULTS.md`.
No model refit or continuous-time clearance follows; all VIX maturities remain.

The cutoff-3 SPX representative panel now checks 24 contracts across 29/57/211/
302 days at 512 steps/day using separate 128-path pilots and two 256-path checking trials.
All 384 estimates (including pilot-half variants) are nonnegative, but all
primary quote-relative sampling screens fail and linear-control full-panel
variance-cost gains are below 1 in all eight checks. Sparse long outer-call
pilots/checks and up to 99.6% residual concentration remain. This does not clear
stock/normalizer bias or market refitting. See
`matteo/model_comparison/option_pricing/SPX_CENTRAL_PANEL_RESULTS.md`.
All four VIX maturities remain active.

On 23 October, the fixed +3 upward-equity sampling comparison increases outer-call
coverage at all four SPX training maturities but does not show a robust general
variance-cost win. All 480 estimates are nonnegative; all primary precision
screens still fail. Original vol7/ridge5 remains reference, cutoff 3 remains
selected and all VIX maturities stay active. See
`matteo/model_comparison/option_pricing/SPX_CALL_TILT_RESULTS.md` for fixed-beta,
pilot-learning, likelihood/martingale checks and all warnings. A conditional-SPX
Gaussian integration route is recorded as a pending structural investigation;
it has not yet been implemented or validated.

A partial conditional-SPX prototype now integrates the equity Brownian direction
independent of all primitive OU drivers. Frozen nu²=0.00636685 limits its raw
benefit: target-Q paired checks show median variance gain 1.025, with no precision
clearance or replacement of original vol7/ridge5. Conditional raw payoffs are
nonnegative, but six zero-hit cases demonstrate persistent rare-volatility-path
sampling uncertainty. Gaussian/source/covariance and Black oracles pass; an
endpoint-rounding oracle failure was corrected and recorded. All 159 tests pass.
See `matteo/model_comparison/option_pricing/SPX_CONDITIONAL_PRICING_RESULTS.md`.
A conditional-control/retained-shock proposal extension remains pending. All
SPX priority quotes and VIX maturities stay included; bias gates stay separate.
