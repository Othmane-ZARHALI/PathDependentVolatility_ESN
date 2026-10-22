# Dated pricing implementation status — 22 October 2026

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
