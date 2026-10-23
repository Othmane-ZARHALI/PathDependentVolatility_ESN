# Pricing and fitting research ledger

Checkpoint: 23 October 2026, branch `codex/m2-option-pricing-benchmark`.
Latest numerical implementation at this checkpoint: `afbf628`.

This is the entry point for the accumulated attempts and their implications for
fitting. Detailed reports and their immutable run directories remain the evidence;
this index does not replace them. Older reports' “next” sections and test counts
describe their historical checkpoints. Use the current queue below when resuming.
The copy-and-paste continuation prompt is [CONTINUE_FITTING_PROMPT.md](CONTINUE_FITTING_PROMPT.md).

## Goal and present conclusion

Fit the 19 February 2026 SPX/VIX snapshot with M2 under the audited common price
objective, eventually comparing fairly with Othmane's preserved ESN engine and
checking physical compatibility. **No reliable joint market fit, final model
ranking or continuous-model pricing convergence has been established.**

The initial exploratory VIX curve fit improved its fixed simulation sample but
reversed on fresh simulations. Subsequent work identifies rare outer/volatility
paths as a major numerical obstacle. Several variance-reduction ideas help
specific diagnostics, but none clears the complete pricing surface. Avoid
mistaking this numerical bottleneck for demonstrated model infeasibility.

## Decisions to preserve

- [active_scope.json](active_scope.json) is authoritative: VIX 27/55/118/153 days,
  all 175 training options and four forwards. An earlier two-maturity experiment
  was explicitly superseded. **Keep all four maturities.**
- SPX priority is provisional maturity-scaled cutoff 3:
  `abs(log(K/F(T))) <= 3*sigma_ATM_proxy(T)*sqrt(T)`, using training-only market
  nearest-ATM IV anchors. This covers **707/980 training quotes (72.1%)** across
  29/57/211/302 days. The membership artifact is
  `runs/spx_scaled_moneyness_003/training_membership.csv`, column `priority_selected`.
  Four old near-ATM contracts and the newer 24-contract panel are representatives,
  not the whole priority population. Wings remain secondary diagnostics.
- Diagnostic baseline: frozen Stage-3 candidate 18/scenario 0, zero projected
  initial state, selected nuisance curve, and original vol7/ridge5 SPX estimator.
  Read actual run configurations before reusing a curve or cache; old initial
  stripped-curve and later selected-curve results are different experiments.
- Recent SPX experiments use 512 steps per trading day. Exact OU transitions do
  not make the nonlinear stock update exact; **512 is not convergence clearance**.
  The old 2/8-step runs are coarse historical diagnostics.
- Preserve signed negative controlled estimates as failures; never clip them into
  acceptable quotes. Positive raw/conditional prices can still be badly sampled.
- Preserve training/validation/final-test boundaries. A fresh simulation is a
  numerical check, not a new held-out market observation.
- Record code/report discrepancies in [AUDIT.md](AUDIT.md); retain unsuccessful
  attempts, source snapshots, seeds, weights, pilot fits and timing.

## Foundations and attempts to fit

| Attempt / record | What was learned | Status for continuation |
| --- | --- | --- |
| [Legacy ESN source/report review](AUDIT.md), [initial benchmark](RESULTS.md), [Othmane replication](../othmane_replication/RESULTS.md) | Actual executable targets, objective and conventions need to be distinguished from report descriptions; original final artifacts are incomplete. | Preserved comparison evidence; no certified common-protocol fit ranking. |
| [Dated M2 SPX engine](M2_SPX_RESULTS.md), [conditional VIX engines](M2_VIX_RESULTS.md) | SPX simulation, full-state conditional VIX Gaussian quadrature and nested-MC checks implemented. Zero-hit payoffs invalidate naive SE-only tail clearance. | Working foundations; accuracy gates remain separate. |
| [Fine stock grids](FINE_GRID_RESULTS.md) | Requested 64/128 and subsequent 128/256 comparisons implemented; fine grids alone do not resolve sampling. | Historical refinement evidence, not global convergence. |
| [Common loss / Gaussian time integration](CALIBRATION_PROGRESS.md) | Frozen-scale bid/ask interval loss and optimizer oracle implemented. Efficient conditional-VIX time integration checked at selected states. | Reuse loss implementation; restricted anchor scores are not full fits. |
| [Outer sampling and curve initializer](SAMPLING_PROGRESS.md) | iid/Sobol comparison and a necessary-moment-constrained nuisance-curve starting point. | Initializer is not a joint-fit solution; curve identification is limited. |
| [Exploratory cached curve fit](CURVE_PROFILE_RESULTS.md) | Apparent VIX improvement on the fitting simulations reverses on fresh scrambles. | Failed validation of improvement; do not adopt as a proven fit. |
| [Reversal and SPX controls diagnosis](NUMERICAL_STABILITY_RESULTS.md) | At the dominant 153-day VIX block, outer-sample changes dwarf tested quadrature changes; SPX hedging controls help but remain insufficient. | Main motivation for numerical work before fit ranking. |
| [Scope experiment](VIX_SCOPE_DECISION.md), [near-only fit](VIX_NEAR_RESULTS.md) | Removing the last maturity does not remove reversal; near-only curve improvement is small relative to uncertainty. | Historical only; full VIX scope subsequently restored. |

## VIX sampling and moment investigation

Each row links the detailed design, run artifacts and limitations. No row clears
all option-price sampling and bias requirements.

| Attempt / record | Result and disposition |
| --- | --- |
| [Defensive importance mixture](VIX_IMPORTANCE_SAMPLING_RESULTS.md) | Removes zero-hit cases at the tested budget, but uncertainty changes are mixed and single-state concentration persists. |
| [Pilot-learned proposals](VIX_PILOT_SAMPLING_RESULTS.md) | Median option SE improves at all expiries; forward uncertainty worsens. Mass-centering can produce negative prices. No broad clearance. |
| [Fixed component allocation](VIX_STRATIFIED_SAMPLING_RESULTS.md) | Mixed precision gains; moment warning motivates separate attribution. No universal adoption. |
| [Second-moment / rare-tail attribution](VIX_MOMENT_DIAGNOSIS.md) | Independent direct-Gaussian moment checks agree with target; tested conditional/normalizer refinements do not explain the original residual. Rare outer tails are the leading explanation, not a proved normalizer bug. |
| [Score-subspace covariance inflation](VIX_SCORE_SAMPLING_RESULTS.md) | Worsens median option and all forward SEs against the learned-mixture comparator. Not adopted. |
| [Payoff-targeted pilots](VIX_TARGETED_PILOT_RESULTS.md) | Some option gains, but other maturities and all forwards worsen; new 118-day moment warning retained. Experimental. |
| [Matched-state refinement / allocation](VIX_ALLOCATION_REFINEMENT_RESULTS.md) | Refining the same offending states barely changes their warning. Fresh allocation checks are mixed. No rescaling of prices or identity enforcement. |
| [Known-mean Gaussian controls](VIX_GAUSSIAN_CONTROL_RESULTS.md) | Forward SE improves at all expiries; long forwards and most option screens remain unresolved. |
| [Larger pilots / count stability](VIX_GAUSSIAN_STABILITY_RESULTS.md) | Pooled forward SE screens pass at all four maturities, but 153-day pass is borderline/pilot-sensitive. Only 45/175 pooled option screens pass; 118-day moment warnings remain. Latest VIX reference checkpoint, not surface clearance. |

The normalized model implies
`E[VIX_T^2] = (10000/delta) * integral_T^(T+delta) xi(u) du`.
Agreement for this identity does not establish VIX forward/option convergence.
The recorded pending VIX extension is payoff-shaped known-mean Gaussian tail
controls, designed on independent pilots and checked at all four maturities.
It has not been launched; the user then prioritized SPX work.

## SPX sampling, coupling and conditional-pricing attempts

| Attempt / record | Result and disposition |
| --- | --- |
| [Equity-path importance sampling](SPX_PATH_IMPORTANCE_RESULTS.md) | Better tail exposure, insufficient precision. Tested same-path GH48→64 changes are small, without a population normalizer-error bound. |
| [Volatility-excursion sampling / score controls](SPX_VOLATILITY_SAMPLING_RESULTS.md) | vol7 improves fresh variance in many cases; extra Gaussian-score controls offer modest benefit. Original vol7 with existing five controls becomes the working reference, not a universal optimum. |
| [MLMC, 64/128/256/512 levels](SPX_MLMC_RESULTS.md) | Loses central measured variance-cost comparison in both fresh trials. Top-level tails and timestep warnings remain. Keep as comparator, not replacement. |
| [Endpoint-preserving bridge reflection](SPX_BRIDGE_COUPLING_RESULTS.md) | Target-law/oracle checks pass, but all four central variance-cost comparisons lose (0.67–0.77 gain). Not adopted. |
| [Nine OU-innovation controls](SPX_INNOVATION_CONTROL_RESULTS.md) | Central variance gain roughly 1.8–1.9 is nearly consumed by runtime; corrections lose on cost and control/tail warnings persist. Experimental. |
| [Three linear innovation controls](SPX_LINEAR_INNOVATION_RESULTS.md) | About 1.14× central measured efficiency in two trials, but far-put negative price and severe pilot sensitivity prevent general adoption. |
| [Moneyness/positivity audit](SPX_MONEYNESS_POSITIVITY_AUDIT.md) | Saved negative estimates are deep-OTM puts; original near-ATM representatives are positive and more stable. This does not generalize to the full central panel. |
| [Maturity-scaled coverage](SPX_SCALED_MONEYNESS_COVERAGE.md) | Training-only coverage checked for multiple cutoffs; user selects 3, covering 707/980 quotes. Membership policy, not a pricing-accuracy result. |
| [24-contract cutoff-3 panel](SPX_CENTRAL_PANEL_RESULTS.md) | Two fresh checks at all four SPX expiries: all 384 estimates nonnegative, but 0/48 primary sampling screens pass. Linear controls lose full-panel efficiency in all eight checks; long outer-call pilots are sparse. |
| [Stronger upward-equity tilt / larger pilots](SPX_CALL_TILT_RESULTS.md) | All outer-call hit counts increase, but only 10/48 primary variances improve and efficiency wins in 4/8 checks. No robust broad-panel replacement for original vol7/ridge5. |
| [Partial conditional Monte Carlo](SPX_CONDITIONAL_PRICING_RESULTS.md) | Positive raw conditional payoff construction works. Only 0.636685% of equity Brownian variance is independent of OU drivers. Median observed variance gain is 1.0253× **against uncontrolled target-Q MC**, not against vol7/ridge5. Six raw zero-hit cases still expose rare retained-path uncertainty; 0/48 conservative sampling screens pass. Experimental benchmark. |

Latest conditional run: `runs/m2_spx_conditional_001/`; 96 prices and 48 paired
differences replay, with covariance/source audits. At `afbf628`, 117 package +
42 study-driver tests pass. These are implementation checks, not numerical
convergence certificates. Small independent equity variance is not an upper
bound on option/control variance gains from conditioning.

## Return to fitting: ordered queue

1. **Bounded conditional-control feasibility check.** Derive conditional stock,
   reference-payoff and adapted-delta controls under the retained shocks U in
   the latest report. A small raw-payoff gain may coexist with useful controlled
   residual reduction. Conditional hedge gains generally are not zero given U.
   Derive valid importance sampling on retained shocks; do not attach old
   full-path likelihoods to an integrated payoff. Validate law/Black/forward
   oracles first, then freeze pilots and a bounded fresh comparison against
   original vol7/ridge5. Measure actual variance × runtime and tail/pilot risk.
   This extension is **pending**, not implemented. If it fails, record the
   result and retain the reference rather than opening an indefinite feature search.
2. **Assess readiness for a restricted exploratory fit.** Summarize independent
   sampling precision, coupled 256/512 stock corrections (refine as needed),
   normalizer sensitivity and VIX conditional/time quadrature separately.
   Use declared budgets and keep unresolved candidates labelled unresolved.
   Do not repeat finished experiments without a concrete new question.
3. **Resume curve/state/Q fitting under the declared model.** Start with the
   existing common loss and restricted curve/state families; distinguish input
   compatibility, weak identification and inadequate model fit from numerical
   noise. Use common randomness for search comparisons and genuinely fresh
   numerical checks for finalists. No proven fit improvement from search loss
   alone. Any explicit provisional subset objective must report its membership
   and full-scope limitations; do not silently redefine the joint objective.
4. **Qualify and compare.** Only after sufficient numerical accuracy, assess
   all scored quotes, physical compatibility and reserved market evaluation.
   Complete the common-protocol ESN comparison before ranking models. This
   phase remains future work, not an achieved deliverable.

The established sampling screen is `2*SE <= scale/4` (equivalently
`SE <= scale/8`); assessed numerical bias is a separate requirement. Zero-hit
empirical SE and single-state dominance cannot establish tail-price precision.
The price objective uses frozen-scale squared distance to bid/ask intervals,
expiry-balanced within SPX and VIX, plus separately scaled VIX-forward errors
and declared penalties. See [CALIBRATION_PROGRESS.md](CALIBRATION_PROGRESS.md)
and `esn_eq/calibration_loss.py`; never normalize market fit by noisy model SEs.

## Navigation and execution

- Overall specification: [CALIBRATION_PLAN.md](../../exogenous_reservoir_eq_m2/documentation/CALIBRATION_PLAN.md).
- Current implementation: [PRICING_STATUS.md](../../exogenous_reservoir_eq_m2/documentation/PRICING_STATUS.md).
- Detailed queue/history: [NEXT_STEPS.md](NEXT_STEPS.md). This ledger's queue
  supersedes older “next hypothesis” wording that later reports completed.
- Reproduction: each report links drivers/verifiers and `runs/` artifacts.
  Use new non-overwriting run directories and source snapshots.
- Environment used: `/private/tmp/pdv-options-venv/bin/python`;
  `OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`;
  package path `matteo/exogenous_reservoir_eq_m2/src` when needed.
  Verify availability in a new session; the temporary environment may disappear.
- A pre-existing unrelated modification to
  `matteo/exogenous_reservoir_eq_m2/experiments/physical_stability_001/RESULTS.md`
  was left untouched at this checkpoint. Check git status before any edits.
