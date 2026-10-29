# Recorded next steps — 23 October 2026

**25 October checkpoint: dedicated stochastic sensitivity trial completed.**
See `STOCHASTIC_VOV_RESULTS.md` and `visualizations/stochastic-vov016`.
The prescribed12-point two-cap grid improves wing IV levels but the selected
added state flattens fresh measured wing slopes under both caps, while body
and SPX fits worsen. All original coefficients/xi/floor were fixed; this was
an exploratory mechanism test, not full recalibration. Keep the selected baseline.
No further experiment is running. A next design should distinguish actual
slope/curvature improvement from lifting the whole smile, and justify its
interaction and cross-market calibration before expensive production searches.
The following older queues remain historical context, not automatic authorization.

For a consolidated attempt index and current ordered return-to-fitting queue,
read [RESEARCH_LEDGER.md](RESEARCH_LEDGER.md). The conversation handoff is
[CONTINUE_FITTING_PROMPT.md](CONTINUE_FITTING_PROMPT.md). Older hypotheses below
remain historical when a later checkpoint completed them.

## Latest user request: recalibrate all retained maturities in sample

Read [ALL_MATURITY_FIT_RESULTS.md](ALL_MATURITY_FIT_RESULTS.md) and its captured
plan. `m2_all_maturity_fit_001` is complete. It fits 1294 cutoff-3 SPX options at
29/57/119/211/302/666/1030 days and all 217 VIX options/four forwards. Exclude
only SPX 1/1401/1765/2129. All old market split labels are metadata in this
feasibility exercise; there is no market out-of-sample result. The user confirms
that actual recalibration is required. Refit the three curve levels at the
selected state, fixed structure/readout/Q, and show previous and recalibrated
parameters on every panel. Numerical checking uses the SAME market targets.
The selected levels are [0.05762055,0.07134085,0.06321057]. Separate numerical
loss improves 351.36→308.64, with substantial long-SPX improvement but only
103/1511 interval hits and persistent short-SPX/VIX shape errors. All four score
records and membership/means/SEs replay. The audit-reader raw-VIX shape error
was fixed without repricing; its failure log is retained. No jobs remain running.
This expanded curve is the new working reference for the Q-scenario next step
below. Keep the full in-sample population; do not revert to four SPX dates or
training-only quotes. Old source labels/data/results are retained.

## Preceding state/five-piece checkpoint; use the expanded baseline above next

[JOINT_STATE_REFINEMENT_RESULTS.md](JOINT_STATE_REFINEMENT_RESULTS.md) records
the complete seven-state and nested five-piece investigations. Runs 003 resume
frozen searches after retained implementation failures, without reoptimization.
State/three-piece refitting improves pooled loss 355.46→316.93 on new shared
draws, with three favorable trials but only 41/882 prices inside bid/ask.
Five pieces improve only to 314.90 (0.64%), with one trial reversal, descriptive
uncertainty exceeding the gain, worse forwards and 38/882 interval hits.
Retain the C0=1/Y0=F0=0 three-piece working reference; do not add more pieces.

Next bounded fit: declared unspanned-risk-price scenarios, feedback pair first,
one at a time, with selected state/readout/architecture fixed and three curve
levels refitted. Rebuild dated normalization/state mapping/stock/VIX exposures,
review controls and sampling-law correspondence, revisit P diagnostics and use
NEW numerical checking streams. This next fit has not been launched. Thirteen
Lean identities and 179 Python tests pass; all SPX precision screens and
deterministic bias qualification remain unresolved.

Model-own same-expiry forward IV is a display diagnostic; it changes no prices
or forward residuals. Historical raw-parity diagnostic exposed 18 validation
and 17 final-test VIX strike keys to a market consistency calculation; no
model scoring/selection used them. Retain the exposure list and label future
evaluation at those keys as quote-exposed. Current diagnostics restrict to
training option keys. See AUDIT.md. No held-out model performance was scored.

## Earlier curve-only checkpoint; its restricted-state step is now complete

The user subsequently asked for actual joint fitting now and Lean checks,
superseding the unlaunched nine-point/24-SPX proposal below.
[JOINT_CURVE_FIT_RESULTS.md](JOINT_CURVE_FIT_RESULTS.md) records optimization
of three curve levels on all 707 cutoff-3 SPX options, 175 VIX options and four
forwards, followed by three fresh numerical checks. Pooled loss improves
597.90→347.10, but only 50/882 options are inside bid/ask and VIX forwards
worsen. Retain a provisional candidate, not a qualified market fit.

Ten Lean identities pass compiler/axiom audits; 170 Python tests pass.
[CORRESPONDENCE.md](formalization/CORRESPONDENCE.md) states exactly what is
proved and what remains outside Lean. Accuracy gates remain separate.

Next test the declared restricted C/Y/F initial-state coordinates, with the
existing prior penalty and Q scenario 0, updating candidate normalization and
state paths. Preserve the priority objective/scales and fresh checking. Do
not reuse the current G/M0/exposure cache after changing state or restart a
standalone variance-reduction feature search. The older curve-profile proposal
below is historical, not a pending task.

The user's subsequent visual review requests undergraduate-level explanations
and a later explicit review of maturity/strike coverage. See the
[complete maturity and variance context](visualizations/joint-curve-context-003/README.md).
The existing four-of-eleven SPX maturity choice is a pilot split, not an
optimized liquidity criterion; primary OTM selection likewise does not verify
executed-volume liquidity. Preserve the current cutoff-3 objective until a new
selection rule is explicitly adopted. The 119-day market curve is now visually
exposed, remains excluded from fitting and has no model score. Additional SPX
expiries still need model repricing before they can support model comparisons.

## Current decision: conditional controls completed, reference retained

[SPX_CONDITIONAL_CONTROL_RESULTS.md](SPX_CONDITIONAL_CONTROL_RESULTS.md) completes
the previously pending extension. Conditional existing controls and a valid
retained-shock proposal pass their Gaussian/Black/control oracles, but lose
whole-panel variance-cost efficiency in both fresh trials (0.7483/0.8541).
Tiny conditional-tail nominal SE passes are not tail clearance. Preserve all
negative reference estimates, frozen pilots, source hashes and failed results.
Keep original vol7/ridge5; do not begin another feature search by default.

[FITTING_REENTRY_DECISION.md](FITTING_REENTRY_DECISION.md) records the separate
accuracy gates and the concrete next fitting experiment: a nine-point local
profile of two curve levels, fixed readout/state/Q, on 24 explicitly declared
SPX representatives plus all 175 VIX options/four forwards. It has not run.
The provisional subset never replaces the 707 cutoff-3 priority population.
No qualified joint fit, feasibility decision or model ranking is ready.

## Historical conditional prototype handoff (extension now completed)

The [partial conditional investigation](SPX_CONDITIONAL_PRICING_RESULTS.md)
now establishes the Gaussian construction and a positive raw-payoff prototype.
Only 0.636685% of frozen-model equity variance is independent of the OU drivers;
the raw target-Q market screen shows a modest median variance gain 1.025 and
no convergence clearance. It does not replace original vol7/ridge5.

Next derive and test conditional expectations of the **existing** stock,
reference-payoff and adapted delta controls. A small raw-payoff gain does not
rule out a larger controlled-residual gain. The conditional gain generally is
not zero when retained shocks include correlated equity; use the Gaussian
expectation identities recorded in the report. Derive a new importance proposal
on retained shocks before pricing with weights; old full-path weights cannot
be reused after integration. Freeze pilots/checking before any market comparison.
This extension is pending, not implemented or validated. Keep cutoff 3/all 707
priority quotes, secondary wings and all four VIX maturities, with independent
stock-timestep and normalizer refinement gates still pending.
Older entries below retain the prior hypothesis and experiment history.

## Latest sampling result and next structural candidate

The [fixed +3 equity-tilt study](SPX_CALL_TILT_RESULTS.md) completes pilot-only
component diagnosis and a fresh, paired comparison with doubled pilot budgets.
Outer-call coverage improves at every maturity, but no robust full-panel
variance-cost gain follows; the original sampler/ridge5 remains reference.
Both proposal and pilot-size comparisons retain failures and mean warnings.
Cutoff 3/all 707 priority quotes and all four VIX maturities remain active.

Next investigate conditional SPX Monte Carlo before adding further mixture
magnitudes or control features: derive the joint-Gaussian equity/OU transition
conditioning in the **full** state, identify which equity noise can be integrated
analytically under the existing discrete stock scheme, and audit numerical rank
and covariance correspondence. Validate Gaussian/constant-volatility and
same-model checks before market experiments. Derive marginal likelihoods if
importance sampling is retained; full-path weights cannot be reused blindly
after integrating noise out. This is a pending candidate, not completed pricing
or a positivity guarantee. Keep stock-timestep/normalizer refinement gates
separate. Older sampling queues below are historical where superseded.

## Latest completed experiment and immediate sampling work

The user-approved [24-contract cutoff-3 panel](SPX_CENTRAL_PANEL_RESULTS.md) now
has independent pilot/checking results at all four SPX training maturities.
No negative estimates were observed, but 0/48 sampling screens pass and the
linear controls lose full-panel variance-cost efficiency in all eight checks.
Long outer-call pilots contain only one positive payoff. Keep the baseline as
reference and cutoff 3/all 707 training quotes as the priority population.

Next inspect saved pilot residual/payoff contributions by side, distance band
and proposal component. Freeze a better-covered pilot and/or horizon-specific
sampling comparison before new checking; no automatic multi-million-path run
is authorized by unstable plug-in budget forecasts. If old checking results
inform a new design, mark them historical/design evidence and use fresh checks.
Keep coupled 256/512 and independent normalizer refinement pending; neither was
performed by this direct 512-step sampling run. All wings/VIX maturities stay recorded.
The earlier queue below is historical where superseded by this checkpoint.

## Current membership rule: maturity-scaled central band

User requested maturity adaptation and meaningful coverage. Prioritize
`abs(log(K/F(T))) <= 3 * sigma_ATM_proxy(T) * sqrt(T)`, with market training-only
nearest-ATM-IV anchors. This supersedes the fixed 5% band for priority membership.
User selected cutoff **3 provisionally**, superseding cutoff 2.
It includes **707/980 training quotes (72.1%)**, over half at every training
expiry (29/57/211/302 days). The four earlier near-ATM contracts remain diagnostic
representatives, not the entire priority set. All quote and VIX scopes remain.
See [coverage and denominator correction](SPX_SCALED_MONEYNESS_COVERAGE.md).

Next broaden representative strike checks within the central band, keeping
sampling, stock timestep and normalizer gates separate. Longer SPX maturities
need their own convergence/cost checks; no new pricing accuracy is established
by the coverage count. Older priority text below is historical where superseded.

## Current priority: near-ATM SPX (user confirmed)

User: **"Yes, I would also focus more on near-ATM options for now."**
This supersedes wing-control stabilization as the immediate research priority.
It does not remove quotes, maturities or historical warnings from the scope.
`active_scope.json` records the four previously declared training contracts:
29-day 6755 put / 7030 call, and 57-day 6775 put / 7050 call (about 2% OTM).
All four VIX maturities remain active; wings remain secondary diagnostics.

Next bounded work:

1. Assess near-ATM residual variance, control-mean uncertainty and pilot stability
   using saved data, then freeze the baseline/linear-control comparison and a
   sample budget before fresh checking. Do not select gains on checking results.
2. Check fresh near-ATM sampling accuracy and actual runtime at 512 steps/day.
   Positive estimates and roughly 1% SE/price are encouraging but do not meet
   the existing quote-relative target (SE <= scale/8).
3. Once sampling permits interpretation, assess coupled 256/512 stock-price
   changes on these contracts and refine further if needed. Independently
   check normalizer quadrature; control variance reduction does not clear bias.

Retain signed negative controlled wing estimates as unusable-price diagnostics,
not clipped prices. No model/curve/state refit or new broad simulation was run
for this priority update. See [the moneyness audit](SPX_MONEYNESS_POSITIVITY_AUDIT.md)
for evidence and limits. The older queue below is retained as research history;
this section determines the current priority.

User direction: record the VIX tail-control next step, move forward for now,
and keep all maturities. The working interpretation is to return to SPX
numerical diagnostics while leaving the VIX extension pending.

The user subsequently confirmed: **"Yes, let's focus now on that for a while."**
SPX sampling, timestep and normalization diagnostics are therefore the current
explicit focus. The first bounded path-sampling experiment is recorded in
[SPX_PATH_IMPORTANCE_RESULTS.md](SPX_PATH_IMPORTANCE_RESULTS.md); the VIX queue
and full maturity scope below remain in force.

SPX follow-up: [SPX_VOLATILITY_SAMPLING_RESULTS.md](SPX_VOLATILITY_SAMPLING_RESULTS.md)
records independently selected volatility-path sampling that improves fresh
variance in 20/24 cases, while no quote-relative sampling screen passes.
Additional Gaussian-score controls help only modestly. Next prioritize paired
time corrections: standalone level costs, independent correction-pilot fits,
and a bounded matched-cost MLMC comparison. The 128->256 warnings, pilot
instability and tail precision remain open; no model refit is cleared.

The [first MLMC prototype](SPX_MLMC_RESULTS.md) now checks independent
64/128/256/512 levels with frozen allocation and correction-specific controls.
It loses on the central variance-cost objective in both fresh trials and
**does not replace direct vol7 simulation**. Sparse top-level tails and
256->512 warnings remain. Before another allocation study, investigate better
coarse/fine Gaussian coupling on tractable oracles; conditional bridge
reflection is a candidate only, with target-law/likelihood/adaptedness checks
still to be established. Do not interpret the completed prototype as a
completed improved coupling or continuous-model accuracy clearance.

The [endpoint reflection study](SPX_BRIDGE_COUPLING_RESULTS.md) has now completed
the algebraic/Gaussian/Black and independent checking step. Endpoint and target
law checks pass, but all four central variance-cost comparisons lose (gains
0.67–0.77). Direct vol7 and ordinary corrections remain the references.
The same-coefficient averaged residual retains 86–96% of controlled variance.

**Next SPX hypothesis, pending:** cheap projected-OU innovation martingale
controls, including centered quadratic innovations with known conditional
covariance, as complements to the existing stock hedges. Preserve full-state
conditional means; the three readout projections alone are not Markov. Prove
the zero-mean/likelihood correspondence, validate Gaussian/Black and forward
oracles, train on separate pilots, then check direct-price and correction
variance/cost separately. This is not an implemented control or a bias fix;
512-step, tail/pilot and normalization gates stay open. VIX scope is unchanged.

The [full OU-innovation study](SPX_INNOVATION_CONTROL_RESULTS.md) now implements
and independently checks that hypothesis. Central direct-price variance gains
around 1.8–1.9 are nearly consumed by a 1.8 runtime multiplier; corrections
lose on cost. Repeated far-call paired mean shifts, contract-control mean
warnings and pilot sensitivity remain unresolved, despite Gaussian/covariance
and bare-innovation checks passing. The full-nine dictionary remains experimental.

**Next candidate:** an actually cheaper three-linear-control kernel. Old-pilot
two-fold ablations retain most direct-price benefit, but are post-hoc design
evidence only. Avoid unused square/mixed calculations; profile standalone
pricing; freeze before fresh checking and retain every representative strike.
Audit contract-specific control tails/paired mean shifts as well as bare
innovation sums. Keep correction performance and timestep bias separate;
baseline controls remain the reference and no market refit is cleared.

## Scope retained

[active_scope.json](active_scope.json) remains authoritative and unchanged:
**27/55/118/153-day VIX maturities, all 175 training options and four forwards**.
No VIX maturity is deferred or removed. The existing SPX maturity/split scope
also remains unchanged; representative pilot contracts are diagnostics, not
an exclusion of the other targets. Held-out market observations stay reserved.

## Pending VIX extension — recorded, not launched

Screen payoff-shaped Gaussian tail controls with analytically known means,
using the existing nine-control/ridge design and pooled 28-pilot fit as the
reference. Learn and compare candidates on design-only pilot replications,
freeze the dictionary and coefficients, then validate on new shared checking
paths at all four maturities. Retain raw/legacy estimates, payoff-hit counts,
negative adjusted estimates and paired sample-count checks. Keep conditional
quadrature, time integration and normalization checks separate.

This extension remains pending for a later bounded experiment. It is not a
completed result, an automatic background job or a convergence clearance.
The current checkpoint is
[VIX_GAUSSIAN_STABILITY_RESULTS.md](VIX_GAUSSIAN_STABILITY_RESULTS.md): all four
pooled forward SE screens pass, but the 153-day pass is borderline and not
supported by every pilot fit; only 45/175 pooled option SE screens pass, and
118-day second-moment warnings remain recorded. No market refit is authorized
by those numerical screens alone.

## SPX work to advance next

Resume the unresolved sampling, stock-integration and normalization work in
[NUMERICAL_STABILITY_RESULTS.md](NUMERICAL_STABILITY_RESULTS.md):

1. Design a small path-level Gaussian importance-sampling pilot against the
   current adapted hedging-control baseline. SPX depends on the whole stock/
   volatility path; the VIX terminal-state proposal cannot simply be reused.
   Preserve the Q model, use exact path likelihood weights and separate pilot
   learning from checking. Check martingale/reference expectations and rare
   option payoffs with tractable Gaussian/Black oracles before market runs.
2. Keep coupled 64/128/256 steps-per-trading-day comparisons. Exact OU state
   transitions do not remove nonlinear stock-integration bias; the old 2/8-step
   diagnostics remain coarse historical evidence.
3. Measure sampling error, paired-grid correction variance/cost and normalizer
   refinement separately. Preserve zero-hit and unsuccessful cases; do not
   infer accurate tail prices from a zero empirical SE or fit a new market
   model before the numerical gates are resolved.

The latest source review confirms that the SPX runner already supplies a
12-contract diagnostic, independent pilot coefficients, adapted hedging controls
and coupled 64/128/256 grids. Its recorded results do not qualify the full
surface. No new SPX simulation has been launched by this queue update.
## Linear-only follow-up completed; next diagnostic

The [direct512 linear-only experiment](SPX_LINEAR_INNOVATION_RESULTS.md) retains
a modest central variance-cost gain (1.14x in both fresh trials), but far-put
outlier leverage and 30.52x half-pilot SE inflation prevent general adoption.
All primary sampling screens still fail. Keep the baseline as reference.
Next examine saved-data coefficient/scaling and tail-contribution diagnostics;
any shrinkage/stability rule is exploratory until frozen using pilot data and
validated on new checking samples. Do not select a rule on existing checking
performance. Stock-timestep and normalization accuracy remain separate open gates.

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
