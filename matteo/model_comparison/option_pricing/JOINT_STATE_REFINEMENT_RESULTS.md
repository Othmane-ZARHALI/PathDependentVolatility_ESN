# Restricted initial-state fitting, five-piece curve and VIX forward convention

Status: **complete exploratory checkpoint**, 23 October 2026. Keep the selected
state with three curve pieces as the working reference. Five pieces give a small,
inconsistent extra improvement and are retained without adoption. Neither is a
qualified market fit. Runs 003 resume the frozen selections from runs 002 after
a checking implementation failure; neither search was repeated.

## Question and scope

Continue actual joint fitting, rather than another variance-reduction feature
search. Preserve the 707 cutoff-3 SPX training options, 175 VIX training options
and all four VIX forwards, frozen price scales, family/expiry weights and signed
failure recording. SPX fitting remains 29/57/211/302 days; 119 stays reserved,
666/1030 remain diagnostics, and 1/1401/1765/2129 are deferred. No held-out model
score is computed. Historical results and unrelated working-tree edits remain.

The reference is the previous **joint** three-level curve at zero state,
`[0.0600941554,0.0748505719,0.0383634764]`. It differs from the earlier VIX-only
anchor. Its earlier pooled loss was 347.10, but only 50/882 options were inside
bid/ask. Neither curve was adopted as a qualified fit.

## What is being changed, in undergraduate terms

The three-level curve specifies average future variance. The latent starting
state determines how the model's random volatility mechanism begins relative
to that baseline. Changing the starting state can affect distribution shape
even when its future average variance is normalized to the same curve.

`m2_joint_state_fit_002` tests seven states: zero and plus/minus one stationary
prior unit along three whitened C/Y/F directions. Each state is the declared
minimum-prior-norm full-state lift, with penalty `0.1*prior_norm_squared/3`.
The three variance levels are optimized at each state, using two starts and
at most 25 evaluations per start, within a 60-minute searching budget. This is
a discrete restricted-state calibration, not continuous/global optimization.
Structure, readout and Q risk-price scenario stay fixed.
The zero-state baseline versus selected-state comparison includes reoptimized
curve levels: it measures a joint state/curve refit, not the isolated effect of
changing the initial state while holding the curve fixed. The subsequent
three-to-five-piece comparison holds that selected state fixed.

The exact OU initial-state change is a deterministic decaying shift of a path
generated with the same random shocks. Candidate G, dated normalization M0,
stock paths and conditional VIX exposures are rebuilt. Original zero-state
anchor controls remain fixed functions of the same primitive Gaussian draws;
they are not recentered as though they belonged to the changed state. Searching
uses old searching streams only. Three NEW checking ensembles use distinct seeds.

After that selection, the prespecified [five-piece refinement](CURVE_REFINEMENT_PLAN.md)
holds its selected state fixed and adds boundaries near days 57 and 151, retaining
the 91.25/273.75-day boundaries. The old curve `[a,b,c]` embeds as `[a,a,b,b,c]`.
Day 151 is rounded to an existing stock-grid event before fitting, using dates
only. Thus stock discretization does not silently change. Bounds and the weak
adjacent-level penalty are retained; two starts have at most 30 evaluations each.

More pieces can fit maturity/window differences. For example, the 27-day VIX
depends on days 27–57, whereas 55-day VIX depends on days 55–85. A flat first
three-month piece ties most of both windows to one variance level. Splitting it
allows those constraints to move more independently. It does not give the model
independent control of every option strike's smile shape.

Curve boundaries divide **future time**. The vertical line in an option smile
divides **strike**, at K equal to the chosen forward, with the displayed OTM puts
and calls on its two sides. Adding time pieces does not directly fix the IV
display mismatch caused by using a forward different from the model's mean.
The [curve/window picture](visualizations/variance-pieces-006/README.md) shows
the two curve families and why the extra dates were chosen. A lower positive
variance level after day 273.75 still adds positive cumulative variance; a drop
in the input level is not a drop in accumulated variance. The last level is
weakly constrained by these fitting expiries and extrapolated beyond them.

The refinement selects only on searching data and reuses the three state-run
checking ensembles for a paired comparison. These are independent of search,
but the state and refinement checks are correlated and are NOT two independent
validation experiments. No checking score enters curve optimization.

## Separate numerical errors

| Source | Treatment | What is not established |
|---|---|---|
| Sampling | 512 SPX paths/expiry, eight independent 64-point VIX nets/expiry; three checking trials; frozen pilots | Tail accuracy, simultaneous confidence, selection coverage |
| Stock timestep | Existing 512 steps per trading day, exact stored events; new knots on that same grid | Continuous-time stock-price convergence |
| Dated-normalizer interpolation | NEW positive log-GH48 spline, adaptive midpoint and independent-node screen at relative 2e-7; zero-state full-grid comparison | A rigorous uniform interpolation bound |
| Normalizer GH | Order 48 fixed, candidate state updated | GH convergence from interpolation agreement |
| VIX conditional GH | Order 32 fixed with full conditional state | Conditional integration convergence |
| VIX time integration | Order 8, split at new curve knots; embedded-old-curve changes reported separately | A time-quadrature change counted as a fitting improvement |

The interpolated normalizer is an additional deterministic approximation needed
for economical candidate-state evaluation. Its screened discrepancy is kept
separate from sampling SE and GH error. It does not reset normalization at expiry
or force a sample moment to a target.

## Forward convention: a useful diagnostic, not a price correction

For a VIX option expiring at T, the model's own forward is **E_Q[VIX_T]**, at the
SAME expiry T, not at T+30 days. Converting a model option price to Black IV using
that model mean is sensible for inspecting its distribution shape. The forward
still has to fit the market quote as a separate residual. No price, strike,
OTM membership, objective weight or loss changes in this diagnostic.

The initial display check, `m2_vix_forward_convention_001`, uses the prior curve-only
pooled checks. Its 27-day selected curve has a nearby call-minus-put IV difference
of 77.71 percentage points under the market-forward convention and 1.58 under
the model-forward convention. The strikes differ, so these numbers are not an
exact same-strike discontinuity or parity test. The unchanged forward miss is
still +1.70 VIX points. Every one of its 175 model mean prices admits native IV.

The [comparison picture](runs/m2_vix_forward_convention_001/vix-forward-conventions.png)
keeps market IVs on their market forward and all horizontal strike coordinates
fixed. Only the model IV conversion changes. Model uncertainty shading is omitted
there because native IV also depends on the estimated forward and its covariance
with price; price-only bands would not capture that combined uncertainty.

An independent option-implied market forward could be obtained from matched
call/put premiums and appropriate discount/settlement conventions. These source
quotes are IVs, and the stored premiums were reconstructed using the supplied
forward. Using those reconstructed prices to infer a supposedly independent
forward would be circular. The saved parity intervals are conditional consistency
diagnostics, not new forward observations.

Scope inconsistency found during review: historical display diagnostic 001
computed raw matched-put/call parity intervals on the broader valid VIX universe,
including 35 of its 176 matched strike keys outside current training strikes
(18 validation and 17 final-test strike keys).
Those market quotes were exposed to that conditional consistency calculation;
they must not be described as wholly unseen. No model price, held-out model
performance score, parameter selection or objective used them. The historical
result is retained. Their key-only exposure list is saved in diagnostic 002's
`historical_parity_exposure.csv`; future evaluation of these keys must be labelled
quote-exposed rather than wholly blind. Diagnostic 002 and the current state/refinement diagnostics
restrict raw pairs to exact training option keys. The OTM-only training set
supplies no matched call/put pair. This is an additional reason that the current
training inputs alone do not identify an independent option-implied forward.

## Formal and implementation evidence

The current [Lean correspondence](formalization/CORRESPONDENCE.md) has thirteen
checked identities under the pinned environment, with final audit attempt 005.
The three additions cover common-forcing affine state translation, exact
normalization and expected-payoff put/call parity. Only standard axioms occur;
theorems do not certify Python, interpolation or market prices. Attempt 004 is a
retained failed proof-development audit. The earlier ten-theorem audit 003 remains.

All **179 Python tests** pass (124 package tests and 55 driver tests), including
nine new implementation tests: exact translated state/shocks, unchanged
importance weights, direct candidate stock repricing, prior coordinates,
normalizer ownership/interpolation, native-forward Black inversion, nested
three-to-five-piece prices/derivatives, five-dimensional derivatives and loss.

Run `m2_joint_state_fit_001` was interrupted before optimization to prime the
known zero-state dense normalizer cache explicitly before fresh anchor hedges.
It has no fitted result; its source snapshot and failure record are retained.
The corrected run is 002. The five-piece startup `m2_curve_refinement_001`
failed before exposure construction because the provenance snapshot requires
absolute input paths; corrected run 002 explicitly resolves its run/output
paths. The failed startup is retained and has no calibration result.
State run 002 completed all seven searches, then failed before its first fresh
score: the final 24-step stock chunk used a 2048-entry normalizer slice. The
correction slices the normalizer to each expiry's local grid before chunking.
State run 003 reuses the frozen selection and searching artifacts, and generates
new checking streams with root seed 2026102403. Five-piece run 002 completed its
search but stopped waiting for the failed parent's checking prices. Run 003
reuses that frozen five-piece selection on the resumed parent's new draws.
Both failed runs, logs, source snapshots and selected parameters remain intact.
This failure does not change either optimizer's recorded searching results.
The corrected sources pass the full 124 package plus 55 driver tests again.
An initial rerun command used unittest directory discovery on a namespace
package and failed before tests; the explicit module loader then passed all 55.
Initial unit fixtures incorrectly counted a terminal
projection and omitted required curve events; fixtures were corrected without
changing the pricing engine. Failed test/Lean development output is not a pass.

## Results and return-to-fitting decision

The seven-state searching stage is complete (2503.34 seconds). All seven
profiles and unsuccessful directions remain in `state_profiles.json`:
zero 333.6764, C-negative 384.5152, C-positive 288.8213, Y/F-negative 412.5854,
Y/F-positive 420.7160, F-negative 631.3122 and F-positive 822.0228. Whitened
Y/F directions reflect the prior covariance, not independent raw projections.
The selected projections are **C0=1, Y0=0, F0=0**, prior squared norm one, with
three levels `[0.0584806825,0.0727920966,0.0357009546]`. These are searching
results only. The selection froze before the new checking draws were generated.

Five-piece corrected run `m2_curve_refinement_002` uses added boundaries
56.999999999999986 and **151.00054327286466 days**, chosen from the stored grid
before optimization. The nominal 151-day rounding is about 47 seconds and
changes no existing stock grid or old three-piece price target.

The five-piece searching loss falls **288.82128 to 284.80479 (1.39%)**;
selected levels are `[0.0579766618,0.0611022826,0.0692341111,0.0761375463,
0.0300158138]`. Their square roots are annualized volatility scales of about
24.08%, 24.72%, 26.31%, 27.59% and 17.33%, respectively. These are curve inputs,
not option implied volatilities and not estimates of average volatility itself.
The stock-grid events and the old three-piece family are preserved. On searching
draws, splitting VIX time quadrature with the old curve held fixed changes only the 55-day
calculation: maximum per-state VIX difference 2.08015e-5 points, weighted mean
forward difference -6.36e-12 points, maximum relative integrated-variance
difference 1.00523e-7. This deterministic integration change is recorded separately
in `embedding_quadrature_changes.csv`; it is not counted as optimization gain.
Across the three fresh draws, the largest per-state change is 6.75054e-5 VIX
points and largest relative integrated-variance change 2.98781e-7. Weighted
forward changes remain at most 1.22221e-11 points in absolute value. These are
same-curve split-rule differences, not general time-quadrature error bounds.

The resumed state/curve run 003 and its independent score audit are complete.
Three paired fresh losses are **418.1567→370.0608**, **465.9733→424.9127** and
**316.9800→288.3994**. These compare the original joint zero-state curve against
the searching-selected C0=1 state and reoptimized three-piece curve. The pooled
price comparison is:

| Component | Previous joint, zero state | Selected state + three-piece refit |
|---|---:|---:|
| SPX | 252.1582 | 235.0226 |
| VIX options | 89.0203 | 70.2905 |
| VIX forwards | 14.2859 | 11.5868 |
| Curve penalty | 0.000155 | 0.000158 |
| State penalty | 0 | 0.033333 |
| Total | **355.4644** | **316.9334** |

The pooled reduction is **10.84%**. The descriptive paired 24-block jackknife
SE of the loss difference -38.5311 is 4.5560; all leave-one-block differences
remain negative (-39.8282 to -35.7750). This is conditional on fixed pilots and
the approximate pricer, with few blocks and unresolved rare tails. It is not a
coverage guarantee. Comparing its absolute total with the older pooled 347.10
would mix simulation ensembles: the relevant baseline here is 355.46 on these
new shared checking draws. The SPX 29-day slice worsens slightly; the other SPX
and all four VIX expiry losses improve. All VIX forwards remain too high:
**+1.7115/+0.6619/+2.2645/+1.7582 points**.

Only **41/882** pooled objective prices are inside bid/ask (baseline 39).
Only 34/882 nominal SE screens pass, all VIX; all 707 SPX screens still fail.
Primary fresh and pooled prices are nonnegative. Secondary-wing negatives are
0/8/3 by fresh trial for both candidates, retained rather than clipped. Signed
unsuccessful search estimates also remain archived. No pooled
VIX squared-moment check exceeds its nominal two-SE screen; this does not certify
conditional GH, time integration or tail convergence. Two mass-warning rows are
the same 211-day trial-2 sample reused for both candidates, not two independent
failures. No raw stock-forward warning occurs. At the estimated selected mean,
the 153-day necessary market-bid moment margin is **-19.7722** (baseline -9.0390).
In plain terms, that curve may supply too little mean-square VIX to support all
the option bids at that mean. Mean uncertainty and numerical bias prevent a
certified infeasibility conclusion.

The current state's 27-day nearby call-minus-put IV difference is 78.4847
percentage points using the market forward and 2.2572 using its own same-expiry
mean. These are neighboring different strikes, not an exact parity error.
All 175 pooled candidate prices admit native IV. The right-hand pictures remain
too flat relative to much of the market's rising VIX IV curve: removing the
forward-convention gap does not solve the distribution-shape mismatch.

The five-piece fresh comparisons are **370.0608→366.1933**,
**424.9127→431.0621** and **288.3994→280.9778**. Two improve and one reverses.
The pooled comparison at fixed selected state is:

| Component | Three pieces, embedded baseline | Five pieces |
|---|---:|---:|
| SPX | 235.0226 | 233.4911 |
| VIX options | 70.2905 | 69.4175 |
| VIX forwards | 11.5868 | 11.9618 |
| Curve penalty | 0.000158 | 0.000225 |
| State penalty | 0.033333 | 0.033333 |
| Total | **316.9334** | **314.9040** |
| Prices inside bid/ask | **41/882** | **38/882** |

The additional pooled gain is **0.64%**, difference -2.0294 with descriptive
jackknife SE **3.1815**. Leave-one-block differences stay negative, -4.3162 to
-1.2411, but this narrow perturbation is not a confidence interval and does not
resolve the trial reversal or uncertainty larger than the observed gain.
Five pieces worsen three of four SPX expiry losses and VIX days 55/153, while
improving SPX day 29 and VIX days 27/118. Its forward errors are
**+1.6173/+1.1158/+1.6849/+2.2940 points**. Only 34 nominal SE screens pass,
again all VIX. The 118/153-day necessary moment margins are -6.0078/-9.5917 at
estimated means; numerical and mean uncertainty remain. No pooled squared-moment
or raw stock-forward warning occurs; the same duplicated mass warning remains.
All 175 five-piece pooled VIX prices admit native IV. No result is clipped.

### Concrete fitting decision

**Keep C0=1, Y0=F0=0 and the three-level curve as the provisional working
reference. Do not add further curve pieces now.** The joint state/curve refit
has a persistent sampled benefit; the extra two time levels do not show a
convincing additional benefit. Preserve the five-piece selection, all profiles
and failed runs for future comparison. This is a decision about the working
reference, not acceptance of a market or physical fit.

The next bounded fitting hypothesis is the **already declared unspanned-risk-price
scenarios**, one at a time, with the selected initial state and three-piece
curve family, refitting curve levels and preserving all objective contracts.
Start with the existing feedback-risk-price pair (+0.25,0) and (-0.25,0), keeping
the equity premium and readout/architecture fixed. Compare them with scenario 0;
only then consider the declared spike-risk-price pair (0,+0.25)/(0,-0.25).
Declare a finite search budget and freeze selection before new numerical checks.
This is **not launched in this checkpoint**. It follows the existing calibration
sequence rather than reopening the completed variance-reduction studies.
Candidate Q mean/state paths, dated normalization, stock paths and conditional
VIX exposures must be rebuilt; recheck primitive-law/control correspondence and
sampling likelihoods. Revisit P diagnostics after a risk-price change because
the Q normalizer also enters the declared P variance convention.

At undergraduate level: extra time pieces adjust **when average variance is
large**. Risk prices can change the pricing distribution and are the next
declared way to test the remaining smile-shape and forward trade-off. There is
no guarantee that those scenarios will solve it. The current pictures show
substantial shape errors even after using the correct model mean for IV.
Continue actual bounded fitting while keeping sampling and deterministic bias
checks separate; reserve qualification until the full numerical gates pass.

All 16 score records across the two runs replay independently, including row
membership, source hashes and archived sample means. The [new plots and reading
guide](visualizations/joint-state-curve-005/README.md) were inspected as PNGs;
SVG exports and input/output hashes are saved. Additional relevant SPX market
slices remain in the seven-maturity focused view, with no new held-out model
scores. The prior broad raw-parity exposure is recorded above and in AUDIT.md.

No new qualified market or physical fit, global optimum, or model ranking is
claimed.
