## 27 October 2026 — continuous measure, forward constraints and multilevel repair

Completed VIX-only continuous M2/stress Beta-measure experiment003, after rejecting002's noisy55-day forward calibration. The model is unchanged: normalized G^(1-u) times the normalized OU stress exponential of strength b*u, integrated against one Beta measure. Four forward equalities use32768-path boundary plus8192-path paired difference; checking uses65536 NEW boundary paths plus16384 reused difference paths. Do not claim wholly fresh validation or exact population equality.

Best bounded fit remains near kappa=.01; a broader measure is retained as a training-selected diagnostic, with b=5 at its bound. No clear improvement established over the two-atom reference. All217 quotes, four maturities, seven free active xi levels and fitted floor preserved. SPX paused; original baseline and previous working VIX candidate unchanged. See CONTINUUM_RESPONSE_RESULTS.md, runs/m2_continuum_response_003 and visualizations/continuum-response023.

# Continuous stochastic vol-of-vol after the cap profiles

**27 October: tighter VIX forwards plus fitted-floor experiment completed.**
See [FORWARD_FLOOR_RESULTS.md](FORWARD_FLOOR_RESULTS.md),
[FORWARD_FLOOR_PLAN.md](FORWARD_FLOOR_PLAN.md), run `runs/m2_forward_floor_001`,
plots `visualizations/forward-floor021`. All16xi intervals retained, seven active
levels independently fitted in both comparisons; SPX remains paused. Four training
forward equality constraints hold within1.9e-10VIXpoint. Fresh32768-point/expiry
comparison: previous mixture loss10.2653; tight forwards/fixed floor12.7780;
tight forwards/free floor1.7830. Largest fresh forward gap falls0.9001->0.02694,
with sampling SEs reported. Releasing the floor chooses12%volatility, at the
search ceiling, b3.1710,w0.66150,half-life104.215days. Best tested, not a global
optimum or settled floor estimate.55-day lowest put is close to market,27-day
left end improves; lowest long-dated puts become too cheap. Wing IV MAE8.605->8.445
on average;27-day far-right level worsens.84/217quotes within bid–ask.
Five small negative fresh native estimates and missing IVs remain, unaltered;
raw nonnegative put-payoff diagnostics are saved separately. No clipping/smoothing.
Six independent price replays,11Jacobian directions,32fresh seed regenerations,
exact budgets,floor support,291source/fiveinput hashes pass. This is numerical
verification, not population-forward equality or complete calibration. Current
working VIX candidate is the free-floor experimental model in this run; original
selected model unchanged. No SPX computation, new component or Lean theorem.
Historical checkpoints below keep their original context.

**26 October: point 2, two-response VIX-only comparison completed.**
See [MIXED_RESPONSE_RESULTS.md](MIXED_RESPONSE_RESULTS.md),
[MIXED_RESPONSE_PLAN.md](MIXED_RESPONSE_PLAN.md), run `runs/m2_mixed_response_001`,
figures `visualizations/mixed-response020`. Both product and additive mean-one
mixtures fitted seven independent VIX-sensitive levels of the original16-piece
curve. The nine off-window entries remain unchanged; all217quotes/four forwards
remain, SPX paused, no clipping/smoothing. Native-price winners and descriptive
shape selections were frozen before new16384-point/expiry numerical checking.
The additive M2/stress price fit improves fresh native loss20.0052->10.2960 versus
the prior fitted single factor (paired change-9.7093, jackknife SE0.2195), wing
IV MAE30.81->8.51points. Long-dated slopes now have enough steepness; short-dated
slopes overshoot.49/217quotes inside bid–ask, residual forward/ATM errors remain.
The additive fit gives99.57% weight to the new component, b at5bound, half-life
6.588days; it largely replaces the old response in this exploratory readout.
Two short variance pieces reach0.20. No acceptable full-smile/joint fit or model
adoption is claimed.72implementation checks and independent14price replays,
32fresh scramble regenerations, exact budgets,285source/fiveinput hashes pass.
Sampled GH/time VIX refinement max0.000118point for additive winner,0.0332 across
all candidates' extreme probes. Same market targets; no market OOS or Lean claim.
Preserve baseline, scope and historical records. Further timescale/readout changes
remain proposals; no second response timescale has been implemented in this run.

## 26 October: VIX-only right-wing investigation, Othmane documents and literature

**Current direction:** focus on all217 VIX options at27/55/118/153 days and their
four forwards. SPX fitting/repricing is deferred until the VIX shape is resolved.
The selected M2 baseline is unchanged. This turn performs research, a static
marginal diagnostic and exact-arithmetic payoff checks; it launches no new
reservoir calibration. The previous experiment checkpoints below are historical.

**Main finding:** freeing the variance budget is necessary before judging new
shape mechanisms. With the current fixed xi curve, no nonnegative VIX marginal
law can put every archived native option inside its bid–ask interval and match
its forward at55/118/153 days. The obstruction persists with forward tolerance
plus/minus0.5 VIX points. This statement concerns the archived quotes and the
frozen second moments, not market arbitrage or impossibility of a richer model.

### 1. A model-independent obstacle in the frozen-curve trials

For the intended M2 model, tower expectation gives

    E_Q[VIX_T^2] = (10000 / Delta) integral_T^(T+Delta) xi(u) du,
    Delta = 30/365.

Every tested mean-one normalized readout/layer has the same second moment when
xi is fixed. Redistribution of variance can change the shape, but it cannot
increase this total budget. Matching the VIX forward only matches E[VIX], not
E[VIX^2]. Matching one does not imply matching the other.

The following are **necessary lower bounds**, not asserted sharp minima. Units
are squared VIX points; the forward is fixed at its observed value.

|VIX expiry days|Current curve E[VIX²]|Certified necessary lower bound|Current curve excluded?|
|---|---:|---:|---|
|27|477.2295|468.9993|No obstruction established|
|55|536.0269|556.6730|Yes|
|118|580.8043|648.6257|Yes|
|153|610.9326|695.2954|Yes|

Allowing the VIX forward within plus/minus0.5 points gives lower bounds
462.8004 /551.3278 /636.8530 /682.3170 respectively. The last three still exceed
the old budgets. At the exact forwards, the necessary increases in the last
three window-average variances are at least about3.85%,11.68%,13.81%.
These are lower bounds on the required budgets, not proposed calibrated curves.

**Method and verification.** A finite-support linear programme first searches
for a cash/forward/native-options portfolio pi(y) below y². Native puts are
converted algebraically to calls plus cash/forward for the payoff check. On
every interval between quoted strikes, y²-pi(y) is a quadratic with leading
coefficient1, so its global minimum over y>=0 is attained at a strike, zero,
or an in-interval stationary point. `fractions.Fraction` checks all these
candidates with exact rational arithmetic after the proposed coefficients are
serialized and read back. The cash term is conservatively reduced to leave a
positive payoff margin. Worst bid/ask prices and the worst allowed forward
then give E[VIX²]>=cost(pi). The audit checks12 such certificates, across four
maturities and forward tolerances0/.1/.5. The result does not rely on truncating
the VIX tail at the LP support, on Monte Carlo, or on the LP solver's claim
of infeasibility. It is exact arithmetic for the recorded decimal inputs;
there is no Lean formalization.

Inputs are the original native quote CSV archived with `m2_postcap_vov_001`
and the selected curve in `m2_maturity_curve_fit_001/run.json`. An independent
canonical-call marginal LP also finds finite-support fits with nonnegative
weights at all four maturities when the second moment is free (within solver
tolerance); with the old moment fixed it fails at the last three. Its support
and spacing refinements are diagnostics only. No unique marginal distribution,
cross-time process or stock/variance coupling is inferred.

**Implication:** the previous frozen-curve mechanism trials could improve wing
shape at some strikes, but could not attain a full within-spread fit on the last
three VIX slices while matching forwards. A new VIX-only calibration should
allow the VIX-window xi levels to move together with the shape parameters.
The minimum volatility/floor policy and off-window levels can be retained.
No new xi values have been fitted or adopted in this research turn.

Evidence: `runs/m2_vix_design_research_001/moment_certificates.json`,
`moment_certificate_audit.json` (PASS), `moment_certificate.py`,
`distribution_feasibility.json`, `distribution_feasibility.py`.
The first certificate search demanded unnecessarily tight proposal convergence;
its failure log/source are retained as `moment_certificate_attempt_001.*`.
The final method accepts a conservative, nonsharp proposal and then verifies
its entire payoff inequality exactly. The adjustment is explicit in the certificate.

### 2. What the email and three new documents actually establish

The prepared packages are readable; the original PDFs need not be requested.
Equations and material figures were cross-checked against prepared page images.
The email reviewed is `/Users/matteogambara/Downloads/ESN.eml`, dated14 October,
subject `Re: ESN`. Six inline images and the body were extracted to the research
run. The body says the final grid-search/historical/joint chain was pending.
Its VIX figure is genuinely much closer visually, but is not the same documented
configuration as the later19 October reports: for example the emailed readout
sketch uses6 rough/readout pieces and3 gain-side pieces, whereas cs_doc uses5
and1. Do not attach the later report's parameters or pricing estimator to the
older email figure as if their correspondence were established.

|Document|What it studies|What is relevant to VIX wings|
|---|---|---|
|gs_doc,21pages|Historical/P architecture selection on five stocks and surrogate stylised facts|Chooses rough bank, modulation and readout complexity. Its modulation ablation is P evidence, not evidence of a VIX option fit.|
|jc_doc,23pages|Historical/P per-asset calibration across three indices and six stocks|Studies roughness, leverage, Taylor, Zumbach, kurtosis; it is not the SPX/VIX option calibration. The report describes an incomplete optimization and some residual validation failures.|
|cs_doc,14pages|Joint Q pricing/calibration on19 February2026|The relevant option report:17 parameters, bounded flexible readouts, fixed P-selected modulation dynamics, four SPX slices, four VIX smiles/forwards.|

**The mechanism differs from our layers.** Its bank satisfies

    dr_i = -lambda_i r_i dt + sqrt(2 lambda_i Gamma_t) dW_i,
    dGamma = -lambda_Gamma Gamma dt + dL,
    sigma = tilde_sigma(t) exp(psi(zeta)) F(P),

where L is a positive compound-Poisson process. Gamma changes the variance of
the noise entering the bank; it does not simply multiply the final volatility.
Given Gamma's path the bank is Gaussian; unconditionally it is a mixture of
Gaussian scales. The endogenous return-feedback bank is also different from
our exogenous M2 projections. The published/source model is not a drop-in
replacement for a normalized M2 readout.

The Gamma jumps do not directly jump the stock price: the bank and stock remain
Brownian-driven continuous processes in this specification. Gamma itself does
jump, so copying the exact mechanism would still change the current no-jump
working design. A continuous Gamma alternative is discussed below, not implemented.

**The later report does not claim the remaining wing is solved.** cs_doc p8,
Table4 gives model slopes1.02/.67/.41/.33 versus market1.21/.86/.54/.45. Its
p12 summary explicitly says the2–5 month wings are too flat. Its p13 recommends
recalibrating the Gamma dynamics under Q rather than keeping their P estimates.
This is evidence for a useful direction, not evidence that our full target is
already attained by that implementation.

**Scope and loss differ materially.** cs_doc p2 selects20 strikes per VIX
expiry, with maximum log(K/F) .91/.98/.95/.93, hence maximum K/F about2.48/2.66/
2.59/2.53. We retain49/58/55/55 quotes, extending to K/F4.98/8.20/9.44/9.21.
Its p6 objective weights VIX quotes by exp(-k²/(2*.2²)), strongly emphasizing
the centre. At k=.9 the unnormalized weight relative to ATM is about4e-5.
Our K/F>2 wing slopes cannot be compared numerically with its slopes over its
whole narrower selection. The overlay `strike-range-comparison.png` shades
its reported range without digitizing or inventing its model curves.

**Estimator correspondence matters.** cs_doc p5 uses a Gaussian base that
replaces the future Gamma path by its conditional mean and freezes feedback
at a decaying conditional proxy. A regression on realized variance/base then
supplies a multiplier clipped to[.1,10]. Gamma mixing is not exactly Gaussian
after averaging over its future path; this base is an approximation. The
report's p13 explicitly leaves a nested Monte Carlo check pending. Its p4
also describes self-normalized weights and sample stock martingale correction;
p13 says its reported jump reference equals the model law (weights1).
These differences are reasons to validate any replication independently, not
claims that the reported plots are invalid. Matching raw code/results listed
in the report (`cs_flexd.py`, `smile_ci.json`, etc.) was not found in this checkout;
no independent repricing of Othmane's model is claimed.

### 3. Why a larger exponential OU multiplier can still flatten

For our latest factor L=exp(gamma Z-gamma² Var(Z)/2), with constant gamma and
OU diffusion sqrt(2a), the new Brownian coefficient of log L is the constant
gamma sqrt(2a). The absolute variance diffusion is state-dependent, but this
adds a Gaussian log-amplitude factor rather than making the existing reservoir
noise scale random. Adding such a factor can broaden VIX outcomes without
creating enough increase in relative responsiveness at high states.

This is a mechanism diagnosis, not a theorem that our finite30-day VIX is
lognormal. In a one-exponential Gaussian benchmark, the short averaging window
leaves VIX close to lognormal, giving an almost flat Black smile. The exact
30-day conditional averaging, old nonlinear readout and floor can change that.
The existing experiments support insufficiency of the tested grid, not a
universal impossibility result for the architecture after recalibration.

Changing only a deterministic mean in a fully mean-normalized lognormal factor
is a redundant way of seeking shape: centering absorbs the mean shift. Othmane's
initial Gamma is different: it controls the future variance of the bank shocks.

### 4. Primary literature: mechanisms and the limits of the evidence

1. **Horvath–Jacquier–Tankov, Volatility options in rough volatility models
   (2020), [paper](https://arxiv.org/html/1802.01641).** Section4 puts an independent
   positive affine Gamma inside the Volterra noise. Its generator permits a
   continuous square-root diffusion as well as jumps. Their VIX calibration
   demonstrates the mechanism; it is not a completed joint SPX/VIX calibration.
   The affine transform uses moment/Riccati conditions; it cannot simply be
   copied to an arbitrary capped M2 readout.

2. **Fouque–Saporito, Heston Stochastic Vol-of-Vol Model for Joint Calibration
   (2018), [paper](https://arxiv.org/html/1706.00873v1).** The coefficient of the
   variance diffusion is eta(Y_fast,Z_slow)*sqrt(V), with fast/slow stochastic
   factors. Joint calibration examples use perturbative pricing approximations.
   This supports continuous diffusion modulation as a separate mechanism.

3. **Bourgey–De Marco–Gobet, Weak approximations and VIX option price expansions
   in forward variance curve models (2022 preprint),
   [paper](https://arxiv.org/html/2202.10413).** Section3 replaces one exponential
   forward-variance response by a convex sum of two with different sensitivities.
   Section3.3 calibrates both mixed standard and mixed rough models to four
   VIX slices, fitting the variance curve to futures and the sensitivities/mix
   to the smiles. Parameters have a term structure. This is direct VIX-only
   calibration evidence with no requirement for roughness or an independent
   stochastic vol-of-vol process. Their weak formulas are approximate; their
   standard-model benchmark also uses direct two-dimensional quadrature.

4. **Jacquier–Muguruza–Pannier, Rough multifactor volatility for SPX and VIX
   options, [paper](https://arxiv.org/html/2112.14310v2).** Adds separate factor
   directions to control ATM level, skew and curvature; different parameter
   choices produce different curvature signs. These are small-time results
   and numerical insights, not a guaranteed finite-wing fit at our maturities.

5. **Romer, Empirical analysis of rough and classical stochastic volatility
   models to the SPX and VIX markets (2022),
   [primary paper](https://www.tandfonline.com/doi/abs/10.1080/14697688.2022.2081592).**
   The two-factor quadratic-hyperbolic transformations obtain useful joint fits
   on the tested dates, even with ordinary exponential kernels. The transformation
   is chosen partly to reduce the extreme tails of exponential mixtures. Thus
   'make the tail as heavy as possible' is not a sound model-selection rule;
   the mapping and component crossover matter. This is another model/dataset.

6. **Gazzani–Guyon, Pricing and calibration in the4-factor path-dependent
   volatility model (2025), [paper](https://arxiv.org/html/2406.02319).** An
   additional one-sided quadratic response to positive return-memory states
   unlocks upside volatility and supports joint calibration. This is an
   endogenous PDV mechanism, not our exogenous M2 quadratic spike term.
   Its neural conditional-VIX approximation speeds pricing; it does not itself
   enlarge the model's attainable smile family. Any adaptation needs independent
   pricing checks and a clear distinction between state and readout changes.

7. **Forde–Smith, Markovian stochastic volatility with stochastic correlation
   (2025 working paper),
   [author PDF](https://martinforde.github.io/SPX-VIX-Calibration-Stochastic-Correlation-New.pdf).**
   Shows how nonlinear variance diffusion/mean reversion shape short-time VIX
   smiles, with stochastic correlation providing separate SPX freedom. Also
   discusses a Bass construction for one VIX marginal. Asymptotic/single-slice
   evidence cannot be treated as an all-maturity calibrated M2 result.

8. **Zaugg–Grzelak, The VIX-Derived Volatility Model: A VIX-first Joint SPX-VIX
   Framework (August2026 preprint), [paper](https://arxiv.org/html/2608.01479v1).**
   Models VIX first with local volatility and derives an SPX variance coupling.
   Useful strategic comparator for the requested order of work. The paper says
   positivity of the general coupling is not guaranteed and gives numerical
   consistency checks. It is an alternative architecture, not an automatic
   reservoir embedding or a recommendation to switch models now.

Further checked: [Alos–Garcia-Lorite–Muguruza](https://arxiv.org/html/1808.03610v1)
for positive VIX skew via mixed responses, and
[Baschetti–Bormetti–Rossi](https://arxiv.org/html/2507.09412) for accelerated
PDV calibration. Faster optimizers or neural surrogates do not repair an
incompatible frozen second moment. Model-free SPX/VIX consistency is studied by
[Guyon–Menegaux–Nutz](https://arxiv.org/abs/1609.05832); our one-dimensional payoff
certificate is a separate, simpler calculation and is not their joint bound.

### 5. Ranked next experiments: VIX only

**First, repair the comparison.** Refit the relevant VIX-window xi levels jointly
with shape, using the exact forward target and all217 original bid/ask quotes.
Keep the minimum volatility/floor policy, original banks/rates and old baseline.
Do not judge future VIX candidates against SPX loss during this phase. A necessary
moment budget is not a sufficient fit criterion, so also check futures and the
entire call/put price curve. Do not set xi equal to the lower bounds and declare
a calibration.

**A. Fair rerun of the current after-cap mechanism with a free VIX curve.**
The prior fixed-curve scan was not a full calibration of the enlarged model.
Let its strength/persistence and the VIX-window levels move together. This is
the cheapest way to learn how much of the failure was the budget restriction.
The analytic conditional factor and old Gaussian reservoir pricing remain useful.
No assertion that this will supply the missing slope is made.

**B. Two response components, rather than one common amplitude. Preferred next
new readout.** For a centered Gaussian bank projection U with variance q(t),

    R_mix(t) = (1-w) exp(b1 U_t-b1² q(t)/2)
             + w exp(b2 U_t-b2² q(t)/2),  0<=w<=1, b2>b1.

Each component has mean1. Ordinary states can be dominated by the less sensitive
component; sufficiently high U shifts the balance toward the stronger one.
The derivative of log R_mix with respect to U is the component-weighted mean
of b1,b2; its derivative is their nonnegative weighted variance. Thus the
relative response increases through the crossover. That is a concrete difference
from a single exponential, whose log response has constant derivative.
This calculation concerns the readout, not a proof of positive Black-IV slope.
The published mixed-standard-Bergomi model is a useful independent benchmark;
an M2 adaptation can mix an old normalized component with a stress component
and retain an explicit old-model boundary. A second correlated Gaussian factor
can then add curvature freedom if one component direction is insufficient.
This is a sum of variance contributions on each path, not a static choice of
separate model regimes or a mixture of completed option prices.

**C. A flexible asymmetric readout below a remote cap.** Replace a few fixed
quadratic/soft-cap slopes by a small positive monotone spline on a stress score
(or two positive components), keeping a bounded or justified tempered tail.
Use a few interpretable knots, not a large black-box fit. Move response where
quoted-wing events live without forcing an equally large ATM change. Othmane's
bounded piecewise-linear readouts and the hyperbolic-transform evidence motivate
this flexibility; neither is a theorem about M2. The old cap permits VIX outcomes
far beyond quoted strikes, so a hard support ceiling has not been established as
the cause. Simply raising the old cap is not the preferred next move.

**D. Actual random reservoir-noise scale. Most direct analogue of Othmane,
but a larger model/pricer change.** A continuous alternative is

    dGamma = kappa(theta-Gamma)dt + nu sqrt(Gamma)dB,
    dH = -a H dt + sqrt(2a Gamma)dW.

Use H in a positive, bounded/tempered variance readout; calibrate Gamma dynamics
and its initial state under Q to VIX. If the original bank must remain untouched,
H is a new experimental bank. This is not the same as exp(gamma Z) multiplying
the final variance. Conditional on Gamma paths, H has Gaussian noise; after
averaging them out it has a random scale. The old unconditional Gaussian
conditional-pricing cache cannot simply be reused for this new bank. A new
conditional integration method needs direct/nested checks. Positive/CIR and
readout integrability conditions must be verified for the chosen specification.
An exact copy of Othmane's positive jump modulation is a separately labelled
comparison, not the default continuous/no-jump next implementation.

**E. Nonlinear local diffusion / VIX-first comparator.** A diffusion coefficient
whose relative response increases over the observed stress range is another
route. A directly fitted VIX local-vol model can serve as a benchmark for the
required marginals. It must not be called an M2 or SPX-consistent solution until
the rolling30-day variance coupling and positivity are established. Single-slice
fits or four independently fitted laws do not establish a coherent time process.

### 6. Acceptance and uncertainty for this phase

- Fit VIX futures and native OTM prices, keeping market-forward IVs for plots.
  Preserve the full strike range; a separate narrow-range reproduction can
  compare fairly with cs_doc but cannot replace the full-wing target.
- Give the body, shoulder and right wing explicit weight within each expiry,
  rather than downweighting the wing with an ATM Gaussian. Keep bid–ask price
  residuals primary; use wing IV slope/curvature to diagnose shape, not as a
  replacement for fitting prices.
- Fit level/variance-window parameters and shape together; keep old results as
  controls. Readout-only changes reuse exogenous states, not completed stock
  paths. No stock simulation is needed for this VIX-only stage.
- Require fresh complete outer scrambles and inner quadrature/refinement checks.
  For a Gamma-driven bank use independent nested checks at both typical and
  wing-dominant states. A positive clipped/regressed estimator is not exact merely
  because it produces smooth smiles.
- Track E[VIX²] against the refitted xi budget and weighted tail contributions.
  Low tail hit counts or concentrated weights can leave wing estimates unstable.
  Numerical checking uses the same market data; no market OOS claim follows.
- Resume SPX only after a credible full-range VIX fit, preserving the selected
  VIX target when assessing the later joint compromise.

### Evidence locations

Research run: `runs/m2_vix_design_research_001`. Figures:
`variance-budget.png/.svg`, `strike-range-comparison.png/.svg`.
Email images/body are local research evidence, not a re-publication or an
independent model reproduction. The prepared PDFs remain byte-preserved.
The run's `research_manifest.json` records sources, reviewed page images,
literature links and exact-certificate checks. All previous run reports,
selected_model.json and numerical results remain historical and unchanged.

---

**26 October checkpoint: independent amplitude AFTER the original cap completed.**
See [POSTCAP_VOV_RESULTS.md](POSTCAP_VOV_RESULTS.md),
[POSTCAP_VOV_PLAN.md](POSTCAP_VOV_PLAN.md), run `runs/m2_postcap_vov_001`,
and figures `visualizations/postcap-vov017`. This is the authorized first
alternative; no mixture was implemented. One independent Q OU exponential
L=exp(gamma*Z-gamma^2*Var(Z)/2) multiplies G_base/M_base after the old cap.
The old normalizer, reservoir, coefficients, all16 xi levels and floor remain fixed.
Nine nonzero strength/persistence cases plus baseline have full VIX diagnostics.
Price selection gamma=.25/7d and shape selection gamma=1/7d were frozen before
full SPX repricing (all seven dates) and32 fresh joint37D VIX scrambles.
The strong case improves fresh measured wing slope at all four dates and lowers
the descriptive slope/curvature error0.44495->0.30875, but remains too flat.
Fresh native VIX+forward loss26.8525->29.7000; SPX loss43.7775->121.1470.
The mild case's training joint-loss reduction70.9777->69.8756 is smaller than
its paired jackknife sensitivity (change -1.1021, SE3.9335); fresh VIX price loss
does not improve. No accepted joint fit or new model adoption follows.
Five unit checks and the saved-price/source/seed/independent-integral audit pass.
Signed estimates and missing IVs are retained. No population convergence,
market OOS, new P calibration or Lean certificate is claimed.
Keep `m2_maturity_curve_fit_001/selected_model.json`; no experiment remains running.
The earlier checkpoints below are historical, not an instruction to rerun them.

**Implementation update, 25 October:** the dedicated scalar OU proposal below
has now been implemented and tested separately in `runs/m2_stochastic_vov_001`.
See `STOCHASTIC_VOV_PLAN.md` and `STOCHASTIC_VOV_RESULTS.md`. The historical
proposal sections below record the earlier state of work. No production baseline
was replaced. The selected added state improves wing IV levels but flattens fresh
wing slopes relative to each cap reference and worsens the tested joint fit.

25 October 2026. This distinguishes published results, proposed M2 adaptations and actual fixed-curve experiments. No jumps, quintic proposal or modification to the existing OU rates/geometric grids. See CAP_TAIL_RESULTS.md for numerical outcomes and CAP_TAIL_MATHEMATICS.md for baseline equations and stock-martingale arguments.

## What stochastic vol-of-vol means here

If v=F(t,X), Ito's rule makes the Brownian part of variance

    dv_t = ... dt + grad_X F(t,X_t) B dW_t.

Its instantaneous random-motion strength is ||B^T grad_X F||. Thus the current nonlinear M2 already has state-dependent instantaneous vol-of-vol. The new question is whether to give that response extra independently adjustable variation, rather than merely making variance bigger. A cluster state can alter the spike/feedback sensitivities; a dedicated latent state can alter their strength and persistence independently of the cluster's current role. Either changes the distribution of conditional future variance, hence potentially the VIX smile.

## Primary literature checked live

[Fouque and Saporito, Heston Stochastic Vol-of-Vol Model for Joint Calibration of VIX and S&P 500 Options, arXiv1706.00873v1](https://arxiv.org/html/1706.00873v1): their variance diffusion coefficient is eta(Y_fast,Z_slow)*sqrt(V). Fast/slow continuous factors govern eta. Under stated regularity and scale assumptions they derive first-order joint pricing approximations around Heston. This motivates modulating response strength with another state. Their perturbation formulas and accuracy theorem do not apply directly to M2's normalized Gaussian readout.

[Alos, Garcia-Lorite and Muguruza, On smile properties of volatility derivatives and exotic products: understanding the VIX skew, arXiv1808.03610v1](https://arxiv.org/html/1808.03610v1): they derive short-time ATM level/skew expressions using Malliavin derivatives. Their mixed generalized rough-volatility example adds two exponential components with different sensitivities; they discuss mixing as a way to overcome flat VIX smiles. This supports examining sensitivity mixtures. It does not prove any particular M2 cap or cluster modulation will fit finite-strike, finite-maturity wings.

[Jacquier, Muguruza and Pannier, Rough multifactor volatility for SPX and VIX options, arXiv2112.14310v2](https://arxiv.org/html/2112.14310v2): they study ATM level, skew and curvature, including a two-factor rough Bergomi model. Their calibration discussion separates the sensitivity/mixing parameters used for VIX from equity correlations used for SPX. The relevant lesson is to assess shape and cross-market coupling separately. We preserve M2 correlations; these asymptotic equations are not M2 pricing formulas.

Additional primary empirical evidence: [Romer, Empirical analysis of rough and classical stochastic volatility models to the SPX and VIX markets, Quantitative Finance2022](https://researchprofiles.ku.dk/en/publications/empirical-analysis-of-rough-and-classical-stochastic-volatility-m/). The author's university record describes a two-OU model with a nonstandard transformation and joint calibration results without jumps or roughness. This supports continuous transformation/multifactor exploration. It is an empirical result for another model and market exercise, not evidence that the present frozen reservoir is sufficient.

## Existing-cluster experiment actually implemented

Let eta0 be the original score and eta1=eta0+b*qJ*(F-betaJ)^2. Set

    w(C)=0.5+0.45*tanh(0.5*C),
    G_b=(1-w(C))*exp(2h(eta0)) + w(C)*exp(2h(eta1)),
    M_b(t)=E_Q[G_b(t)],
    v_b=f+(xi-f)*G_b/M_b.

As C rises, more weight moves to a component with stronger spike sensitivity. b=.5 and1 were profiled, holding curve/floor and all original coefficients fixed. b controls added spike response, the tanh scale controls how quickly cluster states alter that response, and the0.45 amplitude limits mixture weights to[.05,.95]. Those two tanh choices were frozen diagnostic design choices, not calibrated values. b=0 recovers the original model exactly. Two spike components can respond differently before the common cap; after both saturate, modulation cannot restore unlimited tail response. This experiment therefore tests finite-state response redistribution, not elimination of the ceiling.

The mixture remains positive and bounded by exp(2L); Gaussian moments and the bounded-volatility martingale argument remain available. No new random variables or OU parameters are introduced. Existing full outer states, conditional Gaussian laws, stock C/Y/F features, equity increments, likelihood weights and frozen zero-mean controls remain valid primitives. G, Q normalizers, VIX exposures, volatility and stock paths must all be recomputed. Actual cluster experiment reprices all four VIX dates only; it has no claimed SPX tradeoff score or joint calibration. It is separate from the two cap candidates receiving full SPX checks.

The scalar cap convolution is insufficient for state-dependent mixture weights. The experimental kernel adds a second interpolation coordinate for conditional C means (spacing.05), with Hermite score interpolation (spacing.02); direct same-GH checks and independent Python tests assess it. It changes no model equation. Measured full VIX evaluation times include this extra integration cost; compare saved timings with the roughly8-second baseline profile.

## Dedicated continuous latent layer: explicit proposal, not implemented

A minimal additional state under Q is

    dZ_t = -a_Z*Z_t dt + nu_Z dB_Z(t), Z_0=0,
    A(Z)=exp(s*tanh(Z)),
    eta_Z=c_C*C+qY*(Y-betaY)^2+qJ*A(Z)*(F-betaJ)^2,
    G_Z=exp(2h(eta_Z)), M_Z(t)=E_Q[G_Z(t)],
    v_Z=f+(xi-f)*G_Z/M_Z(t).

This adds one state to the unchanged reservoir. a_Z sets persistence, nu_Z controls how far Z moves, s sets the range of spike sensitivities, qJ*A in[qJ*exp(-s),qJ*exp(s)]. To avoid poor identification between nu_Z and s, start with Var_stationary(Z)=1, so nu_Z=sqrt(2a_Z), and fit only a_Z and s. s=0 exactly nests baseline. All states and stock paths remain continuous. A parallel feedback modulation is possible, but adding both simultaneously would obscure which mechanism helped.

For independent B_Z, existing OU/equity draws and their importance likelihoods can be reused unchanged; simulate independent Z values at the SAME stock times and outer VIX times. The new Z likelihood is1 when sampled directly under its Q law. With nonzero correlation, shared driving randomness and covariance/conditional laws must be reconstructed consistently; assigning arbitrary correlation to saved paths afterwards is invalid. Begin independent to isolate the new response mechanism. Correlation is not an authorized baseline change.

Conditional pricing must condition on Z_T as well as the original full OU state. Z_u conditional on Z_T is Gaussian with mean exp(-a_Z*(u-T))*Z_T and variance nu_Z^2*(1-exp(-2a_Z*(u-T)))/(2a_Z). Baseline three-projection marginal/conditional integrals become FOUR-dimensional, although C or Z may be integrated using reusable tables. A direct tensor rule adds approximately one quadrature-order factor of cost (e.g.32x GH nodes conditionally); this is an operation count, not a measured runtime. Outer importance sampling may need to target joint high-response states even with4096 points; existing proposal efficacy cannot be assumed.

Stock cache contains about3.496 billion step/path entries across all seven maturities. Caching one extra float64 Z value per entry adds26.05GiB; streaming Z avoids that storage but requires reproducible seeds/bridges and full stock recomputation for every changed response. New randomness does not require resimulating the unchanged exogenous reservoir, but neither volatility nor terminal stocks are reusable unchanged. Reference-control and hedge/Gaussian zero means survive independent augmentation, although variance reduction may deteriorate. Rebuild candidate-specific VIX exposures/forwards/option payoffs/normalizers and check new sampling, GH, time and stock errors.

G_Z remains bounded by exp(2L), so integrability and the bounded-volatility stock-martingale argument survive this augmentation. Simply using exp(Z)*(F-betaJ)^2 WITHOUT tempering could produce nonintegrable Gaussian/lognormal mixtures; the bounded sensitivity and outer cap avoid that issue. A future log-continuation adaptation can also preserve at-most-quadratic growth when sensitivity A is bounded.

## A second dedicated-layer option

A smooth mixture weight w(Z) can replace w(C) in the implemented two-component mixture. This allows changing mixture persistence independently of the cluster reservoir. It keeps bounded G, nests a constant mixture at zero modulation, and needs the same added state/pricing work. It differs from direct qJ*A(Z): it interpolates between two response functions instead of continuously changing the curvature. Compare these two mechanisms at matched mean variance and forwards, not by assigning the label 'stochastic vol-of-vol' and assuming success.

No dedicated-state simulation or calibration has been launched. The existing-cluster diagnostic is a bounded test, not an adopted architecture. A poor cluster profile does not reject the dedicated layer, whose new persistence and response dependence differ.
