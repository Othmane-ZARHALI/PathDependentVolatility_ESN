# All retained SPX maturities: in-sample joint fit

Status: **complete**, 23 October 2026. This checkpoint implements the user's
explicit instruction to include all retained SPX maturities in recalibration
and pictures, with no market out-of-sample results. The old parameters are only
the baseline; the variance curve is recalibrated on the expanded objective.

## Scope and interpretation

Use SPX **29, 57, 119, 211, 302, 666 and 1030 days**. Exclude only **1, 1401,
1765 and 2129 days**. All eligible original quote splits are calibration inputs
in this feasibility run. Cutoff 3 is unchanged:
`abs(log(K/F))/(sigma_ATM_proxy*sqrt(T)) <= 3`. The proxy is the nearest-to-ATM
midquote Black IV using all eligible quotes at each expiry, frozen before model
pricing. Quotes remain the existing OTM contracts; no call/put pooling occurs.

| SPX maturity, days | Eligible OTM quotes | Cutoff-3 fitting quotes | Current role |
|---:|---:|---:|---|
| 29 | 410 | 277 | In-sample fitting |
| 57 | 399 | 298 | In-sample fitting |
| 119 | 314 | 228 | In-sample fitting |
| 211 | 204 | 141 | In-sample fitting |
| 302 | 209 | 167 | In-sample fitting |
| 666 | 133 | 112 | In-sample fitting |
| 1030 | 82 | 71 | In-sample fitting |
| Total | 1751 | **1294** | All seven enter optimization |

Keep all four VIX maturities 27/55/118/153 days, all **217** eligible VIX OTM
quotes, and all four observed market forwards. There are **1511** option
targets in the expanded joint objective. The earlier 707-SPX/175-VIX totals
describe a different historical objective and are not the current fit scope.

The copied calibration tables explicitly say `split=in_sample_fit` and retain
the old source labels as `original_split` and `original_maturity_role`.
Raw/preprocessed market data and old numerical runs are unchanged. There is no
reserved market maturity or strike subset in this run. The older raw-parity
quote-exposure correction remains historical evidence, but the current run
deliberately uses the whole retained market population, not a blind split.

## What is recalibrated

Start from the prior provisional C0=1/Y0=F0=0 working reference, variance levels
`[0.0584806825,0.0727920966,0.0357009546]`. Reprice it at every requested expiry,
then optimize all three curve levels on **all seven SPX and all four VIX
maturities together**. The new maturities and all their selected quotes enter
the optimizer, not just the display. The curve retains boundaries 91.25 and
273.75 days and the original 2129-day horizon. The long final segment is now
constrained by the 666/1030-day data as well as day 302.

The initial state, readout, structure and Q scenario 0 remain fixed in this
bounded pass. Thus it tests this declared curve/state family on the expanded
data; a poor result cannot prove that every parameterization of M2 fails.
The Q-risk-price comparison was deferred during this scope-expansion pass and
is the next bounded fitting hypothesis on the established expanded baseline.
The unchanged fixed-price-scale loss balances
each expiry within SPX/VIX and includes forward errors separately. More quotes
do not automatically give SPX more total family weight than VIX.

The [prespecified plan](ALL_MATURITY_FIT_PLAN.md) records three starts, at most
30 evaluations each, a 30-minute optimization budget after preparation, and
512 paths per expiry. Compatible old state-specific ratios/equity draws and VIX
exposures are reused. Added strikes use frozen nearest same-side zero-mean
controls. New SPX maturities receive new primitives and the original five-control
family, using frozen 302-day pilot coefficients with forward scaling. Valid
zero means justify transfer; efficient variance reduction at new maturities is
not established. No control-feature search or new beta fitting is performed.

A separate numerical ensemble checks the selected parameters against the SAME
1511 market observations. This is numerical stability checking, **not market
out-of-sample performance**. Prior numerical ensembles are reused at the old
dates; added maturities have distinct new searching/checking streams. The whole
exercise is explicitly in sample. Numerical checks do not create unseen data.

## Numerical and formal limits

Keep 512 stock steps per trading day. Exact OU updates do not certify nonlinear
stock timestep accuracy. New maturity grids preserve exact dates/curve knots
and have their own event placement. Normalizer GH48 and the declared adaptive
interpolation screen, VIX conditional GH32 and time order 8 remain separate
approximations. Sampling bands exclude all those deterministic errors and pilot
uncertainty. Signed estimates and missing IVs remain visible, without clipping.

All **182 Python tests** pass: 124 package and 58 drivers, including three new
checks for seven-maturity/all-split inclusion, exact same-side donor mapping,
and agreement with the production objective under the explicit in-sample label.
An initial new test fixture incorrectly accessed a dictionary as an object;
the fixture was corrected without changing the objective. The failed output
is retained. Current Lean identities apply under the same assumptions; new
pinned compiler/axiom audit **006 passes all thirteen**, warnings as errors,
with only the declared standard axioms and no sorryAx/custom axioms.
No claim of formally verified Python
selection, interpolation, calibration success or price convergence follows.

The new `active_scope.json` deliberately supersedes the historical scope. Its
old bytes and hash are saved in `scope_decisions/2026-10-23-before-all-maturity-fit.json`.
Old strict current-file replay scripts see that intentional scope-version
change; historical source snapshots/audits remain valid for their captured
versions and are not rewritten or silently rehashed.

## Results and pictures

`runs/m2_all_maturity_fit_001` completed in 1131.62 seconds. All three starts
converged in five reported optimizer evaluations each, to essentially the same
curve; total optimization wall time was 350.20 seconds. State/readout/Q stayed
fixed. The selected curve froze before numerical-check scores were evaluated.

| Future interval, days | Previous annual variance | Recalibrated annual variance | Previous square-root scale | New square-root scale |
|---|---:|---:|---:|---:|
| 0–91.25 | 0.05848068 | **0.05762055** | 24.1828% | **24.0043%** |
| 91.25–273.75 | 0.07279210 | **0.07134085** | 26.9800% | **26.7097%** |
| 273.75–2129 | 0.03570095 | **0.06321057** | 18.8947% | **25.1417%** |

The selected exact levels are `[0.057620550549099724,0.07134085119966904,
0.06321057443478936]`. These are variance inputs, not option IVs. The large
last-piece change is the practical consequence of including the long maturities:
the prior curve was only weakly constrained after day 273.75 by day-302 SPX.
That earlier long level should not have been treated as a calibrated prediction
for 666/1030 days. The new fit uses those observations directly, in sample.

| Numerical ensemble, SAME 1511 market observations | Previous parameters | Recalibrated | Loss reduction | Prices inside bid/ask |
|---|---:|---:|---:|---|
| Searching simulations | 300.9412 | **253.4447** | **15.78%** | 63 → 106 |
| Separate numerical simulations | 351.3587 | **308.6356** | **12.16%** | 66 → 103 |

Neither row is market out-of-sample evaluation. Both score exactly the
observations used in calibration. The second numerical ensemble was not used
to choose the curve; its original-date paths were reused from prior numerical
checks, so it is not a new independent validation of the entire research
selection process. The absolute loss is not comparable with the earlier
four-SPX-maturity objective, which used different quotes and expiry weights.

On the separate numerical ensemble the components are:

| Component | Previous parameters | Recalibrated |
|---|---:|---:|
| SPX | 269.4962 | **235.2859** |
| VIX options | 70.4838 | **64.4466** |
| VIX forwards | 11.3452 | **8.8697** |
| Curve penalty | 0.000158 | 0.000025 |
| State penalty, fixed | 0.033333 | 0.033333 |

The descriptive eight-block jackknife SE of the loss change -42.7231 is
14.4435; leave-one-block changes range from -53.2509 to -32.9562. Few blocks,
tail uncertainty, prior ensemble reuse and unresolved deterministic bias mean
this is not a calibrated confidence certificate. The signed change remains
favorable on both numerical ensembles without changing any market observations.

### All seven SPX expiry outcomes

Here the score is mean squared price-interval distance in the frozen quote
scales; lower is better. These are numerical-check estimates on in-sample data.

| SPX days | Fitting quotes | Previous interval MSE | Recalibrated interval MSE | Recalibrated prices inside bid/ask |
|---:|---:|---:|---:|---:|
| 29 | 277 | 202.1838 | 184.8975 | 10 |
| 57 | 298 | 292.5458 | 301.7972 | 6 |
| 119 | 228 | 427.1884 | 457.3828 | 5 |
| 211 | 141 | 269.0271 | 290.4567 | 5 |
| 302 | 167 | 379.3880 | 389.7611 | 3 |
| 666 | 112 | 183.8314 | **20.2998** | 17 |
| 1030 | 71 | 132.3087 | **2.4061** | 34 |

All four VIX option-expiry losses improve on the separate numerical ensemble.
Their recalibrated interval MSEs are 153.2096/69.0089/25.0285/10.5394, with
11/3/5/4 interval hits out of 49/58/55/55 quotes. Forward errors remain
**+1.4437/+0.5354/+2.0525/+1.5120 VIX points**. The market forwards remain
the observations to be matched; they were neither replaced nor adjusted to
improve the picture. They enter as fitting targets, not exact constraints.

### What the pictures say, in undergraduate terms

The [seven-panel SPX picture](visualizations/all-maturity-fit-007/spx-all-maturity-fit.png)
now shows actual model prices at **every requested maturity**. Blue is the
expanded-scope recalibration, orange is the previous working parameters, and
black/gray is the market. The [VIX picture](visualizations/all-maturity-fit-007/vix-all-maturity-fit.png)
and [forward/variance picture](visualizations/all-maturity-fit-007/forwards-and-variance.png)
use the same conventions. Every panel is explicitly in-sample fit. Market and
model IVs both use observed market forwards; no model-own-forward display
transformation is applied here. IV is a display of the priced options, not
the calibrated loss. PNGs were visually inspected; zoomable SVGs, input/output
hashes, CSV scores and a [reading guide](visualizations/all-maturity-fit-007/README.md)
are saved. The searching-sample preview is retained separately and labelled.

The long-maturity SPX level is substantially improved. Shorter maturities still
show visible shape errors: changing three average-variance levels cannot freely
rotate or reshape every smile. The VIX call/put display gaps and right-wing
shape errors also remain on the common observed-forward IV scale. Small kinks
in model curves can come from finite-sample noise and changes in nearest-strike
control reuse; they do not establish discontinuities of the model's expected
option-price function. No cosmetic smoothing or clipping was performed.

Only **103/1511 (6.82%)** numerical-check prices lie inside bid/ask. All 1294
SPX nominal precision screens fail; 26/217 VIX screens pass. All four score
records contain no negative primary estimates. Sampling bands are conditional
on fixed controls and exclude the separate stock-grid/normalizer/GH/time
approximations. New long-maturity same-GH interpolation probes have maximum
relative discrepancies 9.8769e-8 (independent) and 9.9203e-8 (midpoints), below
the declared sampled screen; this does not establish a uniform bound or GH
convergence. Existing preprocessing marks the retained finite-strike quote
interval grids feasible, which is not a proof of joint M2 feasibility.

### Audit correction and return-to-fitting decision

All captured numerical inputs/sources, all seven maturity memberships, all four
score records and sample means/SEs replay independently. The 24 new seed records
are distinct and disjoint from recorded parent seeds; this checks metadata,
not independent regeneration of every random stream.
The 882 contracts common with the old scope reproduce previous searching
prices to maximum absolute discrepancy **2.8422e-13**. The new audit's first
reader incorrectly assumed the archived VIX array contained controlled option
values; it actually contains raw simulated VIX. The failed audit log is retained.
The corrected reader independently reconstructs controlled forward/option
payoffs from the archived VIX values, fixed coefficients and likelihoods and
then checks their means and sampling SE. No simulation, target, fitted parameter
or optimization result was changed or rerun to repair this reader.

The explicit runner capture list did not include every imported Black-conversion,
loader or plot helper. Final validation also archives these dependencies,
records that the tracked ones match HEAD with no working-tree edits, and
records the shared plotting helper as already untracked and unchanged during
this task. Population and exact old-contract replay support the correspondence.
This is post-run dependency capture; it is not represented as an original
pre-run hash, and the run's captured-source manifest is not rewritten.

**Keep this recalibrated three-level curve as the current all-maturity working
reference at C0=1/Y0=F0=0 and scenario 0.** Do not return to four-maturity or
training-only fitting. Every future feasibility comparison should use the same
1294-SPX/217-VIX/four-forward population unless the user changes it. The old
reference remains a historical comparison, not the current calibrated curve.
The requested all-maturity recalibration and pictures are complete.

This is a better restricted-family in-sample fit, **not a satisfactory joint
market fit and not evidence that every M2 parameterization fails**. Full model
capacity remains unresolved because structure/readout/state/risk prices were
held fixed in this pass. The next bounded fitting hypothesis remains the
declared risk-price scenarios, one at a time, now on this full in-sample scope
and starting from this recalibrated curve. Rebuild candidate normalization and
payoff/exposure caches, assess distribution-shape changes and preserve separate
sampling/deterministic errors. That next parameter search is not launched in
this checkpoint; no market out-of-sample requirement is imposed.

### 24 October: parameter scope and strike-jitter diagnosis

The latest optimizer has exactly **three continuous unknowns**, the levels of
the deterministic forward-variance curve. It does not recalibrate the five
readout shape coefficients, their cap, driver correlations, risk prices or
initial state. Calling the plots a full structural calibration would overstate
the work completed; they show a restricted joint SPX/VIX curve calibration.

A cheap decomposition of saved checking prices at 119/211/666 days confirms
that the strike-dependent control correction can introduce price-order errors.
The estimator is `mean(weight * payoff) - mean(controls * beta)`; control
selection/coefficient transfer changes with strike. The positive-weight raw
payoff estimate is monotone across strikes within each option side, but the
corrected estimate need not preserve that finite-sample property.

| SPX days | Adjacent-strike monotonicity violations, corrected | Raw weighted payoff estimate | Largest offending put-strike pair | Corrected price change | Raw price change | Correction change |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| 119 | 3 | 0 | 5750 to 5775 | -12.5104 | +1.7746 | +14.2850 |
| 211 | 3 | 0 | 5375 to 5400 | -4.3695 | +1.5905 | +5.9600 |
| 666 | 2 | 0 | 6150 to 6200 | -17.0773 | +9.8782 | +26.9555 |

Changes are in undiscounted SPX option-price units. A put payoff increases when
its strike increases, path by path. The violations therefore identify a
finite-sample estimator artifact, not a failure of that mathematical property
for M2. Both plotted parameter sets share the correction, which is frozen
across curve candidates. IV inversion can magnify price artifacts where vega
is small. Three pieces in a curve indexed by time do not themselves create
these jumps along the strike axis at a fixed expiry.

Full reconstructed rows, input hashes, the mapping recipe, all eight offending
pairs and limitations are retained in
`runs/m2_all_maturity_fit_001/validation/strike-jitter-decomposition.csv` and
`strike-jitter-decomposition.json`. This recomputes corrections from saved
arrays and reconstructs raw means from saved controlled prices; it does not
independently regenerate paths. No simulation, optimization, fitted price,
original plot or original source manifest was changed. No new Lean theorem or
accuracy certificate is claimed. Sampling, stock-timestep bias and
normalizer/quadrature errors remain separate. A future estimator revision
should preserve continuous behavior across strikes before visual improvements
are interpreted as better model shape; smoothing the plotted IVs is not a
pricing repair. These findings do not by themselves resolve the broader smile
misfit or establish convergence of the raw estimator.
