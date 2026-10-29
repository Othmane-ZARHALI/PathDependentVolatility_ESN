## 28 October 2026 — protected two-speed SPX precision phase

User authorized improving SPX precision with the successful two-speed VIX-only model fixed, then a protected fit of only nine off-window xi levels, with stress/equity correlation a later possible experiment. All model parameters remain unchanged in this phase: no joint candidate adoption, correlation, factor, readout or floor change. All1751eligible SPX quotes repriced,1294primary at all seven retained maturities; all217VIX quotes/four forwards remain protected. See [TWO_SPEED_PRECISION_RESULTS.md](TWO_SPEED_PRECISION_RESULTS.md), [TWO_SPEED_PRECISION_PLAN.md](TWO_SPEED_PRECISION_PLAN.md), `runs/m2_two_speed_precision_002`, `runs/m2_two_speed_precision_conditional_002`, and `visualizations/two-speed-precision027`.

Implemented an independent-stress auxiliary stock, likelihood-weighted differences, terminal/dynamic hedge controls, then an exact discrete one-step conditional Black-value martingale correction and independently scrambled OU-bridge reference integration. Checking median SE/scale improves about8–29x, from4.68/15.03/12.17/12.76/13.08/5.74/5.01 to0.579/0.550/0.903/1.153/1.239/0.435/0.174. Only1/1294primary quotes passes the existing sampling screen. The nine-level price optimizer was therefore NOT launched. Broad frozen-model SPX smile-shape mismatch remains visible; it is not a model-family infeasibility theorem.

Remaining bottleneck is the512-path full-M2 correction. More cheap stress paths alone cannot clear the gate. Pilot equal-component budget projections give median expensive counts about1.8k–93k across maturities and worst quote about1.2million, excluding timestep/normalizer errors; these are not guaranteed budgets. Next implement streaming fresh full-M2 corrections with adaptive/multilevel allocation and coupled stock-grid checks, rather than expanding the full cache blindly. Then free ONLY indices0,5,6,9,11,12,13,14,15; keep floor, both stresses and VIX-sensitive indices1,2,3,4,7,8,10 fixed. The exact VIX parameter map is invariant to those nine changes. Existing forward accuracy is preserved, not upgraded to population equality. Independent target/reference/conditional/RNG/covariance/full-horizon derivative and finite-law control checks PASS. Checking streams reused; not market OOS, full convergence clearance or a Lean proof. Current chat can continue; these notes are the restart handoff if a new context is chosen.

## 28 October 2026 — two-speed joint SPX/VIX refit completed

User resumed SPX and chose the selected two-speed stress model. Restored seven SPX maturities29/57/119/211/302/666/1030, cutoff3(1294primary;1751eligible including diagnostic wings), all217VIXquotes/fourforwards. Freed all16xi levels, floor and six existing stress parameters(23total); original M2 readout/premia/dynamics fixed. See [TWO_SPEED_JOINT_RESULTS.md](TWO_SPEED_JOINT_RESULTS.md), [TWO_SPEED_JOINT_PLAN.md](TWO_SPEED_JOINT_PLAN.md), `runs/m2_two_speed_joint_001`, `visualizations/two-speed-joint026`. Checking SPX loss5467.10→618.74, VIX0.05008→37.9711. Joint stress weight collapses0.92284→0.01943; slow branch effectively inactive. Forward calibration gap1.24e-8, checking max0.01782. Unsatisfactory overall: VIX fit badly worsens, all SPX sampling gates fail. Coupled64/128/512checks retain broad trade-off and unresolved timestep/sampling error. Source/input/56seed/23derivative/NumPy checks pass. This is one bounded local staged fit, both stages iteration-limited, with old readout frozen; no structural impossibility claim. Two-speed VIX-only and joint candidates preserved separately; no production baseline replacement. SPX scope is active for research; do not silently resume the earlier pause. Improve SPX numerical precision before further calibration conclusions.

## 28 October 2026 — fixed geometric four-speed stress bank finalized

Authorized bank experiment completed with half-lives7,28,112,448days fixed before fitting (ratio4, calendar-based endpoints). Four independent normalized stress responses plus frozen M2; globally shared positive weights and amplitudes. All217nativequotes/fourforwards, seven free active xi levels and fitted floor retained. Checking option loss 0.050081→0.062160; forward gaps/SE and all-maturity wings are in [STRESS_BANK_RESULTS.md](STRESS_BANK_RESULTS.md). See [STRESS_BANK_PLAN.md](STRESS_BANK_PLAN.md), `runs/m2_stress_bank_002`, `visualizations/stress-bank025`.16 fitted parameters versus14 previously; this is not a nested comparison. Cached terminal coordinates reused with new additional nodes; checking partly reused and not market OOS. Source/likelihood/scalar audits PASS; no population-forward certificate, SPX work, baseline replacement or Lean theorem.

## 27 October 2026 — shared fast/slow stress, all maturities jointly

The user explicitly requested one common model rather than special treatment of 27 days. The proposed isolated first-expiry diagnostic is superseded as the next fitting step. Completed the joint two-speed pilot: see [TWO_SPEED_RESULTS.md](TWO_SPEED_RESULTS.md), [TWO_SPEED_PLAN.md](TWO_SPEED_PLAN.md), `runs/m2_two_speed_001`, and `visualizations/two-speed024`. Checking option loss 1.725555 → 0.050081; all 217 quotes/four forwards retained. Four training forward equality constraints, seven free active xi levels and fitted floor. Checking gaps and SEs are reported, not declared exact. Reused old/slow checking states plus new fast normals; not wholly fresh. Previous continuous Beta measure not simultaneously refitted; two-speed extension tests the empirically motivated simpler mixture. Numerical audit PASS; SPX paused and both baselines unchanged.

## 27 October 2026 — teaching notes and first VIX maturity diagnosis

Saved CONTINUUM_EMPIRICAL_REDUCTION.md: equations, Bergomi attribution, the Beta-to-two-atom limit, bounded empirical reduction and why the general integral remains useful for future market snapshots. The reduction concerns mixture/readout complexity, not stochastic state dimension; only one stored market snapshot was tested.

Saved VIX_27_DAY_DIAGNOSIS.md with reproducible read-only diagnostics in runs/m2_continuum_reduction_note_001. The27day right wing is too low while central calls are too high. In the matching boundary reference, the strike100 payoff-weighted stress variance share is98.71%; stress dominance and shared response parameters are a leading shape mechanism. Exact-forward necessary second moment468.99933 exceeds selected467.44959. A27day-only constrained calibration is a proposed discriminating experiment, not an authorized scope change or a completed fit. All217quotes/fourforwards, the working candidate and SPX pause are unchanged.

## 27 October 2026 — continuous measure, forward constraints and multilevel repair

Completed VIX-only continuous M2/stress Beta-measure experiment003, after rejecting002's noisy55-day forward calibration. The model is unchanged: normalized G^(1-u) times the normalized OU stress exponential of strength b*u, integrated against one Beta measure. Four forward equalities use32768-path boundary plus8192-path paired difference; checking uses65536 NEW boundary paths plus16384 reused difference paths. Do not claim wholly fresh validation or exact population equality.

Best bounded fit remains near kappa=.01; a broader measure is retained as a training-selected diagnostic, with b=5 at its bound. No clear improvement established over the two-atom reference. All217 quotes, four maturities, seven free active xi levels and fitted floor preserved. SPX paused; original baseline and previous working VIX candidate unchanged. See CONTINUUM_RESPONSE_RESULTS.md, runs/m2_continuum_response_003 and visualizations/continuum-response023.

# Pricing and fitting research ledger

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

**26 October: VIX-only free-curve comparison completed.**
See [VIX_ONLY_CURVE_RESULTS.md](VIX_ONLY_CURVE_RESULTS.md),
[VIX_ONLY_CURVE_PLAN.md](VIX_ONLY_CURVE_PLAN.md),
`runs/m2_vix_only_curve_fit_002`, figures `visualizations/vix-only-curve019`.
The original curve has **16 levels**, already fitted in the earlier maturity-curve
study; recent factor trials had frozen them. Seven pieces affect the four VIX
averaging windows and now vary independently. The nine off-window pieces retain
their old values because these VIX targets have zero sensitivity to them.
No reversion to three levels or four broad scalers. SPX remains paused.
Fresh native VIX+forward losses: baseline27.3491, curve-only20.5245,
fixed strong factor+curve20.9500, fitted factor+curve19.8086. Additional fitted-factor
gain versus curve-only is -0.7159 with paired jackknife SE0.5476. Wings remain too
flat; one short curve piece hits its0.20 bound and fitted half-life hits365days.
All217native quotes and four forwards retained; no smoothing/clipping. Only21quotes
inside bid–ask for the joint fit; no acceptable full-smile calibration or model
adoption claimed. Baseline preserved. Eight candidate/stream price replays,
nine Jacobian directions, exact budgets,282source/fiveinput hashes pass; sampled
GH/time refinement differs by less than4.6e-6VIXpoints. Numerical checks use the
same market quotes. No SPX pricing, new bank/readout or new Lean theorem.
Historical checkpoints below retain their original context.

**26 October: current phase is VIX-only research and calibration planning.**
SPX fitting/repricing is temporarily deferred by explicit user direction.
Read the new section at the top of `STOCHASTIC_VOL_OF_VOL_RESEARCH.md` and
`runs/m2_vix_design_research_001`. All three Othmane prepared PDFs and the14October
email figures are readable; no original-PDF request is needed. The later cs_doc
still reports too-flat2–5month VIX wings, fits a narrower strike range and uses
an approximate VIX estimator with a nested check pending. Its Gamma modulates
bank noise, unlike our final readout multiplier; no matching-code repricing is claimed.

New model-independent finding: exact-rational static-payoff certificates show
that the frozen xi curve's E[VIX²] budget is too low to fit all native VIX
bid–ask quotes and forwards at55/118/153days. This persists with forward tolerance
plus/minus0.5points. `moment_certificate_audit.json` passes12 global payoff/price
checks; it is exact computation for archived decimal data, not Lean or a market
arbitrage claim. Necessary exact-forward lower bounds are556.673/648.626/695.295,
versus old budgets536.027/580.804/610.933. The27day bound does not exclude baseline.

Next comparison should allow VIX-window xi levels and shape parameters to move
together. First re-evaluate the current after-cap layer fairly; then prefer a
mixed-response readout/independent mixed-standard-Bergomi benchmark. A continuous
Gamma-driven new bank is the more direct Othmane analogue but needs a new
conditional pricer. No new reservoir calibration, SPX run, jump adoption or
selected-model replacement has been performed in this research turn.
All217VIX quotes and four forwards remain; old SPX membership/cutoff3 are retained
for the later return. Earlier checkpoints below remain historical.

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

**COMPLETED — dedicated stochastic spike-sensitivity trial, 25 October:**
Read `STOCHASTIC_VOV_RESULTS.md` and `STOCHASTIC_VOV_PLAN.md`.
Run `runs/m2_stochastic_vov_001`; figures `visualizations/stochastic-vov016`.
Selected baseline remains `m2_maturity_curve_fit_001/selected_model.json`.

One added independent Q OU state modulates qJ by exp(s*tanh(Z)); Z0=0,
stationary variance1. Strength ranges 1.5x/2x and half-lives7/30/90 days,
crossed with original cap and logarithmic continuation a=.5. All original
OU rates/grids/state/correlations/premiums/readout coefficients, all16 xi
levels and floor .00264 are exactly fixed. No curve refit, jumps or quintic.
All four VIX dates/four forwards and4096 samples per strike remain. All12
layer cases and both zero-strength references have full VIX diagnostics.

Wing-level training selection picks2x/90d for BOTH caps, before new SPX
results or checking. Full actual SPX repricing covers all seven dates for
both selections and both cap references (four full comparisons,1511 OTM
quotes each, cutoff3/common-call/native model-forward parity preserved).
No full-SPX search over the entire12-point grid is claimed.

Main finding: closer wing IV LEVELS do not mean greater STEEPNESS. On NEW
joint37D scrambles, adding the selected state FLATTENS measured right-wing
slope at ALL FOUR VIX maturities relative to each same-cap no-layer reference.
New-stream wing MAE33.841 ->15.072/11.733, but near-ATM12.670 ->18.994/18.555.
Native VIX+forward loss27.597 ->272.754/482.047. Training SPX loss43.777
->561.511/1376.764; no-layer slow-growth reference112.992. Paired SPX changes
are positive in every leave-one-block comparison for both selected cases.
The frozen-curve grid does not give the desired joint improvement. Keep baseline;
this does not reject a different interaction or a recalibrated stochastic model.

Normalizer M at153d rises2.408x/3.074x for selected layers vs same-cap refs;
stronger raw response plus stronger normalization redistributes variance.
Targets E[VIX^2] stay identical because xi is fixed. New153d baseline moment
606.410 vs610.933, SE2.466; selected original-cap state608.135, SE12.257.
These are numerical diagnostics, not a convergence pass or resolution of the
historical second-moment discrepancy. No fresh SPX ensemble is claimed.

Three tests pass; audit replays22 price/SE/parity/loss records,270 archived
sources and all32 NEW complete joint mixture scrambles. Eight independent
4D oracles have max relative error5.63e-10. Representative normalizer GH48->96
changes <=5.25e-8, extra-state GH and time refinements are recorded separately.
Original GH32/GH48 and time8/stock512 mesh retained; extra Z usesGH128/GH192.
Exact OU stock updates stream reproducible shared innovations without an
extra all-time state cache. OU primitives reused; volatility/stock recomputed.
Nine signed SPX302 wing estimates are negative, all <.82 reported SE below0.
Separate VIX IV failures reflect market-forward inversion; no clipping/smoothing.
Eleven inspected figures have22 PNG/SVG exports. Baseline full price replay
max difference1.37e-9. No new accepted fit, Lean/P/global-optimum certificate.
All numerical jobs finished; finalization includes audit/environment/hashes.

**COMPLETED — cleaner cap comparison, 25 October:**
Read `CAP_CONTROLLED_RESULTS.md` and `CAP_CONTROLLED_PLAN.md`.
Run `runs/m2_cap_controlled_001`; figures `visualizations/cap-controlled015`.
Selected baseline stays `m2_maturity_curve_fit_001/selected_model.json`.
No new accepted joint fit or dedicated latent layer. Jumps/quintic excluded.

Floor .00264 and every off-window curve piece are EXACTLY fixed. Four broad
VIX-window adjustments average overlaps, with maximum15% excess-variance change;
selected changes are only about1--3%. Same OU rates/geometric grids/state/
correlations/Q scenario/score coefficients, beta_J>0,1511 original OTM quotes,
all seven SPX/four VIX maturities/four forwards, cutoff3/common calls and4096
paths per strike. No clipping/smoothing. Nine cap shapes: L increases0/.375/.75
crossed with logarithmic continuation a=0/.25/.5, k=2. Local fits hit the
common-call safeguard45.509602; several also hit +.5-point forward bounds.
Those screens are implementation choices; no family/global-optimum conclusion.

TRAINING selection frozen before NEW draws: shape03 (L+.375,a=0) and shape02
(L unchanged,a=.5). Wing IV MAE33.646->27.717/26.664; near-ATM12.389->14.921/15.180.
Actual full seven-maturity SPX loss43.777476->87.400641/98.795211. Paired SPX
loss-change SE5.959/6.734; all eight leave-one-block comparisons worsen.
Native joint70.977709->114.299387/125.738973. Keep the unchanged baseline.

Genuinely NEW8x512 VIX scrambles per date, seed base2026102511, same proposal/
frozen pilots/SAME in-sample market quotes. Native VIX+forward baseline26.986598;
shape03/02 give26.845187/26.914668: changes-.141411/-.071930, paired eight-net
SE.203733/.228777. Aggregate gains are not precise. New wing MAE33.866->
27.963/26.919; near-ATM12.590->15.183/15.435. No market OOS or fresh SPX check.
New baseline VIX153 moment608.632 vs610.933, SE2.294 (~1SE); historical deficit
is preserved, not cleared. Baseline139/1511 inside bid/ask vs94/87.

Three tests pass. Audit replays15 scores,139 fit trials,270 frozen-source
snapshots and32 complete new scramble/state/likelihood regenerations, with
seed disjointness against recorded historical pilots/production/diagnostics.
Nine PNG/SVG plots cover all agreed maturities, signed errors, forwards and tails.
Original fresh-cache setup failed BEFORE fresh pricing, after fits/SPX/draws
completed. Current driver is corrected. `resume_cap_controlled.py` retained all
fits/SPX/fresh draws and resumed only VIX pricing. Original failure/sources and
separate correction snapshots are preserved. All numerical jobs are finished.
These tested gentler brakes improve wings but do not give the desired joint fit.
Further continuous response-layer work remains proposed, not automatically run.

**COMPLETED — cap/tail profiles and continuous response diagnostics, 25 October:**
Read `CAP_TAIL_RESULTS.md`, `CAP_TAIL_MATHEMATICS.md` and
`STOCHASTIC_VOL_OF_VOL_RESEARCH.md`. New run `runs/m2_cap_tail_001`;
plots `visualizations/cap-tail-014`. Selected baseline remains
`m2_maturity_curve_fit_001/selected_model.json`; no new joint fit is adopted.
Current user priorities override archived queues: cap first, then continuous
stochastic vol-of-vol; exclude jumps and quintic proposals. Existing OU rates,
geometric grids, initial state, correlations and Q scenario remain fixed.
All seven SPX/four VIX dates, original 1511 OTM inputs/four forwards, cutoff 3,
common calls, market/model-forward parity and all 4096 samples remain active.

Nine fixed-curve VIX profiles (baseline, L +/-0.75, k=1/4, algebraic cap,
logarithmic continuations a=.25/.5 and separate channel caps) rebuild Q M.
Baseline prices replay within 1.37e-9. Raising L steepens the right wing but
worsens native VIX+forward loss 27.166865 -> 34.540828; log a=.25 gives
27.153015, algebraic 26.832740. None solves finite-strike shape. Saved checking
maxima exceed quoted strikes; likelihood-weighted probability/payoff tables
show why raw hits do not establish adequate tail prices. Upper cap bounds
about 937--1003 are not themselves attainability proofs. Payoff-relevant
conditional midpoint responses are compressed to about .57--.66.

Two explicitly VIX-focused curve/floor fits improve native VIX+forward training
to 20.044859 / 18.674677, reused separate checking to 20.653505 / 19.124019,
but actual seven-maturity training SPX losses rise 43.777476 -> 1656.376316 /
318.635148; common-call VIX criteria also worsen. NOT accepted joint fits.
Both hit implementation floor/curve bounds; VIX-only data leave off-window
variance pieces weakly identified, so SPX costs cannot be attributed to cap alone.
Next focused VIX diagnostics should freeze off-window pieces and floor rather
than repeat an unconstrained full-curve VIX-only search.

After cap profiles, two existing-cluster bounded spike-response mixtures worsen
VIX loss to 87.177224 / 117.866292. No SPX or joint success claimed for these.
Dedicated continuous latent response state is mathematically specified, with
parameters/randomness/cache/cost/integrability/martingality discussion and primary
literature; it is not implemented. Current nonlinear M2 already has state-dependent
Ito vol-of-vol; the proposal adds control of that response and persistence.

Selected uncapped first Gaussian multiplier moments actually remain finite
(stationary quadratic eigenvalue .897020488), while second moments diverge at
long dates (critical power 1.114801739). Do not claim uncapping necessarily makes
the selected M infinite, or that first-moment existence proves stock martingality.
GH32->48 and time8->16 representative extreme-state changes are at most
2.27e-6 / 3.40e-8 VIX points; no surface convergence certificate. Baseline VIX153
second-moment deficit 606.208 vs 610.933, SE1.635, persists; higher cap does not
clear it. Saved-draw audit passes 17 score records and 259 frozen-source checks.
Independent Python equation/quadrature/stock tests pass (direct runner; pytest
unavailable). Eleven PNGs plus SVGs cover all agreed maturities and diagnostics.
All numerical jobs are finished. Preserve baseline/historical sources, failures,
signed prices/unavailable IVs; no new Lean or P-compatibility certificate.

**COMPLETED — maturity-aligned curve calibration, 25 October:**
`runs/m2_maturity_curve_fit_001` completed after the user pause/resumption.
Report `MATURITY_CURVE_FIT_RESULTS.md`; plots `visualizations/maturity-curve-013`.
Sixteen curve pieces retain old 91.25/273.75-day boundaries and add all retained
SPX expiries and VIX expiry/window ends; previous curve exactly nested. 24 free
parameters across the exploration. Fixed OU dynamics/rates/geometric grids,
Q scenario/correlations/initial state; original OTM scope/cutoff3/seven SPX/four
VIX maturities, all4096 paths/common calls/GH32/GH48/time8/stock512 unchanged.
Hard xi_i>=floor, floor volatility 5--6%, beta_J>0; no xi monotonicity restriction.

Three full proposals, one accepted. Accepted step changes curve/floor only;
subsequent all24-free proposals improve VIX further but fail full joint safeguards.
Selected readout/cap stay unchanged despite being searched. Actual training:
canonical joint95.117702->89.320446; native joint84.982024->70.977709;
SPX49.518173->43.777476; native VIX+forwards35.430506->27.166865 (-23.32%).
Forward component alone worsens0.838006->0.890564; VIX smile still poor at high
strikes, with nonuniform near-ATM effects. Do not call this a satisfactory fit.
Separate reused VIX checking (selection frozen first): native VIX+forwards
35.862036->27.523027 (-23.25%), paired eight-net change-8.339009,
jackknife SE0.191220. Common-call change-0.601058 has SE0.781157, so its small
gain is not precise. Checking streams were used before: NOT pristine holdout or
market OOS. No SPX separate checking repricing in this run.

Selected floor0.00264 =>5.13809% minimum instantaneous volatility;
beta_J0.369729973>0. No negative native prices or unavailable common-call IVs
in training / VIX checking. SPX training still has1 monotonicity/66 convexity
flags; VIX has0/0. Raw signed values remain. Selected VIX153 second moment
606.208 vs target610.933 (SE1.635): rare-tail/numerical issues are not cleared.
Ten focused tests pass; audit replays70 local trials, all full acceptance decisions,
prices/SE/parity/loss/policy and74 sources/inputs. Eight PNG plots visually reviewed,
16 PNG/SVG hashes verified (`validation/plot_review.json`). Source snapshots and
original pause record preserved; continuation restores completed baseline and
reproduces saved proposal0 parameters, without new OU cache generation.
`selected_model.json` is the exploratory selected configuration. Not global
optimality, identified curve, price convergence or a new Lean/P certificate.
All processes are finished. Preserve prior runs and full agreed scope. Do not
restart archived searches automatically; remaining questions concern VIX smile
shape, forward tradeoffs, local-model safeguards and SPX estimator accuracy.

**25 October — explicit forward-variance hard constraint:**
Every new curve level must satisfy xi_i>=variance_floor, with the existing
variance_floor>=0.05^2 policy. `forward_variance_constraints.py` now supplies
this guard for any number of pieces; `calibration_policy.validate` calls it
for the current three pieces. Nonmonotonic levels and equality are allowed.
The fit launcher snapshots the new dependency. Fourteen focused tests pass,
including rejection of a below-floor value in each of14 positions and the
current calibration's relational check independently of its tighter box bounds.
Evidence: `runs/m2_forward_variance_constraint_001/validation.json` and tests.log.
No recalibration or maturity-grid expansion was launched. Historical results
and frozen sources are unchanged. Current policy/launcher source revisions are
intentional; replay historical source hashes against the preserved snapshots.

**Latest completed work — 5% volatility minimum and point1 acceleration:**
User requires minimum annualized volatility around5%, then acceleration point1.
For ALL NEW calibrations use `calibration_policy.py`: variance_floor in
[0.0025,0.0036] (5--6% floor volatility), beta_J>0; all eleven remain free.
The6% upper endpoint is the implementation choice keeping the floor free in
its stated reasonable range. Historical bounds/entrypoints remain reproducible.
Since xi>=floor and G/M>0, instantaneous sqrt(v)>=sqrt(floor)>=5%.
This is not an option-IV bound. Selected floor is0.00254, volatility5.039841%.

Completed pilot `runs/m2_joint_local_fit_001`, report
[JOINT_LOCAL_ACCELERATION_RESULTS.md](JOINT_LOCAL_ACCELERATION_RESULTS.md), plan
`POINT1_ACCELERATION_PLAN.md`, figures `visualizations/joint-local-012`.
`joint_price_local.py` approximates BOTH SPX and VIX prices/forwards at a true
anchor; native puts include MODEL-forward derivatives. `joint_local_search.py`
uses normalized SLSQP, a small radius-squared reserve and the unchanged strict
full-price guards. Accepted full checks refresh all eleven Jacobian columns.
Ten new/existing tests pass;119 local trial scores, saved full means/SE/parity/
losses/acceptance decisions,68 sources/inputs and selected-policy/config replay.
All six PNG figures visually reviewed; PNG/SVG output hashes verify.

Measured119 optimizer evaluations in0.3236s total (about2.54ms/evaluation), versus
historical33.35s nonlinear-VIX trial median. Approximate versus nonlinear work:
NOT an end-to-end speed factor or identical trial accuracy. Full baseline363s,
two full checked proposals377s/379s, whole pilot1143s. Three full pricing calls
still dominate. Both proposals accepted: common joint95.14831 ->95.11770;
native joint85.03863 ->84.98202; native VIX+forward35.51359 ->35.43051. Gains are
small. SPX and forward components individually worsen slightly; all four guarded
joint/VIX criteria improve. Selected beta_J0.36972997 remains positive.

Same1511 OTM inputs/four forwards/cutoff3/maturities, all4096 paths, fixed OU/Q/
initial state/geometric grids and GH32/GH48/time8/stock512. No coarsening or block
optimization. No new independent checking priced: overlays use actual TRAINING
prices. Selected training has0 negative native estimates and2 VIX55 missing IVs;
this does NOT clear the prior fresh-stream SPX negative-price failures or VIX
wing/forward misfit. No convergence/global optimality/P-compatibility claim.

Pilot driver `fit_joint_local.py` preserves its prior003 start for reproducibility.
Future continuation can restore THIS run's selected full-price/Jacobian checkpoint
with `parity_checkpoint.load_result` and build `JointPriceLocal` directly, avoiding
another identical baseline simulation. Preserve old runs and immutable sources.

**Latest completed parameter adjustment, 24 October — modest VIX gain:**
`runs/m2_parity_parameter_fit_003` is complete. Read
[PARITY_PARAMETER_FIT_RESULTS.md](PARITY_PARAMETER_FIT_RESULTS.md) and
[all-maturity plots](visualizations/parity-parameter-fit-011/README.md).
All eleven agreed parameters were free, beta_J>0; same OU/Q/initial state/grids,
1511 original OTM quotes/four forwards, cutoff3/retained maturities, GH32/GH48/
time8/stock512, common-call simulation and uniform4096 paths. One full true step
accepted; two subsequent20-evaluation local SLSQP cycles found no reserve-feasible
proposal. Budget termination is NOT convergence or evidence of model infeasibility.
The final candidate has beta_J0.34972997, cap_sharpness2 and floor0.0025 (the latter
two at the existing bounds). Starting point remains the previous positive-shift fit.

Training common-call joint98.1192 ->95.1483; native joint85.5281 ->85.0386.
Parameters frozen BEFORE checking. NEW independent4096 samples, same frozen
training-pilot control functions: common-call joint96.4690 ->92.2402; native joint
82.8286 ->81.2319. Original-OTM VIX-plus-forward36.1043 ->35.9254 (only0.50%);
paired eight-block loss-change SE0.0454 versus change-0.1789, descriptive only.
High-strike VIX IV errors remain large/slightly worse. Forward errors after:
[+0.3644,-0.5400,+0.6495,+0.0754], much larger than sampling SE at most dates.
This is a small numerical improvement, NOT an acceptable final smile fit.

Fresh checking exposes8 negative selected native SPX prices (57/302 days), hence
8 missing SPX IVs; two VIX55 IVs also unavailable due model/market-forward mismatch.
Do NOT carry the old ensemble's no-negative/fully-defined-SPX-IV result over to
new draws. No clipping/smoothing. Means/SE/parity/loss/acceptance and source/seed/
frozen-control replay passed; five focused fitting tests and production derivative
screens passed. Twelve representative selected-state GH32/GH48 VIX probes differ
by at most7.9e-7 points; NOT whole-surface/time/stock convergence. All3250 frozen input files verify byte-for-byte; all six PNG figures were visually
reviewed and PNG/SVG output hashes verify. Records are saved alongside the run.

Run001's first true canonical-improving proposal worsened native VIX and was
rejected. Run002 supplies the validated full baseline checkpoint and a local warm
point; both were explicitly stopped before checking, not adopted calibrations.
Run003 fully reprices/accepts that warm point; subsequent native-OTM-VIX-priority
local guidance guards BOTH joint/BOTH VIX criteria. Reserve constraints are extra
numerical margins, not a user requirement; failure to meet them does not exclude
smaller improving directions. Review solver scaling/feasible starts and SPX
estimator stability before treating stalled local progress as a family limitation.
Preserve all run/prototype/checking records and fixed-source manifests.

**Latest completed estimator experiment, 24 October — common calls, 4096 paths:**
See [PARITY_UNIFORM_RESULTS.md](PARITY_UNIFORM_RESULTS.md) and
[all-maturity plots](visualizations/parity-uniform-010/README.md).
`runs/m2_parity_uniform_001` implements one call payoff using only the SAME
1511 original OTM inputs/four forwards, and **ALL 4096 paths at EVERY strike**.
This is the user's final choice; adaptive-prefix allocation is superseded.
All seven SPX/four VIX maturities and cutoff 3 remain. VIX uses eight independent
512-point scrambles; SPX uses 4096 iid importance paths per expiry. Dynamics,
OU grids/Q prices/correlations/initial state, GH32/GH48/time8 and stock512 stay
fixed. No new optimization: old and previously selected parameters are repriced
on identical fresh samples. New independently piloted common-call hedge/Gaussian
controls preserve variance reduction. Earlier 512-path losses are not direct
comparators because estimator and samples changed.

The four former SPX57 missing IVs are all available; selected reconstructed OTM
prices have no negatives and no zero raw-payoff-hit contracts. Median SPX SEs
are about 35–38% of the 512-path prefixes on the same new draws. Two VIX55
common-call IVs (K=14.5/15) remain unavailable: positive model calls are below
market-forward intrinsic because E[VIX] misses the observed forward by -0.5167.
Actual forward errors [+0.3929,-0.5167,+0.6504,+0.1243] remain exposed. Converting
market puts uses market F; reconstructing native model puts uses MODEL E[VIX].
Common-call and original-OTM objectives differ for VIX forward mismatches; both
are saved. There is no cosmetic substitution of market F in model parity.

Same-estimator common-call loss old 1941.75 -> selected 98.12; original OTM
1724.67 -> 85.53. Do NOT interpret this as a precise 95% population improvement:
a rare old-reference SPX302 path has normalized stock 56.63 vs selected 5.58
on identical shocks and drives wing prices. Paired eight-block jackknife loss-
change SE is 1687.09 versus observed reduction 1843.63, despite all eight leave-
one-block comparisons favoring selection. Controls are frozen around selected
parameters, so baseline variance can be much larger. Only 100/1511 selected
native OTM prices are in bid/ask; 1 SPX order/72 convexity failures remain (VIX0).
Many SPX SEs remain several quote scales; price convergence is not established.

Fourteen focused tests, independent saved price/SE/parity/both-objective replay,
73 distinct source checks and 1674 immutable input-array hashes pass. Six PNG/SVG
figures cover all maturities, forwards, signed errors and prefix precision.
Reusable cache: `runs/m2_parity_adaptive_cache_001`; frozen independent pilots:
`runs/m2_parity_adaptive_controls_001`. Names retain the initial adaptive proposal;
production is uniformly 4096. Original cache preparation resumed only missing
chunks through a separately snapshotted wrapper after a tuple-vs-JSON-list resume
assertion; original design/source hashes were preserved. Queued runners 001/002
were interrupted before objective evaluation, not additional fits.

`reusable_ou_initial.py` implements exact deterministic starting-state translation
under fixed exogenous OU dynamics: delta X_t = exp(-Lambda t) delta X_0. Full and
projected paths, Gaussian covariances and unchanged equity increments/importance
scores/components/weights pass identical-seed tests. Recompute G, M_Q, stock and
conditional VIX after an initial-state change. Full OU states are required; the
existing minimum-Mahalanobis projection-to-state map is available. Initial state
remains fixed here; it is not silently added to the agreed eleven parameters.
All jobs are complete. Preserve the prior selected positive-beta_J model as the
working exploratory candidate; do not treat this repricing as a new calibration.
Future fitting must use this explicitly selected common-call convention and
all 4096 samples, report native OTM/forward errors too, and keep the agreed eleven
parameters/beta_J>0 unless the user expands scope. Raw tail sensitivity, signed
shape failures and residual VIX forward/wing misfit remain open numerical issues.

**Latest completed exploration, 24 October — restored-precision full readout:**
The user's request to explore beyond the old fit is implemented in
`runs/m2_controlled_readout_fit_002`; see
[CONTROLLED_READOUT_EXPLORATION_RESULTS.md](CONTROLLED_READOUT_EXPLORATION_RESULTS.md)
and [all-maturity plots](visualizations/controlled-readout-fit-009/README.md).
All eleven agreed parameters were free across three prespecified starts, with
49 total curve-prefit/all-parameter evaluations. OU grids/correlations/projection
powers/Q prices/initial state and the full OTM cutoff-3 market scope stay fixed.
Both SPX ensembles are restored to 512 paths and VIX to 8x64 points, conditional
GH32/normalizer GH48/time order 8 and 512 stock steps per trading day. Continuous
zero-mean control draws remove nearest-donor price jumps; signed/shape failures
remain visible. The same GH32 rule is accelerated by a checked C-convolution
interpolant, with final selected checking prices replayed by direct tensor GH32.
The direct preliminary run 001 was interrupted and retained; it is not another
completed calibration. The immutable cache is `runs/m2_controlled_readout_cache_001`.

Searching loss **256.33→83.53**; separate checking **305.21→103.34 (66.14%)**,
with aggregate SPX/VIX-option/forward components all improving. All eight
leave-one-block changes favor the selection; descriptive paired change SE is
51.53 versus reduction 201.87. Checking streams were used in earlier research;
this is not new independent validation of the full research history or market
out-of-sample performance. Selected **beta_J=0.29972997**, three variance levels
**[0.04930160,0.05859463,0.05015265]**, floor **0.00248112**, sharpness **2.02455**.
The floor/cap-sharpness approach upper/lower bounds, and all three runs stop at
their evaluation limits. Only **114/1511** prices are inside bid/ask; four
checking SPX prices are negative, two raw-payoff SPX contracts have no hits,
and finite-sample shape failures persist. SPX 666/1030 and VIX 153 days worsen
(VIX 118 slightly worsens). This is meaningful provisional numerical progress,
not a satisfactory fit, a global optimum or full M2/P compatibility assessment.

Fourteen focused tests pass; all seven score records/SEs replay, 55 run and
86 cache source/input hashes and 260 immutable arrays verify. Direct GH32
price differences are at most 1.4e-9 quote scales; selected-state GH32→48 probes
change raw VIX at most 6.3e-6 points. Same-GH48 normalizer probes pass 2e-7.
VIX second-moment and raw stock/mass diagnostics are retained without enforcing
identities or rescaling prices. No new formal certificate is asserted.
All jobs are complete. `selected_model.json` is a self-contained exploratory
positive-shift working candidate; the Stage-3 loader remains a HISTORICAL
zero-shift diagnostic anchor, not a loader for this new fit. Preserve all old
runs and the full 1511-option/four-forward scope. Next bounded questions are
profiles beyond the newly hit floor/cap-sharpness bounds and candidate-appropriate
continuous controls for residual signed/shape failures. Neither is launched;
do not automatically restart old raw-estimator or Q-scenario searches.

**24 October follow-up — warm start and strike-continuous controls:** the user
requests starting from the preceding fit and eliminating strike jitter.
The optimizer is bounded least squares (TRF), with analytic derivatives.
The proposed first start now preserves all ten other preceding parameters,
changing only beta_J from 0 to 0.01; the tested helper and exact vector are in
`continuous_strike_controls.py` and the diagnostic's `warm_start.json`.
[STRIKE_CONTINUITY_RESULTS.md](STRIKE_CONTINUITY_RESULTS.md) records a saved-draw
experiment on both historical 512-path ensembles, all seven SPX maturities.
Continuous interpolation of frozen zero-mean control DRAWS removes nearest-
donor jumps: checking strike-order failures fall 11→0 and searching 10→0,
without discarding the old controls' precision. It leaves 15/13 convexity
failures, respectively. A global affine comparison preserves convexity but
creates negative prices/order failures and is not adopted. Final retained
diagnostic: `runs/m2_continuous_strike_controls_003` (001/002 retain development
comparisons; these are the same samples, not independent experiments).
Six focused tests pass. No new simulation, calibration or historical price
replacement occurred. A complete pricing repair is not yet established;
do not resume the degraded raw-estimator calibration or claim that continuity
alone establishes positivity, convexity or numerical accuracy.

**Latest completed attempt, 24 October — fixed-reservoir eleven-parameter fit:**
The user confirmed freeing all five score coefficients, two cap parameters,
three variance levels and the variance floor, while keeping the OU grids,
correlations, projection powers, Q prices and C0=1/Y0=F0=0 fixed. See
[FULL_READOUT_FIT_RESULTS.md](FULL_READOUT_FIT_RESULTS.md) and
[all-maturity plots](visualizations/full-readout-fit-008/README.md).
`runs/m2_full_readout_fit_001` completed two starts/63 evaluations, selecting
beta_J=0.01 at its positive lower bound; cap sharpness and variance floor also
hit bounds. The full 1294-SPX/217-VIX/four-forward OTM in-sample scope and
SPX cutoff 3 are unchanged. This run uses positive raw importance payoffs,
without strike-dependent controls; old parameters are repriced on the same
new ensemble, so historical controlled-estimator losses are not comparators.
Search loss 2168.61→586.91; separate checking loss 1044.10→906.61, but VIX
options/forwards worsen and only 61/1511 prices are inside bid/ask. Descriptive
paired loss-change SE 876.70 exceeds the 137.50 reduction, with a leave-one-block
reversal. No reliable improvement or satisfactory joint fit is established.
Some tails have zero sampled payoffs; prices remain scored and their IVs are
marked unavailable. Six focused tests, six score replays, 49 source/input
hashes and 67 distinct seed records pass; within-side strike monotonicity has
zero violations. GH12→24 changes reach 1.92 quote scales; the selected
normalizer interpolation probe error 1.28e-6 exceeds the former 2e-7 screen.
P compatibility and price convergence are unresolved. Keep this positive-shift
candidate as exploratory evidence, not a qualified replacement for prior fits.
All jobs are complete. This broader fit supersedes the older curve-only/Q-scenario
queue as the latest user-requested attempt; do not automatically rerun or treat
the next Q scenario as authorized by old continuation text alone.

**24 October 2026 user decision — strictly positive spike asymmetry:**
Future M2 fitting must enforce `spike_shift = beta_J > 0` in
`eta = c_C*C + q_Z*(Y-beta)^2 + q_J*(F-beta_J)^2`.
With `q_J > 0`, positive beta_J makes negative spike-factor states contribute
more to volatility than equally sized positive states. This is factor-sign
asymmetry, not a separate volatility jump process. The existing working
candidate and completed all-maturity run have beta_J=0: retain them as
historical baselines, not candidates satisfying the new requirement.
The later eleven-parameter attempt selected beta_J=0.01, without establishing
a reliable joint fit (see above). Before any subsequent fit,
explicitly set and record a positive value or a search domain with a strictly
positive lower bound; reject zero/negative candidate values rather than
silently inheriting the loader default. The earlier illustrative [0,0.5]
domain is superseded at its zero endpoint. Recompute the readout-dependent
normalizer, stock paths and VIX exposures for the changed shift. Earlier
next-step instructions to keep the zero-shift readout fixed are conditional
on first establishing a compliant positive-shift baseline. This decision
records a constraint; it does not launch a fit or alter archived results.

24 October explanatory follow-up: the latest all-maturity optimizer fitted
only three variance-curve levels, with structural shape/state/Q parameters
fixed. Saved-price/control decomposition at SPX 119/211/666 days identifies
3/3/2 strike-monotonicity violations in corrected estimates and none in their
raw weighted-payoff counterparts. The frozen strike-dependent correction
explains concrete jagged features shared by the orange and blue curves; this
is not evidence of discontinuous true model prices. See the dated section of
[ALL_MATURITY_FIT_RESULTS.md](ALL_MATURITY_FIT_RESULTS.md) and the retained
`validation/strike-jitter-decomposition.{csv,json}` in the all-maturity run.
No recalibration, price replacement, smoothing or new formal certificate was
performed in this explanatory follow-up.

Checkpoint: 23 October 2026, branch `codex/m2-option-pricing-benchmark`.
Committed handoff base: `df27b4a` (preceding numerical implementation `afbf628`).
The current all-maturity implementation is identified by source snapshots/hashes
in `runs/m2_all_maturity_fit_001/run.json` and its validation manifest. Earlier
curve/state/refinement runs retain their own versions. Scoped Lean checks
are in `formalization/CalibrationChecks.lean` with pinned compiler/axiom audit 006.

This is the entry point for the accumulated attempts and their implications for
fitting. Detailed reports and their immutable run directories remain the evidence;
this index does not replace them. Older reports' “next” sections and test counts
describe their historical checkpoints. Use the current queue below when resuming.
The copy-and-paste continuation prompt is [CONTINUE_FITTING_PROMPT.md](CONTINUE_FITTING_PROMPT.md).

## Goal and present conclusion

**LATEST USER SCOPE CHANGE, 23 October:** fit **all seven SPX maturities
29/57/119/211/302/666/1030 days**, exclude only 1/1401/1765/2129, cutoff 3.
The user explicitly requests **in-sample feasibility, no market out-of-sample
results**, and confirms recalibration is needed. This supersedes older reserved
maturity/strike roles for the current exercise. All eligible original quote
splits now enter fitting: **1294 SPX options and 217 VIX options**, with all
four VIX forwards. Original source labels and historical runs are preserved.
See [ALL_MATURITY_FIT_RESULTS.md](ALL_MATURITY_FIT_RESULTS.md) and the captured
[plan](ALL_MATURITY_FIT_PLAN.md). `runs/m2_all_maturity_fit_001` is **COMPLETE**:
all three curve levels were refitted jointly at fixed state/readout/Q. Selected
levels are **[0.05762055,0.07134085,0.06321057]**. The last square-root variance
scale rises **18.89%→25.14%**, materially improving the 666/1030-day SPX fit.
Searching loss falls **300.94→253.44**; a separate numerical ensemble on the
SAME market quotes gives **351.36→308.64 (12.16%)**, improving both option
families and forwards. Only **103/1511** prices are inside bid/ask, all SPX
nominal precision screens fail, and substantial short-SPX/VIX shape errors
remain. The [final pictures and guide](visualizations/all-maturity-fit-007/README.md)
show model prices on every panel, all marked in sample. All four score records,
membership, source hashes and sample means/SEs replay; 882 prior contracts match
to 2.8422e-13. An audit-reader shape error was corrected without repricing or
reoptimization, retaining its failure log. This expanded curve is now the
working reference; future fitting must keep the full 1511-option in-sample
scope. Next bounded hypothesis: declared Q-risk-price scenarios one at a time,
starting from this curve, not the old four-maturity fit. No numerical jobs remain
running. Numerical checking is not market out-of-sample evaluation. Active scope intentionally changed; older
run source snapshots retain the old scope and old hashes must not be refreshed.

**Preceding completed checkpoint, 23 October:** the user authorized better joint fitting,
same-expiry model-forward IV diagnostics, then asked for more curve pieces.
See [the completed state/refinement report](JOINT_STATE_REFINEMENT_RESULTS.md)
and [the prespecified five-piece design](CURVE_REFINEMENT_PLAN.md).
`runs/m2_joint_state_fit_002` completed the bounded seven-state/three-level search
(333.6764 at zero to 288.8213 at C0=1), then failed before a fresh score on a final
stock-chunk shape mismatch. `runs/m2_curve_refinement_002` completed the nested
five-piece search (288.8213 to 284.8048), then stopped waiting for that parent.
Both failures and selections are retained. Corrected runs
`m2_joint_state_fit_003` and `m2_curve_refinement_003` completed the frozen
selections' new checking streams without repeating search. State plus three-piece
refitting improves all three fresh trials and pooled loss **355.46→316.93**;
SPX/VIX-option/forward components all improve. Only **41/882** prices are inside
bid/ask and all SPX precision screens still fail. Five pieces add only a
**0.64%** pooled gain, **316.93→314.90**, with one trial reversal, uncertainty
larger than the gain, worse forwards and only **38/882** prices inside bid/ask.
Keep the selected C0=1 state and three-piece curve as working reference; do not
add more pieces now. The paired comparisons share checking ensembles and do not
establish numerical qualification. All numerical jobs/audits are complete.
Next fit the already declared Q-risk-price scenarios one at a time, refitting
three curve levels while keeping state/readout/architecture fixed; see the
bounded next hypothesis in the report. There is no new qualified fit. The original
VIX-only anchor and all failed points remain. Thirteen current Lean identities
pass audit 005 and 179 Python tests pass; these do not certify numerical prices.
The same-expiry model-forward IV conversion substantially reduces the display
gap but changes no prices or forward errors; substantial smile-shape errors
remain. Historical raw-parity diagnostic 001 exposed 18 validation and 17
final-test VIX strike keys to a market consistency calculation, not model
performance scoring or selection. The exposure list is in diagnostic 002;
future evaluation on these keys is quote-exposed rather than wholly blind.
New parity diagnostics use exact training option keys, with no matched pairs.

Fit the 19 February 2026 SPX/VIX snapshot with M2 under the audited common price
objective, eventually comparing fairly with Othmane's preserved ESN engine and
checking physical compatibility. **No reliable joint market fit, final model
ranking or continuous-model pricing convergence has been established.**

The initial exploratory VIX curve fit improved its fixed simulation sample but
reversed on fresh simulations. Subsequent work identifies rare outer/volatility
paths as a major numerical obstacle. Several variance-reduction ideas help
specific diagnostics, but none clears the complete pricing surface. Avoid
mistaking this numerical bottleneck for demonstrated model infeasibility.

**Actual joint fitting has now been attempted:** three nuisance-curve levels,
all 707 cutoff-3 SPX priority options, all 175 VIX options and four forwards.
The selected curve improves total loss on all three fresh numerical checks;
pooled loss falls **597.90→347.10**, but SPX improvement sacrifices aggregate
VIX fit, only **50/882** objective options are inside bid/ask, and accuracy
gates remain unresolved. This is a provisional sampled improvement, not a
qualified joint market fit. Ten scoped Lean identities and 170 Python tests
pass; neither certifies numerical tail/bias accuracy.

Visual inspection is available in the [joint-fit plots and reading guide](visualizations/joint-curve-fit-002/README.md):
all eight expiry smiles, scaled price errors, four VIX forwards and the variance
input, using saved pooled fresh checks. Calls/puts and missing IV-band endpoints
remain distinct; bands contain sampling error only. No additional fit was run.
The [expanded maturity/variance context](visualizations/joint-curve-context-003/README.md)
adds all 11 SPX expiries, with the seven additional slices explicitly market-only.
The reserved 119-day market curve has now been visually exposed at user request;
it remains excluded from fitting and has no model score. No new prices or fit
were computed. The annotated variance levels also document weak identification
of the long final segment and the anchor's fixed third coefficient.
The subsequent [focused view](visualizations/joint-curve-focus-004/README.md)
follows the user's reaffirmed cutoff 3 and omits SPX 1/1401/1765/2129 days from
current displays/active analysis. Retain SPX 29/57/119/211/302/666/1030 days,
keeping their existing fitting/reserved/diagnostic roles. Additional market-only
slices use an explicitly recorded display-only cutoff-3 proxy; the 707 fitting
rows and all four VIX maturities remain unchanged. Historical data and plots
are retained. See `active_scope.json` for the current scope.

## Decisions to preserve

- [active_scope.json](active_scope.json) is authoritative: VIX 27/55/118/153 days,
  now all **217 in-sample options** and four forwards. Historical training-only
  runs used 175. An earlier two-maturity experiment
  was explicitly superseded. **Keep all four maturities.**
- SPX priority is provisional maturity-scaled cutoff 3:
  `abs(log(K/F(T))) <= 3*sigma_ATM_proxy(T)*sqrt(T)`, now using all eligible
  nearest-ATM market midquotes at each of the seven retained dates. This covers
  **1294/1751 quotes**. Current membership is in
  `runs/m2_all_maturity_fit_001/spx_membership.csv`. The historical 707/980
  training-only membership remains in `runs/spx_scaled_moneyness_003/`.
  Four old near-ATM contracts and the newer 24-contract panel are representatives,
  not the whole priority population. Wings remain secondary diagnostics.
- Diagnostic baseline: frozen Stage-3 candidate 18/scenario 0, zero projected
  initial state, historical VIX-only selected nuisance curve, and original
  vol7/ridge5 SPX estimator. The new joint curve is a separate provisional
  candidate, not a silently replaced diagnostic baseline.
  Read actual run configurations before reusing a curve or cache; old initial
  stripped-curve and later selected-curve results are different experiments.
- Recent SPX experiments use 512 steps per trading day. Exact OU transitions do
  not make the nonlinear stock update exact; **512 is not convergence clearance**.
  The old 2/8-step runs are coarse historical diagnostics.
- Preserve signed negative controlled estimates as failures; never clip them into
  acceptable quotes. Positive raw/conditional prices can still be badly sampled.
- Current all-maturity feasibility work has no market holdout; preserve old
  split labels only as metadata and label all current scores in sample. A fresh simulation is a
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
| [Actual joint curve optimization + Lean checks](JOINT_CURVE_FIT_RESULTS.md) | Eight starts converge near `[0.0600942,0.0748506,0.0383635]`. All three fresh total-loss comparisons improve; pooled SPX loss 529.20→249.04, VIX option loss 67.88→84.91 and forward loss 0.82→13.14. Only 50/882 objective options fit inside intervals. | Provisional candidate and restricted-family trade-off; no qualified market/P fit or ranking. Next fit restricted state with appropriate normalization updates; do not reopen standalone sampling-feature searches. |

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
| [Conditional existing controls / retained proposal](SPX_CONDITIONAL_CONTROL_RESULTS.md) | Conditional stock, reference payoff and two adapted hedge gains derived and implemented; new retained-shock likelihood validated. Against original vol7/ridge5, 34/48 variances improve but whole-panel variance-cost gains are only 0.7483/0.8541; just 1/8 expiry/trial panels wins. Two nominal small-SE tail cases have effectively empty pilots and >99% single-path concentration. One primary and one half-pilot reference estimate are negative and retained. Do not adopt; bounded extension completed. |

Latest conditional-control run: `runs/m2_spx_conditional_controls_001/`;
288 full/half-pilot prices, 52 source/input hashes, 48 batch seeds, likelihoods,
coefficients and diagnostics replay. **124 package + 42 study-driver tests pass**.
The prior raw conditional run remains at `runs/m2_spx_conditional_001/`.
These are implementation checks, not numerical convergence certificates.
Small independent equity variance is not an upper bound on option/control
variance gains from conditioning. Conditional positive values are not raw
exercise hits; the preliminary audit's misleading screen label and original
outputs are retained and corrected in [AUDIT.md](AUDIT.md).

## Return to fitting: ordered queue

0. **LATEST user priority: all retained maturities, entirely in sample.**
   `m2_all_maturity_fit_001` is complete; do not repeat it.
   Include 1294 cutoff-3 SPX options at all seven dates, 217 VIX options and
   four forwards, with actual curve recalibration and model overlays everywhere.
   No market holdout or out-of-sample interpretation. The recalibrated reference
   is [0.05762055,0.07134085,0.06321057] at C0=1/Y0=F0=0. Separate numerical
   loss 351.36→308.64, 103/1511 inside intervals; long SPX improvement remains a
   poor joint fit. The next Q-scenario comparisons below must use this expanded
   scope/reference, not revert to training-only or four maturity inputs. All 182
   Python tests and thirteen Lean identities (audit 006) pass within stated scope.

1. **Conditional-control extension completed; stop this feature search.** The
   tested implementation loses whole-panel variance-cost in both fresh trials
   and retains serious tail/pilot uncertainty. Retain original vol7/ridge5;
   do not repeat this experiment or substitute another dictionary by default.
2. **Readiness decision completed.** [FITTING_REENTRY_DECISION.md](FITTING_REENTRY_DECISION.md)
   separates sampling, coupled stock-grid, normalizer, VIX conditional and time
   quadrature evidence. No qualified market fit or ranking is ready. This is
   unresolved numerics, not demonstrated M2 infeasibility.
3. **Actual joint optimization and three fresh checks completed.** The user
   superseded the unlaunched nine-point/24-SPX proposal with actual fitting and
   Lean checks. `m2_joint_curve_fit_002` frees all three curve levels and scores
   the complete 707-SPX priority population plus all VIX targets. Other 273 SPX
   training quotes are secondary diagnostics. All three fresh comparisons favor
   the selected curve, but its VIX forward misses and low interval coverage
   expose a poor compromise. Retain the candidate without calling it a validated
   fit. Its next restricted-state stage has now completed, as follows.
4. **State/curve and nested five-piece fitting completed.** Seven initial-state
   profiles select C0=1/Y0=F0=0 under Q scenario 0. All three fresh comparisons
   favor the state/three-piece refit; pooled loss 355.46→316.93. Five pieces
   add only 0.64%, with one reversal and uncertainty larger than the gain.
   Keep the three-piece state-adjusted working reference and all five-piece
   evidence. Do not repeat either search or add more time knots by default.
5. **Next bounded fitting hypothesis: declared Q-risk-price scenarios.** Start
   with the feedback pair (+0.25,0)/(-0.25,0) against scenario 0 at selected state,
   refitting the three curve levels; then consider the declared spike pair.
   Declare budget before launch and freeze before NEW checking draws. Rebuild
   candidate normalization, state mapping, stock and VIX exposures; revisit
   primitive-law/control/likelihood correspondence and P diagnostics. Keep
   state/readout/architecture fixed. This is not launched here. Do not normalize
   market errors by model SEs or conflate sampling and deterministic bias.
6. **Qualify and compare.** Only after sufficient numerical accuracy, assess
   all scored quotes, physical compatibility and reserved market evaluation.
   Complete the common-protocol ESN comparison before ranking models. This
   phase remains future work, not an achieved deliverable.
   Preserve the documented VIX quote-exposure caveat for the 35 raw-parity keys;
   no held-out model performance has been scored, but those quotes are not blind.

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
