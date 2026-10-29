# Controlled cap comparison: fixed floor and four broad window adjustments

Authorized follow-up, 25 October 2026. Explain results at undergraduate level.
Question: can a gentler volatility brake improve high-strike VIX prices without
unrelated curve changes or an unacceptable actual SPX tradeoff?

Baseline remains m2_maturity_curve_fit_001/selected_model.json. Preserve all old
runs and production code. Fixed OU reservoir/rates/geometric grids, initial
state, correlations, Q scenario, positive beta_J, floor .00264. All 1511 original
OTM options, seven SPX/four VIX dates/four forwards, cutoff3, common-call parity
and all4096 samples per strike. No jumps or quintic approach/new latent layer.

Four window indicators cover [27,57], [55,85], [118,148], [153,183] days. On a
curve piece lying in two windows, average their adjustments rather than add.
Let B be these row-normalized indicators on the existing16 curve pieces:
xi_new_i=xi_baseline_i+(xi_baseline_i-f)*(B*s)_i, |s_j|<=.15.
Thus excess variance changes by at most15% per piece; all off-window pieces and
floor are EXACTLY fixed, xi>=f, and no monotonicity is imposed. This15% range
is an experimental choice, not a user requirement. Every window edge is already
a curve boundary. No hidden time-grid change or additional curve pieces.

Nine finite cap settings: L offsets0/.375/.75 crossed with logarithmic
continuation a=0/.25/.5. k=2 fixed. a=0 is the current logistic-cap family;
a>0 adds a*log(1+softplus(k*(eta-L))/k) to h. Rebuild candidate-specific Q M and
conditional exposures. All original score coefficients/floor fixed.

For each shape, optimize only four broad adjustments, starting at zero, at most
60 SLSQP iterations. Objective is existing native-OTM VIX option + forward loss.
Diagnostic safeguards: each VIX forward within .5 point of market, and common-call
VIX+forward loss no greater than the unchanged baseline. These are experimental
screens, not new user requirements or accepted joint calibration rules. Keep
failed/infeasible outcomes. Report all four smile levels/slopes/curvatures,
near-ATM errors, signed prices/unavailable IVs, forwards, tail contributions.
The grid and local fitting are bounded, not a global convergence test.

Freeze two TRAINING representatives before fresh checks: best capped (a=0) and
best continued (a>0), ordered by equal-maturity right-wing IV MAE (K/F>2), among
feasible candidates. Prefer candidates improving native VIX+forward, with
near-ATM IV MAE no worse than baseline+5 percentage points and all right-wing
IVs available. If this preference set is empty in a family, use its feasible
candidate with smallest wing MAE and explicitly mark failed preference tests.
If the entire family is infeasible, preserve the failure and do not force a
full SPX run for it. No baseline adoption is authorized by these screens alone.

At most TWO full all-seven-maturity SPX evaluations, recomputing volatility and
stock on all original4096 paths/expiry and stock512 mesh; same frozen controls.
Restore baseline prices from their checked saved arrays, not an identical full
stock rerun. Record paired eight-block uncertainty and signed shape failures.

After selection freezes, generate genuinely NEW VIX-only outer streams:
base seed2026102511, SeedSequence([base,2,days,net]), eight512-point independent
scrambled nets per date, same existing defensive-mixture proposal and exact p/q
likelihoods. Verify seed disjointness against recorded historical VIX seeds and
prior marginal diagnostic seeds. Reprice baseline plus both frozen candidates on
these NEW states with training-pilot controls unchanged. This checks numerical
sampling on the SAME market observations, not market out-of-sample performance.
No new stock cache, no pilot refitting, no fresh-check-driven reselection.

Separate immutable run m2_cap_controlled_001, source/config hashes, snapshots,
raw arrays and failures. Plots cap-controlled015 cover every agreed maturity,
price errors, forwards, bounded curve changes, weighted tails and fresh checks.
Audit saved prices/SE/parity/loss, constraints, exact frozen coordinates, basis
sensitivities and fresh weights/state mapping. Preserve the unresolved VIX153
second-moment discrepancy; no new Lean/P/continuous-price certificate.
