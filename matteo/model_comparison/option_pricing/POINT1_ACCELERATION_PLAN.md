# Joint local-price acceleration and 5% volatility floor

User instruction, 24 October 2026: enforce an approximately 5% minimum annualized
volatility, then implement acceleration point1. Calibration policy is now
`calibration_policy.py`: floor variance in [0.05^2,0.06^2], i.e. [0.0025,0.0036].
The 6% endpoint keeps the floor free and reflects the user's stated 5--6% range.
All eleven parameters remain free; xi bounds exceed the largest allowed floor,
so v=floor+(xi-floor)G/M >= floor, since G/M>0. This is a pointwise bound for
model instantaneous volatility; it does not directly bound option implied IV.
Historical model configurations/bounds and immutable source snapshots stay intact.

## Controlled implementation experiment

Run `m2_joint_local_fit_001`. Start from parameter fit003's selected vector,
rounding its floor up by about9e-15 to the new exact lower bound. Fully recompute
baseline prices AND all eleven analytic derivative columns once. Subsequent
inner trials use an affine approximation of ALL SPX/VIX call prices and VIX
forwards. Apply the original nonlinear interval loss to these approximate prices.
Reconstruct native model puts using the approximate MODEL forward, including
its derivative. No pricing engine is available to the local-model object.

Use normalized original-OTM VIX-plus-forward merit and four normalized constraints:
both joint and both VIX-plus-forward losses. Small extra reserves shrink as radius
squared; all true acceptance conditions remain strict. Initial normalized box
radius0.02. At most two full proposals. A rejected proposal shrinks the box;
an accepted proposal refreshes the full price/forward Jacobian. Budget20 minutes
after the baseline. No claim of convergence/global optimality/parameter identification.

Scope stays identical: all1511 original OTM quotes/four forwards, SPX cutoff3 and
seven retained maturities, four VIX maturities, fixed OU rates/grids/Q prices/
correlations/initial state; common calls; all4096 paths/strike; GH32/GH48/time8/
stock512 unchanged. No coarser time mesh, parameter-block optimization or new paths.

Record true/approximate scores separately, full call counts, timing and prediction
errors at each checked proposal. Full4096 nonlinear pricing must improve BOTH
joint and BOTH VIX-plus-forward criteria before accepting. No checking data enter
selection. This is a training-ensemble implementation/acceleration pilot; no new
independent checking ensemble is priced. Generate all-maturity overlays from true
training prices and preserve signed estimates, invalid IVs and SE diagnostics.

## Verification and timing interpretation

Tests cover the5--6% floor policy, all eleven free bounds, positivity bound, exact
saved-anchor objective/native-price replay, all44 criterion derivatives, affine
predictions, model-forward parity and no simulator calls inside local optimization.
Replay saved Monte Carlo means/SE/parity/true losses and every acceptance decision.
Freeze executable sources before running; original input caches remain read-only.
Compare approximate-trial latency with recorded old nonlinear-VIX trial latency
only as different computational work. Report full-price work separately; do not
translate a local speed ratio into an end-to-end speed guarantee or pricing accuracy.
