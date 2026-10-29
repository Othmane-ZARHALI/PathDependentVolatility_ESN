# Maturity-aligned forward variance calibration

25 October2026. User requests parameter adjustment toward better VIX fits after
discussing maturity-aligned curves and explicitly approving xi_i>=floor, without
monotonicity. Extend the deterministic curve, keeping OU dynamics/rates/geometric
grids/Q premiums/correlations/initial state fixed. Keep1511 original OTM quotes,
four forwards, SPX cutoff3/all seven retained maturities/four VIX maturities,
common calls/all4096 samples, GH32/GH48/time8/stock512.

## Fixed design before pricing

Curve ends in days:27,29,55,57,85,91.25,118,119,148,153,183,211,273.75,302,666,1030.
These are SPX expiries, VIX expiries and their30-day window ends, plus the two
old91.25/273.75 boundaries. Retaining those makes the previous three-piece curve
an exactly nested model.16 curve levels +7 readout/cap parameters +floor =24.
Start from `m2_joint_local_fit_001` selected model. Each new level is initialized
to its original interval value. Existing independent readout bounds retained;
curve levels in[0.0025,0.2], floor in[0.0025,0.0036], additionally every xi>=floor.
beta_J>0. No nondecreasing restriction. Retain original0.1*sum(diff(xi)^2)
roughness penalty (baseline unchanged because repeated levels have zero jumps).
This is weak regularization; extra pieces are not asserted to be identified.

## Numerical implementation

Generic24-column stock sensitivities, checked against the old three-piece kernel
and finite differences. VIX geometric-lag integration panels are split at EVERY
new curve boundary. Outer states/shocks/importance weights and fixed pilot hedge/
Gaussian controls are reused. Equity uses its unchanged cached left-point512 mesh:
new curve boundaries can lie within a stock step. Record deterministic integrated-
curve left-grid versus exact overlap errors; this is not a stock convergence pass.

Freeze the analytic lognormal SPX reference-control curve at the prior selected
three levels. Its known Q mean remains valid for any new curve, and its baseline
draws/prices are identical. It is no longer changed with candidate xi, so its
parameter derivatives are zero. The fixed hedge controls also retain zero means.
This avoids regenerating112GB-scale primitives or accumulating new curve Brownian
integrals. Actual price derivatives account for this estimator choice.

## Search and checks

Fully reprice the nested baseline and all24 sensitivities. Verify baseline SPX
price correspondence and report changes from VIX panel refinement. Production
finite differences check all16 curve and floor VIX price/forward columns, using
cached readout exposures (no extra expensive Gaussian expectations).

Use the accelerated affine SPX/VIX local price model. First propose a curve/floor
step with readout fixed (radius0.1); subsequent proposals free all24 parameters
(initial radius0.05). Native-OTM VIX+forward merit; guard both joint/both VIX+forward
criteria. Hard xi>=floor constraints. At most3 proposed candidates /30 search
minutes after baseline. Before full stock pricing, check ACTUAL nonlinear VIX
prices/forwards on all4096 samples and reject if either VIX criterion worsens.
Every accepted step must pass FULL SPX/VIX nonlinear pricing with all strict
guards. Rejected/unsuccessful steps remain recorded; no convergence claim.

Freeze selection before numerical checking. Reuse the EXISTING separate seed
ensemble `m2_parity_check_cache_001` to check VIX only (no new simulation), same
fixed pilot coefficients. These streams have been used in earlier research;
this is not a pristine holdout or market out-of-sample test. Report paired
eight-net sensitivity and per-maturity/strike-region diagnostics. SPX overlays
use full TRAINING prices; VIX has training AND separate-checking overlays.
Preserve signed prices, invalid IVs, model/market forwards, shape/SE/tail limits.
