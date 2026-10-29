# Joint SPX/VIX refit of the selected two-speed stress model

28 October 2026. User authorization: "stay for the moment with the 2-speed
stress model and try to refit the SPX options as well". The VIX-only pause is
lifted for this experiment. The four-speed bank remains separate.

Use the selected model in runs/m2_two_speed_001/experimental_model_two_speed.json.
The normalized variance mixture and two independent OU stress factors are
unchanged. Original M2 architecture, initial state, risk premia, capped readout
and positive spike shift remain frozen for this first joint refit.

Restore the seven SPX maturities29,57,119,211,302,666,1030days. Preserve the
previous market-only cutoff3 membership:1294primary quotes and1751eligible
quotes including diagnostic wings. Preserve all217native VIX quotes and four
forwards. No clipping, smoothing, removing signed estimates or expiry-specific
model parameters. Numerical checking uses the same market snapshot, not market
out-of-sample observations.

Free all16 maturity-aligned variance levels (previously only seven VIX-sensitive
levels varied), floor volatility5–12%, and the six existing stress parameters:
two sensitivities0–5, two ordered half-lives1–365days, total stress weight and
fast fraction0–1. There are23 fitted parameters. Active curve levels satisfy
floor<=xi<=0.2. The original M2 readout is not refitted in this first experiment.

Objective: equally weighted SPX and VIX families, with equal maturity weights
and native bid–ask-scaled squared interval distances within each family.
Four VIX-forward equalities are imposed on the calibration estimator. Report
checking gaps and block SEs, not exact population-forward claims. First perform
a curve/floor-only joint fit; then a bounded local joint fit including existing
stress parameters. Select using training loss and forward feasibility. Retain
VIX-only and intermediate curve-only candidates to show the SPX/VIX trade-off.

SPX requires full adapted paths. Reuse the old exogenous M2 projections and
matched equity increments from the existing512-path per-expiry search/check
caches. The corresponding old path importance weights and known-mean controls
remain valid, since the state/equity law is unchanged. New independent stress
innovations are cached with separate seeds. Rebuild stress paths, v, stock,
payoffs and derivatives for every candidate; never reuse a completed stock path.
Exact Gaussian OU transitions are used on each integration grid. Different
candidate rates can share innovation normals as a common-random-number coupling;
state values must change with the rate. Stress factors remain independent of
M2 and equity, as in the selected model.

Use a64-step/trading-day grid for the initial bounded fit, obtained by coarsening
the existing512-step grid while retaining expiry events and both fine-grid points bracketing every
curve knot. Any remaining curve-knot delay is at most one original512-step
interval, and belongs to the stock timestep approximation. Reuse exact fine-state
values at coarse left endpoints and sum equity increments. Check selected
models on the separate checking state/equity stream and coupled finer stock
grids. Sampling, timestep and normalizer approximation remain separate gates.
If accuracy gates fail, report that limitation and do not call the resulting
fit satisfactory or replace a production baseline.

VIX reuses the validated two-speed multilevel training/checking infrastructure.
Map23joint parameters to the14parameter VIX model; the nine off-window curve
pieces correctly have zero VIX derivatives. Reusing VIX checking streams is
explicitly disclosed. SPX checking uses previously stored separate old-state
streams plus new stress innovations, so it is not wholly fresh.

Validation: cached-model provenance and curve coverage; independent NumPy stock
replay; analytic derivatives for all23directions (including rate/OU-state chain
rules); four equality residuals; sampled M2 normalizer checks; signed price and
stock-moment diagnostics; coupled timestep checks; all four VIX smiles and all
seven SPX maturities plotted, with source/input hashes and original candidates
preserved. No new mathematical/Lean certification or global optimum is claimed.
