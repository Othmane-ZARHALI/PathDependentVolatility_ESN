# VIX-only recalibration with the existing 16-level curve

Completed 26 October 2026. Run `runs/m2_vix_only_curve_fit_002`; plan
[VIX_ONLY_CURVE_PLAN.md](VIX_ONLY_CURVE_PLAN.md). SPX remains paused.

The user was right: the earlier maturity-curve calibration already fitted 16
forward-variance levels, expanding the original three-piece curve. Later
stochastic-factor comparisons froze that selected curve to isolate their effects.
This experiment releases the curve again while fitting the existing post-cap
factor. It introduces no new reservoir bank or readout.

A VIX averaging window is [T,T+30 calendar days], not a curve parameter count.
All 16 original intervals remain represented. Seven pieces intersect the four
windows: 27–29, 29–55, 55–57, 57–85, 118–119, 119–148, 153–183 days.
They move independently, without broad-window tying, monotonicity or smoothing.
The other nine have exactly zero effect on these VIX prices and retain their old
values. The 55–57-day piece contributes to both the first and second windows.
The complete levels are in [curve_levels.csv](runs/m2_vix_only_curve_fit_002/curve_levels.csv).

## What was fitted

Same capped old readout, positive beta_J=0.369729973, Q reservoir, premiums and
initial state. Floor fixed at 0.00264. Active variance levels range from floor to
0.20. Normalized independent OU factor L=exp(gamma Z-gamma² Var(Z)/2) remains
outside the old cap; no new cap is added. Its amplitude gamma ranges from 0 to 3
and its half-life from 1 to 365 days. Six joint starts, at most 120 function
calls each; separate curve-only and fixed gamma=1/7-day-half-life fits.

The objective retains all 217 native OTM bid–ask observations and four forwards,
with existing quote scales, equal total weight per expiry and forward scale
0.5 point. It has no SPX term, curve penalty or extra wing weighting. Canonical
calls use model forwards for model parity, market forwards for market parity;
plots use observed market forwards throughout. Signed prices remain saved.

## Results

Lower price loss is better. Training reuses the previous study's checking draws
(4096 per expiry). Fresh numerical checking uses eight new 1024-point scrambles
per expiry, generated after the parameters were frozen, with exact joint37D
importance-mixture weights. No self-normalization or control refitting.

| Candidate | Training loss | Fresh loss | Fresh quotes inside bid–ask |
|---|---:|---:|---:|
| baseline | 26.852498 | 27.349124 | 18/217 |
| curve_only | 20.139371 | 20.524504 | 20/217 |
| fixed_factor_curve | 21.630542 | 20.949991 | 21/217 |
| factor_curve | 19.946095 | 19.808593 | 21/217 |

Freeing the curve yields a clear reduction: fresh loss change -6.82462,
paired eight-net jackknife SE 0.16812. Adding the fitted factor to the free curve
changes loss by only -0.71591, SE 0.54756. These are numerical sensitivity estimates
on the same calibration quotes, not out-of-sample market results or confidence
certificates. The old readout/floor remain fixed throughout.

The joint fit chooses gamma=0.36974 and a 365-day half-life, at its search upper
bound. Several starts instead converge to negligible factor amplitude. This
finite search does not establish a global optimum or that other bounds/readouts
cannot improve the fit.

The wings remain visibly too flat. Measured IV slopes versus log(K/F), using all
retained strikes above 2F:

| Expiry days | Baseline | Curve only | Fixed strong factor + curve | Fitted factor + curve |
|---|---:|---:|---:|---:|
| 27 | 36.649 | 37.253 | 31.600 | 30.594 |
| 55 | 19.636 | 26.749 | 29.866 | 26.589 |
| 118 | 11.493 | 15.846 | 22.829 | 16.650 |
| 153 | 9.992 | 9.651 | 11.026 | 8.711 |

Market slopes are approximately 52.38, 38.31, 31.23, 27.32 respectively. The
lowest overall price loss is therefore not the best wing-shape result. Fresh
wing IV MAE is 33.91 baseline, 29.43 curve-only, 27.72 fixed strong factor,
30.75 fitted factor (volatility percentage points). No prices or IVs were
smoothed or clipped; all four candidates have zero unavailable fresh IVs.

The optimizer concentrates variance into short intervals: the 118–119-day piece
hits the existing 0.20 variance bound (44.7% square-root forward variance); the
55–57-day piece also increases sharply. This is an exploratory irregular curve,
not a production recommendation. More free levels alone do not ensure a sensible
curve or the required smile shape.

The fitted-factor curve implies exact E[VIX²] values 476.64, 559.10, 606.90,
636.69. The last two remain below the archived necessary bounds for fitting every
native quote plus the exact observed forward (648.63, 695.30). This fit still
trades off mismatches; it does not invalidate those necessary bounds. Its forward
errors are about +0.39, +0.30, +0.63, +0.95 VIX points. Only 21/217 quotes lie
inside bid–ask, so no acceptable full-smile fit is claimed.

## Saved plots

- [All four smiles, full strike ranges](visualizations/vix-only-curve019/vix-all-strikes.png)
- [Right wings](visualizations/vix-only-curve019/vix-right-wing.png)
- [All 16 variance pieces](visualizations/vix-only-curve019/variance-curve.png)
- [Forwards and exact second moments](visualizations/vix-only-curve019/forwards-and-moments.png)

PNG and SVG versions are saved. Sampling bands are pointwise and exclude
systematic integration error and independent-pilot uncertainty.

## Verification and limits

Old baseline price replay agrees to numerical precision. Analytic objective
Jacobian checked in all nine fitted directions; maximum relative error below
9e-5. Independent direct nodal accumulation replays all eight candidate/stream
sets, native prices and forwards; maximum price difference below 1.4e-15.
Exact second moments, all 16 levels, unchanged nine off-window levels, fixed old
readout/floor, positive beta_J, 282 source hashes and five input hashes are checked
in `audit.json` (PASS). A fresh C++ library was compiled from captured sources.
Sampled GH32/time8 to GH48/time16 refinements change VIX values by less than
4.6e-6 VIX points over 30–32 selected states per expiry. This is a sampled
refinement, not a uniform integration bound. Likelihood-weighted moments, tail
payoff concentration and ESS are saved per candidate; finite-sample tail error
remains a limitation. No new Lean theorem is asserted.

Run001 is preserved as a preliminary numerically identical attempt. Its sparse
quadrature probes emitted empty-net diagnostic warnings, affecting unused probe
moment SEs, not VIX values or fitted results. Run002 excludes empty probe nets and
compiles its own pricing library; it supersedes001. Repeated fits reproduce the
same selected parameters and fresh prices. No model was adopted or historical
selected_model.json replaced. No SPX paths were read, simulated or priced.

The next model candidate remains the previously discussed mixed-response
alternative, keeping the curve free. It has not been implemented or authorized
by this run. The present bounded fit establishes improvement in the overall
price objective, not a solution of the right-wing problem.
