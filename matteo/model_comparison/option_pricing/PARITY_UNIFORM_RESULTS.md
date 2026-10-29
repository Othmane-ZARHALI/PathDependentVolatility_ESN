# Common calls and 4096 paths at every strike

Completed fixed-parameter estimator experiment, 24 October 2026. No new calibration.

Run: `m2_parity_uniform_001`. All 1294 original OTM SPX quotes, 217 original OTM VIX quotes and four forwards remain. All seven SPX and four VIX maturities are included. The old reference and previous selected parameters are repriced on identical fresh samples.

## What was implemented

- One call payoff at every strike. Market OTM puts become synthetic calls using their observed forward; no ITM observations are added.
- All **4096 paths at every strike**, following the user's update. VIX uses eight independent 512-point scrambles.
- New reusable exogenous OU/equity caches, independent frozen-control pilots and common call controls across the forward.
- A tested deterministic translation for reusing the same innovations from another full OU starting state. The current initial state stays fixed.

Original OU grids, dynamics, Q risk prices, correlations and initial state are unchanged. Keep conditional GH32, normalizer GH48, VIX time order 8 and 512 stock steps per trading day. The earlier adaptive-prefix proposal is superseded and unused; nested-prefix outputs are diagnostics only.

## Fixed-parameter comparison

| Quantity | Old reference | Previous selected parameters |
|---|---:|---:|
| Common-call loss | 1941.75 | 98.1192 |
| Original-OTM loss | 1724.67 | 85.5281 |
| Inside common-call bid/ask | 84 | 93 |
| Inside original-OTM bid/ask | 89 | 100 |
| Unavailable common-call IVs | 0 | 2 |
| Negative reconstructed OTM prices | 0 | 0 |

Both parameter sets use the SAME new estimator and samples. Earlier controlled-estimator losses (including 103.34 at 512 paths) are not direct comparators: controls, samples and, for VIX puts, price convention have changed.

| Instrument | Days | Old common-call MSE | Selected common-call MSE | Inside bid/ask | Missing IVs | Median SE: 4096 / 512 |
|---|---:|---:|---:|---:|---:|---:|
| SPX | 29 | 164.445 | 84.202 | 26/277 | 0 | 0.384 |
| SPX | 57 | 191.942 | 45.729 | 18/298 | 0 | 0.350 |
| SPX | 119 | 334.520 | 77.289 | 6/228 | 0 | 0.363 |
| SPX | 211 | 333.696 | 45.094 | 6/141 | 0 | 0.353 |
| SPX | 302 | 10548.111 | 66.903 | 12/167 | 0 | 0.373 |
| SPX | 666 | 5.856 | 19.219 | 4/112 | 0 | 0.346 |
| SPX | 1030 | 0.957 | 10.184 | 5/71 | 0 | 0.356 |
| VIX | 27 | 452.764 | 75.824 | 5/49 | 0 | 0.404 |
| VIX | 55 | 113.428 | 25.434 | 7/58 | 2 | 0.297 |
| VIX | 118 | 371.110 | 53.506 | 2/55 | 0 | 0.417 |
| VIX | 153 | 177.963 | 34.928 | 2/55 | 0 | 0.254 |

The 512 and 4096 estimates here are nested portions of the SAME fresh ensemble with the SAME controls and parameters. The ratio is descriptive. It does not establish population convergence or precision for every strike. Using eight times as many iid paths would typically reduce SE by about sqrt(8), but observed tail variance can change with additional draws.

The large old-reference score is sensitive to a rare SPX302 path: normalized stock is 56.63 under the old parameters and 5.58 under the selected parameters on the SAME innovations. Its old weighted stock ratio is 48.17; far-wing baseline SE reaches about 450 quote scales. Controls are frozen around the selected model, so baseline variance can be much larger. The paired loss-change jackknife SE is 1687.09, versus observed reduction 1843.63. All eight leave-one-block comparisons favor the selected vector, but the magnitude is highly unstable. Do not interpret the apparent roughly 95% score reduction as a precise population improvement.

Selected sample shape diagnostics: **1 SPX strike-order and 72 convexity failures**, **0 VIX order/convexity failures**. There are no negative reconstructed OTM prices and no zero raw-payoff-hit contracts, but sampling precision remains several quote scales for many SPX contracts. The old SPX57 missing IVs at strikes 8000/8100/8200/8300 are all available. The only two selected unavailable common-call IVs are VIX55 strikes 14.5/15: positive model calls are below intrinsic using the observed forward because the model VIX forward is lower.

## Put–call parity and VIX forwards

All prices are undiscounted: $C_{market}=P_{market}+F_{market}-K$. Shift put bid/mid/ask by the same known amount; retain spread widths and quote scales.

Model puts are $P_{model}=C_{model}-F_{model}+K$. For SPX, the normalized adapted discrete stock has expectation 1, so model and market forwards agree by construction of this price model. For VIX, $F_{model}=E[VIX_T]$ is an actual model output, not the observed forward.

Therefore, for converted VIX puts, $e_{call}=e_{put}+F_{model}-F_{market}$. The common-call objective differs from the former OTM objective if the VIX forward is wrong. Both losses, native OTM prices/SEs and all model-forward errors are retained. A continuous common-call IV display does not establish a correct VIX forward.

| VIX days | Observed forward | Selected model forward | Sampling SE | Error |
|---:|---:|---:|---:|---:|
| 27 | 20.100000 | 20.492885 | 0.020371 | +0.392885 |
| 55 | 20.730000 | 20.213338 | 0.017773 | -0.516662 |
| 118 | 21.190000 | 21.840352 | 0.017186 | +0.650352 |
| 153 | 21.720000 | 21.844290 | 0.011386 | +0.124290 |

## Sampling SE and controls

SPX: for corrected independent observations $Z_i(K)$, use $\widehat C=N^{-1}\sum_iZ_i$ and $\widehat{SE}=s_Z/\sqrt{4096}$. VIX: calculate a complete price in each of eight independently randomized nets, then use SD of those eight prices divided by sqrt(8). Native VIX put SE includes its covariance with the model forward.

The IV bands invert price ±2 sampling SE using the observed market forward and the same call convention everywhere. These are pointwise sampling bands, conditional on frozen pilot controls. They omit pilot/parameter uncertainty and deterministic numerical bias. Prices outside Black bounds have unavailable IVs. Signed prices stay in the tables and loss; no clipping or IV smoothing. A price-band endpoint outside Black bounds also has no IV; missing shading is not zero uncertainty.

SPX uses an analytic deterministic-lognormal-reference call driven by the same equity shocks, plus known-zero importance-mass, normalized-stock, reference-call and two predictable hedging-gain controls. A fresh independent 256-path pilot fits coefficients at nine common call anchors spanning the forward. The hedge holdings change about four times per trading day on the unchanged fine stock path. Continuous interpolation applies to zero-mean correction draws, never prices. VIX controls use Gaussian coordinates and centered squares with coefficients fitted on independent 8x128-point pilots.

The simple reference control alone was screened on OLD draws and was noisier than the old richer controls (57-day strike 8000 SE 1.921 versus 0.485). The independent hedging controls were added before examining fresh production results. Estimator replacement remains a material change; the new results are not a test of sample count alone.

## Reusing paths from a different starting state

For fixed exogenous OU dynamics, $X_t^{new}-X_t^{old}=e^{-\Lambda t}(X_0^{new}-X_0^{old})$ when the driving innovations are shared. Apply its projection to cached C/Y/F paths. Equity Brownian increments, proposal components, scores and importance weights stay unchanged. For VIX, translate the terminal OU state by the same expression at expiry.

Recompute volatility, its date-zero conditional Q normalizer, downstream stock paths and conditional VIX values. The reusable quantities are exogenous OU/equity innovations, not already computed stock returns or option prices. Full initial OU states are required; C/Y/F alone need a chosen map into the full state. `dated.initial_state([C,Y,F])` supplies the project's minimum stationary-Mahalanobis-norm map.

`reusable_ou_initial.py` has three tests: full/projected translations against fresh identical-seed simulations, unchanged Gaussian covariance, and unchanged importance weights/scores/components. States match within 3e-14. The agreed initial state is not added to the eleven calibration parameters in this experiment.

## Verification, limitations and artifacts

Fourteen focused tests pass. Independent saved-price/SE/parity/market-conversion/both-objective replay passes; 73 distinct sources verify. Run inputs are SHA256-frozen before the first objective evaluation and verified unchanged afterwards. The interrupted cache preparation resumed only missing chunks; its original source hashes and JSON-normalized model settings verified. The original resume check confused Python tuples with JSON lists, so the unchanged worker was resumed through a separately snapshotted wrapper. Queued runners 001/002 were interrupted before pricing; neither is another calibration.

| Candidate | Instrument | Days | Contracts with zero raw OTM payoff hits |
|---|---|---:|---:|
| selected | SPX | 29 | 0 |
| selected | SPX | 57 | 0 |
| selected | SPX | 119 | 0 |
| selected | SPX | 211 | 0 |
| selected | SPX | 302 | 0 |
| selected | SPX | 666 | 0 |
| selected | SPX | 1030 | 0 |
| selected | VIX | 27 | 0 |
| selected | VIX | 55 | 0 |
| selected | VIX | 118 | 0 |
| selected | VIX | 153 | 0 |
| baseline | SPX | 29 | 0 |
| baseline | SPX | 57 | 0 |
| baseline | SPX | 119 | 0 |
| baseline | SPX | 211 | 0 |
| baseline | SPX | 302 | 0 |
| baseline | SPX | 666 | 0 |
| baseline | SPX | 1030 | 0 |
| baseline | VIX | 27 | 0 |
| baseline | VIX | 55 | 0 |
| baseline | VIX | 118 | 0 |
| baseline | VIX | 153 | 0 |

An empirical zero SE or zero-hit tail is not evidence that the true tail contribution vanishes. Known-zero controls and one payoff convention do not guarantee finite-sample positivity, monotonicity or convexity. Shape failures are saved in score JSON files. A larger path count is not a price-convergence certificate. No new optimization, P compatibility assessment or Lean certificate is claimed.

Plots: [all seven SPX maturities](visualizations/parity-uniform-010/spx-all-maturity-fit.png), [all four VIX maturities](visualizations/parity-uniform-010/vix-all-maturity-fit.png), [actual VIX forwards](visualizations/parity-uniform-010/forwards-and-path-budgets.png), [sampling precision](visualizations/parity-uniform-010/sampling-precision.png). Signed price-error panels and SVG exports accompany them.

Full outputs: `runs/m2_parity_uniform_001`. Reusable cache: `runs/m2_parity_adaptive_cache_001`; independent frozen controls: `runs/m2_parity_adaptive_controls_001`. Names retain the initial adaptive proposal; active production allocation is uniformly 4096.
