# Broader M2 calibration with restored precision and continuous controls

Completed 24 October 2026. This is a bounded exploration of all eleven agreed parameters; reproducing the old fit is an incumbent, not the research objective.

## What is held fixed and what changes

The complete OU architecture, geometric grids/rates/counts, correlations, projection powers, Q risk prices and C0=1/Y0=F0=0 stay fixed. All five score coefficients, both cap parameters, three forward-variance levels and variance floor are free. Optional orthogonal/echo channels remain disabled. Same 1294 OTM SPX cutoff-3 quotes at 29/57/119/211/302/666/1030 days, all 217 OTM VIX options at 27/55/118/153 days and four VIX forwards enter the objective. All market data are in sample.

Original per-expiry 512-path searching AND checking ensembles and 8x64 VIX states are restored. Equity shocks and importance weights replay exactly; full-state projections are translated by the fixed initial-state shift and cached. Three anchor controls per side are combined with deterministic continuous piecewise-linear weights path by path. These are fixed known-zero random controls, not smoothing of fitted prices or IVs. VIX forward controls remain unchanged. Negative estimates and shape failures remain scored and visible. No nearest-donor discontinuous switching is used.

The continuous controls still have possible convexity failures; the user explicitly requested continuing exploration with these limitations understood. This is not an arbitrage-consistent finite-sample price certificate. The natural-cubic control comparison was rejected because it worsened shape failures.

## Numerical implementation and retained attempts

Stock grid: 512 steps/trading day with exact event dates. Conditional VIX: GH32 in all three Gaussian coordinates, time order 8. Normalizer: GH48, with candidate-specific adaptive interpolation passing midpoint and additional independent same-order probes at 2e-7. These sampled screens do not establish quadrature convergence or a uniform interpolation bound.

The independent linear-Gaussian cluster score permits computing its GH32 convolution once as a one-dimensional function. A checked cubic-Hermite table then reuses that convolution inside the remaining GH32 x GH32 rule. The finite quadrature nodes/weights are unchanged. Interpolation failure or loss of independence falls back to direct tensor GH32. Three start vectors, all four VIX dates and representative lags were compared against direct values and derivatives before adoption. The final selected checking prices were also recomputed using the full direct tensor rule. Standard -O3 compilation is used; fast-math was benchmarked but not adopted.

The direct preliminary attempt runs/m2_controlled_readout_fit_001 was interrupted after preserving baseline/warm results because each readout required roughly 150 seconds. It is not a completed separate calibration. The completed run is m2_controlled_readout_fit_002. No old source manifests or saved prices were rewritten.

## Search and parameters

The first start preserves all ten old entries and moves only beta_J from zero to 0.0001. Two other starts use the previous broad-run candidate and a materially asymmetric beta_J=.5/q_J=.35 configuration. Each has up to six curve-prefit evaluations and twelve ALL-ELEVEN-parameter evaluations, with a total 45-minute / per-start 15-minute wall budget checked between evaluations. The best searching incumbent is retained across every evaluation; selection freezes before separate checking. This is local bounded least squares, not a global-optimum claim.

Three previously hit bounds were expanded once: beta_J lower 0.01 -> 0.0001, cap-sharpness upper 15 -> 30 and floor lower 1e-6 -> 1e-8. Other bounds remain unchanged.

| Parameter | Old reference | Selected | Bounds | Flag |
|---|---:|---:|---|---|
| cluster_loading | 0.2333221411 | 0.110792303 | [0.05, 0.8] |  |
| feedback_curvature | 0.06940074448 | 0.06285028058 | [0.005, 0.4] |  |
| feedback_shift | 2.935896338 | 3.11527412 | [0.05, 5] |  |
| spike_curvature | 0.2495993216 | 0.2053404866 | [0.002, 0.6] |  |
| spike_shift | 0 | 0.2997299726 | [0.0001, 2] |  |
| cap_level | 6.073741979 | 5.802414156 | [2, 8] |  |
| cap_sharpness | 7.956517501 | 2.024553487 | [2, 30] | near lower bound |
| xi_1 | 0.05762055055 | 0.04930159884 | [0.005, 0.2] |  |
| xi_2 | 0.0713408512 | 0.05859462501 | [0.005, 0.2] |  |
| xi_3 | 0.06321057443 | 0.05015265365 | [0.005, 0.2] |  |
| variance_floor | 0.0001 | 0.002481115229 | [1e-08, 0.0025] | near upper bound |

| Ensemble, same in-sample market quotes | Old reference | Positive-shift warm start | Selected |
|---|---:|---:|---:|
| Searching loss | 256.3288 | 256.3957 | 83.5283 |
| Separate checking loss | 305.2138 | 305.3228 | 103.3426 |

Selected checking prices inside bid–ask: **114/1511**, versus 89 for the old parameters repriced with the SAME new estimator. These losses are not direct comparisons with the historical nearest-donor estimator or the low-budget raw run.

| Checking component | Old reference | Selected |
|---|---:|---:|
| SPX | 233.095 | 65.8452 |
| VIX | 63.2158 | 36.6757 |
| forwards | 8.86973 | 0.788438 |
| curve_penalty | 2.54348e-05 | 1.57627e-05 |
| state_penalty | 0.0333333 | 0.0333333 |

## Every agreed maturity

| Instrument | Days | Quotes | Old checking interval MSE | Selected checking interval MSE | Selected inside bid–ask |
|---|---:|---:|---:|---:|---:|
| SPX | 29 | 277 | 209.9719 | 65.8768 | 25 |
| SPX | 57 | 298 | 262.1175 | 70.8852 | 20 |
| SPX | 119 | 228 | 438.9634 | 203.5815 | 7 |
| SPX | 211 | 141 | 310.1085 | 37.5724 | 11 |
| SPX | 302 | 167 | 390.8173 | 11.6574 | 18 |
| SPX | 666 | 112 | 16.6856 | 36.3043 | 5 |
| SPX | 1030 | 71 | 2.9998 | 35.0387 | 2 |
| VIX | 27 | 49 | 150.0014 | 29.5796 | 8 |
| VIX | 55 | 58 | 68.3682 | 48.5627 | 4 |
| VIX | 118 | 55 | 24.3675 | 26.7621 | 8 |
| VIX | 153 | 55 | 10.1262 | 41.7983 | 6 |

Selected VIX-forward errors, points: +0.3866, -0.4571, +0.6283, +0.1877.

## Robustness and what remains unresolved

The estimated checking loss change is -201.8711; descriptive eight-block paired jackknife SE is 51.5327. Leave-one-block changes range -234.9194 to -183.8953. These are few-block diagnostics on reused numerical checking streams, not a calibrated confidence certificate or market out-of-sample evaluation.

Direct tensor GH32 selected checking prices differ from the accelerated prices by at most 1.384e-09 quote scales; maximum forward difference is 4.007e-11 VIX points. Direct checking loss is 103.342624. Agreement checks the implemented same-order acceleration at this candidate, not GH convergence.

| Selected checking estimator shape | Order violations | Convexity violations | Negative prices |
|---|---:|---:|---:|
| SPX | 0 | 16 | 4 |
| VIX | 2 | 1 | 0 |

Saved-price/SE/objective replay passes for all baseline/warm/selected records and the final direct-GH32 check. 55 run source/input hashes, 86 cache source/input hashes and 260 immutable cache array hashes verify. The starting point preserves the ten other old coordinates exactly. Full end-to-end directional derivative check passes. Fourteen focused kernel, convolution and control/warm-start tests pass. These tests are implementation correspondence, not pricing certification.

Physical-measure compatibility has not been re-established. Stock discretization, conditional/time quadrature and rare-tail sampling remain separate accuracy questions. Boundary hits and finite local-search budgets prevent claiming identification or global model capacity. No new Lean theorem, formal audit or full-model certificate is asserted.

## Artifacts

- [Run and configuration](runs/m2_controlled_readout_fit_002/run.json).
- [Independent replay](runs/m2_controlled_readout_fit_002/score_audit.json).
- [Self-contained selected model](runs/m2_controlled_readout_fit_002/selected_model.json).
- [All maturity plots and guide](visualizations/controlled-readout-fit-009/README.md).
- [Prespecified plan](CONTROLLED_READOUT_EXPLORATION_PLAN_V2.md).


## Additional frozen-candidate diagnostics

Paired production GH32/GH48 probes at low/median/high saved VIX states (twelve states total) have maximum absolute change 6.21618e-06 VIX points. This is selected-state evidence, not surface convergence or unseen-tail coverage.

| VIX days | Second-moment target | Raw weighted estimate | Across-net SE | Largest contribution share |
|---|---:|---:|---:|---:|
| 27 | 493.0160 | 491.2401 | 6.3182 | 0.456% |
| 55 | 493.0160 | 499.1708 | 7.9016 | 0.545% |
| 118 | 585.9463 | 588.9976 | 9.0244 | 0.557% |
| 153 | 585.9463 | 582.6310 | 13.5033 | 0.528% |

There are 2 SPX checking contracts with no positive raw payoff; their controlled estimates and empirical SEs do not establish tail accuracy. VIX second moments agree with their normalization targets within the displayed sampling errors, but that identity alone cannot certify option prices. Raw stock forward/mass diagnostics are retained without rescaling weights or prices.

## Interpretation and next bounded questions

The enlarged readout family substantially reduces estimated joint loss in both ensembles, with aggregate SPX, VIX-option and VIX-forward gains. It is not merely reproducing the old curve fit: the warm-start curve prefit changes loss very little, while releasing all eleven parameters produces the large gain. The separate checking comparisons favor the selection in all eight leave-one-block diagnostics. This supports a provisional numerical improvement, not a global optimum or certified model fit.

Important trade-offs remain: checking SPX 666/1030-day losses worsen; VIX 118-day is slightly worse and 153-day is materially worse. Only 114/1511 prices are inside bid–ask, four checking estimates are negative, and finite-sample shape failures remain. The selected floor is near its upper bound and cap sharpness near its lower bound; all three local runs stopped at their evaluation limits. Do not call these parameters uniquely identified or infer global M2 infeasibility. Keep the new positive-shift vector as an exploratory working candidate alongside the old benchmark.

A next bounded hypothesis is to profile the newly hit floor/cap-sharpness boundaries while preserving the same OU cache and quote population, with a variance-excess parameterization if needed to maintain xi > floor. Another issue is candidate-appropriate continuous controls for the remaining signed/shape failures, with any changed estimator applied to both old and new benchmarks. Neither extension is launched in this checkpoint.