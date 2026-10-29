# Fixed-reservoir, eleven-parameter SPX/VIX calibration

Completed 24 October 2026. Exploratory in-sample calibration, with separate numerical checking on the same market quotes.

## Scope and implementation

The user confirmed eleven free parameters. OU grids/counts/rates, driver correlations, projection powers, Q risk prices and C0=1/Y0=F0=0 are fixed. Optional orthogonal/echo channels remain off. Beta_J is strictly positive. All 1294 OTM SPX cutoff-3 quotes at 29/57/119/211/302/666/1030 days, all 217 OTM VIX quotes at 27/55/118/153 days, and four VIX forwards enter the objective. Extremes remain excluded. No market holdout is introduced.

Raw positive importance-weighted payoffs replace the historical strike-dependent control estimator for this run. The old zero-shift reference is repriced on the SAME new paths and quadrature as the new candidate. Its losses here must not be compared directly with historical controlled-estimator loss 308.64. No smoothing or clipping is applied. Some raw estimates are zero because the sampled paths have no positive payoff; their IVs are marked unavailable in the plots, while their prices and errors remain in the scores. These cases do not establish a true zero option price.

The immutable cache stores reservoir projections and equity increments, not candidate G/M_Q. Every candidate recomputes G, its date-zero Q normalizer, stock values and conditional VIX. Stock updates remain at 512 steps per trading day on exact event dates. One common stock grid preserves dependence across SPX maturities.

## Parameters

| Parameter | Old reference | Selected | Search interval | Boundary flag |
|---|---:|---:|---|---|
| cluster_loading | 0.2333221411 | 0.1396709482 | [0.05, 0.8] |  |
| feedback_curvature | 0.06940074448 | 0.2328840917 | [0.005, 0.4] |  |
| feedback_shift | 2.935896338 | 1.166236318 | [0.05, 5] |  |
| spike_curvature | 0.2495993216 | 0.2715626468 | [0.002, 0.6] |  |
| spike_shift | 0 | 0.01 | [0.01, 2] | Near lower bound |
| cap_level | 6.073741979 | 5.267358224 | [2, 8] |  |
| cap_sharpness | 7.956517501 | 15 | [2, 15] | Near upper bound |
| xi_1 | 0.05762055055 | 0.05415407103 | [0.005, 0.2] |  |
| xi_2 | 0.0713408512 | 0.04018991223 | [0.005, 0.2] |  |
| xi_3 | 0.06321057443 | 0.03631150499 | [0.005, 0.2] |  |
| variance_floor | 0.0001 | 1e-06 | [1e-06, 0.0025] | Near lower bound |

Boundaries and finite optimization budget prevent interpreting this vector as a unique, globally identified calibration. The positive beta_J lower bound 0.01 is a declared numerical search choice. A boundary preference does not establish that a particular economic asymmetry is absent.

## Same-estimator comparisons

| Ensemble | Old reference loss | Selected loss | Change | Old / selected inside bid–ask |
|---|---:|---:|---:|---|
| Searching simulations, conditional GH12 | 2168.6127 | 586.9130 | -1581.6997 | 48 / 38 |
| Separate simulations, conditional GH24 | 1044.1042 | 906.6085 | -137.4957 | 72 / 61 |
| Same checking states, conditional GH12 | 1043.5309 | 902.6194 | -140.9115 | 73 / 60 |

The selected candidate has a lower estimated loss on the separate checking ensemble, but the uncertainty diagnostics below do not establish a reliable improvement. These are numerical comparisons on the SAME in-sample market observations, not market out-of-sample results.

| Checking loss component | Old reference | Selected |
|---|---:|---:|
| SPX | 969.37 | 781.98 |
| VIX | 64.3571 | 102.94 |
| forwards | 10.3438 | 21.6554 |
| curve_penalty | 2.54348e-05 | 2.1004e-05 |
| state_penalty | 0.0333333 | 0.0333333 |

## Every agreed maturity

| Instrument | Days | Quotes | Old interval MSE | Selected interval MSE | Selected inside bid–ask |
|---|---:|---:|---:|---:|---:|
| SPX | 29 | 277 | 1069.0393 | 2229.4354 | 6 |
| SPX | 57 | 298 | 655.4947 | 223.4045 | 24 |
| SPX | 119 | 228 | 2590.4659 | 663.0517 | 1 |
| SPX | 211 | 141 | 1383.0759 | 640.4504 | 17 |
| SPX | 302 | 167 | 1057.0735 | 1214.0154 | 0 |
| SPX | 666 | 112 | 25.6139 | 297.2131 | 1 |
| SPX | 1030 | 71 | 4.8266 | 206.2900 | 3 |
| VIX | 27 | 49 | 152.8906 | 52.9824 | 5 |
| VIX | 55 | 58 | 55.7278 | 51.5489 | 4 |
| VIX | 118 | 55 | 39.3167 | 150.0563 | 0 |
| VIX | 153 | 55 | 9.4934 | 157.1711 | 0 |

Selected VIX forward errors (points): +2.0945, +1.0839, -2.9325, -2.7376.

## Numerical checks and limits

Six focused Python tests pass: production readout and Gaussian-expectation correspondence, all eleven stock derivatives against finite differences, and existing all-maturity membership/objective checks. Two end-to-end directional Jacobian checks passed before search. These checks do not certify pricing convergence. No new Lean theorem or compiler audit was requested or claimed; existing formal identities retain their prior conditional scope.

Saved prices and SEs replay, all 49 captured source/input hashes match, and 67 recorded seeds are distinct across searching/checking streams. All six score records reproduce. Raw positive estimates have zero within-side strike monotonicity violations at every maturity. Replay does not independently regenerate the simulations.

The descriptive paired eight-block jackknife SE of the checking loss change is 876.6969. This SE is much larger than the estimated loss reduction, and one leave-one-block comparison reverses sign. The numerical improvement is not robustly established. Few blocks and rare-tail sampling limit interpretation; this is not a calibrated confidence certificate.

On matching checking states, selected GH24 versus GH12 price changes have median 0 and maximum 1.923 quote scales across all options (the median includes unchanged SPX prices); the VIX-only median is 0.1996. Search used GH12; check used GH24, with GH12 checking records retained to distinguish this change from outer sampling. Neither order proves convergence.

Normalizers use GH48 and interpolation on 257 transformed-time knots. Maximum independent same-GH-order interpolation probe errors: search/baseline: 4.46e-09, search/selected: 1.28e-06, check/baseline: 4.46e-09, check/selected: 1.28e-06, check_gh12/baseline: 4.46e-09, check_gh12/selected: 1.28e-06. These are sampled interpolation checks, not a bound on GH bias. The selected candidate exceeds the earlier 2e-7 interpolation screen; normalizer accuracy is not cleared.

Searching stock budget is 128 paths; checking is 256. VIX uses eight scrambles of 16/32 outer points respectively. Normalizer, conditional/time quadrature, stock discretization, rare tails and adaptive parameter selection remain separate sources of uncertainty. Sampling bands exclude deterministic errors. Physical-measure compatibility has not been re-established for the selected readout.

Both starts were bounded to 35 evaluations, within the declared 30-minute optimization budget. Total run wall time: 745.00 seconds. Optimizer outcomes and every candidate remain in the run directory. This is a bounded attempt, not proof of global model feasibility or infeasibility.

## Artifacts

- Run: [m2_full_readout_fit_001](runs/m2_full_readout_fit_001/run.json).
- Replay audit: [score_audit.json](runs/m2_full_readout_fit_001/score_audit.json).
- Search plan: [FULL_READOUT_FIT_PLAN.md](FULL_READOUT_FIT_PLAN.md).
- [Plots and reading guide](visualizations/full-readout-fit-008/README.md).
