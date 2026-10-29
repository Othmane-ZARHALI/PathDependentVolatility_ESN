# Maturity-aligned forward variance calibration

Completed 25 October 2026. The previous three-piece curve is exactly nested in the new 16-piece curve. The fit frees 24 parameters: 16 variance levels, seven readout/cap coefficients and the variance floor. OU rates, reservoir dynamics, random draws, risk prices, correlations, initial state and geometric grids remain fixed.

Original OTM observations only: SPX 29/57/119/211/302/666/1030 days (cutoff 3), VIX 27/55/118/153 days, 1511 option quotes and four VIX forwards. Calls form the common payoff; market puts use market-forward parity and native model puts use model-forward parity. Every strike uses all 4096 paths. All market observations are in sample.

## Actual full-price results

| Training criterion | Previous | Selected |
|---|---:|---:|
| Common-call joint | 95.117702 | 89.320446 |
| Original-OTM joint | 84.982024 | 70.977709 |
| SPX | 49.518173 | 43.777476 |
| Common-call VIX | 44.728178 | 44.619039 |
| Original-OTM VIX | 34.592500 | 26.276301 |
| VIX forwards | 0.838006 | 0.890564 |

Original-OTM VIX-plus-forward training loss changes from 35.430506 to 27.166865 (-23.32%). Accepted full-price steps: 1. This is a bounded exploration, not a convergence or global-optimum claim.

| Proposal | Curve only | Predicted joint | Actual joint | Actual native VIX + forwards | Accepted |
|---:|:---:|---:|---:|---:|:---:|
| 0 | True | 89.305906 | 89.320446 | 27.166865 | True |
| 1 | False | 89.309281 | 89.760750 | 25.581229 | False |
| 2 | False | 89.317655 | 89.392468 | 26.233363 | False |

SLSQP minimizes the native VIX-plus-forward objective on a local affine price/forward model. Both joint and both VIX-plus-forward criteria constrain the inner solve. Nonlinear VIX prices screen proposals before full SPX/VIX pricing. Only full nonlinear scores can accept a step. The first proposal frees curve/floor parameters; subsequent proposals free all 24. Only the first proposal was accepted, so the selected readout/cap coefficients remain unchanged despite having been searched in the later proposals.

## Separate reused VIX numerical checking

| Criterion | Before | After | Change | Paired jackknife SE of change |
|---|---:|---:|---:|---:|
| canonical VIX + forwards | 46.886107 | 46.285049 | -0.601058 | 0.781157 |
| native VIX + forwards | 35.862036 | 27.523027 | -8.339009 | 0.191220 |

Selection froze before these VIX prices were checked. These eight independent scrambled-net streams have been used in previous investigations; this is descriptive numerical checking, not an untouched holdout or market out-of-sample test. SPX checking was not repeated. Leave-one-net changes and second moments are preserved in [uncertainty.json](runs/m2_maturity_curve_fit_001/checking_vix/uncertainty.json).

## Forward variance levels

| Interval in days | Previous annual variance | Selected annual variance | Selected square-root scale |
|---|---:|---:|---:|
| 0–27 | 0.04898962 | 0.04898962 | 22.134% |
| 27–29 | 0.04898962 | 0.04298962 | 20.734% |
| 29–55 | 0.04898962 | 0.04752808 | 21.801% |
| 55–57 | 0.04898962 | 0.05498962 | 23.450% |
| 57–85 | 0.04898962 | 0.05350362 | 23.131% |
| 85–91.25 | 0.04898962 | 0.04898962 | 22.134% |
| 91.25–118 | 0.05827299 | 0.05827299 | 24.140% |
| 118–119 | 0.05827299 | 0.06427299 | 25.352% |
| 119–148 | 0.05827299 | 0.05786689 | 24.056% |
| 148–153 | 0.05827299 | 0.05827299 | 24.140% |
| 153–183 | 0.05827299 | 0.06109326 | 24.717% |
| 183–211 | 0.05827299 | 0.05827299 | 24.140% |
| 211–273.75 | 0.05827299 | 0.05827299 | 24.140% |
| 273.75–302 | 0.05257415 | 0.05257415 | 22.929% |
| 302–666 | 0.05257415 | 0.05257415 | 22.929% |
| 666–1030 | 0.05257415 | 0.05257415 | 22.929% |

Every level satisfies xi_i >= variance_floor; the levels may rise or fall. These square-root scales describe the expected variance curve, not individual option implied volatilities. Close boundaries can make individual levels weakly identified.

| Readout/cap/floor parameter | Previous | Selected |
|---|---:|---:|
| cluster_loading | 0.07579230299 | 0.07579230299 |
| feedback_curvature | 0.04995142993 | 0.04995142993 |
| feedback_shift | 3.39527412 | 3.39527412 |
| spike_curvature | 0.2078592781 | 0.2078592781 |
| spike_shift | 0.3697299725 | 0.3697299725 |
| cap_level | 5.502780057 | 5.502780057 |
| cap_sharpness | 2 | 2 |
| variance_floor | 0.00254 | 0.00264 |

Selected minimum instantaneous volatility is 5.13809%; beta_J is 0.36972997 > 0. Floor bounds remain 5–6% volatility (the 6% upper endpoint is the earlier implementation choice).

## All retained maturities

| Instrument | Days | Native interval MSE before | After | Selected prices in bid/ask |
|---|---:|---:|---:|---:|
| SPX | 29 | 110.54088 | 104.88029 | 28/277 |
| SPX | 57 | 39.94458 | 42.06696 | 20/298 |
| SPX | 119 | 60.63675 | 36.68916 | 31/228 |
| SPX | 211 | 32.48490 | 25.04300 | 22/141 |
| SPX | 302 | 71.62859 | 64.36789 | 8/167 |
| SPX | 666 | 19.52761 | 21.01668 | 6/112 |
| SPX | 1030 | 11.86391 | 12.37835 | 7/71 |
| VIX | 27 | 28.01390 | 26.54672 | 4/49 |
| VIX | 55 | 45.87496 | 20.56882 | 2/58 |
| VIX | 118 | 27.90219 | 28.10265 | 6/55 |
| VIX | 153 | 36.57895 | 29.88701 | 5/55 |

## Verification and numerical limits

Ten focused tests pass. The saved-price audit replays prices, standard errors, parity, losses, acceptance decisions, 70 local trials and 74 frozen sources/inputs. Curve/floor production VIX sensitivities passed finite differences. Baseline SPX prices agree with the preceding selected fit to the recorded tolerance; splitting VIX quadrature at the new boundaries changes baseline prices only by the recorded tiny refinement error.

Selected training estimates: 0 negative native prices; 0 unavailable common-call IVs; 139/1511 native estimates inside bid/ask. Signed estimates, invalid IVs, strike-shape and raw-tail diagnostics are retained. The VIX wing mismatch must be judged from the plots and strike-region diagnostics, not just the aggregate loss.

SPX SE uses the standard deviation across 4096 independent controlled path contributions divided by sqrt(4096). VIX SE uses the standard deviation of the eight independent 512-point scrambled-net means divided by sqrt(8), including model-forward covariance when recovering native puts. Bands are pointwise sampling uncertainty only. They exclude parameter uncertainty, frozen pilot uncertainty, normalizer/quadrature bias and stock timestep error.

GH32 conditional expectations, GH48 normalization, time order 8 and the cached 512-step stock grid remain unchanged. New curve boundaries split VIX integration panels. Stock uses the cached left-point mesh; deterministic curve-integral discrepancies are recorded in grid_baseline.json and grid_selected.json, and are not a full price-convergence bound. No Lean/P-compatibility certificate follows.

The pause interrupted the first proposal, not the completed baseline. The continuation verified frozen original source hashes, restored the baseline and reproduced the first proposal before resuming full pricing. Baseline and original pause/source records are preserved; continuation sources are separately snapshotted. Timing in the final run records continuation wall time; it does not include time spent paused or all interrupted work.

[All seven SPX maturities](visualizations/maturity-curve-013/spx-all-maturity-fit.png) · [All four VIX maturities](visualizations/maturity-curve-013/vix-all-maturity-fit.png) · [Separate VIX checking](visualizations/maturity-curve-013/vix-separate-checking.png) · [Forward variance curve](visualizations/maturity-curve-013/forward-variance-curve.png). Blue is selected, orange is the preceding fit. All PNG/SVG exports and signed-price plots are in the visualization folder.

[Selected model](runs/m2_maturity_curve_fit_001/selected_model.json) · [Audit](runs/m2_maturity_curve_fit_001/audit.json) · [Plan](MATURITY_CURVE_FIT_PLAN.md).

## Actual VIX forwards

| Streams | Days | Market | Previous | Selected |
|---|---:|---:|---:|---:|
| Training | 27 | 20.1000 | 20.4490 | 20.2267 |
| Training | 55 | 20.7300 | 20.1784 | 21.0936 |
| Training | 118 | 21.1900 | 21.8231 | 21.7815 |
| Training | 153 | 21.7200 | 21.8253 | 22.3464 |
| Separate checking | 27 | 20.1000 | 20.4636 | 20.2401 |
| Separate checking | 55 | 20.7300 | 20.1899 | 21.0992 |
| Separate checking | 118 | 21.1900 | 21.8434 | 21.8015 |
| Separate checking | 153 | 21.7200 | 21.8013 | 22.3194 |

## VIX implied-volatility errors on separate checking streams

| Days | Strike region | Previous IV MAE, percentage points | Selected IV MAE | Defined selected IVs |
|---:|---|---:|---:|---:|
| 27 | 1 < K/F <= 2 | 14.446 | 12.942 | 25/25 |
| 27 | K/F <= 1 | 47.771 | 38.185 | 10/10 |
| 27 | K/F > 2 | 26.431 | 29.744 | 14/14 |
| 55 | 1 < K/F <= 2 | 12.437 | 9.923 | 24/24 |
| 55 | K/F <= 1 | 2.055 | 27.564 | 13/13 |
| 55 | K/F > 2 | 41.868 | 37.724 | 21/21 |
| 118 | 1 < K/F <= 2 | 10.350 | 10.545 | 19/19 |
| 118 | K/F <= 1 | 23.230 | 21.944 | 12/12 |
| 118 | K/F > 2 | 35.648 | 35.396 | 24/24 |
| 153 | 1 < K/F <= 2 | 15.194 | 12.349 | 20/20 |
| 153 | K/F <= 1 | 4.407 | 17.911 | 12/12 |
| 153 | K/F > 2 | 34.722 | 33.485 | 23/23 |

IV errors use the observed-market-forward common-call convention. Large far-wing IV errors can remain even when aggregate interval-price loss falls.

Selected VIX 153-day second moment: estimate 606.208, exact forward-curve target 610.933, reported eight-net SE 1.635. The deficit is about 2.89 reported SE; rare-tail/numerical concerns are not cleared by this fit.

All eight PNG figures were visually inspected; PNG/SVG hashes and source hashes are checked in validation/plot_review.json.

Selected training shape diagnostics: SPX: 1 monotonicity and 66 convexity flags; VIX: 0 monotonicity and 0 convexity flags. These finite-estimator flags are not model-level arbitrage proofs.
