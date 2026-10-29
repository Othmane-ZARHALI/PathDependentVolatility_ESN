# Eleven-parameter joint recalibration with 4096 paths

Completed bounded calibration attempt. All eleven agreed parameters were free; OU reservoir, Q prices, initial state and geometric grids remain fixed.

Run: `m2_parity_parameter_fit_003`. Start at the previous fitted positive-shift vector, not the historical zero-shift reference. Keep all 1294 SPX cutoff 3 OTM observations at seven maturities, 217 VIX OTM observations at four maturities and four VIX forwards. Every strike uses 4096 paths; VIX uses eight independent 512-point nets. GH32/GH48/time8/stock512 remain unchanged.

## Results on identical inputs

Loss measures squared price error outside bid/ask, divided by the existing quote scale, then averaged equally over maturities within SPX and VIX. A price inside bid/ask contributes zero. Forward error is scaled by 0.5 VIX points. These are price scores, not percentages or IV errors.

| Ensemble / convention | Before | After |
|---|---:|---:|
| Training: Common-call loss | 98.119205 | 95.148310 |
| Training: Original-OTM loss | 85.528066 | 85.038631 |
| Training: SPX component | 49.802995 | 49.491692 |
| Training: VIX-option component | 47.423157 | 44.791849 |
| Training: Forward component | 0.859704 | 0.831425 |
| Training: original-OTM VIX-option component | 34.832017 | 34.682170 |
| Training: original-OTM VIX plus forwards | 35.691722 | 35.513595 |
| Independent checking: Common-call loss | 96.469009 | 92.240237 |
| Independent checking: Original-OTM loss | 82.828650 | 81.231878 |
| Independent checking: SPX component | 46.691042 | 45.273163 |
| Independent checking: VIX-option component | 48.856835 | 46.081785 |
| Independent checking: Forward component | 0.887783 | 0.851944 |
| Independent checking: original-OTM VIX-option component | 35.216476 | 35.073426 |
| Independent checking: original-OTM VIX plus forwards | 36.104259 | 35.925370 |

The checking streams are NEW and independent of training, with 4096 paths per SPX expiry and 8x512 VIX points. Parameters were frozen before checking prices were read. Both candidates use identical checking draws and fixed controls. All market observations still enter calibration: this is numerical checking, not market out-of-sample performance.

The first prototype (run001) was stopped after its first true proposal: common-call joint loss fell from98.119 to79.216, but native OTM VIX-plus-forward loss worsened by9.060. The guard rejected it. Run003 aligns its inner constraints with that safeguard; the rejected vector is not an adopted calibration.

## Parameters

| Parameter | Before | After | Search bounds | Near bound |
|---|---:|---:|---|---|
| cluster_loading | 0.110792303 | 0.085792303 | [0.05, 0.8] |  |
| feedback_curvature | 0.0628502806 | 0.0528502806 | [0.005, 0.4] |  |
| feedback_shift | 3.11527412 | 3.31527412 | [0.05, 5] |  |
| spike_curvature | 0.205340487 | 0.206964524 | [0.002, 0.6] |  |
| spike_shift | 0.299729973 | 0.349729973 | [0.0001, 2] |  |
| cap_level | 5.80241416 | 5.56118274 | [2, 8] |  |
| cap_sharpness | 2.02455349 | 2 | [2, 30] | lower |
| xi_1 | 0.0493015988 | 0.0489543719 | [0.005, 0.2] |  |
| xi_2 | 0.058594625 | 0.0582136627 | [0.005, 0.2] |  |
| xi_3 | 0.0501526536 | 0.0526139733 | [0.005, 0.2] |  |
| variance_floor | 0.00248111523 | 0.0025 | [1e-08, 0.0025] | upper |

## Independent VIX forwards

| Days | Market | Before | After | After SE | Error after |
|---:|---:|---:|---:|---:|---:|
| 27 | 20.100000 | 20.506898 | 20.464411 | 0.016905 | +0.364411 |
| 55 | 20.730000 | 20.225480 | 20.190003 | 0.017261 | -0.539997 |
| 118 | 21.190000 | 21.866977 | 21.839516 | 0.019146 | +0.649516 |
| 153 | 21.720000 | 21.816847 | 21.795378 | 0.010669 | +0.075378 |

## Every retained maturity

| Instrument | Days | Checking common-call MSE before | After | Native OTM MSE before | After | Native quotes in spread after |
|---|---:|---:|---:|---:|---:|---:|
| SPX | 29 | 61.0602 | 85.1705 | 61.0602 | 85.1705 | 31/277 |
| SPX | 57 | 91.8712 | 69.2958 | 91.8712 | 69.2958 | 21/298 |
| SPX | 119 | 44.6707 | 24.8386 | 44.6707 | 24.8386 | 17/228 |
| SPX | 211 | 29.5362 | 36.1854 | 29.5362 | 36.1854 | 12/141 |
| SPX | 302 | 65.2132 | 63.6358 | 65.2132 | 63.6358 | 4/167 |
| SPX | 666 | 16.4475 | 18.4564 | 16.4475 | 18.4564 | 8/112 |
| SPX | 1030 | 18.0384 | 19.3297 | 18.0384 | 19.3297 | 4/71 |
| VIX | 27 | 77.6088 | 65.5020 | 32.0427 | 27.7116 | 5/49 |
| VIX | 55 | 25.2114 | 27.2976 | 46.6515 | 46.0370 | 3/58 |
| VIX | 118 | 55.8020 | 52.5527 | 25.4354 | 27.3921 | 7/55 |
| VIX | 153 | 36.8051 | 38.9748 | 36.7363 | 39.1531 | 6/55 |

## High-strike VIX smile still misses

| Days | IV MAE before, percentage points | After |
|---:|---:|---:|
| 27 | 25.8819 | 26.7260 |
| 55 | 41.1559 | 41.9558 |
| 118 | 35.0739 | 35.6915 |
| 153 | 34.1683 | 34.7392 |

The price objective is not an IV-error objective. Large IV discrepancies can correspond to comparatively small prices in the wings. The lower aggregate price loss does not clear the visible smile misfit, and the high-strike errors above must not be hidden.

## Independent paired VIX loss sensitivity

| Criterion | Before | After | Change | Eight-block jackknife SE of change |
|---|---:|---:|---:|---:|
| canonical_VIX_plus_forwards | 49.744618 | 46.933729 | -2.810889 | 0.142027 |
| native_VIX_plus_forwards | 36.104259 | 35.925370 | -0.178888 | 0.045424 |

These paired changes use the same checking paths for both candidates and frozen controls. The eight-block diagnostic is descriptive; it excludes parameter and deterministic quadrature uncertainty.

## What the optimizer did

The objective keeps the previous family weights, quote scales, interval loss and penalties. It prices one call payoff at all strikes, converting only the existing OTM market inputs. Native model puts use the MODEL forward; no observed-forward substitution is used to hide a VIX mismatch.

At each outer anchor, calculate full joint 4096-path prices and analytic derivatives. Inside a bounded parameter box, use first-order SPX prices and actual nonlinear 4096-point VIX prices/forwards. A VIX-focused SLSQP solve proposes the next vector, prioritizing original-OTM VIX-option loss plus forwards. The warm proposal comes from the preceding average-VIX local solve. Both local VIX and joint criteria are constrained to improve, with margins for SPX approximation error. A normalized 1e-4 tolerance applies only to these extra reserve margins; the actual-price improvement tests remain strict. Reprice the proposal with the FULL SPX/VIX simulator before accepting it. An acceptance requires lower true joint and VIX-plus-forward losses in BOTH common-call and original-OTM conventions. This checked local approximation reduces expensive SPX evaluations; no surrogate score is represented as an actual fit.

1 full outer proposal(s), 1 accepted. Two subsequent local cycles reached their 20-evaluation limits without a reserve-feasible proposal. Maximum 3 proposals and 75 optimization minutes were prespecified. All inner/outer evaluations and rejection reasons are retained. Termination at the declared budget is not convergence, global optimality or parameter identification.

## Numerical checking and limitations

Saved prices/SEs/parity/both objectives and every acceptance decision replay; 91 source/input files verify. Production baseline replay and derivative correspondence screens pass. Five new tests check cached VIX equivalence, local objective/native criteria derivatives and true-price acceptance rules. The unchanged true pricer passed production VIX/full-price directional checks in retained run001; run002 adds a production native/canonical-criteria derivative screen. Run003 restores the verified full baseline Jacobian checkpoint and reuses those screens.

SPX SE is sample SD/sqrt(4096). VIX SE is SD of eight complete independent scrambled-net prices/sqrt(8). They are pointwise sampling estimates, conditional on frozen pilot controls, excluding pilot/parameter uncertainty and deterministic bias. No clipping or post-price/IV smoothing is applied. Invalid IVs and signed prices are retained in the diagnostic tables.

Independent checking paired loss change: -4.228772; descriptive eight-block jackknife SE: 2.242997. Inspect leave-one-block changes and largest stock contributions; this is not a price-convergence certificate.

| Checking candidate | Inside native OTM bid/ask | Negative native prices | Unavailable common-call IVs | SPX order failures | SPX convexity failures | VIX order / convexity |
|---|---:|---:|---:|---:|---:|---|
| baseline | 93/1511 | 16 | 18 | 11 | 61 | 0 / 0 |
| selected | 118/1511 | 8 | 10 | 9 | 56 | 0 / 0 |

Independent selected IV gaps by instrument/expiry: {('SPX', 57): 6, ('SPX', 302): 2, ('VIX', 55): 2}. Eight selected SPX signed estimates are negative at 57/302 days; two VIX55 positive calls remain below observed-forward intrinsic. These fresh-stream failures show that common-call parity and4096 paths do not establish estimator convergence. The previous fixed-parameter ensemble had no negative selected native prices; that fact did not carry over to independent draws. No clipping is applied.

Checking zero-raw-payoff-hit contract counts: {'baseline': 0, 'selected': 0}. Zero empirical SE/hits cannot clear unseen tail contributions. The model fit can remain outside bid/ask even when sampling error decreases. No new P compatibility assessment or Lean certificate is asserted. The selected model remains exploratory until its actual residuals and numerical limitations are accepted.

Plots: [all seven SPX maturities](visualizations/parity-parameter-fit-011/spx-all-maturity-fit.png), [all four VIX maturities](visualizations/parity-parameter-fit-011/vix-all-maturity-fit.png), [VIX forwards](visualizations/parity-parameter-fit-011/forwards-and-path-budgets.png). All plots compare old/new parameters on identical NEW checking draws. Signed errors and SVG exports are included.

Selected configuration: [selected_model.json](runs/m2_parity_parameter_fit_003/selected_model.json); full records and checking data under `runs/m2_parity_parameter_fit_003`. New checking cache/control directories are separate from training; previous runs and plots remain untouched.

Final integrity checks: all 1674 training and 1576 independent-checking frozen input files verify byte-for-byte. The selected model configuration replays the saved parameter vector. All six PNG figures were visually reviewed; PNG/SVG output hashes verify. See `validation/plot_review.json` and `immutable_inputs_verified.json` in the run. Representative selected-state conditional GH32/GH48 probes differ by at most 7.9e-7 VIX points across twelve tested states; this is not a whole-surface quadrature or stock-timestep convergence result.
