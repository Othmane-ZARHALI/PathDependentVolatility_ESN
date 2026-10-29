# Point1 acceleration and minimum volatility policy

The current calibration policy requires minimum instantaneous volatility of5%, with its floor calibratable between5% and6%. The6% upper endpoint is the implementation choice used to keep the floor free within the stated reasonable range. All eleven parameters remain free. Historical bounds/configurations are preserved.

With $v_t=f+(\xi(t)-f)G_t/M(t)$, positive $G_t/M(t)$ and $\xi(t)\geq f$, $f\geq0.05^2=0.0025$ ensures $\sqrt{v_t}\geq5\%$. This is an instantaneous model-volatility bound, not an option-IV constraint.

## What changed

At each fully priced anchor, store call prices, VIX forwards and all eleven derivative columns. Inner optimizer trials now approximate BOTH SPX and VIX prices/forwards by an affine Taylor model. The original interval loss is evaluated on those approximate prices. Native model puts use the approximate MODEL forward and its derivative. No simulator is held by the local-model object. Full nonlinear4096-path repricing checks every proposal before acceptance; accepted anchors refresh the Jacobian.

Pilot `m2_joint_local_fit_001` starts from the previous selected vector with its variance floor rounded upward by about9e-15 to the exact new lower bound. Same1511 original OTM observations, four forwards, SPX cutoff3/seven retained maturities/four VIX maturities, OU dynamics/rates/grids/Q prices/initial state, common-call payoff, all4096 paths, GH32/GH48/time8/stock512. No coarse mesh or parameter-block optimization was introduced.

## Measured work

| Operation | Seconds |
|---|---:|
| Old nonlinear-VIX local trial, historical median | 33.345576 |
| New approximate joint-price evaluation, measured median | 0.002540 |
| All new inner solves (119 optimizer evaluations) | 0.323587 |
| Full baseline, including eleven sensitivities | 362.988 |
| All full checked proposals | 756.597 |
| Whole pilot | 1143.436 |

Full nonlinear pricer calls: 3; local optimizer simulator calls:0. The old/new local rows perform different work (nonlinear calculation versus approximation). This establishes local computational acceleration, not equal accuracy at every trial or an end-to-end speed factor. Whole-pilot timing also differs in validation scope from the earlier run: no new independent checking is priced here.

## Actual training-ensemble fit

| Criterion | Before | After |
|---|---:|---:|
| Common-call joint | 95.148310 | 95.117702 |
| Original-OTM joint | 85.038631 | 84.982024 |
| SPX | 49.491692 | 49.518173 |
| Common-call VIX options | 44.791849 | 44.728178 |
| Original-OTM VIX options | 34.682170 | 34.592500 |
| VIX forwards | 0.831425 | 0.838006 |

2 fully checked proposal(s), 2 accepted. Acceptance requires true improvements in both joint and both VIX-plus-forward criteria. These are actual full-pricer scores. Local approximations are used only to propose steps; rejected points and prediction errors remain in the run records. A two-proposal pilot does not establish optimization convergence, global optimality, identification or a satisfactory smile fit.

| Proposal | Predicted common-call joint | Actual joint | Predicted native VIX + forwards | Actual | Accepted |
|---:|---:|---:|---:|---:|---|
| 0 | 95.102639 | 95.128448 | 35.464852 | 35.463530 | True |
| 1 | 95.082787 | 95.117702 | 35.433485 | 35.430506 | True |

## Parameters

| Parameter | Before | After | New bounds |
|---|---:|---:|---|
| cluster_loading | 0.08579230299 | 0.07579230299 | [0.05, 0.8] |
| feedback_curvature | 0.05285028059 | 0.04995142993 | [0.005, 0.4] |
| feedback_shift | 3.31527412 | 3.39527412 | [0.05, 5] |
| spike_curvature | 0.2069645236 | 0.2078592781 | [0.002, 0.6] |
| spike_shift | 0.3497299725 | 0.3697299725 | [0.0001, 2] |
| cap_level | 5.561182737 | 5.502780057 | [2, 8] |
| cap_sharpness | 2 | 2 | [2, 30] |
| xi_1 | 0.04895437191 | 0.04898961623 | [0.005, 0.2] |
| xi_2 | 0.05821366266 | 0.05827299307 | [0.005, 0.2] |
| xi_3 | 0.05261397332 | 0.05257415274 | [0.005, 0.2] |
| variance_floor | 0.0025 | 0.00254 | [0.0025, 0.0036] |

Selected floor volatility: 5.039841%; beta_J=0.369729973>0.

## Every retained maturity, actual original-OTM prices

| Instrument | Days | MSE before | After | Quotes in bid/ask after |
|---|---:|---:|---:|---:|
| SPX | 29 | 106.02360 | 110.54088 | 27/277 |
| SPX | 57 | 39.36155 | 39.94458 | 18/298 |
| SPX | 119 | 62.31250 | 60.63675 | 5/228 |
| SPX | 211 | 33.48936 | 32.48490 | 27/141 |
| SPX | 302 | 73.61559 | 71.62859 | 20/167 |
| SPX | 666 | 19.99475 | 19.52761 | 7/112 |
| SPX | 1030 | 11.64448 | 11.86391 | 7/71 |
| VIX | 27 | 27.42891 | 28.01390 | 5/49 |
| VIX | 55 | 46.07289 | 45.87496 | 0/58 |
| VIX | 118 | 28.24451 | 27.90219 | 6/55 |
| VIX | 153 | 36.98237 | 36.57895 | 6/55 |

## Verification and limits

Ten focused new/existing tests pass. Saved actual means/SE/parity/both losses, 119 local trial scores, all true acceptance decisions and 68 frozen sources/inputs replay. Bounds and selected configuration verify. Core datasets are the existing read-only caches, whose byte identities were verified in parameter-fit003; no new dataset was generated.

Selected actual prices: 0 negative native estimates; 2 unavailable common-call IVs; 128/1511 native prices inside bid/ask. Signed estimates, invalid IVs and shape/tail diagnostics are retained without clipping or smoothing. SPX pointwise SE uses iid4096 paths; VIX uses eight independent512-point nets. Bands exclude parameter/pilot uncertainty and deterministic error.

No new independent checking or price-convergence claim is made for this pilot. The prior VIX wing mismatch and SPX estimator stability concerns remain relevant. A local linear model can extrapolate poorly; the saved full checks measure this directly. No P-compatibility or Lean certificate is asserted.

[All seven SPX maturities](visualizations/joint-local-012/spx-all-maturity-fit.png) · [All four VIX maturities](visualizations/joint-local-012/vix-all-maturity-fit.png) · [Forwards and path allocation](visualizations/joint-local-012/forwards-and-path-budgets.png). Blue is the selected vector and orange the previous fitted vector, on the same TRAINING draws. PNG/SVG and signed-error/precision figures are included.

[Selected model](runs/m2_joint_local_fit_001/selected_model.json). Current bounds live in `calibration_policy.py`; use these for new fits. Historical scripts keep their old bounds for reproducibility.

All six PNG figures have been visually reviewed and PNG/SVG hashes verify; see `validation/plot_review.json`. Exact command/environment are recorded in `command.json` in the run. A future continuation can restore the selected full-price/Jacobian checkpoint with `parity_checkpoint.load_result` rather than repeat an identical baseline simulation.
