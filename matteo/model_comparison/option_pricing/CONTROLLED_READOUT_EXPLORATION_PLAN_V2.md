# Broader fixed-reservoir exploration — 24 October 2026

The user explicitly asks to continue exploring beyond the unsatisfactory old
fit. A warm start is an incumbent, not the scientific objective or a stopping
criterion. Keep all eleven agreed free parameters; do not freeze shape simply
because the old vector can be reproduced.

Fixed: original OU grids/counts/rates, correlations, projection powers, Q prices,
C0=1/Y0=F0=0, optional channels disabled. Same 1294 OTM SPX cutoff-3 quotes at
29/57/119/211/302/666/1030 days, 217 OTM VIX quotes at 27/55/118/153 days and four
forwards. All market observations remain in sample.

Restore original per-expiry 512-path importance samples for BOTH searching and
separate checking, all 8x64 VIX states, GH32 conditional VIX, GH48 normalizer,
time order 8 and 512 stock steps per trading day. Reuse original equity shocks
and weights with exact replay. Existing complete checking streams are reused;
this is not a new independent validation of the whole research history.
Recompute all readout-dependent quantities. Reuse original independently fitted
zero-mean controls, continuously interpolated as RANDOM CONTROL DRAWS at three
market-only anchors per option side. SPX uses scaled moneyness .5/1.5/2.5;
VIX uses first/middle/last sorted strikes. Forward controls remain unchanged.
This removes nearest-donor discontinuities without smoothing prices or IVs.

The previous piecewise-linear test has residual convexity failures. A natural
cubic/linear-extrapolation comparison was also screened: it restored several
order failures and worsened convexity (5 order/93 convexity failures checking),
so it is NOT adopted. Continue with the more stable continuous piecewise-linear
controls as an explicitly exploratory estimator, report all signed/shape failures,
and do not make pricing-convergence or shape-qualification claims. This run is
not blocked on a universal finite-sample shape guarantee. The user has explicitly
requested further model exploration after those limitations were explained.

Free: five score coefficients, cap level/sharpness, three variance levels, floor.
Expand only previous boundary-hit domains once: beta_J lower 0.01 -> 0.0001
(still strictly positive), cap sharpness upper 15 -> 30, floor lower 1e-6 -> 1e-8.
Keep all other previous bounds. Record boundary hits; do not infer identification.

Three prespecified starts: prior all-maturity reference with only beta_J moved
to the new positive lower bound; prior eleven-parameter raw-run selection;
and a materially asymmetric start with beta_J=.5, q_J=.35, q_Z=.12, beta=2,
L=6.5, sharpness=6 and the old curve. Each gets at most 6 curve-prefit evaluations
(three nuisance levels only), then at most 12 ALL-ELEVEN-parameter evaluations.
Use the same interval loss, expiry weights, forward scale and penalties as before.
Keep the best searching incumbent across every evaluation. Analytic Jacobians
and identical deterministic quadrature protocol at all candidates.

Total optimization wall budget 45 minutes, shared with per-start 15-minute caps
(checked between evaluations). Preparation/checking are separate. Freeze final
selection before separate checking. No clipping or post-price smoothing. Candidate
controls are fixed functions of the original primitives and remain zero-mean
under unchanged target law. Adaptive candidate normalizer interpolation must pass
midpoint and independent same-GH48 probes at 2e-7; this is not GH convergence.
Native parallelization must match the serial production formula and derivatives.

Retain sources, config, seeds/provenance, every objective evaluation, failures,
all quote prices/SEs, terminal/payoff arrays, shapes and all eleven maturity plots.
Report P compatibility as unresolved; no new formal certificate is asserted.


## Version 2: deterministic GH32 convolution acceleration

The direct tensor implementation reproduced the original baseline but cost
about 150 seconds per readout. The first run was explicitly interrupted before
completed optimization/selection; its source, baseline/warm results and log remain.
Do not count it as an independently completed calibration.

The cluster feature is linear in the score and Gaussian-independent of Y/F.
For each conditional lag, precompute its GH32 convolution as a function of
z=c_C*mean(C)+q_Z*(Y-beta)^2+q_J*(F-beta_J)^2, with cubic-Hermite interpolation
of values and derivatives. Keep the same GH32 nodes and weights in ALL three
dimensions. Adaptive midpoint/quarter probes screen value interpolation at
1e-9 relative (small absolute floor) and derivative interpolation at 1e-7 scaled.
If independence structure or interpolation screens fail, fall back to the full
GH32 tensor. Saturated-score treatment starts beyond the lowest C-node score
L+40/k, at floating-point negligible soft-cap differences. Record comparison
against direct GH32 over representative lags, all four maturities and each of
the three starts. This is a numerical interpolation screen, not a rigorous
uniform error bound or new model approximation. Use standard -O3 compilation;
the separately benchmarked fast-math build is NOT adopted.

Re-evaluate selected checking prices using the direct full GH32 tensor after
selection freezes. Retain both accelerated and direct outputs and differences.
No cache or old price record is overwritten. All other scope, parameter bounds,
starts and optimization budgets above remain unchanged.
