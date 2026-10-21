# Coupled 64/128-step pricing checks — 21 October 2026

## Fixed design

The user requested the finer resolution required by the existing M2 time-step
audit. The model remains the frozen Stage-3 anchor, candidate 18/scenario 0;
the projected date-zero state remains (C,Y,F)=0. There is no recalibration or
substitution of the later Stage-4 candidate. The diagnostic training-mid curve
and its existing cross-market conflict flag remain unchanged.

`fine_grid.py` adds bounded-memory batching and a deterministic normalizer
cache keyed by times, initial state, readout and GH order. Two new tests check
cached versus direct normalization and exact coupled constant-volatility paths
across independently seeded batches. Captured sources and input hashes accompany
both run directories. Original experiment records are preserved.

## SPX: `runs/m2_spx_64_128_001`

- 4096 evaluation paths, 2048 independent pilot paths, batches of 512.
- Normalizer order 48 fixed across both resolutions.
- Same exact fine-grid OU states and equity shocks for 64/128 stock calculations;
  coarse increments are sums of the fine increments.
- 2520 coarse and 5040 fine intervals over the two exact calendar expiries.
- Maximum lambda*dt is 0.309875 at 64 and 0.154938 at 128 steps per trading day.
- Common bid/ask training targets: 328 at 29d and 320 at 57d.

| Expiry | Largest absolute controlled fine-minus-coarse price change | Paired SE at that strike | Strike | Largest controlled fine-price SE |
| --- | ---: | ---: | ---: | ---: |
| 29d | 0.71955 | 0.42794 | 5930 | 2.14626 |
| 57d | 1.35498 | 0.57438 | 6910 | 2.69454 |

All price quantities are undiscounted SPX index points. The largest plain
paired changes are 0.72432 and 1.83675. The maximum absolute controlled changes
in frozen price-scale units are 3.83 and 3.07. The strike maximizing that ratio
need not maximize the absolute change. These maxima are selected across many
quotes; their pointwise SEs are not simultaneous confidence bounds.

The fine-grid normalized stock means are 0.999077 (SE 0.000818) and 0.999694
(SE 0.001079), compatible with the forward-one check at this sampling budget.
None of these checks establishes convergence within the predeclared price
tolerances. Further coupled refinement and improved sampling remain necessary.

## VIX: `runs/m2_vix_64_128_states_001`

This is a fixed-state conditional integration check at the 55-day expiry,
not a new outer-distribution market fit. Full OU states are obtained from the
minimum-prior-norm map of the three declared (C,Y,F) projection vectors.
The global normalizer uses GH order 48 and the conditional expectation uses
GH order 32. There are 1327 and 2652 midpoint integration nodes at 64/128.

| Projected state (C,Y,F) | Quadrature VIX, 64 | Quadrature VIX, 128 | Fine minus coarse |
| --- | ---: | ---: | ---: |
| (0,0,0) | 16.862582 | 16.862579 | -0.00000289 |
| (0,-1,-1) | 25.383440 | 25.383432 | -0.00000793 |
| (1,-2,-2) | 51.554177 | 51.554107 | -0.00007001 |

These deterministic differences isolate time-integration changes at fixed
quadrature orders and states. They do not measure quadrature-order error or
establish convergence for every possible outer state.

Nested full-state OU continuations were also run at both grids with 512 inner
paths per state. The VIX estimates were [18.2315,27.9310,52.4962] at 64 and
[16.1831,23.6902,52.1527] at 128. The integration grids have different node
counts and hence different random streams: those changes cannot be interpreted
as paired time-discretization errors. The inner integrated-variance SEs are in
`run.json`; at this inner budget the nested estimator remains visibly noisy.
The small quadrature grid differences do not resolve that sampling uncertainty.

## Pricing interpretation

The conditional engines evaluate

\[
\mathrm{VIX}_T^2(x)=\frac{10000}{\Delta}\int_T^{T+\Delta}
\left[\underline v+(\xi_0(u)-\underline v)
\frac{E^Q[G(\eta_u)\mid X_T=x]}{M_Q^0(u)}\right]du.
\]

The date-zero normalizer M_Q^0 remains fixed across states. Quadrature uses the
exact three-dimensional Gaussian projection law conditional on the **full**
OU state x. Nested MC simulates full conditional OU continuations, averages
their integrated variances and only then takes the square root. It does not
simulate inner SPX returns or re-normalize separately at expiry. Outer Monte
Carlo is still needed for VIX futures/options when evaluating the distribution
of X_T; these three diagnostic states do not replace that outer distribution.

## Next and defaults

The standalone M2 SPX runner now defaults to 64 stock steps/day, and the VIX
runner to 64 window-integration steps/day. Explicit 2/8-step controls remain
available for reproducing earlier coarse diagnostics. The batched driver is
preferred for high-resolution SPX work. The older commands in dated results
files reflect their original defaults; source snapshots and recorded arguments
define those historical runs.

Reproduce with fresh output directories:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -m matteo.model_comparison.option_pricing.fine_grid --mode spx \
  --output /tmp/m2_spx_64_128_new
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -m matteo.model_comparison.option_pricing.fine_grid --mode vix \
  --output /tmp/m2_vix_64_128_new
```

Next numerical work is a coupled 128/256 SPX check with frozen parameters and
independent confirmation seeds, plus sufficient variance reduction to assess
quote-relative errors. No final fit or model ranking is claimed here.
