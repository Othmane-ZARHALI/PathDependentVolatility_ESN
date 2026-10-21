# Cached curve profiling — 21 October 2026

This checkpoint implements the curve-exposure cache and an exploratory
VIX-only nuisance-curve calibration. It is not a joint SPX/VIX fit, a full
readout/state calibration, a physical-feasibility decision or a model ranking.

## Mathematical and implementation correspondence

With fixed readout, full dated state and Q dynamics, the date-zero normalizer
M0(u) is independent of the deterministic variance curve. For each sampled
full outer state x and curve segment j, cache

\[
 A_j(x)=\int_{[T,T+\Delta]\cap I_j}
 \frac{E^Q[G(\eta_u)\mid X_T=x]}{M_0(u)}du.
\]

For segment variance levels xi_j, the cached conditional integral is

\[
 \mathcal V(x;\xi)=\underline v\Delta+
 \sum_j(\xi_j-\underline v)A_j(x),\qquad
 VIX_T(x;\xi)=100\sqrt{\mathcal V(x;\xi)/\Delta}.
\]

Thus the expensive conditional quadrature is reusable across curve candidates.
The affine relation concerns integrated variance, not VIX or option prices.
The implementation uses floor times the sum of time weights in place of
floor times Delta, preserving the direct numerical integration exactly up to
floating-point summation. It does not rescale the stock, reset normalization
at expiry, linearize the square root or remove the variance floor.

`esn_eq/vix_curve_profile.py` stores read-only segment exposures with the
owning date-zero model configuration. Changed parameters, state, risk prices
or template configuration invalidate reuse. Time grids split at the template
curve knots. Candidate levels change through the cache's evaluate method.
Tests cover direct repricing with two different curves, affine variance with
a nonzero floor, constant-variance segment lengths, inactive columns, invalid
levels and stale parameters. These are implementation/numerical checks, not
formal proofs or error certificates.

## Exploratory optimization design

`vix_curve_fit.py` fixes the Stage-3 candidate-18/scenario-0 architecture,
readout, premium and projected initial state zero. It uses the declared
three-segment nuisance family. Only the first two variance levels are free.
The third segment starts at 0.75 years, beyond every fitted VIX window, so
its exposure column is identically zero. That coefficient is held at the
training diagnostic curve's terminal variance, approximately 0.0642081;
it is not estimated from VIX prices.

All 175 VIX training options and all four quoted forward levels enter the
shared interval objective, with VIX and forward family weights one, SPX weight
zero and the declared curve penalty reported separately. There is no midpoint
tie-breaker or auxiliary stripped-variance residual. Eight optimization starts
use the same cached training exposures. Selection uses training loss only.

Per maturity, four independently scrambled Sobol nets of 64 points build
training exposures; four fresh independent nets check the selected candidate.
The check uses the same training contracts: it is numerical repricing, not
held-out market validation. This preliminary curve-only run is deliberately
below the plan's 4096-path screening budget and cannot qualify final prices
or fit rankings. Conditional GH order is 32, date-zero normalizer order 48,
Gaussian time order 8 and Delta=30/365 under ACT/365. Outer SE is calculated
across four complete scrambles, with only three variance degrees of freedom.

The initial point uses the first two coefficients from the necessary-moment
initializer, with the third reset to the fixed value above. Original audit
flags are preserved. After fitting, necessary moment bounds are recomputed
using sampled model VIX means. Eight ordered check-sample VIX positions per
maturity separately assess conditional GH 32/48 and normalizer 48/64 changes.
Those selected states do not certify full-distribution quadrature accuracy.

The source and inputs are frozen and hashed. Outputs retain all optimizer
starts, train/check contract and forward residuals, payoff-hit counts,
second moments, necessary moment margins, numerical budgets, seeds and raw
segment-exposure arrays for further curve evaluation.

## VIX results: `runs/m2_vix_curve_fit_001`

The selected variance levels are [0.05424961, 0.06003190, 0.06420810], with
the third fixed. All eight starts terminate successfully at essentially the
same two free coefficients. There were 217 residual callbacks across starts;
optimization took 1.46 seconds after the exposures were built. The complete
fit/check/quadrature experiment took 210.83 seconds. This establishes the
computational benefit of caching at the frozen model, not market-fit accuracy.

| Evaluation sample | Starting curve loss | Selected curve loss | Selected VIX option component | Selected forward component | Quotes inside interval |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fitting scrambles | 72.6724 | 62.3668 | 61.7137 | 0.6531 | 22/175 |
| Fresh checking scrambles | 66.8887 | 72.2576 | 71.7047 | 0.5528 | 8/175 |

The curve penalty is approximately 0.00000509 in both selected-curve scores.
**The apparent improvement reverses under independent repricing.** The run
does not establish that the selected curve improves market fit. Numerical
selection sensitivity is material even though the optimizer consistently
solves the fixed-sample problem. Do not replace the initializer as a preferred
fit on this evidence, tune against the checking scrambles or call the check
out-of-sample market performance.

| VIX maturity | Fresh model forward | Quoted forward | Replication SE | Zero-hit training-contract payoffs in check sample |
| --- | ---: | ---: | ---: | ---: |
| 27d | 20.58685 | 20.10 | 0.11059 | 4 |
| 55d | 20.64463 | 20.73 | 0.15615 | 7 |
| 118d | 20.98583 | 21.19 | 0.08541 | 10 |
| 153d | 21.20346 | 21.72 | 0.17659 | 7 |

Outer uncertainty remains too large relative to the declared tolerances.
Zero-hit tails cannot be classified accurately from their zero empirical SE.
At the selected checking states, conditional GH 32-to-48 changes VIX by up to
0.03212/0.03015/0.04526/0.06406 points across the four maturities. Holding
conditional order 48, normalizer 48-to-64 changes it by up to
0.000151/0.000363/0.000788/0.000969 points. These are selected-state changes,
not full-distribution price-error bounds.

Necessary second-moment margins at the sampled model means remain negative
at 118/153 days: approximately -17.97/-48.56 squared VIX points on the check
sample, versus positive margins at 27/55 days. Model-mean uncertainty and
quadrature error prevent treating these as certified feasibility decisions.
They are retained as input/fit diagnostics; the original conflict is not cleared.

## SPX consequence: `runs/m2_curve_spx_check_001`

The VIX-selected curve was repriced on the 29/57-day SPX training contracts,
without optimizing any SPX residual. The checker verifies that its fixed
model/state, family and target snapshot match the fitted run. It uses 4096
evaluation paths, 2048 independent pilot paths, batches of 512, normalizer
order 48 and coupled 128/256 stock steps per trading day. There are 5038/10076
intervals, with maximum lambda*dt approximately 0.1550/0.0775. No negative
controlled option-price estimates occurred.

| Curve | Descriptive SPX interval loss | Quotes inside interval |
| --- | ---: | ---: |
| Earlier training-mid stripped curve, archived 256/day run | 963.9281 | 31/648 |
| VIX-selected nuisance curve, new 256/day run | 384.1083 | 19/648 |

The lower squared-distance loss accompanies **lower interval coverage**;
this is not a uniform improvement. Both are point diagnostics at unresolved
numerical accuracy. The comparison concerns the two short SPX maturities only.
It is not a qualified joint-model score or held-out market evaluation.

Both short maturities lie before the first nuisance knot at 0.25 years, so
only xi0 affects their prices. The initializer already used xi0=0.0545936,
close to the selected 0.0542496. Much of the difference from the stripped
curve is therefore attributable to the higher near variance in this restricted
family, not to the VIX optimizer. Selected-curve total variance is 41.70%/25.53%
above the stripped diagnostic at 29/57 days; these are variance differences,
not percentage errors in option prices.

| Expiry | Largest controlled 256-minus-128 price change | Paired SE at that strike | Largest controlled fine-price SE |
| --- | ---: | ---: | ---: |
| 29d | 0.68126 | 0.36231 | 2.35805 |
| 57d | 0.93661 | 0.54853 | 2.80206 |

Price quantities are undiscounted SPX points. Maxima select across many quotes;
pointwise paired SEs are not simultaneous error bounds. Stock-forward ratios
at 256/day are 0.999828 (SE 0.000880) and 0.998809 (SE 0.001156).
Quote-relative sampling and time-grid gates remain unresolved. Longer SPX
maturities, including identification of the third curve coefficient, remain open.

## Interpretation and next work

The cache and optimization machinery work in the tested scope. A reliable
market improvement has not been established. The next numerical priority is
to stabilize candidate comparisons with larger independent scramble budgets,
better tail sampling/controls and conditional-quadrature refinement. Preserve
the current train/check reversal rather than selecting anew on checking data.
Further readout/state fitting and the joint SPX/VIX objective remain open.

81 package tests and eight legacy/driver checks pass. No held-out quote or
final-test maturity performance has been evaluated.

## Reproduction

From repository root, use the project Python environment:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -m matteo.model_comparison.option_pricing.vix_curve_fit \
  --output /tmp/m2_vix_curve_new
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -m matteo.model_comparison.option_pricing.curve_spx_check \
  --fit /tmp/m2_vix_curve_new/run.json --output /tmp/m2_curve_spx_new
```

The environment used here is `/private/tmp/pdv-options-venv/bin/python`.
Commands refuse existing output directories. The archived curve initializer
is needed for the optimizer's first starting point; its hash is recorded.
