# Results: nested M2 time-step convergence audit

## Executive conclusion

The exact OU transition prevents numerical instability, but two integration steps per trading
day are not sufficient for the nonlinear variance readout and stock-return integral of the
selected M2 architecture. The coarse grid preserves average variance almost exactly while
materially changing the realised variance path and several stylised-fact statistics.

This is a numerical convergence failure, not a volatility explosion. Cap activity remains near
zero and maximum volatility is stable. The problem is unresolved intraday variation and temporal
alignment between the fast OU factors, variance, and the equity shock.

## Design

- Frozen selected Stage-3 M2 architecture and readout; no recalibration by resolution.
- Five new seeds: `910003`, `930011`, `950009`, `970019`, `990001`.
- 12 paths per seed, six years, one-year burn-in, 252 daily outputs per year.
- Internal resolutions: 2, 4, 8, 16, 32, and 64 steps per trading day.
- One exact 64-step OU/equity path per seed. Coarser OU states are exact endpoint subsamples and
  coarser equity Brownian increments are sums of the same fine increments.
- The 64-step result is the numerical reference, not an assertion of an exact continuous-time
  limit. At that resolution the fastest feedback mode still has
  \(\lambda_{\max}\Delta t=0.31\).

## Main results

Across-seed means are:

| Steps/day | Max \(\lambda\Delta t\) | Return corr. vs 64 | Variance corr. vs 64 | Return RMSE / daily SD | Variance relative RMSE |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 9.921 | 0.883 | 0.861 | 0.484 | 0.672 |
| 4 | 4.960 | 0.908 | 0.946 | 0.430 | 0.396 |
| 8 | 2.480 | 0.938 | 0.983 | 0.353 | 0.215 |
| 16 | 1.240 | 0.965 | 0.996 | 0.265 | 0.107 |
| 32 | 0.620 | 0.986 | 0.999 | 0.167 | 0.047 |
| 64 | 0.310 | 1.000 | 1.000 | 0.000 | 0.000 |

The mean variance bias relative to the 64-step result is below 0.06% at every resolution. That
agreement is misleading if considered alone: at two steps/day the pointwise daily-variance RMSE
is 67% of mean variance and the pooled path correlation is only 0.86.

| Steps/day | Rough H | Rank leverage | Rank Zumbach | Excess kurtosis | Log-var. ACF |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 0.0923 | -0.0200 | 0.0524 | 6.43 | 0.2037 |
| 4 | 0.1144 | -0.0146 | 0.0485 | 6.73 | 0.2347 |
| 8 | 0.1273 | -0.0120 | 0.0425 | 7.13 | 0.2523 |
| 16 | 0.1340 | -0.0088 | 0.0361 | 7.41 | 0.2605 |
| 32 | 0.1364 | -0.0073 | 0.0317 | 7.58 | 0.2635 |
| 64 | 0.1370 | -0.0058 | 0.0270 | 7.63 | 0.2644 |

Paired seed differences between two and 64 steps/day are statistically separated from zero:

| Statistic | Mean difference, 2 minus 64 | Paired 95% CI |
| --- | ---: | ---: |
| Rough H | -0.0447 | [-0.0489, -0.0406] |
| Rank leverage | -0.0142 | [-0.0174, -0.0110] |
| Rank Zumbach | +0.0254 | [+0.0078, +0.0430] |
| Excess kurtosis | -1.20 | [-2.16, -0.24] |
| Log-variance ACF | -0.0607 | [-0.0696, -0.0517] |

Thus the coarse grid makes volatility look less rough and less persistent, strengthens measured
negative leverage and positive rank Zumbach, and suppresses return kurtosis.

## Acceptance-rate warning

The declared joint proxy acceptance rate falls from 100% at two steps/day to 20% at 64
steps/day. This should not be read as evidence that the coarse model is better. It shows that
the acceptance decision is resolution-dependent. In particular, the fine-grid mean rank
leverage is -0.0058, weaker than the declared -0.008 boundary, whereas the two-step estimate is
-0.0200. Coarse discretisation can therefore manufacture an apparent pass.

Positive Zumbach asymmetry survives refinement, but its rank magnitude falls from 0.0524 to
0.0270. The mechanism remains present; its previously reported magnitude was overstated by the
coarse grid.

## Explosion and cap diagnostics

No numerical explosion is visible:

- the mean cap-exceedance fraction stays around \(3\times10^{-5}\);
- maximum annualised volatility is approximately 0.77 at every resolution;
- the 99.9% volatility quantile stays near 0.56;
- mean-variance bias is negligible.

The dominant issue is therefore missed and mistimed intraday variation, not unstable OU paths
or cap failure.

## Resolution decision

Two, four, and eight steps/day should not be used for final M2 evidence with the current rate
grid. Sixteen steps/day is adequate only for exploratory work. Thirty-two steps/day gives close
variance paths and nearly converged roughness, clustering, and tails, but leverage and Zumbach
still move between 32 and 64.

Use 64 steps/day as the new working baseline for the current architecture. Before treating the
asymmetry magnitudes as final, run one targeted 64-versus-128 confirmation. If 128 steps/day is
too expensive, reduce the fastest feedback rates and recalibrate the architecture on a rate grid
whose half-lives are resolvable at 32--64 steps/day.

All future comparisons must keep continuous-time rates fixed when changing the integration
grid. Recalibrating separately at each resolution would conceal rather than measure numerical
error.

## Files

- `per_seed.csv`: all seed-resolution rows and pathwise-reference errors.
- `summary.csv`: across-seed means and standard deviations.
- `manifest.json`: frozen model, rates, seeds, and exact nesting convention.
- `run.py`: deterministic reproduction script.
