# Active two-maturity VIX pilot — 21 October 2026

## Explicit user-authorized scope assumption

The active exploratory VIX objective includes **27 and 55 days only**:
87 training options and their two associated forward levels. The 118/153-day
options and forwards are deferred. SPX scope is unchanged. This supersedes the
earlier decision to keep all four VIX maturities active. The authoritative
decision and reintroduction gate are in `active_scope.json` and
[VIX_SCOPE_DECISION.md](VIX_SCOPE_DECISION.md).

This is a research simplification, not evidence that the deferred quotes are
wrong or that M2 cannot fit them. Deferred training maturities are not pristine
new held-out tests. Raw/processed data, historical results and global conflict
flags remain preserved. A two-maturity result does not establish full-maturity
joint feasibility.

## What is identifiable in this pilot

The two continuous VIX windows are 27-to-57 days and 55-to-85 days. Both end
before the first declared nuisance-curve knot at 0.25 years. Consequently
their prices depend only on the first curve coefficient. No independent
post-0.25-year coefficient is identifiable from these prices.

`vix_near_fit.py` therefore fits one near variance level while keeping all
reservoir/readout/state/risk-price parameters fixed. Inactive curve segments
are given the same level solely as computational placeholders; their levels
are not inferred or forecast. This makes the declared adjacent-level penalty
zero and avoids feeding a longer-maturity estimate into the near fit through
that penalty. In-scope price comparisons use the same placeholder convention
for every candidate. The old historical candidates are fixed references only.

The full objective at the preceding four-maturity fit was essentially a sum
of a near block depending on xi0 and a far block depending on xi1, with only
a tiny smoothness coupling. It is therefore expected that removing the far
block need not materially change the near coefficient.

## Experiment and results

`runs/m2_vix_near_fit_001` reuses only the original fitting exposures at the
two retained maturities (four scrambles x 64 states each). Eight one-dimensional
optimization starts use the common bid/ask interval objective, VIX/forward
family weights one and SPX weight zero. No deferred quotes or forwards enter
the objective. The old checking draws are not used for this selection.

After freezing the fitted level, eight newly seeded Sobol scrambles of 128
states per maturity are generated. Their first 64 points provide a nested
sampling-budget comparison. These are independent numerical checks on the
same training contracts, not held-out market performance. Each full fresh
estimate uses 1024 states per maturity. Conditional GH order 32, normalizer
order 48 and Gaussian window-time order 8 are fixed for those price estimates.

The selected near variance is **0.0542496100** (about 23.29% annualized
volatility). It differs from the earlier four-maturity fit's near coefficient
by only about 1.2e-9. The full experiment took 209.59 seconds.

| Numerical sample | Previous starting near level, 0.0545936 | New near fit, 0.0542496 | Selected minus starting |
| --- | ---: | ---: | ---: |
| Fitting: 4 x 64 states | 107.796558 | 107.649592 | -0.146966 |
| Fresh: 8 x 64 states | 122.671479 | 122.419119 | -0.252360 |
| Fresh: 8 x 128 states | 109.855757 | 109.820535 | -0.035222 |

The large earlier full-scope ranking reversal is absent in these point
comparisons, but this does **not** establish a reliable improvement. A paired
leave-one-whole-scramble-out calculation gives descriptive jackknife SEs
0.8512 and 0.2760 for the fresh 64/128-point loss differences, respectively.
At 128 points, leave-one-out differences range from -0.1813 to +0.1958.
The tiny apparent gain is unresolved. Jackknife uncertainty with eight
scrambles and unseen tails is not a tail-error certificate.

At the full fresh budget, 14/87 options lie inside bid/ask, and three 55-day
payoffs have no positive observed payoff. Median option SE divided by fixed
market scale is 1.039, versus the target of at most 0.125. Neither adequate
market fit nor numerical convergence has been established.

| Maturity | Model forward | Market forward | SE across fresh scrambles |
| --- | ---: | ---: | ---: |
| 27d | 20.76664 | 20.10 | 0.09199 |
| 55d | 20.39612 | 20.73 | 0.05201 |

Only one curve level is free, so the two forward levels cannot be adjusted
independently. These residuals belong to this restricted, fixed-readout/state
experiment; they do not establish a limitation of unrestricted M2.

## Separate numerical-error checks

At eight ordered fresh-sample VIX positions per maturity, each integration
setting is varied separately while the others are held fixed:

| Maturity | Max conditional GH 32-to-48 VIX change | Max normalizer 48-to-64 VIX change | Max time-rule 8-to-16 VIX change |
| --- | ---: | ---: | ---: |
| 27d | 0.024000 | 0.00007982 | 0.000001419 |
| 55d | 0.016689 | 0.00050251 | 0.000000169 |

These are selected-state level changes, not bounds on full-distribution option
price bias. The conditional rule, normalizer and time rule remain separate
approximations. Original global preprocessing flags include deferred blocks;
they are not automatically either cleared or evidence against the near scope.

The next studies should keep separate:

1. Sampling and tail error: vary path/state budgets and independent seeds at
   fixed integration settings; test likelihood-weighted rare-event sampling
   in pilots as well as evaluation.
2. SPX stock time bias: use coupled 64/128/256 grids and finer levels as needed;
   more paths alone cannot remove this bias.
3. VIX window-time integration: refine Gaussian time rules independently of
   the conditional Gaussian expectation rule.
4. Normalizer error: refine date-zero M0 integration independently. A biased
   deterministic normalizer affects every simulation and is not removed by
   increasing the outer count.

Importance sampling and multilevel Monte Carlo remain next implementations;
this run uses the existing Sobol/Gaussian machinery. SPX's adapted hedging
controls are already implemented and tested. No final numerical gate has been
relaxed. There are 84 package tests and 11 legacy/driver checks passing.

## Reproduction and handoff

Use `vix_near_fit` for the active pilot. The earlier `vix_curve_fit` driver and
its frozen sources remain historical four-maturity reproductions.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -m matteo.model_comparison.option_pricing.vix_near_fit \
  --output /tmp/m2_near_new
```

The environment used here is `/private/tmp/pdv-options-venv/bin/python`.
Sources and the active-scope JSON are captured before calculation and checked
unchanged at completion. Cached training input hashes, fresh seeds, candidate
ledger, prices/forwards, SEs, payoff hits, exposures and refinement checks are
saved. Reintroduce the deferred maturities only under the recorded return gate
before any full-maturity claim.
