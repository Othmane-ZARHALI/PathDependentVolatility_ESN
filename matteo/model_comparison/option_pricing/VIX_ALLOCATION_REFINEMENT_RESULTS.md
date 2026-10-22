# Targeted VIX allocation and warning attribution — 22 October 2026

## Retained scope and frozen objects

All four VIX maturities (27/55/118/153 days), 175 training options and four
forwards remain active. The Stage-3 model, risk premium, initial state, selected
variance curve and 17-component targeted proposal from
[VIX_TARGETED_PILOT_RESULTS.md](VIX_TARGETED_PILOT_RESULTS.md) are frozen. No
new model fit, proposal learning, maturity removal or price rescaling occurs.
SPX is unchanged. The lean-numerics-bridge workflow separates the model
identity, implementation correspondence, sampling evidence and numerical bias.
No Lean certificate or rigorous price-error bound is claimed.

Two bounded questions are studied independently: whether quadrature explains
the saved 118-day moment warning, and whether fixed component allocations help
fresh randomized-Sobol pricing at the same targeted density.

## Exact saved-state 118-day attribution

`m2_vix_targeted_refinement_001` loads **all 8 x 64 original targeted 118-day
states**, their exact weights and GH32/M48/time8 VIX values. Compute GH48 node
expectations once, then integrate them with normalizer orders 48 and 64 on
the same time nodes. No new outer sample enters. The normalizer remains
conditioned at date zero; it is not reset at expiry.

| Conditional GH / normalizer | E[VIX squared] ± scramble SE | Residual from 600.31899 | Paired change from original ± SE |
| --- | ---: | ---: | ---: |
| 32 / 48 | 588.51408 ± 4.30771 | -11.80491 | 0 |
| 48 / 48 | 588.51875 ± 4.30421 | -11.80024 | +0.00467 ± 0.02205 |
| 48 / 64 | 588.50805 ± 4.30413 | -11.81094 | -0.00603 ± 0.02206 |

The changes do not explain the warning on this exact sample. The corresponding
forward change is -0.000205 ± 0.000233 VIX points. Across all 44 118-day options,
the largest absolute paired price change divided by quote scale is 0.00464.
These results isolate these rule changes on the observed states; they do not
bound bias over unseen tail states or the population finite-GH rule.

Limited additional checks, relative to GH48/M64/time8:

- Time8 to time16 on the first eight states per original scramble: weighted
  moment change +6.88e-7 ± 8.34e-7; maximum VIX-level change 2.99e-6.
- Conditional GH48 to GH64 on the first four states per original scramble:
  weighted moment change +0.01722 ± 0.01010; maximum VIX-level change 0.01386.

Those subsets contain 64 and 32 states, respectively. They cannot replace the
full-sample results or certify global time/GH convergence. Other maturities
remain active but are not newly refined in this attribution study. Together
with previous direct Gaussian moment diagnostics, the evidence continues to
favor outer-tail/small-scramble sensitivity; it does not establish every bias
as negligible or imply that a normalizer patch is needed.

The verifier matches original states/weights, reconstructs node integration,
and replays 138 summary rows and both subset checks from saved artifacts.
It does not recompute Gaussian node expectations. All 38 captured source/input
hashes match. The run took 416.8 seconds while the allocation run also used CPU.

## Fresh same-density allocation comparison

`m2_vix_targeted_allocation_001` compares two designs with identical density

    g(z) = sum_k pi_k phi(z-m_k),  f(z) = phi(z).

Random selection estimates mean[(f/g) H]. Fixed selection takes n_k draws from
each component and estimates

    sum_k pi_k mean_k[(f/g) H]
      = (1/N) sum_samples (N pi_k/n_k) (f/g) H.

The extra factor N*pi_k/n_k is essential. It corrects deterministic allocation,
not the model or a measured normalization error. The density g and its exact
likelihood ratios remain unchanged. No sample-average normalization is used.

At the coarse count N=64, allocate 32 origin draws and two to each of the 16
shifted components. At N=128, double every count. The origin has pi=0.5, the
six inherited shifts have pi=1/24, and ten new shifts have pi=0.025. Allocation
factors are therefore 1, 4/3 and 0.8 respectively. Likelihood weights f/g remain
bounded by 2; **effective allocation weights are bounded by 8/3**, not 2.

Each component gets its own independent scrambled 36-dimensional Gaussian
net, with seed derived from (replication seed, component index). Every count
is a power of two; at the fine count there are still only four draws in each
shifted component. A coarse estimate uses the first n_k(coarse) points inside
every component, not the first 64 rows of the component-blocked fine array.
The complete set of component scrambles is one uncertainty unit.

Use 12 complete independent replications per method/maturity, with nested
64/128 points per replication: **1536 expensive conditional values per method
per maturity, 12288 total**, plus 12 selected-state GH48 checks for the fixed
method. This increases both replication count and points per net from the
preceding R8 x 64 study. The two current methods have equal budgets and
independent streams; count changes are paired within replications. Reported
± values are one scramble SE and exclude numerical bias.

The baseline here is random allocation under the **same targeted density**.
The older seven-component proposal is a historical reference; this study does
not freshly compare its relative performance. It isolates the complete fixed
versus random allocation designs, including their different Sobol coordinates,
not a pure variance effect of component counts alone.

## Allocation results

The same-density comparison is mixed:

| Days | Median option SE: fixed / random | Contracts with lower fixed SE | Random SE screen | Fixed SE screen |
| --- | ---: | ---: | ---: | ---: |
| 27 | 0.797 | 30/40 | 4/40 | 6/40 |
| 55 | 1.124 | 15/47 | 9/47 | 9/47 |
| 118 | 0.931 | 28/44 | 12/44 | 11/44 |
| 153 | 0.883 | 31/44 | 11/44 | 9/44 |

These are medians of individual-contract SE ratios. Totals are **36/175 random
versus 35/175 fixed** observed SE screens, requiring SE <= quote scale / 8
and at least two positive payoffs. This is not a total error or convergence
screen. All options have positive payoffs in every replication for both methods;
that does not establish tail accuracy. Maximum single-observation shares across
option estimates are 4.80% random and 5.48% fixed. Those smaller shares relative
to the preceding experiment also reflect the threefold increase in total
observations and cannot be attributed to allocation alone.

| Days | Random forward ± SE | Fixed forward ± SE |
| --- | ---: | ---: |
| 27 | 20.88844 ± 0.08478 | 20.85717 ± 0.08500 |
| 55 | 20.28596 ± 0.10821 | 20.38502 ± 0.08774 |
| 118 | 21.20466 ± 0.11159 | 21.44278 ± 0.12734 |
| 153 | 21.12794 ± 0.08632 | 21.19731 ± 0.09201 |

Fixed allocation improves forward SE at 55 days, is essentially unchanged at
27 days, and worsens it at 118/153 days. Neither method clears the planned
forward SE tolerance 0.5/8 = 0.0625 at any maturity. Removing random component
selection is therefore not a uniform remedy for the forward bottleneck.

### Fresh moments and count checks

| Days | Exact normalized-model target | Random second moment ± SE | Fixed second moment ± SE |
| --- | ---: | ---: | ---: |
| 27 | 542.49611 | 541.27760 ± 2.96467 | 547.35144 ± 2.89395 |
| 55 | 542.49611 | 546.45377 ± 3.99778 | 549.54750 ± 4.18300 |
| 118 | 600.31899 | 596.45131 ± 5.04099 | 603.33812 ± 4.24759 |
| 153 | 600.31899 | 596.34376 ± 5.91338 | 592.37098 ± 5.06279 |

All eight residuals are within about 1.7 reported SE. In particular, the old
118-day -11.8 residual is not reproduced on these new streams. Coupled with
its exact-state quadrature attribution, this supports sampling sensitivity
as the leading explanation for that sample's warning. It does not prove
reliable confidence coverage or negligible finite-GH population bias.

Largest absolute paired 64-to-128 option changes divided by quote scale are
1.082/2.148/1.261/1.452 (random) and 1.531/2.211/2.791/0.450 (fixed), at
27/55/118/153 days. Paired SEs are saved. These unresolved changes prevent
calling the full surface converged even when moment estimates appear plausible.
Counts here are outer states per replication, not daily stock time steps.

The 12 selected fixed-allocation states have maximum GH32-to-GH48 VIX-level
changes 0.00285/0.00147/0.00553/0.01163 by maturity. These checks do not certify
complete-price or population bias and are distinct from the full original-state
118-day attribution above.

### Replication sensitivity

`replication_splits.csv` compares the first six with the last six checking
replications. These subsets are independent of each other within a frozen
method, but are reused checking data, not new validation or calibrated tests.
Across option contracts, second-half/first-half SE ratios range from
0.323 to 2.508 (random) and 0.229 to 3.062 (fixed). The fixed 153-day moment
SE changes from 4.33 to 9.50 between halves; the fixed 55-day moment means
are 541.63 and 557.47. The descriptive splits show how variable small-group
uncertainty estimates remain; they do not prove a bias or supply simultaneous
confidence coverage. Do not promote a favorable SE from one subset as stable.

## Decision and next bounded work

Do not adopt fixed allocation as a universal sampling fix or select a different
method at each expiry and call that fresh validation. Retain random targeted
allocation as the comparison reference, with fixed allocation as an experimental
alternative. Current option/forward accuracy gates remain unresolved under both.
All four maturities and forwards stay active; the original 118-day warning is
preserved with this additive attribution, not deleted or corrected away.

The completed studies narrow the problem: tested quadrature changes do not
explain the saved warning; fresh moments fluctuate around their targets;
allocation affects precision but does not consistently remove it. Outer tail
coverage and stability across scrambles remain priorities. No optimizer is
launched and no joint model ranking is accepted.

Next, test simple **known-mean Gaussian controls** at a fixed proposal, starting
with the forward bottleneck. For an orthonormal Gaussian direction u, E[u'z]=0
and E[(u'z)^2-1]=0 under the target f. Importance-weighted controls therefore
have known mean zero without assuming that the implemented VIX-squared
quadrature has the exact model mean. Fit coefficients on independent pilot
scrambles, freeze them, and evaluate fresh price/count changes and uncertainty
across complete replications. Preserve negative adjusted estimates and pilot
variation rather than clipping them. This is a proposed next experiment,
not a result of the present study. Keep conditional/time/normalizer gates
separate; controls would not by themselves certify those biases.

Allocation runtime was 1421.2 seconds. Correspondence replay passes for 36
captured source/input hashes, 864 unique actual scramble seeds disjoint from
prior studies (48 random-method nets plus 816 component nets), all 12288
full states and weights, 732 summary rows, per-replication estimates, payoff
hits/concentration, component counts and paired prefix differences. These
checks establish saved-artifact correspondence, not numerical accuracy.

## Reproduction and checks

Use the local pricing Python environment, single-thread BLAS and fresh output
directories:

```sh
env OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /private/tmp/pdv-options-venv/bin/python \
  -m matteo.model_comparison.option_pricing.vix_targeted_refinement \
  --output matteo/model_comparison/option_pricing/runs/NEW_REFINEMENT

env OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /private/tmp/pdv-options-venv/bin/python \
  -m matteo.model_comparison.option_pricing.vix_targeted_allocation \
  --output matteo/model_comparison/option_pricing/runs/NEW_ALLOCATION
```

Correspondence replay uses `verify_vix_targeted_refinement` and
`verify_vix_targeted_allocation`; descriptive method comparisons use
`audit_vix_targeted_allocation`. Full states, weights, allocation factors,
component identities, prefix changes and contribution/hit diagnostics are saved.
96 package tests and 25 driver tests pass (**121 total**), including the
17-component exact-prefix, effective-weight and analytic Gaussian tail-payoff
oracles. A successful artifact replay does not certify price accuracy.
