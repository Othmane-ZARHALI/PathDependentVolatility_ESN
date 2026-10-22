# Gaussian-control pilot and count stability — 22 October 2026

## Scope and frozen design

Continue [VIX_GAUSSIAN_CONTROL_RESULTS.md](VIX_GAUSSIAN_CONTROL_RESULTS.md)
without changing the model, risk premium, state, variance curve, targeted
17-component density, four score directions, nine control functions or ridge
parameter (lambda=1). All four maturities (27/55/118/153 days), 175 options
and four forwards remain included. No market calibration, maturity exclusion,
price rescaling or forced VIX moment is performed. SPX is unchanged.

The prior 12-pilot fit improved fresh forward SE but varied across six-pilot
halves. This study asks whether more independent pilot replications stabilize
that gain and whether it survives a larger Sobol net. The lean-numerics-bridge
workflow keeps exact Gaussian control identities, code correspondence, pilot
selection, fresh checking and numerical bias separate. No Lean certificate or
rigorous price-error bound is claimed.

## Larger pilots, explicit reuse and frozen alternatives

Promote raw observations from two compatible completed studies to pilot data:
12 random-targeted 128-point replications from
`m2_vix_targeted_allocation_001`, followed by 16 raw 128-point replications from
`m2_vix_gaussian_controls_001`. Reconstruct all states and exact density weights;
verify model, curve, target, rule and proposal correspondence before pooling.
These historical observations cease to be fresh validation for the new fit.

There are 28 independent complete pilot replications per maturity. The primary
`gaussian` method fits all 28. Two diagnostic fits use disjoint, balanced
14-replication batches: six from each half of the original batch plus eight
from each half of the recent batch. Thus each fit contains both historical
sources. The original 12-replication coefficient matrix is retained unchanged as
`legacy12`, and checked against re-fitting its original pilot.

Six estimators share every fresh conditional valuation: raw, mass-only fitted
on all 28 replications, primary Gaussian fitted on all 28, the two balanced
14-pilot Gaussian fits, and the unchanged legacy12 Gaussian fit. No alternative
is selected after checking. Same zero-mean control rows as before:

    C = (w-1, w U'z, w((U'z)^2-1)),   w=f/q.
    estimator = mean(w H) - mean(C)' beta.

Pilot regression uses complete-replication covariance, standardized controls
and the same fixed ridge lambda=1. Checking controls are not centered by an
observed pilot or checking mean. Exact Gaussian control expectations do not
force the VIX-squared numerical integral to its exact-model target. Negative
adjusted estimates, if any, are retained.

## Fresh nested-count check

`m2_vix_gaussian_stability_001` uses **12 new independent scrambles per
maturity with nested 128/256 points**, or 3072 expensive conditional values
per maturity and 12288 total. Six estimators share those paths, so estimator
and count changes are paired. The complete scramble is the uncertainty unit.

Compared with the previous R16 x 128 study, the total point budget is 1.5 times
larger and per-net resolution doubles, but checking SE has 11 rather than 15
variance degrees of freedom. Historical cross-run SE differences are therefore
not a clean proposal/control comparison. Within the current run, all methods
have equal budgets and identical paths.

All pilot replications have 128 points. Keep every coefficient fixed at both
128 and 256 checking counts, so paired changes isolate the complete fixed
estimator's count sensitivity. Its coefficients need not be optimal for a
256-point net; this transfer is part of the test. The existing GH32/M48/time8
rules remain unchanged. Twelve minimum/median/maximum first-scramble states
receive a GH48 sensitivity check; no full bias gate is cleared by that subset.

## Results

### All four primary forward SE screens pass, with a narrow final pass

| Days | Raw forward ± SE | Unchanged legacy12 forward ± SE | Pooled 28-pilot forward ± SE |
| --- | ---: | ---: | ---: |
| 27 | 20.77313 ± 0.06684 | 20.79905 ± 0.04432 | 20.78816 ± 0.04230 |
| 55 | 20.29204 ± 0.06675 | 20.29824 ± 0.04127 | 20.30596 ± 0.03750 |
| 118 | 21.15207 ± 0.07295 | 21.16889 ± 0.06057 | 21.16054 ± 0.03961 |
| 153 | 21.10439 ± 0.08276 | 21.10376 ± 0.06728 | 21.13570 ± 0.06157 |

The pooled fit passes the observed forward SE threshold 0.5/8 = 0.0625
at all four maturities; legacy12 passes three, mass-only two and raw none.
The 153-day primary pass is narrow: SE 0.06157 is only about 1.5% below the
threshold. It is not a robust accuracy clearance or a reason to ignore pilot
sensitivity, sample-count changes or numerical bias.

Primary/legacy12 forward-SE ratios on these identical checking paths are
0.954/0.909/0.654/0.915. The strongest improvement from the larger pilot is at
118 days. Primary/raw ratios are 0.633/0.562/0.543/0.744. These are observed
conditional checking gains, not guaranteed population efficiency or calibrated
confidence statements.

### Balanced-pilot sensitivity

| Days | Batch A forward SE | Batch B forward SE | Both pass forward SE screen? |
| --- | ---: | ---: | --- |
| 27 | 0.04775 | 0.04796 | Yes |
| 55 | 0.04814 | 0.04065 | Yes |
| 118 | 0.04309 | 0.04325 | Yes |
| 153 | 0.06075 | 0.06666 | No |

Both 14-pilot fits support the first three forward sampling passes. One
still fails at 153 days, as does the unchanged legacy12 fit there. The pooled
fit shares both pilot batches, so the three estimates are not independent
training experiments; all also share checking paths. This supports the first three passes under two pilot fits, without
establishing uniform stability or allowing the best fit
to be selected per maturity after checking.

### Option-tail accuracy remains the main unresolved sampling target

| Days | Median option SE: pooled / raw | Raw SE screens | Pooled SE screens |
| --- | ---: | ---: | ---: |
| 27 | 1.019 | 8/40 | 7/40 |
| 55 | 0.887 | 11/47 | 11/47 |
| 118 | 0.955 | 12/44 | 12/44 |
| 153 | 0.918 | 16/44 | 15/44 |

Only **45/175 primary option estimates pass**, versus **47/175 raw** and
44/175 legacy12. The screen requires positive adjusted estimates, at least
two raw payoff hits and SE <= quote scale / 8; it is not full accuracy.
Median SE decreases at three maturities but screen counts do not improve.
The prior 27-day median option gain did not persist on these new larger nets.

All option estimates are positive for every variant, and every raw option
has positive payoffs in every scramble. No negative value is clipped. Observing
hits everywhere still does not establish tail precision. The forward-oriented
linear/quadratic controls therefore remain useful, but do not solve the full
option-surface sampling problem.

### Retained moment warnings

| Days | Exact normalized-model target | Raw second moment ± SE | Pooled second moment ± SE |
| --- | ---: | ---: | ---: |
| 27 | 542.49611 | 542.51042 ± 2.86213 | 543.53483 ± 2.36037 |
| 55 | 542.49611 | 539.08578 ± 2.66329 | 539.73402 ± 2.09350 |
| 118 | 600.31899 | 593.77160 ± 3.43867 | 595.76038 ± 2.24254 |
| 153 | 600.31899 | 598.60876 ± 4.34392 | 598.08407 ± 3.03565 |

**118-day warning:** pooled residual -4.55861 ± 2.24254 is about -2.03 reported
SE. Batch A gives residual -4.21342 ± 1.98663, about -2.12 SE. These are the
two entries in `moment_warnings.csv`, which flags all fine-count variants
beyond two reported SE. Fits and maturities are dependent/multiple diagnostics;
these are warning flags, not calibrated tests of model bias or a normalization
bug. The unchanged legacy12 residual is -4.38674 ± 2.80996. No moment target is
imposed and no price rescaling removes the warnings.

These data do not clear the population finite-GH moment question, unseen-tail
coverage or the reliability of small-replication SEs. Earlier exact-state
refinements remain separate configuration-specific evidence, not a uniform
bound for the current new states or all model parameters.

### Paired count and conditional-rule checks

Primary 128-to-256 forward changes ± paired SE are +0.04842 ± 0.03944,
+0.06131 ± 0.03672, -0.01881 ± 0.07156 and -0.06073 ± 0.06967 at
27/55/118/153 days. The changes reflect finite-sample fluctuation of two
estimators with the same expectation; they are not direct estimates of a
time-discretization bias. Their uncertainty still prevents treating a single
fine-count SE screen as full price convergence.

Largest absolute paired option changes divided by quote scale are
0.433/0.551/1.459/0.563 (raw), 0.387/0.961/1.234/0.579 (pooled) and
0.270/1.174/1.202/0.746 (legacy12). Full paired SEs are saved. These are outer
sample counts, not time steps per trading day. Stock integration is unchanged.

The 12 selected-state GH32-to-GH48 VIX changes have maxima
0.00055/0.00041/0.00326/0.00120 by maturity. These subset checks do not bound
global conditional, time-rule or normalizer bias. No new full-state bias
refinement was performed in this study.

## Decision and next bounded work

The larger pilot strengthens forward variance reduction, especially at
118 days, and the first three forward SE passes are supported by both
balanced fits. Retain the pooled 28-pilot method as an experimental reference,
with raw and legacy estimates alongside it. The 153-day pass is borderline,
not stable across pilot batches. No market refit, joint feasibility claim,
held-out result or model ranking follows.

Prioritize the remaining option-tail problem rather than assuming more pilot
replications alone will solve it. Next run a bounded pilot screen of
payoff-shaped Gaussian controls with analytically known means, retaining the
current dictionary/ridge as the baseline. For example, centered Gaussian tail
hinges can target rare payoff regions without imposing a VIX moment identity.
Any new dictionary or coefficient selection must use design-only data and
fresh subsequent checking. Keep all four maturities and all accuracy/bias
gates; do not drop the 153-day case or increase tolerances to declare success.

Validation: **99 package + 27 driver tests = 126 passing checks**. Pilot/state/
control/statistics replay passes for 45 captured source/input hashes, 48 unique
new checking seeds disjoint from earlier studies, all 12288 full states and
control rows, every pilot fit/split/legacy coefficient, 2196 summary rows and
all paired estimator/count changes. This is saved-artifact correspondence,
not pricing certification. The run took 1393.8 seconds.

## Reproduction and correspondence

Use a fresh output directory and the local single-thread pricing environment:

```sh
env OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /private/tmp/pdv-options-venv/bin/python \
  -m matteo.model_comparison.option_pricing.vix_gaussian_stability \
  --output matteo/model_comparison/option_pricing/runs/NEW_STABILITY_RUN
```

Inputs are explicitly the two historical source directories above. Source
snapshots are captured before pricing and mid-run changes are rejected. Pilot
arrays, split indices, coefficient/proposal hashes, full checking states and
control rows, raw VIX values, all replication estimates and paired changes are
saved. `verify_vix_gaussian_stability` reconstructs the two pilot sources,
balanced fits, legacy fit and checking statistics; `audit_vix_gaussian_stability`
produces descriptive screens and pilot-fit sensitivity. Neither certifies price
accuracy, tail coverage or the finite conditional quadrature approximation.
