# Payoff-targeted VIX pilot study — 22 October 2026

## Scope

Keep the frozen Stage-3 model, state, risk premium, selected variance curve,
GH32 conditional rule, GH48 date-zero normalizer and Gauss time8. All four
maturities (27/55/118/153 days), 175 training options and four forwards remain
included. This study fits numerical sampling proposals, not model parameters
or market residuals. No new market calibration, SPX work or maturity exclusion.

The previous covariance-inflation proposal worsened actual pricing uncertainty.
Here test a narrower mechanism: payoff-learned **mean shifts** in a deterministic
score-sensitive subspace, retaining unit covariance in every component and the
full 36-dimensional state. Keep the old learned mixture as the baseline.

## Design, screening and fresh checking have different roles

1. Promote all saved observations in `m2_vix_score_sampling_001` to **design
   data**. They were fresh checks for the preceding design; they cease to be
   fresh validation for a proposal fitted to them. Regenerate their full
   whitened Gaussian states from saved seeds and proposals, checking exact
   likelihood weights. Per maturity, four disjoint batches contain 512 states
   each: four replications from each of the old mean and inflated proposals.
   Weights correct each observation's actual generating law before pooling.
2. Fit each candidate independently on those four batches. Also fit a final
   pooled proposal on all 2048 states; this is not an independent fifth pilot.
3. Screen all fits on the older `m2_vix_stratified_sampling_001` **random-mixture**
   observations only (8 x 128 states/maturity). These seeds are disjoint from
   the four pilot batches. This set is independent of the new centroid fits,
   conditional on the old frozen proposal, but was seen historically and now
   informs selection: it is design screening, not fresh validation. The four
   fits share this screening set, so their assessed performance is dependent.
4. Freeze one complete proposal and test against the baseline on new independent
   streams in `m2_vix_targeted_check_001`: eight scrambles/method/maturity with
   nested 32/64 points. That is 512 costly valuations per method/maturity,
   4096 in total, plus 24 selected-state GH48 checks. This is a bounded
   exploratory check, with fewer replications and half the total pricing
   budget of the preceding R16 x 64 comparison. It is not the larger final
   count-convergence study. Current methods have equal budgets; historical
   cross-run comparisons do not isolate proposal effects.

## Candidate proposals

Six predeclared variants combine geometry, target families and a shift bound:

- Eight leading terminal PCA coordinates, three option targets (scaled puts,
  calls and upper-strike-quartile calls), shift norm capped at 3 or 4.
- The same eight coordinates with two additional targets, VIX and VIX squared,
  shift norm capped at 3 or 4.
- Four deterministic model-score directions U from the preceding experiment,
  with all five targets and shift norm capped at 3 or 4.

Each target learns two centers from its weighted pilot mean/covariance using
`fit_payoff_mixture`. For the score geometry, fit on U'z and lift each center
back as U m. Only shifts are restricted to this subspace; all 36 Gaussian
coordinates are sampled. Pilot target scaling cancels in its normalized
weighted law. This internal centroid calculation does not self-normalize any
pricing estimate.

Blend the newly fitted proposal equally with the old learned mixture:

    g_new = 0.5 g_old + 0.5 g_learned.

Both have half their probability at the standard Gaussian. After merging the
two identical origin components, the new design retains that half-Gaussian
component and exact weights f/g_new <= 2. The five-target design has 17 mixture
components versus seven in the baseline. This changes the full Sobol integrand,
including component selection; unit component covariances do not imply that
its randomized-Sobol variance must improve.

## What the pilot screen estimates

For a nonnegative payoff h and candidate proposal g, the IID importance
estimator's **second moment** is

    J_g(h) = integral f(z)^2 h(z)^2 / g(z) dz.

If a saved observation came from source density q, estimate it with

    (f(z)/q(z)) * (f(z)/g(z)) * h(z)^2.

Both density ratios are essential. This is not the squared price estimate,
not its variance (which subtracts the squared mean), and not a prediction of
scrambled-Sobol SE. Ratios to the baseline are formed after pooling complete
screening scrambles, then averaged over contracts within five families
(puts, calls, far calls, forward, second moment). Far calls also occur in the
call family, intentionally. The overall score equally averages those five
families and four maturities. It is a proposal-design proxy, not the market
calibration objective or an accuracy gate.

For the score-four/cap-four candidate, the overall ratio across the four
separate fitted pilots is 0.88705/0.88907/0.88936/0.88930; the pooled fit gives
0.88251. This is more stable across pilot fits than the PCA alternatives.
Cap-three score fits are close (0.89160–0.89393); selection within this small
grid can still be optimistic. The pooled cap-four proposal has these ratios:

| Days | Puts | Calls | Far calls | Forward | VIX-squared target |
| --- | ---: | ---: | ---: | ---: | ---: |
| 27 | 0.9999 | 0.8865 | 0.9272 | 0.9193 | 0.8949 |
| 55 | 0.9945 | 0.8856 | 0.8547 | 0.9245 | 0.9044 |
| 118 | 1.0193 | 0.7892 | 0.6095 | 0.9353 | 0.9058 |
| 153 | 1.0288 | 0.7562 | 0.6187 | 0.9221 | 0.8738 |

The proxy favors tail calls but slightly worsens long-maturity puts. Select
this single pooled proposal for the fresh exploratory check; do not report
these screen ratios as achieved price-variance reductions.

## Fresh results

The fresh equal-budget check is mixed:

| Days | Median option SE: targeted / baseline | Contracts with lower targeted SE | Baseline SE screen | Targeted SE screen |
| --- | ---: | ---: | ---: | ---: |
| 27 | 1.165 | 12/40 | 0/40 | 2/40 |
| 55 | 0.592 | 43/47 | 4/47 | 9/47 |
| 118 | 1.066 | 16/44 | 8/44 | 8/44 |
| 153 | 0.500 | 34/44 | 7/44 | 5/44 |

The ratios are medians of per-contract SE ratios. The SE screen requires
observed scramble SE <= quote scale / 8 and at least two positive payoffs;
it is necessary but does not establish accuracy. Totals are **19/175 baseline
versus 24/175 targeted**. Although 105/175 contracts have lower targeted SE,
no uniform gain or full-surface accuracy is established. At 153 days the
median improves while the screen count worsens: these summaries measure
different aspects of the contract distribution.

Far-call median ratios (upper strike quartile) are 0.800/0.625/0.861/1.031.
Thus even the strong 153-day overall call benefit is not a uniform far-wing
gain. All 175 pooled option estimates have positive payoffs under both designs;
two targeted 153-day quotes nevertheless have at least one zero-hit scramble.
Maximum pooled single-observation shares across option prices are 21.31%
for baseline and 14.05% for targeted sampling. Better concentration does not
clear price-error gates.

| Days | Baseline forward ± SE | Targeted forward ± SE |
| --- | ---: | ---: |
| 27 | 20.74465 ± 0.16486 | 20.91660 ± 0.26994 |
| 55 | 20.78199 ± 0.11850 | 20.30932 ± 0.20791 |
| 118 | 21.17184 ± 0.21865 | 21.04728 ± 0.22854 |
| 153 | 21.44688 ± 0.26628 | 21.42343 ± 0.32535 |

Every forward SE worsens, despite the favorable IID second-moment screen.
No forward passes the planned observed SE tolerance 0.5/8 = 0.0625. This is
an explicit discrepancy between a favorable design proxy and actual randomized
Sobol pricing precision, not an implementation identity failure. The screen
cannot substitute for fresh scrambling experiments.

### Retain the new moment warning

| Days | Exact normalized-model target | Baseline second moment ± SE | Targeted second moment ± SE |
| --- | ---: | ---: | ---: |
| 27 | 542.49611 | 536.88009 ± 7.88328 | 545.95883 ± 12.33406 |
| 55 | 542.49611 | 558.22397 ± 7.94602 | 542.36699 ± 7.94601 |
| 118 | 600.31899 | 587.65337 ± 7.57303 | 588.51408 ± 4.30771 |
| 153 | 600.31899 | 589.12354 ± 18.55850 | 609.87187 ± 10.22932 |

The **118-day targeted residual is -11.80491 with reported SE 4.30771**,
roughly -2.74 SE. Baseline at 55 days also has a +15.72785 residual with
SE 7.94602. These warnings must remain in the record; smaller SE does not
mean greater accuracy when tails or quadrature error are unresolved. The
preceding fresh moment agreement was one sample-specific result, not a durable
convergence clearance. No normalizer change, forced moment control or price
rescaling is applied.

Actual sampling error, finite-scramble uncertainty estimation and conditional
quadrature bias remain possible contributors. The earlier direct Gaussian
moment reference addresses a different integral and cannot certify this
finite-conditional-GH outer estimator. This run does not resolve that
attribution. An empirical SE from eight scrambles is not a calibrated tail
confidence bound.

### Count and conditional checks

Largest absolute paired 32-to-64 state option changes, divided by quote scale,
are 1.936/3.465/1.103/1.468 (baseline) and 3.571/0.874/2.173/2.064 (targeted).
Paired scramble SEs are saved; the changes leave count-convergence open.
These counts are outer states, not time steps per day.

The 24 selected-state GH32-to-GH48 checks have maximum VIX-level change
0.00358, at the 153-day baseline median state. Targeted maxima by maturity
are 0.00081/0.00191/0.00298/0.00313. They do not bound aggregate moment or
option bias. No new time-rule or normalizer refinement is run here; keep
those gates separate from the observed sampling SE.

## Decision and next bounded step

Retain the targeted mean-shift design as an **experimental alternative**, not
a replacement for the baseline or a qualified calibration engine. It improves
some option families, but every forward remains less precise and the new
118-day moment warning blocks a convergence claim. Keep all four maturities.

Next isolate the sampling mechanism before fitting another proposal: use the
same frozen 17-component density to compare random component selection with
fixed allocations and correct allocation weights. This can test whether
component selection interacts badly with Sobol integration, but that explanation
is only a hypothesis. Include the baseline, enough independent scrambles and
larger nested counts; retain moment/forward gates and replay the offending
118-day states for paired conditional/normalizer refinement. Do not infer a
quadrature bug from this warning alone or tune to force its sign away.

Current validation: **96 package + 23 driver tests = 119 passing checks**.
Pilot replay passes for seven files/24 array blocks; fresh-check replay passes
for 40 captured source/input hashes, 64 unique independent checking seeds,
all 4096 full states and weights, and 732 summary rows including paired
prefix differences. Full states reproduce exactly; maximum weight replay
error is 2.22e-16. These are artifact/implementation checks, not price
certificates. Fresh pricing took 440.7 seconds.

## Reproduction and evidence boundary

Use single-thread BLAS and the local pricing environment. Fresh output
paths are mandatory:

```sh
env OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /private/tmp/pdv-options-venv/bin/python \
  -m matteo.model_comparison.option_pricing.vix_targeted_pilots \
  --output matteo/model_comparison/option_pricing/runs/NEW_PILOT_DIRECTORY

env OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /private/tmp/pdv-options-venv/bin/python \
  -m matteo.model_comparison.option_pricing.vix_targeted_check \
  --output matteo/model_comparison/option_pricing/runs/NEW_CHECK_DIRECTORY
```

The checking driver explicitly reads `m2_vix_targeted_pilots_001`; its selected
candidate is `score4_balanced_cap4`, pooled pilot index 4. An alternative design
requires an explicit source revision and new frozen run, not overwriting inputs.
`verify_vix_targeted_pilots` reproduces all seven design files and 24 array
blocks exactly, including 42 source/input snapshot hashes. This reuses the
tested deterministic design implementation; it is an artifact replay, not an
independent performance assessment. `verify_vix_targeted_check` replays fresh
states, weights and estimates. `audit_vix_targeted_check` compares methods.

New tests independently check the blended density, the two-density-ratio
second-moment formula against a Gaussian exponential oracle, and subspace
lifting/full-state retention. The initial synthetic test accidentally used a
put below the minimum synthetic VIX level and correctly triggered the existing
zero-target rejection; the test fixture was corrected, not the production
pricing code. No Lean theorem or numerical error certificate is claimed.
