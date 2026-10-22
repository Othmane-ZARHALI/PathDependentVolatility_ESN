# Model-score VIX sampling — 22 October 2026

## Scope and question

Follow the rare-tail attribution in [VIX_MOMENT_DIAGNOSIS.md](VIX_MOMENT_DIAGNOSIS.md):
can an outer proposal chosen from the frozen model's score geometry improve
conditional VIX sampling? Keep **all four maturities (27/55/118/153 days),
175 training options and all four forwards**. The model, state, risk premium,
normalizer and previously selected nuisance curve remain frozen. No optimization,
new market fit, maturity removal, or price/moment rescaling is performed.

## Proposal fixed before evaluation

Write the complete 36-dimensional terminal Gaussian state as X_T=mu_T+L_T z,
where the target law of z is N(0,I). At future lag s, the conditional projection
mean changes by A_s z, with A_s=P exp(-Lambda s) L_T. For the frozen readout

    score(C,Y,F) = a C + b (Y-y0)^2 + c (F-f0)^2,

let D=diag(0,b,c), q_s be the score gradient at the unconditional future mean,
and v_s=A_s' q_s. Construct the positive sensitivity matrix

    H = sum_s (time_weight_s / delta) [A_s' D A_s + v_s v_s'].

Its leading four orthonormal eigenvectors form U. This is a deterministic
heuristic for directions that change the score: no pilot price, checking state,
option residual or observed checking error selects U, rank or inflation scale.
The four directions capture 99.74%–99.77% of trace(H) across the maturities.
That fraction concerns this sensitivity matrix, not captured option variance
or a mathematical bound on price error. The remaining 32 directions are
sampled too; the state is not truncated to four dimensions.

The proposal is

    g(z) = 0.5 phi_I(z) + 0.5 phi_Sigma(z),
    Sigma = I + (2.5^2-1) U U',
    w(z) = phi_I(z)/g(z)
         = 1 / [0.5 + 0.5 * 2.5^(-4)
                  * exp((1-2.5^(-2)) ||U'z||^2 / 2)].

A scrambled Sobol coordinate selects the Gaussian component; remaining
coordinates generate a full Gaussian vector in an orthogonal basis beginning
with U. Only the four leading coordinates are inflated in the broad component.
Every forward, second moment and option payoff uses its exact weight w;
0<w<=2. We never divide by the observed mean weight. Broader Gaussian tails
therefore change the numerical integration distribution and preserve the target
pricing measure.

For comparison, replay the preceding pilot-learned mean-shift proposal on new
independent streams. Both complete designs are frozen before evaluation. The
new design changes covariance and the Sobol rotation, so the comparison cannot
isolate the effect of inflation alone. It also compares against one fixed old
pilot, not the distribution of outcomes over independently learned pilots.

## Budget, uncertainty and reproducibility

`runs/m2_vix_score_sampling_001` uses 16 independent scrambles per method per
maturity, with nested 32/64 points per scramble: **1024 expensive conditional
valuations per method per maturity** at the fine count. This matches the total
number in the earlier R8 x 128 studies, but changes per-net resolution and
uncertainty degrees of freedom. Historical estimates are not paired controls.
The two current methods have independent evaluation streams. Count differences
are paired within a scramble; the scramble, not a point or mixture component,
is the uncertainty unit. Reported values after ± are one standard error, not
confidence guarantees, and exclude quadrature bias.

Conditional GH32, date-zero normalizer GH48 and geometric-panel Gauss time8
match the preceding full-state design. At the first scramble of each method
and maturity, re-evaluate the minimum, median and maximum sampled VIX states
with GH48 (24 selected states total). These are sensitivity diagnostics,
not a uniform conditional-quadrature accuracy check. The independent direct
Gaussian moment reference is retained separately; no moment identity is
imposed on option pricing.

The run saves source snapshots captured before evaluation, input/proposal
hashes, seed records, all 8192 full terminal states, their exact weights,
conditional VIX values, replication estimates, contribution concentration,
payoff hits and paired count changes. Source mutation during evaluation is
rejected. `verify_vix_score_sampling` independently regenerates proposal draws,
weights and full states and recomputes saved statistics; it does not recompute
the costly conditional functional or certify convergence.

Reproduce with a fresh output directory and single-thread BLAS:

```sh
env OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /private/tmp/pdv-options-venv/bin/python \
  -m matteo.model_comparison.option_pricing.vix_score_sampling \
  --output matteo/model_comparison/option_pricing/runs/NEW_DIRECTORY
```

## Numerical results

### The score proposal did not improve pricing precision

At the equal fine budget (16 x 64 observations), median **per-contract** SE
ratios and nominal SE screens are:

| Days | Median option SE: score / learned | Quotes with lower score SE | Learned SE screen | Score SE screen |
| --- | ---: | ---: | ---: | ---: |
| 27 | 2.242 | 0/40 | 3/40 | 0/40 |
| 55 | 1.858 | 7/47 | 7/47 | 5/47 |
| 118 | 1.358 | 3/44 | 9/44 | 2/44 |
| 153 | 1.806 | 11/44 | 9/44 | 4/44 |

The screen is observed scramble SE <= quote scale / 8, with at least two
positive payoffs; it is necessary, not sufficient for accuracy. Total passes
are **28/175 versus 11/175**. Every option has positive payoffs in every
scramble for both methods, yet most estimates fail even this sampling screen.
The maximum single-observation contribution across option prices is 18.25%
for the learned design and 36.80% for inflation. Broad exploration of extreme
states does not itself give good weighted-payoff precision.

| Days | Learned forward ± SE | Score forward ± SE |
| --- | ---: | ---: |
| 27 | 20.85556 ± 0.15259 | 20.64398 ± 0.21062 |
| 55 | 20.17249 ± 0.18921 | 20.18098 ± 0.24190 |
| 118 | 21.14720 ± 0.16872 | 20.89871 ± 0.18598 |
| 153 | 21.23163 ± 0.19760 | 21.17468 ± 0.30404 |

Forward SE also increases at all four maturities. These are numerical estimates
at a frozen exploratory curve, not a fit acceptance or market-forward claim.

### Fresh moment checks favor sampling sensitivity, not a normalizer patch

| Days | Exact normalized-model target | Learned second moment ± SE | Score second moment ± SE |
| --- | ---: | ---: | ---: |
| 27 | 542.49611 | 543.72377 ± 5.74914 | 523.70744 ± 14.00625 |
| 55 | 542.49611 | 540.53907 ± 6.89610 | 528.50917 ± 11.93495 |
| 118 | 600.31899 | 600.97753 ± 6.89445 | 580.64744 ± 9.26286 |
| 153 | 600.31899 | 593.86271 ± 6.18110 | 614.21284 ± 20.54145 |

All four fresh learned-mixture residuals are within 1.1 reported SE of their
targets. In particular, the old 27-day +13.65 residual is not reproduced on
these fresh draws. Together with the independent direct-Gaussian reference
and old-state refinement checks, this supports sampling sensitivity as the
leading explanation. It does not prove absence of conditional quadrature
bias, reliable tail confidence coverage, or convergence of options/forwards.
The score proposal is less precise for the moment at every maturity and
still has sizable residuals. No identity is imposed on either method.

Average weights range from 0.98919 to 1.00645 for learned sampling and from
0.98450 to 1.00061 for inflation, with SEs roughly 0.0063–0.0081. These are
consistent with target mass one at this diagnostic resolution; they are not
used to rescale prices. Individual inflated weights can be extremely small:
retaining exact weights matters even when a simulated VIX level is very large.

### Count and conditional-rule checks remain separate

The largest absolute paired 32-to-64 point option changes, divided by quote
scale, are 3.147/3.051/0.729/1.529 for the learned design and
2.348/3.244/2.372/1.834 for inflation at 27/55/118/153 days. The saved paired
SEs are available in `prefix_changes.csv`; these observed changes and finite
scramble SEs do not establish convergence. Here 32/64 are **outer states per
scramble**, not SPX time steps per trading day. SPX is unchanged by this run.

The largest selected-state GH32-to-GH48 VIX-level change is 0.01135, at the
118-day inflated maximum-VIX state. Learned-design maxima by maturity are
0.00070/0.00158/0.00034/0.00679; inflated-design maxima are
0.00179/0.00057/0.01135/0.00150. These 24-state checks give no global bound.
This experiment does not repeat the separate time/normalizer refinements;
those earlier results remain configuration-specific evidence.

## Decision and next bounded work

**Do not adopt this rank-four, scale-2.5 score proposal.** It has poorer median
option and forward precision at every maturity at the tested equal budget.
This is a negative result for this complete design, not a rejection of all
covariance inflation or all score-informed importance sampling. The 99.7%
sensitivity-trace fraction must not be reported as captured payoff variance.
Keep the pilot-learned mixture as a benchmark, without calling it converged.

Next, compare more targeted tail proposals using independent design pilots
and measure pilot-to-pilot variation. Include weighted conditional variance,
forwards and actual option payoffs in pilot diagnostics, rather than optimizing
a raw-score surrogate alone. Freeze the chosen design before checking it on
new draws with both more points per scramble and enough independent scrambles.
Compare the learned baseline at the same budget; keep paired count refinements
and conditional/time/normalizer error checks separate. Do not tune on this
checking run and present the same observations as fresh validation. All four
maturities, 175 quotes and four forwards remain mandatory. No new fit or model
ranking is justified yet.

## Validation and retained artifacts

96 package tests and 20 driver/legacy tests pass (**116 total**), including
analytic Gaussian mass, first/second moments and a three-sigma call oracle
under oblique low-rank inflation. Full-dimension and exact-prefix checks pass.
The independent artifact verifier checks 38 source/input snapshot hashes,
128 unique checking seeds disjoint from prior outer-study seeds, all 8192 full
states and weights, 732 summary rows, per-scramble estimates, contribution
concentration/hit counts and paired prefix differences. These are implementation
and correspondence checks, not formal Lean certification or price-error bounds.

The verifier's first invocation stopped on a Python-tuple versus JSON-list
representation comparison for risk-premium metadata. Canonical JSON conversion
fixes that comparison; no model parameter, saved state or price was changed.
The successful verification and the failure's nature are recorded here.

Run `verify_vix_score_sampling RUN_DIRECTORY` and `audit_vix_score_sampling
RUN_DIRECTORY` through the same Python module prefix as the pricing command.
See `verification.json`, `comparison.json`, `method_comparison.csv` and
`option_comparison.csv` in the run directory. The frozen pricing run took
1197.4 seconds. No optimization or held-out market evaluation was run.
