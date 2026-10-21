# VIX moment-discrepancy diagnosis — 21 October 2026

## Question and retained scope

The preceding 27-day estimates of E[VIX^2] were 554–556 against the exact
normalized-model target 542.49611, with reported scramble SE around 3. The
mass-centered residual remained nonzero. Investigate that warning without
rescaling prices or imposing a moment constraint.

All four maturities, 27/55/118/153 days, and associated forwards remain in
scope. The present diagnostics examine a model identity at every maturity;
they do not reprice or refit the 175 market options or the forwards. Preserve
the frozen Stage-3 candidate-18/scenario-0 model, zero projected initial state,
risk premium and selected variance curve. No held-out market observations enter.

## Separate the mechanisms

For the exact date-zero normalization M0(t)=E[G(eta_t)|X0],

    V_t = floor + (xi(t)-floor) G(eta_t)/M0(t),
    E[VIX_T^2] = (10000/delta) integral_T^(T+delta) E[V_t] dt
               = (10000/delta) integral_T^(T+delta) xi(t) dt.

The identity uses the tower property and conditional averaging before the
square root. It concerns mean VIX squared, not the mean VIX level.

1. **Gaussian-law algebra.** At every integration node, check that propagating
   the full terminal OU mean/covariance and adding conditional innovations
   reproduces the direct future projection law. Check the terminal covariance
   factor independently. This tests the law used in outer/inner composition.
2. **Direct marginal expectation.** Integrate G under the unconditional future
   three-dimensional projection law, bypassing full-state outer sampling and
   the conditional GH engine. Use the same date-zero normalizer at every state.
   Compare normalizer orders 48/64 and time orders 8/16 on paired random inputs.
   Check sparse node references at GH96/128 independently.
3. **Population conditional-GH rule.** Integrate the implemented finite GH32,
   GH48 and GH64 conditional rules over the exact Gaussian distribution of
   conditional means. This bypasses the old full-state outer sample while
   retaining its conditional quadrature approximation. Compare each rule with
   a continuous conditional innovation using paired inputs.
4. **Multiplier concentration.** At the middle future date in each VIX window,
   estimate probabilities and contributions from G>100, 1000, 10000, 100000.
   These are fixed-model diagnostics, not empirical market tail probabilities.

These numerical checks are not formal Lean proofs or universal error bounds.
The implementation/equation correspondence and numerical outcomes are distinct.

## Runs and numerical designs

- `m2_vix_collapsed_moment_001`: plain three-dimensional Sobol, eight scrambles,
  nested 2^12/2^14/2^16 points, time orders 8/16 and normalizer orders 48/64.
- `m2_vix_collapsed_inflated_001`: direct marginal calculation with eight new
  scrambles, nested 2^14/2^16/2^18 points and a defensive Gaussian proposal.
- `m2_vix_composed_moment_002`: population conditional-GH diagnostic, eight new
  scrambles and nested 2^12/2^14/2^16 points, conditional orders 32/48/64,
  normalizer orders 48/64 and time orders 8/16.
- `m2_vix_multiplier_tail_001`: independent middle-date tail audit at all four
  maturities, eight scrambles of 2^18 points, with GH48/64/96/128 references.

The first composed-run directory, `m2_vix_composed_moment_001`, is a preserved
startup failure: a relative source-input path was rejected while saving its
snapshot, before numerical computation. The subsequent driver resolves the
path explicitly. This is unrelated to the pricing moment discrepancy.

### Direct expectation and importance weights

For future date t, draw z standard normal in three dimensions and evaluate
G(mu_t + sqrt(Sigma_t) z). Reusing the same z at different dates need not give
the correct path law, but suffices for this linear expectation of an integral.
**The resulting integrands are not VIX draws; taking their square roots cannot
price VIX futures or options.**

The defensive proposal is g(z)=0.5 phi(z)+0.5 phi_{2.5}(z), where the latter
Gaussian has covariance 2.5^2 I. Its likelihood weight is

    w(z) = 1 / [0.5 + 0.5 * 2.5^(-3)
                         exp((1-2.5^(-2)) ||z||^2 / 2)].

Apply exact weights to the direct marginal integrands. Weights are bounded by
2 and are never divided by their observed mean. This changes the sampling
proposal, not the model distribution or normalizer. The plain run is retained,
including its unexpectedly large sampling uncertainty.

### Conditional composition without full-state outer pricing

At lag s=t-T, the Gaussian conditional mean of eta_t is distributed as
mu_t + sqrt(B_s) z_outer, where B_s=P D_s Sigma_T D_s P'. Conditional
innovations have covariance C_s. Independently sample z_outer and z_inner;
then mu_t+sqrt(B_s)z_outer+sqrt(C_s)z_inner has the correct future marginal law.

For the finite GH rule, replace z_inner by independent discrete GH nodes with
their exact product probabilities. Averaging over these nodes and the outer
Gaussian then reproduces the population mean of the existing conditional GH
approximation. This is a six-dimensional numerical integral, not a claim that
the full OU state is three-dimensional or that the projections are Markov.

The variance-inflated proposal broadens the outer Gaussian and tilts inner GH
atom probabilities by exp((1-2.5^(-2)) node^2/2). A half-original defensive
mixture retains support; exact continuous/discrete joint density ratios are
used for every rule. Finite GH total-mass roundoff is retained explicitly.
The constant variance floor is integrated analytically in this diagnostic;
this is not an imposed G normalization or a forced VIX second moment.

Standard errors use complete independent scrambles. Time/normalizer/conditional
rule changes reuse random inputs and have paired uncertainty. Comparisons
between independent diagnostic runs are not treated as paired.

## Results and attribution

### Gaussian laws agree algebraically at the tested nodes

Across all four maturities and both time rules, maximum absolute mean,
covariance and terminal-root reconstruction errors are 4.16e-17, 3.33e-16
and 2.00e-15. These checks give no indication of a mismatch between the
terminal sampler and future conditional Gaussian laws. They are floating-
point correspondence checks, not formal proofs.

### Direct independent moment reference

The variance-inflated direct-marginal calculation gives:

| Days | Exact normalized-model target | Direct estimate ± scramble SE | Residual |
| --- | ---: | ---: | ---: |
| 27 | 542.49611 | 542.49136 ± 0.09509 | -0.00476 |
| 55 | 542.49611 | 542.67142 ± 0.13291 | +0.17531 |
| 118 | 600.31899 | 600.41035 ± 0.14271 | +0.09136 |
| 153 | 600.31899 | 600.26792 ± 0.12181 | -0.05107 |

These results agree with the identity at the reported numerical resolution.
They sharply reduce the concern that a large date-zero normalization error
explains the old 27-day discrepancy. They do not certify the normalizer
uniformly over other parameters, dates or price tolerances.

Changing normalizer48→64 changes these moments by -0.00101, -0.00626,
-0.01091 and -0.00998 at 27/55/118/153 days. At the sparse reference nodes,
maximum relative differences from GH128 are 2.60e-5 (GH48), 1.46e-5 (GH64)
and 1.19e-6 (GH96). The paired direct-marginal time8→16 changes are below
2e-13 at this configuration; this smooth expectation calculation does not
by itself validate state-specific conditional integration or option prices.

The plain 3D run is retained: even 8x65536 samples give moment SEs of roughly
5.5–8.2. More points alone do not guarantee informative precision for this
readout. The larger weighted reference uses 8x262144 samples; its performance
should not be attributed solely to a proposal change without matching budgets.

### Exact original-state quadrature attribution

`m2_vix_matched_moment_refinement_001` reuses **all 8x128 original 27-day
random-mixture states**, exact original likelihood weights and the frozen
curve. It changes conditional GH and normalizer orders separately:

| Conditional GH / normalizer | Moment | Residual from 542.49611 | Paired change from original |
| --- | ---: | ---: | ---: |
| 32 / 48 | 556.16551 | +13.66939 | 0 |
| 48 / 48 | 556.15012 | +13.65401 | -0.01539 ± 0.01193 |
| 48 / 64 | 556.14907 | +13.65296 | -0.01643 ± 0.01193 |

Conditional GH32→48 and normalizer48→64 therefore do not explain or remove
the observed 27-day residual on that sample. On the first 16 states of each
replication, time8→16 changes the weighted moment by 8.44e-7 ± 4.48e-7, with
maximum VIX-level change 3.56e-5. On the first four states per replication,
conditional48→64 changes the moment by -0.00959 ± 0.01311, with maximum
VIX-level change 0.01012. These subset checks do not clear global price bias.

At 55/118/153 days the matched diagnostic uses only the first eight states
per replication, not the full 128. Combined GH32/M48→GH48/M64 changes are
-0.10036 ± 0.03593, -0.07212 ± 0.03726 and -0.02653 ± 0.04637 squared VIX
points. Their subset means have large outer uncertainty and must not be
substituted for the preceding full-sample price/moment results.

### Independent composition diagnostic remains inconclusive

The 6D/7D population-composition method is implemented and checked against
known Gaussian/discrete-GH moments and a linear-exponential oracle. However,
at 8x65536 its actual-model moment SE remains about 3–7 squared VIX points.
At 27 days and M48/time8 its continuous, GH32, GH48 and GH64 residuals are
-2.84, -4.04, +2.80 and -1.92. Its GH48-minus-continuous paired difference
is +5.64 ± 1.16, whereas other order comparisons differ in sign/magnitude.
Retain these results as an unresolved, noisy diagnostic; they do not certify
the population conditional quadrature error or override the deterministic
rule changes on the original states. More scrambles/variance reduction would
be needed before interpreting those paired SEs as reliable evidence of bias.

### Rare multiplier states dominate mean variance

The frozen readout is exponential quadratic with a smooth cap:

    L = 0.23332 C + 0.06940 (Y-2.93590)^2 + 0.24960 F^2,
    G = exp(2 h(L)),   h(L) <= 6.07374.

Thus G is bounded, but its upper bound is about **188619**. The cap ensures
finiteness; it does not make ordinary sampling efficient. At the representative
middle date in each VIX window:

| Days | P(G>10000) | Fraction of E[G] from G>10000 | P(G>100000) | Fraction from G>100000 |
| --- | ---: | ---: | ---: | ---: |
| 27 | 0.0556% | 48.71% | 0.0074% | 23.03% |
| 55 | 0.0657% | 51.74% | 0.0092% | 25.15% |
| 118 | 0.0711% | 52.89% | 0.0100% | 25.90% |
| 153 | 0.0717% | 52.96% | 0.0102% | 26.02% |

For example, about 0.056% of states contribute about 49% of the 27-day
middle-date mean multiplier, and about 0.0074% contribute 23%. These fractions
also describe the mean variance **above the floor** at that date. They explain
why even apparently large Gaussian samples can have poor precision.

These are future marginal G-state probabilities, not frequencies of VIX_T
option payoffs. Conditional averaging changes the distribution. Do not infer
a VIX outer-path budget by taking their reciprocals. The middle-date audit
does not establish a statement about every date in each window.

## Conclusion and next bounded step

The evidence favors **outer-sample/tail sensitivity** as the leading explanation
of the original 27-day warning at this frozen configuration. The Gaussian laws
agree, an independent accurately sampled direct expectation reproduces the
moment identity, and tested conditional/normalizer/time refinements leave the
original discrepancy essentially unchanged. This is empirical attribution,
not a proof that every quadrature bias is negligible or that an SE from eight
scrambles has reliable confidence coverage.

The concrete numerical difficulty is the concentration of expected variance in
rare multiplier states. The existing full-state proposal and small scramble
budget have not resolved it. This conclusion updates the interpretation of
the preceding warning; its original data and audit remain unchanged. No model
normalizer, parameter, market input, option price or forward is corrected here.

Next, use the independent moment reference as a diagnostic while designing
and validating an outer proposal that targets the conditional variance tails
across all four maturities. Price options/forwards on fresh independent draws,
report moment residuals alongside price errors, and assess pilot-to-pilot
variation or increase the number of complete scrambles. Do not force the
moment by rescaling prices or assume it is an exact known-mean control for an
unresolved finite-quadrature implementation. Retain quadrature and price-relative
bias gates; no new optimizer or joint-fit ranking is justified yet.

Validation: 93 package tests plus 20 driver/legacy checks pass (113 total).
Source snapshots and saved statistics are verified for all five completed
diagnostic runs. The startup failure is preserved separately. The shared
lean-numerics-bridge workflow was applied; no Lean-certified theorem or price
certificate is claimed.
