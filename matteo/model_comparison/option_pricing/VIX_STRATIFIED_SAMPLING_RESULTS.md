# Fixed component allocations for VIX — 21 October 2026

## Predeclared scope and design

Keep 27/55/118/153-day VIX options and all four forwards: 175 training option
quotes. SPX scope is unchanged. The model, initial state, risk premium, selected
variance curve, conditional GH order 32, normalizer order 48 and time order 8
match the preceding full-scope pilot study. No fitting or proposal learning
occurs here. Reuse its frozen seven-component Gaussian mixture byte for byte.
The study changes only allocation of checking draws across these components.

Run: `runs/m2_vix_stratified_sampling_001`. Compare fresh random component
selection against fixed allocations, with eight independent complete
replications at 128 points and nested 64-point prefixes at every expiry.
Both methods have identical conditional-pricing counts and the same proposal
centers and probabilities. There are 8192 base conditional evaluations plus
12 selected-state conditional quadrature refinement evaluations. Previous
ordinary-Sobol and mixture outcomes remain historical design evidence.

## Exact estimator with integer allocations

Let f be standard Gaussian density in the full whitened OU coordinates,
phi_k(z)=f(z-m_k), and g(z)=sum_k pi_k phi_k(z), with pi_0=1/2 and the other
six probabilities 1/12. Random component selection uses mean[(f/g) H].

For fixed component counts n_k and independent scrambled draws in each
component, use

    estimate = sum_k pi_k mean_{i in k} [(f(z_i)/g(z_i)) H(z_i)].

This preserves the original expectation because sum_k pi_k phi_k=g. It does
not require n_k/N=pi_k. Saved raw data express it as a plain mean with

    effective_weight_i = (N pi_k/n_k) (f(z_i)/g(z_i)).

The allocation factor is deterministic. Neither the density ratio nor the
effective weights are divided by their observed sample sum. No market forward
is substituted or used to rescale option payoffs. All estimator weights are
positive. The exact density ratio remains bounded by 2; the allocation-adjusted
plain-mean coefficient is bounded by 32/15 (about 2.1333).

The 64-point allocation is (32,6,6,5,5,5,5); the 128-point allocation doubles
every count. This deterministic rounding rule is declared before evaluation,
without tuning counts to checking payoffs. The proposal probabilities stay
unchanged. Independent component scrambles have deterministic seeds derived
from the complete-replication seed and component index. Where a component
count is not a power of two, it uses a prefix of the next full scrambled Sobol
net; no full-net balance property is claimed for that prefix.

Refinement takes the first n_k(coarse) draws *within each component*. Taking
the global first N_coarse rows is incorrect and would discard later components.
Tests explicitly protect this distinction and verify exact allocation correction
against known Gaussian moments under strongly unequal allocations.

## Error and claim boundaries

Standard errors use eight independent complete replications, including every
component scramble. Points within a component net are not treated as iid.
Fine/coarse changes are paired across the same replication. Cross-method draws
are independent, so cross-method uncertainty is not treated as paired.

Only one previously learned pilot proposal is tested. No pilot-to-pilot
robustness, final candidate ranking, joint SPX/VIX fit, physical feasibility or
held-out performance is established. Zero hits/observed SE, payoff contribution
ESS and moment agreement are diagnostics, not acceptance certificates.
The failed simple mass control from the previous study is not used here.

Conditional 32/48 checks at low/middle/high states in the first stratified
replication at each maturity are separate diagnostics; outer SE excludes
conditional, time and normalizer biases. Existing input conflicts stay flagged.
Pricing uses the exact terminal OU Gaussian law, with no SPX stock time step.
The separate SPX time-grid and sampling work remains unresolved.

Validation before pricing: 93 package tests and 13 driver/legacy tests pass
(106 total). Frozen sources, proposal bytes, seeds and raw weighted observations
are saved. The verifier regenerates density/allocation weights, checks component
prefix counts and independently recomputes saved means, replication SEs and
paired prefix changes. Those checks establish artifact correspondence only.

After completion, the existing fixed-curve ranking audit will compare the same
original starting and selected curves using the fresh weighted observations.
It uses every fine-grid point, so its full-sample estimator remains valid with
the saved allocation-adjusted weights. It performs no fitting or selection.
Its nonlinear-loss jackknife is descriptive; full numerical gates remain open.

The random-mixture sampler uses one (dimension+1)-dimensional scrambled net,
with its first coordinate selecting a component. The fixed-allocation sampler
uses a separate dimension-dimensional scrambled net per component. Thus this
comparison changes both component allocation and the associated QMC layout.
It does not isolate count fluctuation from Sobol-coordinate/discontinuity
effects or prove the iid stratification variance formula for these nets.
Conditional-pricing counts match; random-number generation overhead can differ.

## Separate known-Gaussian diagnostic

`runs/m2_gaussian_allocation_oracle_001` tests the actual four frozen proposals
against known expectations of 1, each of the first eight coordinates, their
squares, and a call on the first standard normal with strike 3. It uses separate
seeds, eight replications, 128-point prefixes and 16,384-point fine samples.
No conditional VIX pricing is involved. At the fine budget the maximum absolute
errors across maturities are:

| Integrand family | Random mixture | Stratified |
| --- | ---: | ---: |
| Total Gaussian mass | 0.000134 | 0.000108 |
| First moments | 0.000377 | 0.000697 |
| Second moments | 0.001393 | 0.001705 |
| Strike-3 Gaussian call | 0.0000121 | 0.00000843 |

These finite checks support the sampling/weight correspondence for the actual
proposals. They do not prove accurate VIX prices or general superiority of one
method. Total likelihood-weight mass is distinct from the model's M_Q(t)
normalizer; no conclusion about its quadrature accuracy follows here.

Reproduce with a new output directory from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m matteo.model_comparison.option_pricing.vix_stratified_sampling --output /path/to/new/run
python -m matteo.model_comparison.option_pricing.verify_vix_stratified_sampling /path/to/new/run
```

## Completed results

Run wall time: 848.2 seconds (14.1 minutes). Verification passed:
35 frozen source hashes, 256 actual scramble seeds disjoint from prior pilot/mixture
streams, 8192 saved observations, 732 recomputed summary rows, regenerated
likelihood/allocation weights, component counts and paired prefixes.

| Days | Random forward ± SE | Stratified forward ± SE | Median option SE / scale: random → stratified |
| --- | ---: | ---: | ---: |
| 27 | 21.05871 ± 0.12199 | 21.01088 ± 0.12726 | 0.605 → 0.737 |
| 55 | 20.43087 ± 0.18668 | 20.61894 ± 0.09130 | 0.722 → 0.720 |
| 118 | 21.10757 ± 0.11558 | 21.17682 ± 0.20293 | 0.497 → 0.436 |
| 153 | 21.41430 ± 0.24845 | 21.06347 ± 0.12958 | 0.411 → 0.412 |

The forward SE falls at 55/153 days, is similar at 27 days and rises at 118
days. Median option SE is worse at 27 days, similar at 55/153 and somewhat
better at 118 days. Neither sampler dominates. Median per-option variance
ratios (random / stratified) are 0.67, 2.13, 1.42 and 0.89; these differ from
ratios of median SEs because each quote has a different market scale and
sampling response. Eight-replication comparisons remain noisy.

Both methods have positive observed payoffs for all 175 options at the fine
budget. Only 31 random-mixture and 26 stratified option estimates meet the
observed SE / scale <= 0.125 threshold. No full sampling gate passes, and
meeting this observed threshold does not clear tails or deterministic biases.

Stratification does not uniformly improve concentration either. At 153 days,
the largest pooled single-observation payoff share is 21.36%, versus 11.96%
for random allocation; minimum contribution ESS is 11.40 versus 14.42. These
are diagnostic concentrations, not iid equivalent sample sizes.

The 27-day stratified forward shifts by +0.53718 ± 0.14982 (paired SE) from
64 to 128 points per replication. Its small-looking fine-budget SE therefore
does not establish count convergence. Prefix changes and all other maturities
are retained in `prefix_changes.csv`.

### Curve ranking is still unresolved

`runs/m2_vix_stratified_rank_audit_001` compares the same original starting
and selected curves, using exact allocation-adjusted weights and the original
four-maturity loss, without refitting:

| Method | Selected minus starting loss | Descriptive jackknife SE | Leave-one-scramble-out range |
| --- | ---: | ---: | ---: |
| Random mixture | -0.21042 | 1.44738 | -1.23710 to +0.45901 |
| Stratified | +0.93316 | 1.11633 | +0.41400 to +1.82213 |

Negative favors the selected curve. These fresh comparisons do not establish
a reliable gain and differ from the preceding run, where both methods had
positive changes around +2.7 to +2.9. Preserve both records: the ranking is
sensitive to fresh draws and does not support a final fitted-model judgment.

### Conditional quadrature: retain the larger tail discrepancy

At the first stratified replication, selected-state conditional GH 32→48
maximum changes are 0.000540, **0.022312**, 0.001601 and 0.002325 VIX points.
The 55-day maximum occurs at a tail state near VIX 262. This is larger than
the preceding selected-state checks (around 0.0034 or less), so those older
checks must not be generalized to all new proposal states. A pointwise change
is not a full-price bias bound. No new global-normalizer or time-order gate
is cleared.

### New moment-identity warning — recorded without correction

For the exact normalized model, the tower property gives

    E[VIX_T^2] = (10000/delta) integral_T^(T+delta) xi(t) dt.

The target is 542.49611 at 27/55 days and 600.31899 at 118/153 days. At 27
days the raw estimates are 556.16551 ± 3.23481 (random) and 554.03172 ± 3.10259
(stratified). Their departures are several reported SEs. This is a numerical
identity discrepancy to investigate, not proof of a normalization bug. Only
eight scrambles, multiple inspected diagnostics and rare tails limit inference.

`runs/m2_vix_stratified_moment_audit_001` additionally removes the random
constant-mass contribution by measuring w*(VIX^2-target), using the exact
identity E_g[w]=1. The 27-day residuals remain +15.66458 ± 2.75395 (random)
and +12.60904 ± 3.46012 (stratified). This diagnostic neither self-normalizes
weights nor forces model prices or moments to match their targets. Its design
was prompted by the observed discrepancy and is explicitly exploratory.

The prior study found noisy second-moment agreement; that historical result
remains preserved. The new warning requires renewed checks. Possible sources
include outer-tail sampling/underestimated SE, conditional-GH composition,
date-zero normalizer quadrature and time integration. The separate known-
Gaussian mass checks validate none of these model-normalization approximations.

## Decision and next work

Keep all four VIX maturities and forwards and retain both sampler references.
Do not adopt fixed component allocation as a universal improvement. The method
is implemented and useful for controlled comparisons, but this experiment does
not solve the original forward/tail trade-off or establish convergence.

Prioritize the normalized-model moment discrepancy before further optimizer
runs: freeze matching states and integration nodes, vary conditional GH and
global-normalizer orders separately, and compare the moment residual across
independent outer scrambles and larger nested counts. Use 27 days as the first
diagnostic focus while retaining all four maturities in the full analysis scope.
Record computational budgets against quote scales before launching larger
runs. The separate SPX sampling and time-grid work remains open. No joint-fit,
physical-feasibility or held-out market claim has been made.
