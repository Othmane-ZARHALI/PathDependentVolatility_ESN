# SPX path sampling, timestep and normalization pilot — 22 October 2026

The user confirmed SPX numerical diagnostics as the current focus. The model,
initial state and selected nuisance variance curve remain frozen; no market
optimization or held-out evaluation occurs. All four VIX maturities, 175 training
options and four forwards remain active. The VIX tail-control extension stays
pending as recorded in [NEXT_STEPS.md](NEXT_STEPS.md).

## Sampling law and correspondence

Let H=57/365 and Z=W_H/sqrt(H). Independently for each path choose a component
m in {0,+2,-2} with probabilities {1/2,1/4,1/4}. Simulate the complete equity
Brownian path with constant drift m/sqrt(H). The corresponding correlated OU
innovations shift too: on an interval dt, their mean shift is
Cov(OU innovation, dW) * m/sqrt(H). The conditional Gaussian OU residuals and
the model's Q drift remain unchanged. This changes the proposal, not the model.

The exact full-path target/proposal weight is

w(Z) = 1 / [1/2 + (1/4) exp(2Z-2) + (1/4) exp(-2Z-2)].

It is positive and at most two. Conditional residual densities cancel in the
likelihood ratio. The same full-horizon weight applies to earlier expiry
payoffs and coupled coarsenings, though it can be inefficient for earlier
expiries. Weighted estimates use ordinary means of w*payoff, without dividing
by the observed mean weight. These are iid path draws, not Sobol draws.

The four previous target-zero controls (stock-minus-one, reference Black payoff
minus its exact expectation, two adapted stock trading gains) become w*control.
An additional control w-1 has known mean zero. Separate pilot paths fit the
coefficients; coefficients remain fixed across checking paths and refinements.
The target sampler uses the same construction with w=1 (a zero mass column).
Negative controlled estimates and zero-payoff counts are retained.

## Experiment

[spx_path_importance.py](spx_path_importance.py) compares independent target and
importance samples on 12 predeclared actual training contracts at 29/57 days.
Each of two trials uses 512 pilot and 1,024 checking paths per method, in
batches of 256. The predeclared shift is two; checking data do not select it.
The representative contracts do not replace the full SPX calibration scope.

Each checking batch shares fine OU/equity paths across coupled 64/128/256 steps
per trading day at normalizer GH48. A separate same-path GH48-to-GH64 comparison
uses the 256-step grid. Thus stock time integration and normalizer refinement
are measured separately, conditional on the same frozen model/curve/state.
All stock prices use the adapted exponential-Euler scheme. Exact OU transitions
do not imply exact continuous-time stock integration.

Sampling uncertainty is the standard error across iid checking paths,
conditional on each independent pilot. Refinement uncertainty uses paired path
differences. The two methods use independent samples; their price differences
are not paired. Small samples, observed zero SEs and two trials do not establish
population error bounds or stable variance ratios. Many simultaneous diagnostic
screens are descriptive, without familywise statistical guarantees.

## Validation and reproducibility

101 package tests and 28 comparison-driver tests pass. New independent oracles
check Gaussian mass/first/second moments, full OU means and equity cross-moments,
Black calls and puts including rare payoffs, the weighted stock/reference means,
zero-shift backward correspondence, and likelihood validity after coarsening
and for earlier expiry payoffs. A quadratic oracle checks independently trained
weighted controls against its exact target expectation. These are numerical
and code-correspondence checks; no Lean theorem or pricing certificate is claimed.

Run (single-thread BLAS, existing Python environment):

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.spx_path_importance --output matteo/model_comparison/option_pricing/runs/m2_spx_path_importance_001 --paths 1024 --pilot 512 --batch-size 256 --trials 2
/private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.verify_spx_path_importance matteo/model_comparison/option_pricing/runs/m2_spx_path_importance_001
```

The runner snapshots sources and the curve-fit input before computation,
records model/target hashes and seeds, rejects mid-run source edits and saves
all checking price/control draws needed to replay estimates and paired changes.
An existing output directory is never overwritten.

## Results and retained warnings

The completed [run](runs/m2_spx_path_importance_001/run.json) took 347.8 seconds.
Grids contain 2,520/5,040/10,080 intervals; max(lambda*dt) is approximately
0.30988/0.15494/0.07747. The [audit](runs/m2_spx_path_importance_001/audit.json)
and raw CSVs retain every method/trial/contract. The artifact verifier reproduces
**384 estimates, 288 paired changes and 20 mass/stock/reference monitors** from
saved checking draws and validates captured source hashes and seed separation.

At the 256-step/GH48 primary checkpoint (12 contracts x two trials):

| Diagnostic | Target sampler + controls | Importance sampler + weighted controls |
|---|---:|---:|
| Positive-hit SE screens passing SE <= scale/8 | 0/24 | 0/24 |
| Zero raw-payoff contract/trial cases | 2/24 | 0/24 |
| Negative controlled price estimates | 2/24 | 0/24 |
| Median observed controlled SE / quote scale | 4.7033 | 2.3512 |
| Median observed raw SE / quote scale | 10.2845 | 7.0107 |

The paired-by-contract descriptive ratio of target to importance **controlled
sample variances** has median **1.2944**, range **0.2538 to 71.6783**; it favors
importance sampling in **13/24** cases. These compare independent samples and
independently trained pilots. The two separate median SE/scale values above do
not imply a fourfold median variance gain. No method was selected by checking
results, and these small-sample ratios are not population guarantees.

For example, the 29-day 3400 put has target payoff-hit counts 0/1 across the two
trials, versus 5/7 under importance sampling. Its controlled importance price
estimates are 0.26581 ± 0.12052 and 0.41420 ± 0.17160 (mean ± one observed SE).
The quote scale is 0.10434, so even this improvement is far from the scale/8
sampling screen. At 57 days the 8800 call records 2/3 target hits versus 14/13
importance hits, yet its controlled SE is worse under importance sampling in
one trial. More payoff hits are useful diagnostic evidence, not sufficient
price precision.

**Preserved negative estimates:** target-sampler trial 1 gives -0.48866 ± 1.14178
for the 29-day 7900 call and -0.98739 ± 1.47133 for the 7950 call. These are noisy
control-adjusted estimators, not negative payoff values or an arbitrage finding.
They remain in the record without clipping. Both raw zero-payoff cases occur
in target trial 0 (29-day 3400 put and 7950 call); their adjusted estimators
have nonzero values and SEs. Neither raw zero SE nor a nonzero controlled
estimate establishes an accurately priced rare option.

Mass estimates under importance sampling are 0.99248 ± 0.01934 and
1.01791 ± 0.01945. All weighted stock/reference monitor discrepancies are below
two observed SEs. This finite diagnostic is consistent with correct weighting;
it does not demonstrate resolution of rare stock/volatility tails.

| Paired controlled refinement | Target: median / max absolute change in quote scales | Importance: median / max |
|---|---:|---:|
| 64 -> 128 steps/day | 0.87910 / 8.13025 | 0.57354 / 2.92499 |
| 128 -> 256 steps/day | 0.75054 / 6.46213 | 0.89337 / 3.55711 |
| GH48 -> GH64 at 256 steps/day | 0.000482 / 0.001263 | 0.000336 / 0.000959 |

**Time integration is unresolved:** observed 128-to-256 changes exceed scale/8
in 20/24 target and 21/24 importance cases, but most paired changes are noisy;
only one case per method exceeds two paired SEs. For the 57-day 5550 put,
trial 1 changes are +1.29399 ± 0.92166 under target sampling and
-0.71228 ± 0.47191 under importance sampling, against scale 0.20024. The
opposite noisy signs do not establish opposite population biases or an
integration order. Do not treat 256 steps/day as convergence clearance.

**Normalization looks smaller in this bounded check:** all observed GH48-to-64
changes are below 0.0013 quote scales, versus a sampling screen of 0.125
scales. Some are statistically distinguishable from zero but economically
small relative to the frozen tolerance. This is same-path, sample-conditional
sensitivity, not a population normalizer error bound or a global GH certificate.

## What to advance next

This first equity-only path mixture helps expose tails but is insufficient.
Keep it as a comparator. On separate design pilots, investigate path shifts
that also target volatility-factor excursions, and check pilot coefficient
stability before fresh evaluation. Better tail hits alone are not the objective:
the criterion remains weighted payoff/paired-correction variance at measured
cost, with valid mass and reference expectations.

For time integration, retain the coupled 64/128/256 corrections and measure
variance *and cost* by level before choosing an MLMC allocation or a further
512-step check. This experiment records correction variances but has not built
an MLMC estimator or established a decay rate. Sampling and grid errors remain
open; the small normalizer sensitivity is a useful prioritization result, not
permission to refit. All VIX maturities remain active, with its pending
extension unchanged.
