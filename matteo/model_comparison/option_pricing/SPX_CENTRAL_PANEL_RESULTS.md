# Cutoff-3 SPX representative sampling panel

## Design declared before simulation

User approved representative-panel selection and independent sampling checks.
Freeze 24 contracts from `spx_scaled_moneyness_003/training_membership.csv`:
for each 29/57/211/302-day expiry and put/call side, choose the nearest z to
0.5, 1.5 and 2.5 within the corresponding [0,1), [1,2), [2,3] band. Ties sort by
strike. No model output or held-out prices enter selection. All 707 priority
quotes and all four VIX maturities remain in scope; 24 are representatives only.

Run `m2_spx_central_panel_001`: 128-path independent pilot, two 256-path fresh
checking trials per expiry, batches of 128, seed base 20261101 plus expiry days.
Keep all Stage3 parameters, premium, zero-projection initial state and selected
nuisance curve fixed. Compare ridge8 (five baseline plus three linear OU gains)
against ridge5, both penalty 1; retain legacy OLS5 and both half-pilot fits.
Freeze each expiry's coefficients before its checking samples. No tuning or
clipping based on checking results. Record negative estimates as unusable prices.

Use direct512 pricing on a common nested event grid, including all expiries,
7/21-day anchors and curve jumps. Each expiry uses its exact prefix. Normalizer
GH48 is computed pointwise once on the longest grid in bounded chunks, then
reused by exact prefix correspondence. This is not a lower-order interpolation
or a normalizer-convergence check. Per-expiry vol7 path proposals use their own
horizon and retain the exact full-path target/proposal likelihood without mass
rescaling. Extending that proposal to long maturities is experimental.

Small budgets screen instability and estimate conditional sampling cost; they
do not establish rare-tail coverage. Report price/SE, SE/quote-scale, positive
payoff counts, pilot sensitivity, residual leverage, contract-control means,
paired estimator shifts, and warm-normalizer runtime by expiry and distance.
A prospective budget `N*(SE/(scale/8))^2` is only a plug-in planning calculation:
it assumes variance is representative and stable and cannot clear numerical bias.
No timestep refinement or market fit is run in this step.

## Results

Run `runs/m2_spx_central_panel_001` completed in **1218.55 seconds** (20.31
minutes), including **629.97 seconds** of common normalizer setup. Per-expiry
wall times, including basis construction, pilot, checking and output, were
28.5 / 55.4 / 206.2 / 298.4 seconds. No checking-based alteration was made.

The auditor reproduces **384 estimates**, **48 paired differences**, **384
contract-control means** and **48 nonduplicated bare monitors**. Frozen source
hashes, deterministic panel selection, all pilot coefficients, exact defensive
weights, mass controls and pathwise adjusted arrays replay. All **39 study-driver
tests pass**; the end-to-end synthetic audit fixture and a small nonconstant
GH48 chunk oracle pass. The first **20,168** normalizer times and values match
the previous independent short-horizon setup **exactly** (see
`normalizer_prefix_check.json`); that is correspondence evidence, not refinement.

### Prices and sampling

**No negative estimate occurs among any of the 384 entries**, including raw,
legacy, primary ridge and both half-pilot variants. This is an observed finite
sample result; linear controls still do not guarantee positivity. Every contract
has at least one positive raw payoff in each checking trial.

All raw/legacy/base/augmented methods pass **0/48** nonnegative positive-hit
`SE <= quote_scale/8` screens. Median SE/scale is 27.94 raw, 8.22 legacy,
8.15 baseline and 6.64 augmented. For the augmented estimator:

| Maturity | Inner-band median SE/price | Inner-band median SE/scale | Outer-band minimum hits / 256 | Largest outer residual-path share |
| --- | ---: | ---: | ---: | ---: |
| 29 days | 2.99% | 3.27 | 3 | 90.57% |
| 57 days | 3.68% | 4.08 | 3 | 50.47% |
| 211 days | 5.16% | 6.95 | 2 | 99.50% |
| 302 days | 4.87% | 8.22 | 1 | 99.59% |

Inner means z in [0,1), outer means [2,3]. Medians combine two contracts and two
checking trials. Residual-path share is the largest share of the sum of centered
squared adjusted contributions, not a population-variance or paired-difference
claim. The middle band also shows concentration, reaching 99.23% at 302 days.

The 302-day 10400 call is **0.7910 +/- 0.3399** with two raw hits in trial 0,
and **9.5494 +/- 9.3213** with one hit in trial 1. The latter path accounts for
99.59% of centered squared residuals. Neither is a reliable tail-price validation.
Both long-expiry outer calls have **one positive pilot payoff out of 128**;
one half-pilot has none. See `pilot_payoff_hits.csv` for all counts.

### Variance and actual cost

Linear controls reduce variance in **46/48** primary contract/trial comparisons.
The two exceptions are the 211-day 9600 call, trial 1, and 302-day 10400 call,
trial 0. Improvement is not enough for an all-band computational win:

| Maturity | Normalized variance gain, trials 0 / 1 | Pipeline cost multiplier, trials 0 / 1 | Variance-cost gain, trials 0 / 1 |
| --- | --- | --- | --- |
| 29 days | 1.247 / 1.194 | 1.406 / 1.428 | 0.887 / 0.836 |
| 57 days | 1.313 / 1.411 | 1.430 / 1.439 | 0.918 / 0.981 |
| 211 days | 1.252 / 1.424 | 1.433 / 1.441 | 0.874 / 0.989 |
| 302 days | 1.236 / 1.200 | 1.450 / 1.459 | 0.853 / 0.823 |

Variance objective averages variance/quote_scale² equally across the six
representatives of an expiry, not across all market quotes. Gain >1 is better.
Cost is the separately timed six-contract standalone512 pipeline, including
simulation and innovation recording, excluding deterministic setup, basis
construction and pilot fitting. Pooled checking costs across two trials are
9.30 / 18.14 / 67.52 / 96.88 seconds baseline, versus 13.18 / 26.02 / 97.03 /
140.88 seconds augmented. All eight full-panel variance-cost comparisons lose.

Inner-band variance gains remain stronger (1.35–2.29), especially at 29 days.
Band-specific efficiencies in `comparisons.csv` use the **six-contract** runtime;
they are not standalone inner-only timing measurements or grounds for selecting
a band-dependent production method on this checking data.

### Stability and expectation warnings retained

- Half-pilot SE inflation reaches **4.20x**, on the 302-day middle-band 9000
  call; baseline reaches 4.13x there. The 302-day inner put reaches 3.09x for
  augmented-A (2.84x baseline-A). Inner-band stability is not assured with this
  small pilot. Prior <3% sensitivity concerned different contracts/seeds and a
  512-path pilot; it was not a universal near-ATM result.
- **Two of 48** augmented-minus-base paired means exceed two empirical SEs:
  57-day 7600 call, trial 0, **+1.34748 +/- 0.50868** (2.65 SE), and 8100 call,
  **+0.38265 +/- 0.15957** (2.40 SE). Their second-trial differences do not repeat
  a >2-SE warning. Retain both; no centering bug is established by these flags.
- **7/384** correlated contract-control means exceed two SEs, including existing
  deterministic/spot deltas and one OU-linear-C gain. **3/48** bare monitors flag:
  29-day linear-C, 57-day linear-Y, and trial 1's 211-day stock mean
  **0.92749 +/- 0.03343** against 1. All mass and deterministic-reference monitors
  are within two empirical SEs. Passes do not prove adequate tail coverage.
- Independent checking-trial price differences also exceed two combined SEs in
  five primary comparisons (three augmented inner puts at 29/57/302 days and
  two baseline 57-day calls). See `independent_trial_differences.csv`. These
  descriptive, correlated warnings are neither independent significance tests
  nor proof that Gaussian conditional-mean identities fail.

### Conditional sample-budget forecasts

The plug-in relation `N_required = 256*(SE/(scale/8))²` gives median **inner-band**
augmented path counts of approximately **175k / 274k / 813k / 1.13m** by maturity.
Outer-band medians are approximately **1.88m / 1.13m / 3.35m / 2.30m**; the largest
single forecast is about **39.8m**. These huge, unstable forecasts are **not**
approved allocations. They assume representative fixed variance, frozen fits
and N^-1/2 scaling despite observed sparse/high-leverage tails. They also omit
stock timestep/normalizer bias and future pilot costs. Do not launch those
budgets or infer convergence from this arithmetic.

### Decision and next bounded work

Keep user-selected cutoff 3, all maturity/quote scopes and the baseline estimator
as reference. The linear dictionary remains experimental; positive output and
46/48 variance improvements do not justify adopting it across the band. Do not
refit the market model or launch all-707 pricing yet.

The next sampling-design work should diagnose saved **pilot** residual/payoff
contributions by side, band and proposal component, then freeze a larger-pilot
and/or better horizon-specific proposal comparison before **new** checking
samples. In particular, a seven-day OU-excursion proposal was selected on short
maturities and is not validated for long-maturity tail efficiency. Reusing this
checking set to choose a new rule would turn it into design evidence, requiring
fresh confirmation. Avoid adding features before rare-event sampling/pilot
coverage improves. Keep planned coupled256/512 and independent normalizer
refinement pending as separate accuracy gates; this run only assesses direct512
sampling. No failures, negative values or earlier warnings have been deleted.


## Frozen contract panel

| Days | Put strikes (z near 0.5 / 1.5 / 2.5) | Call strikes (z near 0.5 / 1.5 / 2.5) |
| --- | --- | --- |
| 29 | 6740 / 6440 / 6160 | 7050 / 7380 / 7725 |
| 57 | 6700 / 6290 / 5890 | 7135 / 7600 / 8100 |
| 211 | 6550 / 5750 / 5050 | 7475 / 8500 / 9600 |
| 302 | 6500 / 5550 / 4725 | 7625 / 9000 / 10400 |

Actual z values are retained in contracts.csv. Long-expiry outer calls reach
z=2.4386 (211 days) and 2.4372 (302 days), the closest available training
quotes in their bands; no unobserved strike is invented. Both sides and all
three distance bands remain represented at every maturity.
