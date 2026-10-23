# Fixed upward-equity tilt and larger pilots: cutoff-3 SPX panel

## Pilot-only diagnosis and design

Continuation authorized after the 24-contract central-panel sampling screen.
Use only `m2_spx_central_panel_001/expiry_*/pilot_draws.npz`, frozen coefficients
and contracts for design. Do not use its checking draws to tune a proposal.
`spx_sampling_design_001` records per-component pilot contributions and
cross-fold off-policy second-moment proxies at positive equity shifts
2/2.5/3/3.5/4. All outer-call pilot hits came from the positive-equity component;
each long outer call had one hit, so fold information is incomplete.

The stronger tilt is **not** a robust aggregate pilot optimum: at +3, the
quote-scale objective ratios in the two folds are 1.004/1.013 (29d),
0.954/1.006 (57d), 1.016/0.998 (211d), and 1.003/1.075 (302d).
Increasing to +3.5/+4 often worsens these proxies further. Select a **fixed +3**
as a moderate, targeted call-coverage hypothesis; no all-band win is assumed.

Math: with mass control `w-1`, the old estimator is `Z_old=w_old*h+beta_mass`.
For fixed beta, its nonconstant second moment under a new proposal can be
estimated from old proposal samples as
`mean((Z_old-beta_mass)^2 * w_new/w_old)`. It is incorrect to omit the constant
or invert the density ratio. Cross-fold coefficients use one pilot half and
score the other; the result is still noisy pilot evidence, not fresh validation.

## Frozen fresh experiment

`m2_spx_call_tilt_001`: same 24 contracts and four SPX maturities, 512 steps/day,
Stage3 parameters/premium/initial state/nuisance curve unchanged. Compare the
original vol7 mixture (+2 equity) against +3 equity. Negative-equity and all
volatility components, probabilities and orthonormal path basis remain fixed.
Half the mixture is target Gaussian; exact full-path weights remain bounded by 2.
The +3 is a shift of the normalized Brownian path score W_H/sqrt(H), distinct
from the independently selected market-moneyness cutoff 3.
Same seeds pair component labels and underlying normals, not identical paths.

Use **256 new pilot paths per proposal per expiry** (twice the old pilot),
two fresh **256-path** checking trials, batches of 128, seed 20261122. Each proposal
gets ridge5 penalty 1 full and half-pilot fits, frozen before checking. Retain
raw prices and the identical original128-path beta under **both** proposals.
This fixed-beta comparison isolates sampling law; new256 vs old128 under the
original sampler isolates pilot learning on identical checking paths. No new
innovation features, clipping, model fitting, timestep or quadrature refinement.

Reuse GH48 normalizers only after exact master-grid, parameters, premium,
initial state, target hashes, curve and all package-source correspondence checks.
Freeze cached bytes and provenance. This avoids repeating the 10-minute setup;
it does not make the approximate quadrature exact. All 707 cutoff-3 quotes,
secondary wings and all four VIX maturities remain in scope.

## Results

Run `m2_spx_call_tilt_001` completed in **618.82 seconds**, including **1.00
second** of cache validation/setup. All **480 estimates** replay, including
fresh independent pilots, exact likelihood/mass controls, paired component
labels and Gaussian-score shifts. All **42 study-driver tests pass**, including
Gaussian/Black target-law and off-policy mass-constant oracles. Cached-normalizer
hashes match the independently recorded previous provenance.

### Coverage improves, but not reliably variance

Outer-call counts below are out of 256 paths; pilot and checking samples are
separate. The same component assignments/underlying normals pair the proposals.

| Days | New pilot hits: original -> call3 | Checking trial 0 | Checking trial 1 |
| --- | --- | --- | --- |
| 29 | 6 -> 13 | 1 -> 10 | 7 -> 18 |
| 57 | 3 -> 10 | 6 -> 12 | 3 -> 17 |
| 211 | 2 -> 7 | 5 -> 9 | 1 -> 4 |
| 302 | 1 -> 9 | 1 -> 7 | 3 -> 8 |

Thus all eight checking outer-call hit counts increase. However, the primary
own-pilot controlled variance improves in only **10/48** contract/trial cases;
with identical old128 coefficients, only **9/48** improve. All put variances
increase slightly with the fixed coefficients (median gains 0.973 / 0.983 /
0.985 by distance band). Increasing payoff frequency does not guarantee that
the weighted residual variance falls; proposal density and fitted controls matter.
Low reported variance in a one-hit original sample is also unreliable evidence
of an accurate tail price.

### Aggregate cost and variance tradeoff

The objective is equal mean variance/quote_scale² over the six representatives
of each expiry. The table reports original variance / call3 variance. Efficiency
divides that ratio by the measured call3/original pipeline time.

| Days | Own256-pilot variance gain, trials 0 /1 | Fixed old128 variance gain, trials 0 /1 | Own256 variance-cost gain, trials 0 /1 |
| --- | --- | --- | --- |
| 29 | 0.963 / 0.969 | 0.989 / 1.012 | 0.945 / 1.036 |
| 57 | 0.917 / 1.168 | 0.896 / 0.951 | 0.917 / 1.168 |
| 211 | 1.410 / 0.968 | 1.166 / 0.970 | 1.410 / 0.964 |
| 302 | 0.775 / 1.008 | 0.984 / 1.061 | 0.817 / 1.021 |

Only **4/8** primary aggregate variance-cost comparisons win. Cost is roughly
unchanged (paired trial multipliers 0.935–1.019); small timing differences are
not a structural speedup. Checking pipeline totals at 29/57/211/302 days are
11.59/18.08/71.31/106.81 seconds original versus
11.38/18.09/71.44/103.19 call3. These exclude pilot/basis/setup costs; all band
comparisons reuse the six-contract timing and are not standalone band timings.
No robust broad-panel efficiency win is established.

### Learning and precision remain unresolved

- **No negative estimate among all 480**, including both half-pilot variants.
  Linear controls still do not guarantee positivity in future samples.
- Original and call3 each pass **0/48** primary quote-scale sampling screens.
  Median SE/scale is **7.665 original**, **7.482 call3**, against target 0.125.
  Fixed-old128 medians are 8.171/8.174. No zero observed raw checking payoff
  case occurs here, but several counts remain very small.
- On identical original checking samples, new256 vs old128 coefficients give a
  median variance benefit **1.073**; under call3 it is **1.024**. These are not
  pure causal sample-size effects: pilots use different data and fixed ridge 1
  has different regularization per observation at different N.
- Larger-pilot learning is not uniformly stable. At 302 days, trial 0, call3,
  strike10400: new256 gives **2.99674 +/- 2.25281**, old128 gives
  **0.44628 +/- 0.31909** on the **same checking paths**. New-fit variance is about
  **49.85x** higher. Do not infer that the smaller pilot is generally superior
  from this sparse case; it demonstrates fitted-control sensitivity.
- Half-pilot/full-pilot SE ratios reach 1.77 original and 1.68 call3, lower than
  the previous panel's worst 4.20 but on different pilots/checks; no universal
  stability improvement is proved. Nine within-sampler new-vs-old learning
  differences exceed two empirical paired SEs (retained in learning_differences.csv).
- The 302-day outer call under original sampling has one hit in trial 0:
  **0.07528 +/- 0.07520**, with **99.61%** of centered squared residuals on one
  path. Call3 has seven hits, but **2.99674 +/- 2.25281** and **78.99%** maximum
  residual share. More hits improve coverage without yet establishing accuracy.

### Mean and likelihood diagnostics

One of 48 primary proposal-paired mean comparisons exceeds two SEs: 57-day
7135 call, trial 0, call3-minus-original **+2.58537 +/- 1.04038** (2.49 SE).
Two of 48 fixed-beta comparisons flag: 29-day 6740 put (trial 1, -2.31 SE) and
211-day 5750 put (trial 1, +2.34 SE). All are retained without clipping/recentering.

There are **30/480 correlated contract-control mean flags**, of which 18 are
repeated mass-control flags shared across six contracts; the others are
reference payoff and deterministic/spot delta controls. **18/176** correlated
bare forward/Gaussian-moment checks exceed two SEs. These include weights:
29-day trial 1 call3 mass **0.92137 +/- 0.03594**, and 211-day trial 1 original
mass **1.07224 +/- 0.02992**. Never divide weights by these sample means to make
them pass. Exact Gaussian/Black oracles and paired-score identities pass, but
those implementation checks do not make small-sample means exact or certify
rare-event convergence. Flags are descriptive, not independent tests.

## Decision and next direction

Keep **original vol7/ridge5 as reference**, cutoff 3/all 707 priority quotes and
all four VIX maturities. Keep call3 as an experimental coverage option; do not
adopt it globally, refit the market model, or launch multi-million-path budgets
on the strength of this screen. Both stronger tilt and larger pilots have been
tested, including an unchanged-coefficient comparison, without a robust general
sampling solution. All prior warnings remain historical evidence.

Before another round of mixture magnitudes/features, investigate a structural
conditional-Monte-Carlo route for SPX: whether the Gaussian equity noise that
remains after conditioning on the full exogenous OU state path can be integrated
analytically in the current discrete stock scheme. First establish the actual
joint-Gaussian decomposition and source correspondence, including numerical
rank/covariance issues, then Gaussian/constant-volatility and same-model checks.
A nonnegative conditional payoff average is a candidate, not implemented or
validated here; do not reuse full-path importance weights after marginalizing
noise without deriving the correct marginal likelihood. Stock timestep and
normalizer refinement remain separate pending gates.
