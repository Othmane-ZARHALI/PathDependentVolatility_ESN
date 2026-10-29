# Conditional existing SPX controls — 23 October 2026

## Question and fixed scope

Can analytical integration of independent equity noise improve the **controlled**
estimator enough to justify replacing original vol7/ridge5? The preceding
[raw conditional experiment](SPX_CONDITIONAL_PRICING_RESULTS.md) did not answer
this question. This is the bounded experiment in the research ledger; it does
not start an open-ended search over controls.

Run: `runs/m2_spx_conditional_controls_001`. Frozen Stage-3 candidate 18,
scenario 0, zero projected state, selected nuisance curve, GH48 normalization,
512 stock steps per trading day. The unchanged 24 market-selected contracts
cover put/call sides and three distance bands at 29/57/211/302 days. They are
representatives of **707/980 cutoff-3 priority training quotes**, not the entire
priority population. All **27/55/118/153-day VIX maturities, 175 training options
and four forwards remain active**. No held-out market observation is opened.

## Conditional expectations of the existing controls

Write normalized stock as `Y_k=exp(L_k)`. Given all retained Gaussian shocks U,
the OU path and start-of-step variance V_k are known. The eliminated independent
equity normals are still independent across steps and independent of U. Thus

`L_k | U ~ N(m_k,v_k)`,

`d_k=-V_k*dt_k/2+sqrt(V_k)*a_k`, `q_k=nu²*V_k*dt_k`,

`m_next=m_k+d_k`, `v_next=v_k+q_k`.

The frozen independent fraction is `nu²=0.0063668537395`. This is a small
Brownian variance fraction, not an upper bound on variance-reduction benefit.
The conditional target payoff is Black/lognormal with log mean m_T and log
variance v_T; conditional stock is `exp(m_T+v_T/2)`.

For the deterministic-volatility reference stock, replace V_k by xi_k in
the same two recurrences. Its conditional reference payoff is another Black
expectation. Subtract its **unconditional** Black price with total variance
`sum xi_k*dt_k`. The stock control similarly subtracts its unconditional mean 1.
Neither conditional stock nor conditional reference payoff is forced to its
unconditional mean path by path.

The existing two adapted deltas use remaining variance W_k, either the
deterministic curve integral or that integral multiplied by `V_k/xi_k`.
Given U, both choices of W_k are known positive numbers. Gaussian integration
gives, for a call,

`E[Y_k*delta(Y_k;K,W_k)|U]`

` = exp(m_k+v_k/2)*Phi((m_k+v_k-log(K)+W_k/2)/sqrt(W_k+v_k))`.

For a put replace Phi(d) by `-Phi(-d)`. This follows by multiplying a Gaussian
density by exp(L), which shifts its mean by v, and then integrating the normal
CDF. It is **not** delta evaluated at the expected stock.

The conditional hedge gain for each step is

`E[delta_k*(Y_next-Y_k)|U]`

` = expm1(d_k+q_k/2)*E[Y_k*delta_k|U]`.

Sum this over the original hedge steps. It is generally **not zero** given U.
Its unconditional expectation is zero by the original adapted discrete
martingale argument. Dropping it after conditioning would change the estimator.

In undergraduate terms, averaging out one source of randomness means averaging
the payoff **and every subtracted hedge** over that same randomness. A hedge
that gains nothing on average across all scenarios can still have a gain or
loss when the volatility-driving scenario has been fixed.

The four conditional controls are therefore the conditional versions of stock,
reference payoff and two hedge gains. The fifth control below is importance
mass. Signed control-adjusted estimates can be negative, although the
unadjusted conditional payoff is nonnegative. No negative estimate is clipped.

These are mathematical derivations for the discrete exponential-Euler stock
scheme under the implemented Gaussian law. No Lean formalization, rigorous
floating-point payoff bound or continuous-time convergence theorem is claimed.

## Retained-shock proposal and likelihood

For each step, use the previous covariance-audited retained root R_k. Draw
standard retained coordinates z_k and set `(a_k,eps_X,k)=R_k*z_k`. The integrated
independent equity shock does not occur in this proposal.

Construct four orthonormal directions B over **all retained-normal coordinates**:
retained equity W at the horizon, then orthogonalized projected C at the horizon,
Y at day 7 and F at day 7. Gram-Schmidt retains all state modes in the sampler;
it only defines the four tilted score directions. With `s=B' z`, sample the
defensive mixture `q(z)=sum pi_j*phi(z-B*mu_j)`. Its exact target/proposal ratio is

`w(z)=1 / sum pi_j*exp(mu_j' s - ||mu_j||²/2)`.

There is an unshifted component of mass 1/2; equity shifts +/-2 have mass 1/8
each, and shifts +/-2 on each other axis have mass 1/24 each. Hence `w<=2`.
Multiply the conditional target and all four conditional controls by w;
append `w-1` as the fifth known-zero control. Never self-normalize the weights.

This is a **new retained-coordinate proposal analogous to vol7**. It is not
the marginal density of the original full-path mixture. The comparator uses
the unchanged original full-path vol7 proposal. Thus measured differences
compare whole estimators: conditioning and proposal geometry both differ.
No old full-path likelihood multiplies an already-integrated payoff, and no
pure Rao-Blackwell efficiency claim is inferred across these two proposals.
For fixed coefficients under the same retained proposal, integrating all
controls does satisfy the ordinary conditional-variance inequality; the
independently fitted comparison need not preserve that inequality.

## Implementation checks and frozen numerical design

`esn_eq/conditional_controls.py` streams states and analytical controls, without
storing the full path. The legacy sampler and pricing code are unchanged.
Seven new tests cover independent one- and two-dimensional Gaussian integration,
nonzero conditional hedge gains, piecewise-constant independent Black oracles,
nonlinear correlated raw/conditional mean correspondence, target forward/control
means, retained-proposal mass and Gaussian-score means/second moments, weighted
Black prices, perfect-correlation limits and invalid inputs.

Each method/expiry uses **256 independently sampled pilot paths**, followed by
**two fresh 256-path checks**, in batches of 128. Both pilots freeze before
checking at an expiry. Standardized ridge penalty is fixed at 1. Full-pilot and
two disjoint half-pilot coefficients are all retained; none is selected using
checking results. Seeds include method, expiry, role, trial and batch.

Timing measures the standalone six-contract pipelines on warmed common GH48
normalizers. Basis construction and pilots are separately retained. Ordering
alternates by checking trial. Primary efficiency is the ratio of
`mean_contract(variance/scale²) * measured checking runtime`, by expiry/trial;
gain above one favors the conditional method. The quote scales and contract
membership are frozen. Cross-method price changes use **independent combined
SEs**; half-pilot changes within a method are paired.

The declared adoption condition is improvement in both fresh trials across the
whole representative panel, without worsened tail/pilot risk. An adoption would
still leave stock timestep, normalizer and all-price sampling gates open.
`design.json` and all executable source snapshots were written before pricing.
The run refuses an existing output directory and checks for mid-run edits.

## Results

The run completes in **904.99 seconds**. All **288 prices** (96 full-pilot
primary estimates and 192 half-pilot estimates), 52 captured source/input
hashes, 48 distinct pilot/checking batch seeds, weights, fitted coefficients,
known-mean diagnostics and adjusted arrays replay. The maximum standardized
transition covariance discrepancy against the legacy sampler is **6.37e-11**
(maximum absolute error **2.08e-16**). These finite-precision covariance checks
do not bound nonlinear payoff error or prove bitwise equivalence of samplers.

**124 package tests + 42 existing study-driver tests = 166 passing tests.**
The seven additions are the conditional-control oracles above. Test success
and saved-array replay are implementation evidence, not price convergence.

### Variance reduction does not pay for its cost

Individual primary variance improves in **34/48** contract/trial comparisons.
The following ratios aggregate the six representatives at an expiry using
their frozen quote scales. Efficiency above one favors the conditional method.

| SPX days | Trial | Reference / conditional scaled variance | Reference seconds | Conditional seconds | Variance-cost efficiency gain |
| --- | --- | ---: | ---: | ---: | ---: |
| 29 | 0 | 1.245 | 5.09 | 11.59 | 0.547 |
| 29 | 1 | 1.027 | 4.56 | 9.74 | 0.482 |
| 57 | 0 | 2.078 | 8.95 | 18.96 | 0.981 |
| 57 | 1 | 1.191 | 8.96 | 18.98 | 0.562 |
| 211 | 0 | 1.647 | 34.66 | 72.19 | 0.791 |
| 211 | 1 | 3.234 | 33.38 | 70.85 | 1.524 |
| 302 | 0 | 1.307 | 48.32 | 100.99 | 0.625 |
| 302 | 1 | 0.805 | 48.10 | 100.52 | 0.385 |

Only **1/8** expiry/trial whole-panel comparisons wins on cost. Inner-band
efficiency loses in **all 8/8** comparisons (range **0.430–0.888**).
Across all four equally weighted expiry panels, variance-cost gains are
**0.7483** and **0.8541** in the two fresh trials. The conditional pipeline
takes about 2.1 times the checking runtime. Pilot costs are also retained:
95.09 seconds reference versus 209.64 seconds conditional, with per-expiry
basis costs in `run.json`. No production speedup is established.

### Tiny conditional tails can pass a nominal SE screen misleadingly

Median primary SE/quote-scale is **7.6167** reference and **6.5122** conditional,
far above the required **0.125**. Reference has **0/48 nominal SE passes**;
conditional has **2/48**, both the **302-day 10400 call**. Neither conditional
pass clears tail accuracy:

| Quantity at that call | Observation |
| --- | --- |
| Conditional pilot maximum weighted payoff | **3.08e-27 SPX points** across 256 paths, despite 201 positive values |
| Largest full-pilot conditional coefficient magnitude | **8.24e-31**; essentially no useful tail regression was learned |
| Fresh conditional estimate, trial 0 | **0.0219090 +/- 0.0218790**; 205 positive conditional values |
| Fresh conditional estimate, trial 1 | **7.13e-11 +/- 6.88e-11**; 198 positive conditional values |
| Largest centered-residual-square share | **99.609% / 99.470%** |
| Fresh reference estimates | **1.63269 +/- 1.35970** and **0.22755 +/- 0.28085**, with only 3/2 actual payoff hits |

The two fresh reference/conditional differences do not exceed two independent
combined SEs, but neither estimator resolves this tail. The same is true of all
48 cross-method mean comparisons: absence of a warning does not certify missing
tails. Counting tiny positive conditional probabilities as independent raw
exercise hits would be wrong. The preliminary audit used a misleading
`sampling_screen` label for that positive-count proxy. Its original JSON/CSV
are retained as `audit_preliminary_positive_count.json` and
`screens_preliminary_positive_count.csv`; the corrected audit reports
`nominal_se_screen`, distinguishes raw-hit counts, and gives **no numerical
quote clearance**. No pricing simulation was changed or reselected.

In plain language, the pilot saw many astronomically small conditional prices,
but none of the volatility scenarios that matter for that call. Integrating a
small independent equity noise does not supply those missing volatility paths.
A small observed standard error can therefore accompany a severely unresolved
expected price.

### Signed prices and pilot/control warnings retained

Reference trial-0 **211-day 9600 call** is **-0.000516 +/- 0.001140** with the
full pilot, and **-0.001531 +/- 0.002699** with half A. These two negatives are
saved without clipping. This new observation does not change the historical
central-panel run's positive outputs; positivity was specific to those draws.
All 144 conditional full/half-pilot estimates are nonnegative in this run,
which is an empirical observation for the controlled estimator.

Maximum half-pilot SE inflation is **2.415x conditional** (57-day outer put)
and **2.368x reference** (211-day outer call). Maximum primary residual-square
concentration is **99.61% conditional** versus **98.02% reference**. Conditional
positive counts are not comparable to reference exercise-hit counts. There are
**11 conditional and 14 reference contract-control mean warnings** beyond two
empirical SEs; none of the 48 bare mass/stock/reference monitors flags. These
are multiple correlated diagnostics, not calibrated independent tests or proof
of a centering bug. The exact identities and empirical uncertainty remain
separate.

## Disposition and return to fitting

**Do not adopt this implementation as the production replacement. Retain
original vol7/ridge5.** The bounded condition fails on whole-panel efficiency
in both trials, and conditional pilot/tail risk is not removed. This rejects
the tested implementation/proposal/budget as a general replacement; it does
not prove that every conditional-control construction must lose.

Keep the derivation, streaming implementation, failed cost comparison, tiny-tail
cases and signed reference failures as evidence. Do not open a new feature or
proposal search by default. Sampling precision, stock timestep, normalization,
VIX conditional quadrature and VIX time integration remain separate gates.

The concrete decision is in
[FITTING_REENTRY_DECISION.md](FITTING_REENTRY_DECISION.md): **no qualified fit or
ranking is ready**; a bounded curve-only joint local profile is permissible as
explicitly unresolved exploration. It specifies the two free coefficients,
24-SPX/all-VIX provisional objective, frozen scales, nine candidate points,
budgets, fresh checking and full-scope limitations. That market profile has
not been launched here. No search loss is presented as a validated fit gain.

## Reproduction

From repository root, with the existing pricing environment:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
PYTHONPATH=matteo/exogenous_reservoir_eq_m2/src \
/private/tmp/pdv-options-venv/bin/python \
  -m matteo.model_comparison.option_pricing.spx_conditional_control_study \
  --output /private/tmp/NEW_CONDITIONAL_CONTROL_RUN

OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
PYTHONPATH=matteo/exogenous_reservoir_eq_m2/src \
/private/tmp/pdv-options-venv/bin/python \
  -m matteo.model_comparison.option_pricing.audit_spx_conditional_controls \
  --run /private/tmp/NEW_CONDITIONAL_CONTROL_RUN
```

`validation.json` in the recorded run retains test commands/counts, validation
source hashes and the original negative-result/screen-label corrections.
