# SPX volatility-excursion path sampling — 22 October 2026

The user asked to continue improving SPX numerical pricing. This experiment
extends [the equity-only path pilot](SPX_PATH_IMPORTANCE_RESULTS.md), with the
same frozen Stage-3 model, zero projected initial state and selected nuisance
curve. It changes sampling and pilot controls only. All four VIX maturities
(27/55/118/153 days), 175 training options and four forwards remain active;
VIX tail controls remain pending. Held-out quotes are not used.

## Exact path proposal

Each joint OU/equity transition is generated from independent standard normals
u_j: one equity normal and the OU innovation normals conditional on that equity
increment. A deterministic matrix B_j defines four orthonormal full-path scores:

Z = sum_j B_j' u_j,    sum_j B_j' B_j = I_4.

The first score is W_H/sqrt(H), with H=57/365. The other three are orthonormalized
linear combinations of conditional OU residual innovations: cluster projection
C at H, feedback projection Y at an excursion anchor, and spike projection F
at that anchor. The exact OU semigroup maps each innovation to the endpoint
projection. Gram-Schmidt normalizes and decorrelates these scores. No state or
path dimension is discarded: all remaining innovation coordinates are drawn.

Given mixture component m, normals become u_j + B_j*m. Thus Z has component
mean m and identity covariance, and the exact full-path likelihood is

w(Z) = 1 / sum_k pi_k exp(m_k' Z - ||m_k||²/2).

The original Gaussian law retains probability 1/2, so 0 < w <= 2. The weights
are never divided by their sample mean. Both stock and OU innovations shift
consistently in the same transition representation; the model's Q drift,
readout and date-zero normalizer remain fixed. These are iid path draws.

Three proposals are predeclared:

- **equity:** half original Gaussian, quarter each +/-2 equity-score shifts;
  this reproduces the previous proposal up to floating-point arithmetic.
- **vol7:** half original, 1/8 each equity +/-2, and 1/24 each +/-2 shift in
  the three OU score coordinates; feedback/spike excursion anchor 7 days.
- **vol21:** identical mixture probabilities/magnitudes, 21-day excursion anchor.

The OU cluster endpoint remains H for both candidates. These are targeted
Gaussian directions, not a claim that every volatility excursion or nonlinear
rare event is covered. Their path weights remain valid for earlier expiries
and coupled coarsenings of the finest cache.

## Design versus checking

Each candidate has **512 independent training paths** and **512 separate
selection paths**, all on the finest 256-step/day grid. Weighted stock,
reference-payoff and two adapted hedging controls, plus the mass control w-1,
are fit on the training sample. Three frozen fits use all 512, the first 256,
and the second 256 paths. Selection uses only the pooled fit.

The predeclared selection objective is mean log(controlled payoff variance /
quote_scale²) over contracts **2/3/8/9**: the two nearest-ATM representatives
at each of 29 and 57 days. It intentionally does not select on unobserved tail
payoffs. The selected proposal and all design scores are saved before any
checking draw. No checking-based proposal or coefficient reselection occurs.

All three candidates are retained in **two fresh 1,024-path checking trials**
each, with coupled 64/128/256-step/day grids. Every one of the 12 representative
training contracts is evaluated, including tails. All fits share checking paths
within a candidate. Candidates use independent draws; checking trials share
frozen pilot coefficients. Uncertainty is across iid checking paths conditional
on those fits and the selection. Pilot-half comparisons are dependent checks.

Normalizer GH48 is held fixed: the prior same-path GH48/64 result was small,
not a population convergence certificate. No market calibration is performed.
The 12 representatives do not replace the full SPX maturity/split scope.

Deterministic quadrature is warmed separately before timed sampling. Cost rows
separate fine-state generation from stock/control computation at each grid.
Coarse stock/control timings use the same finest state cache and do not measure
standalone coarse-state simulation cost; no optimal MLMC allocation follows
from them alone.

## Validation and launch record

**133 package/driver tests pass** (103 package and 30 driver tests). New tests
cover asymmetric volatility shifts against joint score/OU Gaussian moments,
constant-volatility stock and Black payoffs, backward correspondence with the
old equity-only sampler, coarsening, invalid/duplicate tilt use, independent
quadratic payoff/control oracles for all pilot fits, and small end-to-end
pilot/three-grid checking collection. This is numerical/code correspondence,
not Lean verification or a full pricing certificate.

Run **001 failed before any market path draw** after deterministic setup:
seven empty lists were assigned to six collection variables. The code was
corrected, an end-to-end collection test added, and a new directory used.
Its [failure status](runs/m2_spx_volatility_sampling_001/status.json) and frozen
sources remain intact. No failed run is reported as numerical evidence.

Reproduce with the existing Python environment and single-thread BLAS:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.spx_volatility_sampling --output matteo/model_comparison/option_pricing/runs/m2_spx_volatility_sampling_002
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.verify_spx_volatility_sampling matteo/model_comparison/option_pricing/runs/m2_spx_volatility_sampling_002
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.audit_spx_volatility_sampling matteo/model_comparison/option_pricing/runs/m2_spx_volatility_sampling_002
```

Sources, the full active scope and curve-fit input are captured before execution.
Training, selection and checking seeds are disjoint. The verifier reconstructs
weights from all path scores, all full/half pilot coefficients, the selection
objective and winner, adjusted payoff draws, all estimates and paired grid
changes. Source snapshots and negative/zero-hit cases remain inspectable.

## Fresh checking results

The completed [run 002](runs/m2_spx_volatility_sampling_002/run.json) takes
275.6 seconds, including 86.9 seconds of setup. The
[audit](runs/m2_spx_volatility_sampling_002/audit.json) retains all comparisons.
The verifier reproduces likelihoods, all pilot fits and the design-selected
winner, **864 estimates, 576 paired corrections and 30 expectation monitors**.

The pilot selected **vol7**, with design geometric variance gain 1.1554 relative
to equity on the four central contracts. This is about 13% lower estimated
variance, not a 15% standard-error reduction. On fresh 256-step/GH48 checks:

| Pooled-fit diagnostic | Equity reference | vol7 (selected) | vol21 (retained comparator) |
|---|---:|---:|---:|
| Improved contract/trial variances versus equity | reference | **20/24** | 17/24 |
| Median equity/candidate variance ratio | 1 | **1.4192** | 1.2929 |
| Improved central-contract variances | reference | **8/8** | 8/8 |
| Median central variance ratio | 1 | **1.2957** | 1.3380 |
| Median controlled SE / quote scale | 3.5215 | **2.4683** | 3.3219 |
| Positive-hit SE <= scale/8 screens | 0/24 | **0/24** | 0/24 |
| Zero raw-payoff cases | 0 | 0 | 0 |
| Negative pooled price estimates | 0 | 0 | 0 |

The selected vol7 improves 10/12 cases in each checking trial; trial-wise median
variance gains are 1.3520 and 1.5925. Its range across all 24 cases is
0.5432 to 12.7840: four cases worsen. These are observed sample-variance ratios
from independent candidate samples, conditional on the frozen design. They are
encouraging replication evidence, not confidence bounds on population gains.
Do not retrospectively switch the primary proposal to vol21 on central quotes.
Its all-contract geometric gain is only 0.9618 because a rare put worsens sharply.

Checking-only generation plus all three stock/control grids costs approximately
43.82 seconds for equity and 47.82 for vol7. Combining these measured costs with
the median observed variance ratio gives about **1.30 times descriptive
variance-cost efficiency** for vol7. This is a single warmed timing comparison,
not a standalone finest-grid benchmark or a certified speedup. The original
normalizer setup cost is separate and shared.

### Retained expectation and pilot warnings

**Equity-reference trial 1 mass warning:** E_q[w] is estimated as
**1.04935 ± 0.01945**, or +2.54 observed SE from one. Its weighted stock and
reference expectations also exceed two SEs, moving with the same sample mass.
These are correlated diagnostics, not five independent failures. No weights
are rescaled to force mass or stock expectations to one. The independent
Gaussian likelihood oracles pass; this finite warning remains a reason to
retain fresh replication and uncertainty checks rather than declare convergence.

**Pilot instability:** vol21 half-B produces **-4.77501 ± 7.84805** for the
29-day 3400 put in checking trial 1, with four raw payoff hits. Its SE is
13.60 times the pooled-fit SE in that case. The pooled vol21 estimate there is
2.86491 ± 2.49888; neither estimator is accurate relative to scale 0.10434.
Both estimates are preserved. Even selected vol7 has a half-fit SE up to
2.14 times its pooled SE. Better central performance does not settle rare
payoffs or training-pilot stability. A zero half-fit difference in one vol7
contract is recorded as zero difference-SE, not price convergence.

### Coupled time integration remains a separate obstacle

For selected vol7, median absolute grid changes in quote scales are
**0.5794** for 64->128 and **0.9522** for 128->256. At 128->256,
20/24 observed changes exceed scale/8 and five exceed two paired SEs. Those
multiple dependent flags are diagnostic, not simultaneous statistical bounds.
For the 57-day 5550 put in trial 0, the change is **-1.62367 ± 0.67967** against
scale **0.20024**. The 29-day 6755 put changes -0.88172 ± 0.43381 and
-0.98754 ± 0.45936 across the two trials. These persistent warnings preclude
treating 256 steps/day as an integration clearance.

The median correction/fine-payoff variance ratio for vol7 falls from
**0.2254** at 64->128 to **0.1182** at 128->256. This supports investigating
multilevel allocation, but does not establish a decay order. Those ratios are
higher than the equity comparator's 0.1377/0.0770: optimizing fine-price
sampling has not automatically optimized the level corrections. Finest-cache
costs also do not identify standalone coarse-level costs. No MLMC estimator
has been completed by this experiment.

## Additional Gaussian score controls: frozen before checking results

While the checking simulations were being generated, before their completed
archives or numerical results existed, a separate
[score-control record](runs/m2_spx_score_controls_001/frozen.json) froze an
additional dictionary and ridge coefficients using only the saved design data.
The freeze command explicitly rejects an existing checking archive or completed
run record. It did not read partial checking paths or change proposal selection.
This is a separately recorded exploratory extension to the original experiment.

The original five target-zero controls are augmented with **w*Z_i** and
**w*(Z_i²-1)** for all four path scores: 13 controls in total. Their means are
zero under the proposal because weighted expectations recover the target
standard Gaussian moments. Ridge penalty one on centered, pilot-standardized
columns stabilizes fitting; full/half fits are frozen independently of checking.
The exact quadratic target oracle passes. Ridge and dictionary change together,
so their individual effects are not identified by this comparison.

The [replay audit](runs/m2_spx_score_controls_001/audit.json) validates all
**648** additional estimates. Against original pooled controls on identical
checking paths, pooled vol7 score controls improve 18/24 observed variances,
but the median gain is only **1.0178** (about 1.8% as a variance ratio).
Median SE/scale changes from 2.4683 to 2.4212; **0/24 sampling screens still
pass**. Equity/vol21 median gains are only 1.0006/1.0063. Half-pilot variants
remain less stable and vol21 half-B retains a negative estimate. This modest
add-on does not warrant claiming the sampling problem solved or selecting a
favorable half after checking.

Replay this extension (the freeze itself must precede completed checking data):

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.verify_spx_path_score_controls matteo/model_comparison/option_pricing/runs/m2_spx_volatility_sampling_002 matteo/model_comparison/option_pricing/runs/m2_spx_score_controls_001
```

There are now **134 passing tests: 103 package and 31 driver tests**, including
the score-control oracle. This supersedes the 133-test count at the first launch.

## Next experiment

Retain vol7 as a promising comparator, not a universally superior sampler.
Next prioritize **paired timestep corrections**: measure standalone level costs,
train any correction-specific controls/proposals on new design data, and compare
a bounded MLMC prototype against a single fine-grid estimator at matched cost.
Use Gaussian/Black oracles and independent evaluation; retain 128/256 warnings
and consider a small coupled 512-step check before any continuous-model claim.
Pilot stability and remaining tail precision stay explicit gates. The small
Gaussian score-control gain is recorded rather than expanded blindly. No market
refit, VIX maturity exclusion, forced normalization or price clipping occurs.
