# Gaussian controls for conditional VIX — 22 October 2026

## Frozen scope and question

Continue the forward-precision bottleneck from
[VIX_ALLOCATION_REFINEMENT_RESULTS.md](VIX_ALLOCATION_REFINEMENT_RESULTS.md).
Keep all four maturities (27/55/118/153 days), 175 training options and four
forwards. The model, initial state, risk premium, selected variance curve,
17-component targeted proposal, conditional GH32, normalizer GH48 and Gauss
time8 remain fixed. This study learns sampling-control coefficients, not
model parameters or a market fit. SPX is unchanged.

The Gaussian-control expectations do not depend on the VIX quadrature being
exact. The study does not force E[VIX squared] to its exact-model target or
rescale prices to market forwards. The lean-numerics-bridge evidence boundary
remains: algebraic identities, implementation checks, artifact correspondence
and numerical evidence are distinct. No Lean proof or rigorous price-error
certificate is supplied.

## Controls and estimator

Under the target standard Gaussian latent law f, write the full terminal state
as X_T=mu_T+L_T z. The same fixed four score-sensitive orthonormal directions
U from the preceding studies give y=U'z. Under f, E[y_j]=0 and E[y_j^2-1]=0.
For draws from the frozen proposal q with exact density weights w=f/q, define
nine control rows

    C(z) = (w-1, w*y_1, ..., w*y_4,
                  w*(y_1^2-1), ..., w*(y_4^2-1)).

Every component has expectation zero. The full 36-dimensional state is still
sampled; controls do not truncate it or change the proposal. For a VIX forward,
second moment or option payoff H, estimate

    mean(w H) - mean(C)' beta.

The coefficient beta is frozen independently of checking paths. Hence the
subtraction preserves the expectation of the implemented payoff integral,
conditional on any pilot fit. There is no division by the observed mean w and
no substitution of a market or model VIX moment. This formula applies to the
random-mixture draws used here, not directly to component-allocation weights.
Control values and adjusted estimates may be negative; nothing is clipped.

## Pilot training at the right sampling unit

Promote **only** the prior `m2_vix_targeted_allocation_001` random-targeted
12 x 128 observations per maturity to pilot data. They were checks for that
previous study and are not fresh checks for the new control fit. Regenerate
and verify their states and weights before forming controls. All new checking
seeds are independent of those pilot seeds.

Train on complete-replication means, not individual paths. Let Y be the
12 x number-of-payoffs matrix of raw pilot estimates and C the 12 x 9 matrix
of pilot mean controls. Center these pilot arrays for covariance estimation;
standardize controls by their sample standard deviations. With standardized
pilot matrix Z, solve

    A = [Z'Z/(R-1) + lambda I]^{-1} Z'(Y-mean(Y))/(R-1),
    beta_j = A_j / pilot_sd_j,     lambda = 1.

The dimensionless ridge is fixed in advance, not selected on checking results.
The checking controls are **not** centered by the pilot mean: their known
expectation is zero. Constant pilot controls get zero slopes. Twelve pilot
replications for nine controls are few; regularization improves conditioning
but does not guarantee useful out-of-pilot variance reduction.

Five estimators use each identical checking path:

- `raw`: no correction.
- `mass`: ridge fit of the mass control w-1 alone, from all 12 pilot replications.
- `gaussian`: all nine controls, all 12 pilot replications; **primary method**.
- `gaussian_half_a` and `gaussian_half_b`: the same fixed ridge/dictionary fitted
  on the first six or last six pilot replications. These disjoint half fits
  diagnose pilot sensitivity; the pooled fit shares their data. They are not
  candidates chosen after inspecting checking results.

Fit separate coefficients for every payoff and maturity. Use the exact same
coefficients at both checking counts; the paired 64-to-128 change therefore
does not change the fitted control. The pilot trains at the fine count 128,
so coarse performance is not optimized separately.

## Fresh checking design

`m2_vix_gaussian_controls_001` uses 16 independent scrambles per maturity,
nested 64/128 points per scramble: **2048 costly conditional valuations per
maturity, 8192 total**. Every estimator uses these same states, weights and
VIX values. Estimator differences and count changes are paired; the complete
scramble is the uncertainty unit. Five estimators do not require five times
the conditional-pricing cost or constitute five independent checking samples.

At the first scramble of each maturity, re-evaluate minimum/median/maximum
sampled VIX states with GH48 (12 states total). These checks are limited
sensitivity diagnostics, not global quadrature clearance. The original time,
normalizer, finite-conditional-rule and tail gates remain open unless separately
resolved. Reported ± values are one checking SE conditional on a frozen pilot
fit; repeated-pilot efficiency requires a separate assessment.

## Results

### Forward precision improves, especially at short maturities

| Days | Raw forward ± SE | Mass-only forward ± SE | Primary Gaussian forward ± SE | Gaussian / raw SE |
| --- | ---: | ---: | ---: | ---: |
| 27 | 20.84947 ± 0.09392 | 20.82834 ± 0.07667 | 20.81306 ± 0.05735 | 0.611 |
| 55 | 20.40928 ± 0.08802 | 20.45402 ± 0.07474 | 20.38003 ± 0.04649 | 0.528 |
| 118 | 21.24809 ± 0.09301 | 21.25109 ± 0.08852 | 21.18557 ± 0.08702 | 0.936 |
| 153 | 21.04119 ± 0.13883 | 21.06468 ± 0.11934 | 21.05835 ± 0.11230 | 0.809 |

All primary forward SEs decrease, with variance gains about 2.68/3.58/1.14/1.53
at 27/55/118/153 days on these checking paths. The **27- and 55-day primary
forward SEs pass the observed threshold 0.5/8 = 0.0625**; raw and mass-only
pass none. The two longer forward maturities remain unresolved. Passing this
sampling screen does not establish total price accuracy, bias control or
market-forward fit.

The estimates differ because controls correct observed sampling fluctuations,
not because the market target or model mean was changed. For example, the
118-day forward correction is -0.06252 ± 0.02951 paired SE. These are
many dependent diagnostic comparisons; neither a realized nonzero correction
nor a small paired SE is a calibrated test of a model/normalizer defect.

### Option gains are modest and not uniform

| Days | Median option SE: Gaussian / raw | Raw SE screens | Gaussian SE screens |
| --- | ---: | ---: | ---: |
| 27 | 0.827 | 6/40 | 7/40 |
| 55 | 0.891 | 9/47 | 9/47 |
| 118 | 1.000 | 10/44 | 11/44 |
| 153 | 0.959 | 12/44 | 12/44 |

Total primary screens are **39/175 versus 37/175 raw**; the mass-only method
also has 37/175. The screen requires a positive adjusted estimate, at least
two raw positive payoffs and observed SE <= quote scale / 8. All pooled option
estimates are positive for every reported estimator; no negative result was
clipped. Every raw option has positive payoffs in every checking replication,
but most quotes still fail the SE screen. These controls are not a solution
to the complete option-wing accuracy problem.

### Moment checks and pilot sensitivity remain explicit

| Days | Exact normalized-model target | Raw second moment ± SE | Primary Gaussian second moment ± SE |
| --- | ---: | ---: | ---: |
| 27 | 542.49611 | 539.67449 ± 2.94311 | 539.80318 ± 2.23683 |
| 55 | 542.49611 | 549.44773 ± 4.34092 | 546.03824 ± 2.92058 |
| 118 | 600.31899 | 604.56687 ± 4.43201 | 602.79969 ± 4.20151 |
| 153 | 600.31899 | 598.69410 ± 3.81934 | 599.78679 ± 3.75925 |

Primary residuals are within about 1.22 reported SE. This is not a forced
identity, a confidence guarantee or a population conditional-quadrature check.
The earlier saved 118-day warning and its higher-order attribution remain
preserved independently.

Pilot variation matters. The half-A/half-B forward SEs are respectively:

| Days | First-six pilot fit SE | Last-six pilot fit SE |
| --- | ---: | ---: |
| 27 | 0.06809 | 0.05660 |
| 55 | 0.06029 | 0.04218 |
| 118 | 0.09954 | 0.09160 |
| 153 | 0.12142 | 0.12206 |

One half fit fails the 27-day forward screen and the half-A fit actually
worsens 118-day forward SE relative to raw. Its median 118-day option SE ratio
is 1.394, versus approximately 1.000 for the pooled primary fit. This is a
limitation of the small pilot and fixed dictionary; do not select the most
favorable half after checking and label it an independently validated method.

**Retained diagnostic warning:** the half-B 153-day controlled second moment
is 605.93740 ± 2.69630, residual +5.61841 from 600.31899 (about 2.08 reported
SE). Half-A gives 595.65633 ± 4.05269, while the primary pooled fit gives
599.78679 ± 3.75925. These alternatives share checking paths and are not
independent replications of fitted-model performance. The warning is retained,
not interpreted as a normalizer bug or concealed by reporting the pooled fit
alone. Small-group conditional SEs and unresolved tails/bias still limit inference.

### Count and quadrature checks remain open

Largest absolute paired 64-to-128 state option changes divided by quote scale
are 0.825/0.715/1.087/1.089 (raw) and 0.394/0.928/0.459/0.704 (primary), at
27/55/118/153 days. The full paired changes/SEs are saved. Reduced variance
at the fine count is not complete count-convergence clearance.

The 12 selected-state GH32-to-GH48 VIX-level changes have maxima
0.00076/0.00082/0.00256/0.00118 by maturity. Those checks do not certify global
conditional-GH, normalizer or time-integration bias. No new time/normalizer
refinement was performed in this control experiment.

## Decision and next bounded work

Retain the frozen Gaussian-control method as a useful experimental variance
reduction for subsequent convergence tests, with raw estimates always reported.
It improves the primary forward SE at all four maturities in this fresh check,
but its gains are weakest at 118 days, pilot fits vary, and the full option
surface remains far from the required sampling accuracy. No new calibration,
joint feasibility, held-out result or model ranking is accepted.

Next assess the stability of replication-trained coefficients with larger,
independent pilot batches and count refinements, retaining the same control
identities as a baseline. Prioritize 118/153-day forwards and option tails while
keeping all four maturities in the study. New dictionaries, ridge settings or
tail controls must be designed on pilot data and tested on new independent
paths; the present checking set would then become design data. Do not keep
adding controls on the same checking observations until an SE screen passes.
Keep the conditional/time/normalizer gates separate from sampling improvement.

Artifact replay passes for 40 captured source/input hashes, 64 unique checking
seeds disjoint from prior studies, all 8192 full states and control rows, every
pilot coefficient, 1830 summary rows, per-replication estimates, raw hit counts,
negative-estimate flags and paired estimator/count changes. Gaussian node
expectations are not recomputed by this verifier. Runtime was 795.7 seconds;
these are correspondence/implementation results, not price certification.

## Validation and reproduction

99 package tests plus 25 driver tests pass (**124 total**). New tests check
weighted Gaussian means, an exact polynomial control oracle under shifted
importance sampling, regularized slopes/constant controls, and the possibility of
negative adjusted values. Source review confirms that the driver does not clip
its adjusted estimates. These implementation checks do not establish
market-price accuracy.

Reproduce with a new output directory:

```sh
env OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /private/tmp/pdv-options-venv/bin/python \
  -m matteo.model_comparison.option_pricing.vix_gaussian_controls \
  --output matteo/model_comparison/option_pricing/runs/NEW_CONTROL_RUN
```

The driver explicitly uses `m2_vix_targeted_allocation_001` as the pilot and
`m2_vix_score_sampling_001/proposals.json` for the frozen directions. It saves
sources before pricing, rejects mid-run changes, and freezes coefficient and
proposal hashes. Pilot arrays, coefficients, full checking states/control rows,
raw values, replication estimates and paired count/estimator changes are saved.
`verify_vix_gaussian_controls` regenerates pilot fits and checking arrays and
recomputes statistics. `audit_vix_gaussian_controls` reports descriptive SE
ratios, negative-price counts, pilot-fit sensitivity and necessary SE screens.
Those screens require a positive adjusted option estimate and at least two
raw payoff hits; a small SE alone is not tail or bias certification.
