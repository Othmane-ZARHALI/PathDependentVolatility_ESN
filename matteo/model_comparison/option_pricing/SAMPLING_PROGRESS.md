# Outer sampling and nuisance-input checkpoint — 21 October 2026

## Matched-cost VIX outer sampling pilot

`runs/m2_vix_outer_sampling_001` compares iid Gaussian sampling with
independently scrambled Sobol sampling at the 55-day training maturity. Model,
dated state and diagnostic curve remain the frozen Stage-3 candidate 18,
scenario 0, projected state zero and `training_mid_piecewise_v1`. No parameters
are fitted in this experiment. All 47 training option contracts are retained.

`esn_eq/outer_sampling.py` constructs the exact full terminal OU mean and
covariance from dated X0. Its PCA factor retains every column and orders them
by contributions to deterministic decayed C/Y/F projections over the VIX window.
No dimensional truncation or state substitution is applied. Both methods use
the same covariance factor. Sobol uses scrambled 52-bit digital points and an
inverse-normal transform; endpoints fail explicitly rather than being clipped.

There are eight independent complete replications per method, each with 128
outer states. Their first 64 points also produce a nested lower-budget estimate.
For both methods the displayed estimate averages replication means, and
SE = sample_SD(replication means)/sqrt(8). Sobol points inside one net are
dependent; a pathwise iid SE is not used. Across methods, seeds are independent.

Conditional GH order 32, date-zero normalizer order 48 and Gaussian time order
8 are fixed. Captured sources, target hashes, seeds, raw VIX values, replication
estimates, payoff-hit counts, paired count changes and variance ratios are saved.

| Method | Points per replication | Total points | Forward estimate | SE across replications |
| --- | ---: | ---: | ---: | ---: |
| iid | 64 | 512 | 20.99200 | 0.48262 |
| Sobol | 64 | 512 | 21.49042 | 0.17430 |
| iid | 128 | 1024 | 21.23435 | 0.49890 |
| Sobol | 128 | 1024 | 21.42255 | 0.06927 |

At 128 points per replication, the estimated iid/Sobol forward variance ratio
is 51.87. The median option variance ratio is 13.33 over the 46 contracts with
nonzero Sobol variance. Restricting descriptively to the 31 contracts with a
positive payoff in every replication under both methods gives median 15.24.
These ratios use only eight replications: they are exploratory, noisy and
parameter/maturity specific. They are not precise efficiency guarantees.

The largest option SE at that budget is 0.48954 for iid versus 0.09797 for
Sobol. The iid forward SE increased from 64 to 128 points in this finite study;
rare large observations make estimated variances unstable. This does not
establish deterioration of the underlying estimator with sample count.

At 128 points, Sobol still observes no positive payoff for the strike-170 call.
The strike-160 call has only one positive payoff across all eight replications.
Those are the two options passing a naive empirical-SE screen; neither is
evidence of accurately resolved tail pricing. The runner reports raw counts
and no acceptance classification. Zero sampled values do not bound missing tails.

The Sobol second-moment estimate is 599.12 (SE 19.86), compared with the curve
target 601.08. This is a numerical tower-moment diagnostic, not a fit to market
second moments. The forward's 64-to-128 paired change is -0.06787 (SE 0.13053),
still too uncertain to establish the planned forward tolerance.

At eight selected ordered VIX quantile positions from the first Sobol sample,
conditional GH refinement 32 to 48 changes VIX by up to 0.01637 points.
Holding conditional order 48, normalizer refinement 48 to 64 changes it by up
to 0.000304 points. This exposes remaining quadrature error at sampled tail
states; it is not a full-distribution error bound. Outer SE excludes both
effects. Total experiment runtime was 285.85 seconds with single-thread BLAS.

## Declared curve, state and risk-price families

`esn_eq/market_calibration.py` now implements:

- The three-piece nuisance curve already declared in `forward_variance.csv`,
  with boundaries 0, 0.25, 0.75 and the finite dataset horizon; each annual
  variance coefficient is bounded by 0.00010001 and 0.25. Its separate penalty
  is 0.1 times the sum of squared adjacent coefficient differences. Metadata
  expressions are checked against the declared contract, never executed.
- The minimum stationary-P-prior-norm full state for explicit `(C,Y,F)`
  projections, with penalty 0.1 times prior squared norm divided by three.
  This is a restricted latent state, not a history-filtered state.
- Five predeclared primitive risk-price scenarios, retaining equity risk price
  0.35 while feedback/spike unspanned prices are zero or individually +/-0.25.
  Scenario zero retains the rho*equity contribution to the OU Q drift.

These construct inputs; they do not establish physical compatibility or clear
the original preprocessing flags. The state-order mismatch in the calibration
plan has been corrected and recorded in `AUDIT.md`. Existing zero-state runs
are unaffected.

## A constrained curve starting point

`runs/m2_curve_initializer_001` constructs a starting point for later search.
It minimizes mean squared relative deviations from diagnostic stripped SPX
total variances at 29/57/211/302 days, plus the declared curve penalty. It
enforces the existing conservative VIX option moment inequalities and Jensen
inequalities **assuming model means equal quoted forwards**. These are linear
constraints in piecewise variance levels; the auxiliary objective is quadratic.
The SLSQP solve succeeded with constraint residuals within floating-point error.

The resulting annual variance levels are [0.0545936, 0.0671315, 0.0231561].
Necessary second-moment constraints bind at 55 and 153 days. Relative changes
in stripped SPX total variance are +42.60%, +26.33%, +9.88% and +1.10% at the
four training maturities. This exposes the trade-off imposed by this restricted
family and these moment conditions; it does not establish option-price errors
of those percentages or inadequacy of M2.

The initializer does not evaluate M2 option prices, fit readout/state parameters
or establish a feasible joint market model. Its auxiliary strip-deviation
objective must never be added as independent information to the common quote
loss. Candidate moment checks must be recomputed at actual model VIX means;
passing necessary inequalities is insufficient to match option intervals.
Original diagnostic curves and compatibility flags are unchanged.

## Validation and next step

77 package tests and eight legacy/driver checks pass. New tests cover full
covariance reconstruction, Q mean/projection correspondence, prefix preservation,
independent scrambles, zero-expiry behavior, replication-based SEs, curve bounds,
penalties, coordinate ordering and preservation of the equity-induced drift.

Next, use independently scrambled outer sampling in the exploratory fitter,
with quadrature refinement and payoff-hit diagnostics retained. At fixed
readout, state and risk prices, conditional integrated variance is linear in
piecewise curve coefficients: precomputing its segment exposures should make
curve profiling substantially cheaper than rerunning quadrature for each trial.
This exposure cache and a common-loss market optimization remain to implement.
The third curve coefficient starts after 0.75 years, beyond all four VIX
windows; VIX prices alone cannot identify it. SPX prices and the full joint
objective are needed. Keep all held-out quotes untouched.

Reproduce from repository root in the project Python environment:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -m matteo.model_comparison.option_pricing.vix_outer_sampling \
  --output /tmp/m2_outer_new
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -m matteo.model_comparison.option_pricing.curve_initializer \
  --output /tmp/m2_curve_new
```

Use `/private/tmp/pdv-options-venv/bin/python` for the environment used here.
No model-comparison ranking or final accuracy claim is supported yet.

The subsequent completed checkpoint is [CURVE_PROFILE_RESULTS.md](CURVE_PROFILE_RESULTS.md):
cached curve exposures and an exploratory VIX-only fit. Its fitting improvement
reverses under fresh numerical repricing; the failure to establish improvement
is retained explicitly.
