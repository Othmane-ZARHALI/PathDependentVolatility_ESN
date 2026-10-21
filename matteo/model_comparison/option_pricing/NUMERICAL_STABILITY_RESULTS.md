# Numerical stability diagnosis — 21 October 2026

This checkpoint answers two questions: why the first VIX curve fitting gain
reversed on fresh draws, and how to reduce SPX sampling error while separately
assessing time bias. Fixed model parameters and both curve candidates are
preserved. No further market fitting, candidate reselection or held-out quote
performance evaluation is performed.

## VIX: paired comparison of the same two candidates

`runs/m2_fit_reversal_audit_001` reconstructs prices from the archived curve
exposures. Its selected-curve prices and loss reproduce the original records.
For each sample, both curves use exactly the same outer states and numerical
rules. Define D = loss(selected) - loss(starting). Lower D favors selection.

| Sample | D | Leave-one-scramble-out range | Descriptive jackknife SE |
| --- | ---: | ---: | ---: |
| Fitting | -10.3056 | -17.5355 to -5.3255 | 8.7974 |
| Independent checking | +5.3688 | +3.9822 to +5.9753 | 1.3643 |

The fitting sample was used for selection, so its uncertainty diagnostic is
not an independent test. Four scrambles, the nonlinear interval loss and
poorly sampled tails also prevent treating these jackknife SEs as reliable
confidence certificates. The checking ranking is insensitive to deleting any
one of its four observed scrambles, while the fitting ranking is much more
variable in magnitude. Neither sample resolves unseen tail expectations.

The selected-curve median price SE divided by frozen market scale is 3.52
in fitting and 1.58 in checking. The target SE/scale is at most 1/8. The
numerical error is large enough to influence the loss materially. It must
not be hidden by changing quote denominators or dropping difficult options.

The squared interval loss is convex in price. In the simpler case where a
noisy price remains outside the same interval edge c,

\[
 E[(\widehat p-c)^2]=(p-c)^2+\operatorname{Var}(\widehat p)
\]

for an unbiased price estimator. Thus the loss of noisy prices has its own
sampling effect. Searching a fixed sample can additionally select parameters
that compensate for its particular price errors. Convergence of an optimizer
on that sample does not establish improvement in the underlying expected-price
objective. Common random numbers help compare candidates, but do not remove
finite-sample selection sensitivity or missing tails.

## Tail concentration and quadrature attribution

The 153-day block, including its forward residual, contributes D=-8.0288 in
fitting but D=+2.7761 in checking. The 118-day option block also changes from
favoring selection to opposing it. The largest selected-curve simulated VIX
values at 118/153 days are 161.80/212.45 in fitting, versus 85.29/114.73 in
checking. These are risk-neutral model draws, not observed market VIX values.

For the actual strike-47.5 training calls, the largest payoff supplies the
following shares of the simulated price sums:

| Expiry | Fitting hits / 256 | Largest payoff share | Checking hits / 256 | Largest payoff share |
| --- | ---: | ---: | ---: | ---: |
| 118d | 8 | 53.50% | 8 | 38.94% |
| 153d | 11 | 42.78% | 8 | 42.70% |

These shares follow directly from the archived exposure arrays and the formula
max(VIX-47.5,0), using the frozen selected levels. A few scenarios dominate
these expectations. Different observed tails can change both quoted-price
errors and the curve levels that appear optimal.

`runs/m2_vix_reversal_quadrature_001` isolates quadrature at 153 days. It
regenerates all 256 fitting and 256 checking states with their original seeds,
then changes conditional GH 32 to 48 and date-zero normalization 48 to 64.
The conditional G expectations are reused between normalizer comparisons,
so changing normalization introduces no sampling confound.

| Quadrature rule | D for fitting 153-day block | D for checking 153-day block |
| --- | ---: | ---: |
| Conditional 32, normalizer 48 | -8.028768 | +2.776144 |
| Conditional 48, normalizer 48 | -8.029117 | +2.778208 |
| Conditional 48, normalizer 64 | -8.027493 | +2.779250 |

Changing the outer sample moves this block's comparison by about 10.8 loss
units; these quadrature refinements move it by at most about 0.0031. This
strongly supports outer-sample/tail sensitivity as the main cause at this
important block. It does not prove full quadrature convergence, explain every
other maturity or establish the globally preferable curve. No candidate is
refitted or reselected using these diagnostics.

## SPX: predictable hedging controls

`esn_eq/hedging_controls.py` adds two Black-inspired adapted controls to the
existing stock-martingale and deterministic-volatility reference controls.
For normalized stock Y, a control has the form

\[
 Z=\sum_k\delta_k(Y_{t_{k+1}}-Y_{t_k}).
\]

The exponential-Euler scheme satisfies
E[Y_next-Y_current | F_current]=0. The control coefficient delta_k uses only
start-of-step information and is bounded. The tower property therefore gives
E[Z]=0 for the discrete pricing scheme. This is an ordinary mathematical
identity, not a Lean-certified theorem. Tests check its implementation.

One delta uses remaining deterministic integrated variance
W(t,T)=integral_t^T xi(u)du in the Black d1 formula. The other replaces W
with W*V_t/xi(t). The second is an adapted variance proxy, not a claim about
the true conditional variance or true option delta. Multiplication by the
actual discrete stock increment, rather than an approximate diffusion gain,
preserves zero expectation for this scheme.

For each payoff H, the adjusted estimate averages

\[
 H(Y_T)-\beta_1(Y_T-1)-\beta_2(H(Z_T^{ref})-Black)
       -\beta_3 Z^{hedge,1}-\beta_4 Z^{hedge,2}.
\]

Coefficients are fitted on independent pilot paths, then frozen. Evaluation
controls are not centered using their own sample means. No model parameter,
price expectation or stock forward is changed, and negative estimates would
remain visible rather than being clipped.

## Matched-path SPX experiment

`runs/m2_spx_hedging_001` uses the same fixed model and the previously selected
nuisance curve. Twelve actual training contracts cover puts/calls at 29/57
days. Targets are nearest available distinct strikes to predeclared moneyness
levels; the available upper wing at 29 days is limited, so its two far calls
are nearby. Selection precedes payoff evaluation. This is not a full-surface
acceptance test.

Two independent trials each use 1024 pilot paths and 2048 evaluation paths.
All three 64/128/256 grids share one finest exact-OU/equity path sample within
each trial. Counts are 2520/5040/10080 intervals and maximum lambda*dt is
0.30988/0.15494/0.07747. Normalizer order stays at 48. Fine-grid pilot
coefficients remain fixed across every evaluation grid. Nested count prefixes
512/1024/2048 are retained, but are not independent repetitions.

At 256 steps/day and 2048 evaluation paths per trial:

- Median estimated variance reduction versus the existing static controls is
  5.20 across the 24 contract/trial pairs. Its range is 1.00 to 24.98.
- The near-the-money 29-day put improves by factors 21.67/24.98 across trials;
  the corresponding 57-day put improves by 12.80/10.34.
- The largest SE across these pairs is 8.596 for plain pricing, 3.863 with
  static controls, and 1.733 with the augmented hedging controls.
- No augmented estimate meets the quote-relative sampling tolerance. Variance
  ratios are finite-sample diagnostics, not guaranteed efficiency factors.
  Extra control-computation cost was not timed separately; there is no measured
  full-surface speedup claim.

Some tails show no improvement because the independent pilot observes no
positive payoff: its regression target is identically zero, giving zero
control coefficients. This occurs for the 57-day strike-8800 call in trial 0
and the 29-day strike-3400 put in trial 1. In trial 1 the strike-8800 evaluation
sample also has no positive payoff. Its raw SE is zero; its adjusted estimate
is 0.3234 with SE 0.1602. Neither observation certifies tail accuracy.
Tail-aware sampling is needed in pilots as well as pricing evaluations.

`runs/m2_hedge_pilot_audit_001` independently regenerates the saved pilot
streams with the existing stock pricer and confirms exactly zero positive
payoffs in both cases. Their nearest normalized terminal-stock distances
from the strikes are 0.0742 and 0.0630, far beyond floating-point differences
between the two algebraically equivalent stock-update implementations.

The entire matched-path study took 213.49 seconds. Constant-variance oracle
tests show agreement with the existing stock pricer and decreasing hedging
error under refinement. A nonconstant-variance test checks control means
against zero and stock-path correspondence. The selector has a regression
test for the case where desired upper strikes exceed available contracts.

## Sampling and timestep convergence are separate

Using adjusted payoffs with variance sigma_adj^2, ordinary iid sampling has
SE approximately sigma_adj/sqrt(N), conditional on the independent pilot.
The agreed gate 2*SE <= s_i/4 requires SE <= s_i/8. Reducing variance helps
reach it, but reducing dt alone does not remove this sampling error.

For the four central contracts, naive fixed-variance extrapolation from this
pilot suggests roughly 140,000 to 430,000 paths per contract estimate to meet
that SE target. This is planning arithmetic, not a validated path budget;
tail variance estimates are much less reliable and can suggest millions.
Simply scaling current Monte Carlo is therefore not the preferred first step.

Time bias also remains open. Pooling the two independent trials, the 29-day
strike-6755 put has a 128-to-256 controlled change of -0.4913 SPX points
(paired SE 0.1959), against frozen price scale 0.5898. A near-zero change in
another contract is likewise not enough when its paired uncertainty exceeds
the desired bias budget. These are pointwise diagnostics, not simultaneous
confidence bounds. No asymptotic convergence order or Richardson correction
is inferred from the three levels. SPX normalizer error was held fixed, not
independently resolved by this experiment.

## Next implementation priorities

1. Test Gaussian mixture importance sampling for tails, with an unshifted
   component retained, exact likelihood weights and separate pilot/evaluation
   samples. This changes the simulation proposal, not the model Q dynamics or
   risk premium. Weighted martingale/reference controls must retain their
   known target expectations. Test both central prices and rare payoffs.
2. Amortize adapted controls across the surface using a modest shared strike
   basis, rather than computing a separate hedge at every quoted strike.
   Validate variance gains and cost on independent evaluation paths.
3. Measure level-correction variance and cost for multilevel Monte Carlo.
   The identity E[P_256]=E[P_64]+E[P_128-P_64]+E[P_256-P_128] targets the fine
   scheme even though its base term is coarse. Coupled differences and cost
   allocation can reduce fine-path work. Coarse terms do not certify coarse
   pricing accuracy; final continuous-time bias still needs finer checks.
4. Recheck the two frozen VIX candidates with larger independent scramble
   budgets and tail-aware sampling before further parameter selection. Keep
   quote scales fixed and retain every difficult contract. Do not use existing
   checking results as a new pristine validation set.
5. Require sampling, paired-grid and normalizer/quadrature tolerances separately
   before a fit ranking. The held-out market contracts remain reserved.

84 package tests and nine legacy/driver checks pass. Source snapshots and
input references accompany all new runs. Reproduce with fresh directories:

```sh
python -m matteo.model_comparison.option_pricing.fit_reversal_audit --output /tmp/reversal_new
python -m matteo.model_comparison.option_pricing.vix_reversal_quadrature --output /tmp/quad_new
python -m matteo.model_comparison.option_pricing.spx_hedging_study --output /tmp/hedging_new
python -m matteo.model_comparison.option_pricing.audit_hedge_pilots --output /tmp/pilot_audit_new
```

Use the existing project environment and single-thread BLAS settings recorded
in earlier checkpoints. No qualified joint fit or model ranking is claimed.
