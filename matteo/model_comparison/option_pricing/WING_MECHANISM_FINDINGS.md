# M2 VIX right-wing mechanism investigation — 26 October 2026

The combined formal/numerical approach can exclude some explanations and test
specific mechanisms. It does not identify a unique missing process from four
option slices. Current selected M2 and active VIX-only scope remain unchanged.

## Findings from the current checkout

1. The intended normalized model has
   E[VIX_T²] = 10000/Delta * integral_[T,T+Delta] xi(u) du.
   This is a tower-expectation identity under the usual integrability and
   conditional-variance representation. That continuous-model bridge has not
   been formalized here. Existing exact-rational portfolio certificates exclude
   the frozen budgets at 55/118/153 days if every archived quote and forward must
   fit. Those certificate instances remain Python exact arithmetic, not Lean.
2. The latest free-curve study already released seven VIX-sensitive pieces of the
   original 16-piece curve. Fresh numerical loss improves 27.3491 to 20.5245
   with curve alone, or 19.8086 with the fitted single post-cap factor. Only
   21/217 quotes fit for the latter; its 118/153-day budgets still violate the
   necessary exact-forward lower bounds. Neither a full-smile fit nor structural
   impossibility of the model family follows from the bounded optimization.
3. More stochastic amplitude is not a sufficient shape prescription. The earlier
   spike-sensitivity layer flattened measured slopes; the strong after-cap
   exponential layer steepened some slopes but did not solve the fit. A single
   normalized exponential has constant log response to its Gaussian state.

Sources: VIX_ONLY_CURVE_RESULTS.md, POSTCAP_VOV_RESULTS.md,
STOCHASTIC_VOV_RESULTS.md, STOCHASTIC_VOL_OF_VOL_RESEARCH.md, active_scope.json.

## New controlled mechanism experiment

For Z standard normal set

    R(Z) = sqrt(.9 exp(.18 Z - .18²/2) + .1 exp(1.8 Z - 1.8²/2)),
    X = R(Z)/E[R(Z)],
    Y = exp(s Z - s²/2), s² = log(E[X²]).

Both have mean 1 and second moment approximately 1.03184937. Here E[R²]=1
exactly by Gaussian exponential centering. Numerical E[R]=0.984445884579356.
Y has an exactly flat Black smile; X has two responses to the same Gaussian
state. No new randomness, jumps or increase in the second moment is needed
for this illustrative comparison.

At an illustrative 153-day expiry, quadrature prices imply:

| K/F | Single-response annual IV | Mixed-response annual IV |
| --- | ---: | ---: |
| 1 | 27.35% | 21.14% |
| 2 | 27.35% | 48.03% |
| 4 | 27.35% | 66.29% |
| 8 | 27.35% | 78.95% |

This isolates a shape mechanism, not an across-strike improvement: the mixed
case lowers the ATM price/IV. These are artificial marginal distributions,
not market-fitted numbers or a simulation of M2's 30-day conditional variance.
There is no continuous-time consistency or SPX compatibility claim.

Direct scalar quadrature is compared with 16 independent scrambled Sobol nets,
8192 points each, sampled under a four-component full-support Gaussian mixture.
Exact density ratios are used without self-normalization. Every mixed call
estimate is within 1.37 reported scramble SE of its quadrature oracle. The
single-response 8F call is about 1.59e-33: the simulation observes no payoff,
and its zero observed SE is explicitly NOT an accuracy certificate.

The integration interval is [-12,12]. The analytic bound from
sqrt(x+y)<=sqrt(x)+sqrt(y) bounds omitted mean and call tails by about 1.35e-29
after normalization. QUADPACK errors and floating-point evaluation of that tail
bound are numerical evidence, not interval-arithmetic certificates. Exact mean
and moment matching refer to exact normalization; numerical implementation
uses its quadrature approximation. Preliminary exploration is retained and
excluded from final evidence.

Reproduce with /private/tmp/m2-postcap312/bin/python
runs/m2_wing_mechanism_001/experiment.py from this directory.
Full definitions, seeds, proposal, per-scramble values, errors and source hash:
runs/m2_wing_mechanism_001/results.json.

## What Lean establishes

formalization/WingMechanisms.lean proves four restricted statements:

- The curvature numerator identity
  (a² A+b² B)(A+B)-(a A+b B)² = A B(a-b)².
- Its ratio to (A+B)² is strictly positive if A,B>0 and a!=b.
- The global payoff inequality 4K max(y-K,0)<=y² for real y,K.
- An integrable portfolio dominated by X² has expected cost no greater than
  E[X²]. The payoff domination and portfolio cost are explicit premises.

For A=c1 exp(a z), B=c2 exp(b z), ordinary differentiation identifies the first
ratio with d² log(A+B)/dz². Thus the log sensitivity increases with z: the
stronger branch takes over in high states. That derivative identification is
reviewed mathematics, not encoded in these algebraic Lean lemmas. None of the
four statements proves positive implied-volatility skew, the numerical table,
the archived market certificate instances or an M2 fit.

check_wing_lean.py reuses the existing pinned read-only Mathlib environment,
checks tracked dependency cleanliness and revision, compiles with warnings as
errors, rejects admitted/custom-axiom source and audits all four dependencies.
Audit attempts are retained under formalization/wing_audits.

## Concrete next M2 candidate and discriminating experiment

Keep the original capped readout B_t=G_base(X_t)/M_base(t), floor f and positive
beta_J. Extend the existing independent OU factor, rather than altering the
reservoir law immediately:

    L_i(t)=exp(gamma_i Z_t - gamma_i² q(t)/2),
    L_mix(t)=(1-w)L_1(t)+w L_2(t),
    v_t=f+(xi(t)-f) B_t L_mix(t),
    0<=w<=1, 0<=gamma_1<gamma_2, q(t)=Var(Z_t).

Independence and E[L_i]=1 preserve E[v_t]=xi(t). The original single factor is
nested at equal sensitivities or a zero branch weight. The exact conditional
multiplier is the same convex sum of the two existing analytic conditional
exponentials. This offers a small implementation change and reuses exogenous
base conditional nodes; affected VIX values and option payoffs must be rebuilt.
If later returning to SPX, variance and stock paths must also be rebuilt.

This is a proposed M2 adaptation; no production mixed kernel or market fit was
implemented in this turn. A fair comparison frees the same seven active curve
pieces in all cases, retains all 217 native quotes/four forwards and compares
curve-only, single-response and mixed-response families. Track ATM fit, forward,
second-moment feasibility, high-strike slopes/curvature, weighted payoff ESS and
concentration, paired price differences and independent numerical checking.
Choose/declare the price and wing-shape objectives before selection. A total
loss reduction alone need not solve wing steepness. Check several plausible
parameter bounds and starts before inferring a structural limitation.

If a mixed readout fails, continuous stochastic noise amplitude in the reservoir
is a distinct next hypothesis: dX_i=-lambda_i X_i dt + sigma_i sqrt(Gamma_t)dW_i,
with a positive continuous stochastic Gamma. It changes the transition law and
needs a new dated normalizer/conditional pricer, positivity and integrability
analysis, and new paths. The current evidence does not establish that this
larger extension is necessary. Jumps and quintic readouts remain deferred.

Mixed-exponential VIX mechanisms and calibration benchmarks appear in Bourgey,
De Marco and Gobet, section 3, https://arxiv.org/html/2202.10413 . This supports
the choice of experiment, not validation of the proposed M2 adaptation.
