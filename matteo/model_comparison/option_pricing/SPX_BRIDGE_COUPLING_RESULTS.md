# SPX endpoint-preserving Gaussian bridge coupling — 22 October 2026

This follows the negative [first MLMC benchmark](SPX_MLMC_RESULTS.md). The
question is whether changing the fine/coarse coupling can reduce correction
variance enough to pay for its additional work. The frozen Stage-3 M2 model,
selected nuisance curve, zero projected state and vol7 proposal are retained.
No market parameters are fitted. All four VIX maturities, 175 training options,
four forwards and original SPX scope remain active; held-out quotes are unused.

## Gaussian reflection and target-law correspondence

Let u1,u2 be the independent standard normals for two successive fine joint
OU/equity transitions. Joint roots L1,L2 map them to equity increments and OU
innovations. With A=diag(1,exp(-lambda*dt2)), the coarse innovation is

G = M u,    M = [A L1, L2],    u = (u1,u2).

A singular-value decomposition constructs an orthonormal row basis V for the
represented endpoint map, completed with null rows when necessary. Define

T = 2 V' V - I.

Then T'=T, T²=I and T'T=I, and M*T=M to numerical precision. Thus u and T*u
have the same target standard Gaussian law and the same coarse joint equity/OU
endpoint. All SVD rows are retained; at numerical rank deficiency this fixes
some additional null directions, so it need not be the full conditional-mean
bridge reflection or a variance-optimal choice. No covariance inverse is used.

The original and reflected fine OU paths are simulated separately. Their
coarse stock shocks and volatility-state endpoints coincide to monitored
floating-point tolerance. Their fine stock paths each use their **own adapted
start-of-step variance** in the existing exponential-Euler update. No future
variance is inserted into a stock step. The payoff average is taken afterwards.
No states are reset or forced to agree after simulation.

For fine payoffs P(u), P(T*u) and the common coarse payoff Pc, compare

ordinary correction: w(u) [P(u)-Pc],
reflected correction: w(u) [(P(u)+P(T*u))/2-Pc].

Here the original full-path proposal is q and w=f/q. Since target f is invariant
under T, E_q[w(u)*P(T*u)] = E_f[P(T*u)] = E_f[P(u)]. Therefore both corrections
estimate the same discretized fine-minus-coarse expectation. In density terms,
the reflected proposal is q_T(v)=q(T*v), so at v=T*u its likelihood weight is
f(T*u)/q(u)=w(u). **Do not evaluate the original f/q formula at reflected scores
and use that as its transported weight.** No sample-mass normalization occurs.

The fine/coarse target laws and option definitions remain unchanged. This is
an estimator/coupling change, not a new stock model. These algebraic statements
are accompanied by floating-point checks, not a Lean formalization.

## Independent pilots, pairing and cost accounting

Two corrections are studied: **128->256** and **256->512** steps per trading day.
Each has 256 independent training pairs. Ordinary and reflected corrections
fit their own weighted control dictionaries, using fine-minus-coarse controls
or the average of fine/reflected controls minus coarse controls. The original
four controls plus mass control are retained. Full-pilot and two disjoint
128-pair half-pilot fits are frozen before checking.

Two fresh checking trials use **512 pairs per correction**. One original and
one reflected path form **one iid observation**, not two independent samples.
Methods and pilot-fit variants share original paths; checking trials share
frozen fits but have independent samples. The difference between the two
estimators uses a paired SE. No method/fit is selected on checking results.

The ordinary full-path sampler is rerun with identical underlying seeds for
actual end-to-end ordinary timing and numerical path correspondence. Ordinary
and bridge execution order alternates across batches. The bridge timing pays
for reflection, the extra fine path and its pricing; coarse pricing is shared
within the bridge contribution. Normalizers/kernels are prepared separately.

Comparisons use equal counts of iid contributions and measured variance-cost
products. They are **not** another allocation experiment at exactly equal
runtime. A variance gain below the runtime multiplier does not establish an
efficiency gain. Kernel/normalizer and pilot overheads are reported separately.

## Validation and reproducibility

**141 tests pass: 108 package and 33 comparison-driver tests.** New checks cover:
reflection orthogonality/involution and endpoint preservation at unequal steps;
represented Gaussian covariance composition, including singular roots; original
sampler correspondence; weighted original/reflected Gaussian moments; reflected
intermediate OU projections; constant-volatility Black prices and vanishing
corrections; stochastic-volatility stock-forward expectations; and end-to-end
ordinary/bridge correction and alternating timing execution.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.spx_bridge_coupling --output matteo/model_comparison/option_pricing/runs/m2_spx_bridge_coupling_001
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.verify_spx_bridge_coupling matteo/model_comparison/option_pricing/runs/m2_spx_bridge_coupling_001
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.audit_spx_bridge_coupling matteo/model_comparison/option_pricing/runs/m2_spx_bridge_coupling_001
```

Sources, model/target hashes and full scope are captured before computation.
Raw fine/reflected/coarse payoffs, full-path Gaussian scores, weights, controls,
all fits and checking draws are retained. The verifier reconstructs likelihood
bookkeeping, coefficients, estimates, paired differences and expectation checks.
Signed corrections are retained; they are not negative option-price findings.
512 steps/day and normalizer GH48 remain numerical approximations, with no
continuous-model convergence or joint-market fitting clearance.

## Results: valid coupling, insufficient savings at measured cost

The [run](runs/m2_spx_bridge_coupling_001/run.json) completes in **364.0 seconds**,
including **181.8 seconds** of setup. The
[audit](runs/m2_spx_bridge_coupling_001/audit.json) retains all raw/controlled
comparisons. The verifier reproduces **384 correction estimates, 192 paired
method differences and 28 expectation monitors**, plus likelihoods, fits and
kernel diagnostics.

Maximum observed coarse OU endpoint mismatch is **6.31e-14**, and maximum
coarse equity-increment mismatch is **5.38e-17**. Kernel endpoint-map residuals
are at most 8.33e-16; orthogonality/involution residuals at most 7.78e-15.
The independently rerun ordinary sampler agrees with the streaming original
path contribution to the declared numerical tolerance. These establish code
correspondence at floating-point precision, not a population pricing bound.

For the pooled controls and fixed central contracts 2/3/8/9:

| Correction | Trial | Ordinary/bridge variance ratio | Bridge/ordinary runtime | Variance-cost efficiency gain |
|---|---:|---:|---:|---:|
| 128->256 | 0 | 1.1630 | 1.5169 | **0.7667** |
| 128->256 | 1 | 1.0481 | 1.5037 | **0.6970** |
| 256->512 | 0 | 1.0103 | 1.5103 | **0.6689** |
| 256->512 | 1 | 1.0994 | 1.5011 | **0.7324** |

A gain below one favors the ordinary correction. Thus **none of the four
central cost comparisons supports adopting this reflection**. Ordinary checking
takes about 9.62/9.65 seconds for 128->256 and 19.21/19.46 seconds for
256->512. Bridge checking takes 14.59/14.51 and 29.01/29.21 seconds respectively.
These are measured costs for the same 512 iid contributions, not predicted
costs or an exactly matched-runtime allocation.

Observed variance improves in 7/12 and 8/12 contracts at 128->256, and 6/12 and
11/12 at 256->512, across trials 0/1. Some observed variance reductions repeat on fresh samples, but these
descriptive sample ratios are not population guarantees. Raw corrections also fail the central cost test (efficiencies
0.64–0.83). The variance savings do not cover the extra fine-path work, with or
without the existing controls. No favorable estimator, contract subset or
pilot half is promoted after checking.

### What the reflection removes

An explanatory replay applies the **same frozen ordinary pooled coefficients**
to both correction constructions, without refitting on checking data. Write

ordinary residual = averaged residual + difference residual.

The central-contract variance decomposition is:

| Correction / trial | Averaged variance / ordinary variance | Difference variance / ordinary variance | 2 covariance / ordinary variance |
|---|---:|---:|---:|
| 128->256 / 0 | 0.86369 | 0.12450 | +0.01181 |
| 128->256 / 1 | 0.92705 | 0.13047 | -0.05753 |
| 256->512 / 0 | 0.96230 | 0.14289 | -0.10520 |
| 256->512 / 1 | 0.89736 | 0.13550 | -0.03287 |

Each row sums to one up to roundoff; the audit checks that identity. The
averaged residual retains **86–96%** of the original controlled variance.
Removing the difference component therefore cannot pay a roughly 50% runtime
premium in this implementation. The covariance term matters: there is no
claim that these components are orthogonal under the weighted proposal. This
is a post-hoc explanation of this coupling, not selection of new coefficients,
a universal attribution of all pricing noise, or evidence of a model defect.

### Retained checking and timestep warnings

All 28 mass/stock/reference discrepancies are within two observed SEs in this
run. No identity is forced. There are no contract/trial/level cases with all
original, reflected and coarse payoff-hit counts zero; nevertheless some
individual hit counts are as small as **two**, so tail corrections are not
accurately established by nonzero hits alone.

Three pooled paired method differences exceed two SEs and remain recorded:

- Trial 1, 128->256, **29-day 6755 put:** bridge minus ordinary
  **-0.43498 ± 0.19928** (-2.18 paired SE).
- Trial 1, 128->256, **57-day 5550 put:** **+0.42483 ± 0.21037**
  (+2.02 paired SE).
- Trial 1, 256->512, **29-day 7950 call:** **-0.00853 ± 0.00413**
  (-2.06 paired SE).

The two constructions have the same population correction target by the
Gaussian reflection/weight argument. These finite, dependent multiple-check
warnings neither establish different target expectations nor justify changing
weights to recover agreement. Independent algebraic and Gaussian/Black oracle
checks pass, while the numerical warnings are preserved.

Pilot-half fits remain imperfect: largest SE inflation relative to pooled is
**2.55** for bridge half-B and **3.59** for ordinary half-B at 256->512.
Signed negative corrections are permitted; no correction estimate is clipped.

**No timestep clearance:** pooled reflected 256->512 changes exceed scale/8
in **20/24** cases, with median absolute change **0.3648 quote scales**.
None exceeds two observed correction SEs in this particular check; this does
not erase the earlier independently recorded 512-step warnings or bound the
remaining continuous-time bias. At 128->256, reflected corrections for the
57-day 8200/8800 calls in trial 1 are +0.03730 ± 0.01578 and
+0.02955 ± 0.01472, exceeding two SEs. This experiment estimates corrections
only; it supplies no new complete MLMC/fine option-price accuracy screen.

## Decision and next hypothesis

Keep the verified reflection implementation as an exploratory comparator.
**Do not replace ordinary corrections or direct vol7 pricing with it**, and do
not rerun the failed MLMC allocation using assumed reflection savings.
The bounded experiment rejects this implementation as a cost improvement at
the tested grids and proposal; it does not reject every antithetic coupling.

The next candidate is **cheap volatility-factor martingale controls**, to
address residual noise without pricing an entire additional fine path. This
is a pending hypothesis, not implemented or established to help. For the exact
full-state update

X_next = D X + b + eta,    E[eta | past]=0,

projected innovation eps=P eta has known covariance R=P Q P'. Predictable
linear controls a' eps and quadratic controls tr(B*(eps eps'-R)) have zero
target mean. They can complement stock-delta controls and then be multiplied
by the exact path likelihood under proposal sampling. Any vega/readout-based
coefficient must use start-of-step information.

Implementation must retain the full-state conditional projection mean or
projected innovations explicitly: **(C,Y,F) alone is not a Markov state for the
multi-rate OU bank**, so its next mean cannot simply be guessed from those
three coordinates. Exact covariance/conditional-mean checks, pilot stability,
Black/forward oracles and fresh cost comparisons should precede adoption.
Apply candidate controls to direct prices and paired corrections separately;
variance reduction improves diagnostic power but does not remove timestep bias.
All VIX scope and separate quadrature/normalization gates remain unchanged.
