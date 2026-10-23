# Partial conditional SPX pricing investigation — 23 October 2026

## Mathematical construction and source correspondence

The Gaussian M2 OU states are exogenous: their transitions depend on Gaussian
factor drivers and deterministic Q drift, not on the simulated stock level.
The equity Brownian driver can be projected onto the primitive driver rowspace
and an independent orthogonal direction. For the current bank architecture,
let nu² denote the squared equity loading in that orthogonal direction.

At each dated step, draw retained joint shocks `(a_k, eps_X,k)` with covariance
`K_k - nu²*dt_k*e_equity*e_equity'`, then add independent `nu*sqrt(dt_k)*Z_k`
only to equity. The full OU state uses eps_X,k, including its exact decay/Q drift;
no projected-state Markov approximation or fine-step covariance inverse is used.
Given all retained shocks U, all V_k are known and the existing exponential-Euler
stock has

`log(S_T/F_normalization) = m_T(U) + sum sqrt(V_k(U))*nu*sqrt(dt_k)*Z_k`,

`m_T = sum [-V_k*dt_k/2 + sqrt(V_k)*a_k]`,

`v_T = nu² * sum V_k*dt_k`.

Thus log stock is conditionally Gaussian and the conditional option payoff is
a closed lognormal/Black expectation, evaluated with nonnegative log-domain tail
differences and put-call parity. Average these payoffs under target Q. No fitted
control subtraction or clipping of a negative final estimate is involved.

Conditioning is on the retained driving shocks (which determine the full OU
path, and may contain more information), not necessarily on the state path alone.
Under the exact matching Gaussian law, this is Rao-Blackwellization:
`Var(E[payoff|U]) <= Var(payoff)`. This does not guarantee a large gain, low Monte
Carlo error or continuous-model accuracy. Stock discretization/GH48 remain fixed.

## Crucial frozen-model limitation

Feedback correlation is **0.99625996**, spike correlation **0.97905483**. Their
shared primitive drivers leave only **nu²=0.0063668537395**, or **0.636685%** of
equity Brownian variance independent. The model can therefore integrate only a
small equity component with this simple construction. This is not a bound on
option-payoff variance reduction; near a payoff kink, a small component can
still matter. No large speedup is assumed.

The primitive-driver nullspace is computed on the constant loading matrix,
not by thresholding a nearly singular fine-step OU covariance. No OU mode is
removed. The prototype starts from the implemented legacy joint covariance
(reconstructed from its conditional root), subtracts the independent equity
variance, and refactorizes the retained joint shocks. Finite-precision PSD
repair errors are checked and recorded. Exact Gaussian algebra and numerical
covariance correspondence are distinct assessments; no formal global payoff
error bound or bitwise equivalence of different RNG factorizations is claimed.

## Oracle checks and recorded discrepancy

Five package tests cover driver fractions/covariances, stable nonnegative Black
tails and parity/zero variance, exact independent constant-volatility prices,
correlated nonlinear readout conditional means and the unchanged legacy stock
rule, and a perfectly equity-driven limiting case.

The first independent-Black oracle failed at roughly 1e-9 price error: SVD
returned an independence fraction just below 1, causing a sqrt(epsilon) spurious
retained equity shock. This was corrected by handling exactly uncorrelated and
perfectly equity-driven structural endpoints explicitly, rather than loosening
the oracle tolerance. All five tests now pass. The frozen model is not such an
endpoint, so its nu² remains unchanged.

A separate frozen-grid check covers all 30 distinct step sizes: maximum joint
covariance discrepancy versus the legacy numerical sampler is 2.08e-16 in
absolute units and 6.37e-11 after standardizing each covariance entry. This is
numerical/source evidence, not a proof of identical nonlinear payoff errors.

## Frozen bounded market experiment

`m2_spx_conditional_001`: same 24 representatives across 29/57/211/302 days, cutoff 3
and all 707 priority quotes retained; all four VIX maturities remain active.
Two independent 256-path target-Q trials, batches of 128, seed 20261123. Raw and
partial-conditional payoffs share the same retained U shocks. No pilot fitting,
importance sampling, old full-path weights, new features or model/curve/state
refitting. GH48 normalizers are reused only after exact grid/model/state/target/
curve/source/hash checks. Price mean/SE, paired differences, positivity, residual
concentration and transition covariance errors are retained.

The comparator is **uncontrolled target-Q Monte Carlo**, not the currently
controlled/importance-sampled reference pricer. Costs cover a coupled prototype
that computes both methods; no standalone production variance-cost claim will
be made. Zero empirical raw payoffs and underflowed conditional tails do not
establish accurate far-option prices, even if a sampling SE is tiny.

## Results

Run `m2_spx_conditional_001` completes in **71.04 seconds**. All **96 prices**,
**48 paired differences**, raw payoffs, lognormal conditional payoffs (also via
an independent direct-CDF Black expression), source hashes and all 30 transition
covariance audits replay. All **117 package tests + 42 study-driver tests** pass.
No estimate is negative. Conditional payoff values are nonnegative by construction;
this positivity claim applies to the raw conditional average, not to a future
linear-control subtraction.

### Observed sampling benefit

Observed variance improves in **28/48** contract/trial comparisons. Median
raw/conditional variance gain is **1.0253** over all 48, or **1.0343** over the 42
cases with positive raw payoff observations. The remaining six zero-hit raw
cases have zero empirical raw variance; that is not evidence that raw pricing
is more accurate than conditional pricing.

| Maturity | Median observed raw/conditional variance gain | Raw median SE/scale | Conditional median SE/scale |
| --- | ---: | ---: | ---: |
| 29 days | 1.0822 | 17.72 | 17.41 |
| 57 days | 1.0605 | 25.53 | 25.01 |
| 211 days | 1.0156 | 41.39 | 40.58 |
| 302 days | 0.9846 | 43.30 | 43.93 |

These are noisy finite-sample ratios, not a contradiction of the population
conditional-variance inequality. Both methods share retained shocks, all strikes
and all maturity checkpoints; two independent 256-path trials do not establish
convergence. No separate standalone cost or comparison with ridge5/vol7 has
been measured. The 71-second coupled prototype timing cannot be compared directly
to earlier studies with different proposals, control loops and expiry sampling.

### Zero-hit and false-precision warnings

Six outer-call contract/trial cases have **no positive raw payoff**. Conditional
integration gives positive estimates, but four would pass the old SE-only
threshold while still lacking outer-path evidence. They are explicitly excluded
from conservative sampling-pass counts (0/48 passes for either method):

| Trial | Days / strike | Conditional estimate +/- SE | Raw hits |
| --- | --- | --- | ---: |
| 0 | 29 /7725 call | 2.406e-10 +/-2.223e-10 | 0 |
| 0 | 57 /8100 call | 0.006152 +/-0.006152 | 0 |
| 1 | 211 /9600 call | 1.793e-10 +/-1.793e-10 | 0 |
| 1 | 302 /10400 call | 0.001645 +/-0.001619 | 0 |

For the same 29-day 7725 call, the second conditional trial gives
**0.310452 +/-0.310452**, rather than 2.4e-10. Many paths have tiny positive
conditional payoffs, but almost all meaningful contribution can still come
from a single retained volatility path. Maximum centered residual shares reach
99.61%. Conditional positivity removes the negative-control-price mechanism;
it does not remove rare-volatility sampling problems or underflow limitations.

### Remaining numerical warnings

Four of 48 paired conditional-minus-raw price means exceed two empirical paired
SEs: trial 0, 302-day 5550 put (-2.33 SE) and 6500 put (-2.85 SE); trial 1, 211-day 7475
call (+2.04 SE) and 302-day 6500 put (+2.14 SE). A 302-day trial 0 paired stock mean
also flags 2.72SE, while all individual raw/conditional stock-forward means are
within two empirical SEs. Related strikes/maturities share paths; flags are
correlated descriptive diagnostics, not independent tests or proof of a bias.
The maximum covariance errors remain **2.08e-16 absolute / 6.37e-11 standardized**.
No covariance or payoff failure is hidden, and no centering or clipping is applied.

## Assessment and pending extension

The Gaussian construction works and supplies a **nonnegative partial-conditional
benchmark** under the existing intended discrete model, with numerical transition
correspondence audited. With feedback/spike almost fully revealing equity noise,
its raw-payoff variance benefit is modest in the frozen model. Keep it
experimental; **original vol7/ridge5 remains reference**. Do not refit the market
model, interpret tiny SEs as accurate tail prices, or remove maturities/quotes.

A modest raw-payoff gain does **not** rule out a larger gain from conditioning
the existing controlled residual. Next derive/test conditional expectations of
the stock, deterministic reference payoff and adapted delta gains, then derive
an importance proposal on the **retained** Gaussian shocks. Old full-path weights
cannot simply multiply an already-integrated payoff. This extension is pending,
not implemented or validated by the present market screen.

For example, under retained U let log S_t ~ N(m_t,v_t), and let A_t denote the
conditional expected stock growth over the next step. For the existing Black
proxy with known remaining variance W_t,

`E[S_t * delta_t | U] = exp(m_t+v_t/2) * [Phi((m_t+v_t-log K+W_t/2)/sqrt(v_t+W_t)) - 1_put]`,

`E[delta_t*(S_next-S_t) | U] = (A_t-1) * E[S_t*delta_t | U]`.

Here `A_t=exp(-V_t*dt/2 + sqrt(V_t)*a_t + V_t*nu²*dt/2)`; it is generally not 1
because U includes correlated equity shocks. Setting the conditional gain to
zero would be an implementation error despite its zero **unconditional** mean.
These identities are algebraic design candidates; they have not been added to
the prototype or assigned performance results. A controlled estimator may again
be negative even when its conditional payoff component is nonnegative. Keep
stock-timestep and GH-normalizer refinements as independent pending gates.

## Cache-runner maintenance

Adding an unrelated package module exposed an overly strict assumption in the
older call-tilt runner's cache validation: it required a prior hash for every
currently present module. It now checks every package source recorded when the
cache was produced, allowing additional unrelated modules while still rejecting
changes to all original normalizer dependencies. Archived sources/results remain
unchanged. The conditional run also validates the independently recorded cache
file hash and all original package dependencies before reuse.
