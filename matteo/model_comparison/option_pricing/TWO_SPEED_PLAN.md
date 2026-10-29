# A shared fast and slow stress response

Date: 27 October 2026. Authorized joint VIX-only experiment; no maturity-specific rules.

## Model and interpretation

The empirical endpoint reduction of the continuous response experiment motivates
starting from its simpler mixture. This experiment adds a second stress state,
rather than reopening the Beta spectrum simultaneously. The old M2 state and
readout are frozen. Let $R_t=G_t/E_Q[G_t]$ be that positive normalized readout.
For $i\in\{f,s\}$ introduce independent OU states, also independent of M2:

$$
dZ_t^i=-a_iZ_t^i\,dt+\sqrt{2a_i}\,dB_t^i,\quad Z_0^i=0,
\qquad a_i=\log(2)/h_i,\quad q_i(t)=1-e^{-2a_it}.
$$

The normalized stress responses and instantaneous variance are

$$
A_t^i=\exp(b_i Z_t^i-\tfrac12b_i^2q_i(t)),\qquad
v_t=f+(\xi(t)-f)\left[(1-w)R_t+w\{(1-r)A_t^s+rA_t^f\}\right].
$$

Here $h_f\le h_s$ are globally shared half-lives, $r$ is the fraction of stress
variance assigned to the fast response, and $w$ is the total stress fraction.
All weights are nonnegative and sum to one. Each response has mean one, hence
$E_Q[v_t]=\xi(t)$ and $v_t\ge f$. At $r=0$ the previous single-speed model is
recovered exactly. This is a mixture of variance responses before taking the
VIX square root, not a mixture of completed option prices. Independent stress
innovations are a parsimonious assumption to test, not a market finding.

Every maturity uses these same dynamics and parameters. Different expiries see
different state distributions and different overlapping future 30-day windows.
There is no special treatment of 27 days and no parameters indexed by option expiry.

## Conditional calculation

For $s\ge T$,

$$
E_Q[A_s^i\mid\mathcal F_T]=
\exp\left(b_i e^{-a_i(s-T)}Z_T^i
-\tfrac12b_i^2 e^{-2a_i(s-T)}q_i(T)\right).
$$

Combine this with the existing conditional M2 calculation, integrate the
conditional variance over $[T,T+30/365]$, then take the square root:

$$
VIX_T=100\sqrt{\frac1\Delta\int_T^{T+\Delta}E_Q[v_s\mid\mathcal F_T]\,ds}.
$$

The two exponential terms permit both different sensitivities and different
persistence. They do not guarantee steep right wings: a dominant exponential
can still produce an approximately lognormal tail. The experiment must decide.

## Fit and numerical evidence

Fit all 217 native quotes and all four forwards jointly, preserving native price
signs and the existing interval-based objective. Keep seven free active levels
of the existing 16-piece variance curve; nine inactive levels stay unchanged.
Bounds: $0\le b_i\le5$, $1\le h_f\le h_s\le365$ calendar days,
$0\le w,r\le1$, floor volatility 5–12%, and $f\le\xi_j\le0.2$.
The slow half-life parameterization depends on the fast half-life to enforce
ordering; its full chain-rule derivatives are included.

Use four equalities on the calibration estimator, then report checking forward
errors and block standard errors separately. Equality of sampled means does not
establish exact population means. No forward adjustment or quote clipping.

Use the multilevel identity

$$
E[P_{two}]=E[P_{slow}]+E[P_{two}-P_{slow}],
$$

with coupled old/slow states for the difference. Training reuses 32,768 old/slow
states per expiry, checking reuses 65,536 previously checked old/slow states.
Fast states are newly drawn independently with defensive Gaussian importance
sampling and an independently randomized within-block permutation. Multiply the
old likelihood by the fast likelihood in BOTH paired difference terms. The
large reference retains its original likelihood. Legacy zero-mean controls are
recomputed under the corresponding likelihood. Eight randomized blocks supply
SE estimates; this is not a wholly fresh or market out-of-sample study.

Select among four bounded starts using training loss and forward feasibility.
Include a refitted single-speed boundary. Check all 14 price/forward derivative
directions, exact boundary nesting, likelihood-weighted tail concentration,
conditional-response formula, and time-integration refinement before interpreting
results. Preserve optimizer failures and numerical problems. No new Lean theorem,
SPX pricing, baseline replacement, or claim of global optimality is implied.
