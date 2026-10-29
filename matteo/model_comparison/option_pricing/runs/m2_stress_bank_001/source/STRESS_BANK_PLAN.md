# Fixed geometric stress bank

27 October 2026. Joint VIX-only experiment authorized after the two-speed study.

Choose four half-lives before fitting: **7, 28, 112, 448 calendar days**,
$h_j=7\,4^j$, $j=0,1,2,3$. The endpoints represent a week through approximately
15 months. These are calendar-based scales rather than the previously fitted
71/365-day speeds. Geometric spacing gives equal resolution on a logarithmic
time axis. Four nodes keep the test modest. This is a design choice, not a
claim of uniquely optimal endpoints. A later, separately recorded grid shift
could assess sensitivity; do not select this grid using checking prices.

The first node can represent stress shorter than the 30-day VIX averaging
window; its contribution is naturally damped by that window. The longest node
allows persistence beyond both the longest option expiry and the earlier
365-day calibration ceiling. A half-life is a shock decay time, not an option
expiry. All nodes act on all maturities.

For independent OU states, also independent of the old exogenous M2 bank,

$$dZ_t^j=-a_jZ_t^jdt+\sqrt{2a_j}dB_t^j,\quad
Z_0^j=0,\quad a_j=\log2/h_j,$$

$$A_t^j=\exp(b_jZ_t^j-\tfrac12b_j^2(1-e^{-2a_jt})).$$

Use globally shared nonnegative weights $\alpha_j$, with $\sum_j\alpha_j\le1$:

$$v_t=f+(\xi(t)-f)
\left[(1-\sum_j\alpha_j)R_t+\sum_j\alpha_j A_t^j\right].$$

Hence $E[v_t]=\xi(t)$ and $v_t\ge f$ when $\xi\ge f$. Total stress weight is
$w=\sum_j\alpha_j$, and the earlier conditional weights are $\pi_j=\alpha_j/w$
when $w>0$. The speed measure is atomic on the fixed grid; a continuous measure
is not simultaneously fitted.

Fit four amplitudes $0\le b_j\le5$, four weights, seven independent active
levels of the same 16-piece variance curve, and floor volatility 5–12%:
The four weight units use stick-breaking internally to preserve nonnegativity
at every trial (including the M2 remainder), without altering the admissible simplex.
There are 16 parameters in total, versus 14 for the fitted-speed two-node experiment.
Nine inactive levels remain unchanged; active levels lie between the floor
and 0.2. The original M2 bank, premiums, readout, positive spike shift and cap
are frozen. Keep all 217 native quotes and four forwards, and preserve signed
estimates. Three starts, at most 100 SLSQP iterations each; select only by
training option loss among candidates with maximum forward gap <=1e-5 points.
Four equality constraints apply to the training estimator, not certified
population means. SPX remains paused; neither baseline is replaced.

## Reuse and estimator

For European VIX prices, cache terminal normal coordinates at the four expiries
and evaluate future conditional variance analytically. This implements one
common OU model; it does not simulate complete joint time paths, which these
marginal prices do not require. Amplitude and weight changes rebuild responses,
variance integrals, VIX and payoffs, while terminal coordinates are reused.

Reuse the original old/slow coordinates as node3 and the previous independent
fast standard-normal coordinates as node2. Rescale them using the newly fixed
448/112-day variances. Standard normal coordinates are reusable; previous
OU state values with different speeds are not silently reused unchanged.
Generate nodes0/1 jointly using a new two-dimensional defensive Gaussian
mixture and an independent within-net permutation. Its target is two independent
standard normals. Multiply all three component likelihoods, without sample
normalization. The proposal dependence does not change target independence.
Training reuses 32,768 old/slow/fast coordinates per expiry; checking reuses
65,536 and draws new additional coordinates. Checking is not wholly fresh.

Use

$$E[P_{bank}]=E[P_{reference}]+E[P_{bank}-P_{reference}],$$

where the reference replaces nodes0/1/2 by their unconditional response mean
one, retaining node3 and M2. It is a full normalized variance process before
taking the square root. Evaluate its large term using only old/node3 likelihoods;
the paired difference uses all four target normal coordinates and matching
product likelihoods. Recompute legacy controls under each likelihood; paired
controls cancel. At weights0/1/2=0 the difference is exactly zero. This identity
changes the estimator, not the bank model.

Compare with the previously selected two-speed model on its stored checking
blocks, disclosing shared old-state streams. This is not a nested comparison:
the fixed grid does not contain its fitted speeds. Report option loss, all
maturity wing errors, quote counts, checking forward gaps/SE, weighted tails,
and bound hits. Check all 16 derivatives, response normalization, nested
reference replay, scalar conditional/time-integral reconstruction, likelihoods
and source/input hashes. No new Lean theorem or global optimum is claimed.
