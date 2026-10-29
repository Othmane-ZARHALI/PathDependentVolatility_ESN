# More stress half-lives: grid resolution versus additional random factors

27 October 2026. Subsequent user steering: four nodes are a first test, not a
structural limit; more half-lives may be added if justified. The current
four-node independent-bank experiment remains in progress as specified.

There is no positivity or Gaussian-normalization obstruction to a larger finite
bank. With nonnegative weights summing to one, each normalized response has
mean one, and the common variance curve remains $E_Q[v_t]=\xi(t)$. More nodes
increase pricing cost and fitted parameters; neighbouring speeds can be hard
to distinguish on the available horizon. This is an identification issue, not
a mathematical prohibition.

However, increasing the number of independent drivers is not simply increasing
the quadrature resolution of one fixed continuous stochastic spectrum.

Let $S_t=\sum_j\alpha_j A_t^j$. For independent states,

$$\operatorname{Var}(S_t)=\sum_j\alpha_j^2\operatorname{Var}(A_t^j),$$

and for the normalized Gaussian exponentials,

$$\operatorname{Var}(A_t^j)=e^{b_j^2q_j(t)}-1.$$

With equal weights $1/N$ and equal response variances $V$, this gives $V/N$.
The individual shocks average out. More generally, with fixed total stress mass
$w$ and uniformly bounded response variances,

$$\operatorname{Var}(S_t)\le
\left(\sup_j\operatorname{Var}(A_t^j)\right)
\left(\max_j\alpha_j\right)w.$$

Thus if a finer grid spreads mass so its largest atom tends to zero, the
independent stress contribution converges in mean square to its deterministic
mean $w$. The old M2 contribution can remain random. This is not the intended
limit of a persistent random stress spectrum. Uniformly bounded amplitudes,
such as the current bound5, imply the required variance bound because $q_j\le1$.
Finite calibrations can instead concentrate weights or raise sensitivities, but
those compensations can create unstable or weakly identified parameters.

## A way to refine a spectrum while retaining common shocks

One possible future construction uses the same Brownian shock at different
speeds:

$$Z_t(a)=\sqrt{2a}\int_0^t e^{-a(t-s)}\,dB_s.$$

This gives OU states with the same marginal variance as before, but correlated
across speeds:

$$\operatorname{Cov}(Z_t(a_i),Z_t(a_j))=
\frac{2\sqrt{a_i a_j}}{a_i+a_j}
\left(1-e^{-(a_i+a_j)t}\right).$$

Normalize each exponential as before. If $C_{ij}(t)$ denotes this covariance,
then

$$\operatorname{Cov}(A_t^i,A_t^j)=e^{b_i b_j C_{ij}(t)}-1.$$

The cross terms in $\operatorname{Var}(\sum_j\alpha_jA_t^j)$ then preserve shared
stress as the grid becomes finer. Splitting one component into two identical
copies driven by the same shock leaves the model unchanged; splitting it into
independent copies generally reduces its variance.

A single common Brownian motion is only one proposal. A small number of common
shock families, each filtered at multiple speeds, would retain several sources
of stress while permitting a smoother persistence spectrum. This would change
the dependence structure of the current independent bank and must be tested as
a separate model. It should not be silently substituted during an ongoing fit.

The present four-node test addresses a modest finite independent bank. Its
results cannot by themselves establish that ever denser independent banks are
a suitable approximation to a continuous random measure over speeds. No new
experiment, correlation choice, or Lean theorem is claimed in this note.
