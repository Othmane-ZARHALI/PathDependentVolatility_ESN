# Why the 27 day VIX right wing remains too flat

The 27-day fit has too much option value near the centre of the smile and too little far out on the right. The evidence points to a restriction in how the model distributes risk between moderate moves and rare large moves, together with limited control of shape across maturities. It does not identify a unique missing factor.

This diagnosis uses the selected continuous experiment and its matching two-component reference. It adds read-only calculations; it does not recalibrate parameters, change the current model or alter the original market scope.

## 1 What the first maturity means

The stored pricing date implied by the quote dates is 19 February 2026. The first VIX option expires 27 days later, on 18 March 2026. At that expiry, the model VIX measures expected variance over the following 30 days, days 27–57:

$$
\operatorname{VIX}_{27/365}^2
=\frac{10000}{\Delta}
\mathbb E\!\left[\int_{27/365}^{57/365}v_s\,ds
\mid\mathcal F_{27/365}\right],\qquad\Delta=30/365.
$$

There are therefore two different time intervals: randomness accumulates until option expiry, and the resulting VIX looks forward over a future variance window. A factor that fluctuates strongly before expiry need not affect the whole future window strongly.

## 2 The error is present in option prices

The following figures use canonical undiscounted call prices. Where the original quote was a put, the existing parity conversion is retained. Implied volatilities use the observed market forward for every curve.

| Strike in VIX points | Market call mid | Model call | Market Black IV | Model Black IV |
|---|---:|---:|---:|---:|
| 20 | 2.22636 | 2.45686 | 100.33% | 111.04% |
| 40 | 0.52673 | 0.35984 | 183.06% | 167.20% |
| 60 | 0.21182 | 0.10249 | 211.27% | 188.55% |
| 100 | 0.06232 | 0.01651 | 239.30% | 207.99% |

These are implied volatilities of **VIX options**, not stock volatilities. Values exceeding 100% are possible.

At strike 100 the market bid–ask is approximately $[0.02854,0.09610]$. The model is below even the bid. Its reported checking standard error is approximately $0.0000479$, much smaller than the price shortfall. Numerical quadrature refinement at sampled stressed states changes VIX by only a few millionths of a point. These diagnostics make a pricing-noise explanation unlikely at these parameters, although they are not uniform error guarantees.

Over the 14 quotes with $K/F>2$, the average absolute IV error is 22.65 percentage points. The fitted slope against log strike-moneyness is

$$
\frac{d\,\mathrm{IV}}{d\log(K/F)}
\approx52.38\quad\text{for the market},\qquad
43.28\quad\text{for the model},
$$

measured by linear regression and expressed in IV percentage points per unit log-moneyness. It is a descriptive slope over observed strikes, not an asymptotic theorem.

Using the model forward instead of the market forward changes right-wing IVs by at most about 0.262 percentage points. This is far smaller than the 22.65-point average error. The remaining forward discrepancy therefore does not explain the displayed right-wing shortage through the IV-inversion convention. No price or forward was shifted to obtain this comparison.

![The smile and the actual right-wing call prices](runs/m2_continuum_reduction_note_001/27day_shape_and_prices.png)

The 27-day slice contributes approximately 0.9311 of the total 1.7307 native option loss, or 53.8%. It is already a major part of the calibration objective. It would be misleading to say it is simply ignored because three later maturities are also fitted.

## 3 A concrete mechanism in the fitted model

The best continuous measure is close to the two-component boundary. That boundary gives a useful exact decomposition. Conditional on the old reservoir state at expiry, write

$$
\operatorname{VIX}_T(Z_T)
=100\sqrt{H_T(Z_T)},\qquad
H_T(z)=B_T+\int_T^{T+\Delta}A_T(s)e^{\lambda_T(s)z}\,ds,
$$

where

$$
\begin{aligned}
B_T&=f+\frac{1-w}{\Delta}\int_T^{T+\Delta}
(\xi(s)-f)\,
\mathbb E[R_s^{\mathrm{M2}}\mid X_T]\,ds,\\
\lambda_T(s)&=b e^{-a(s-T)},\\
A_T(s)&=\frac{w(\xi(s)-f)}{\Delta}
\exp\!\left(-\frac{\lambda_T(s)^2q(T)}2\right).
\end{aligned}
$$

The constant term $B_T$, conditional on $X_T$, contains the floor and old M2 contribution. The integral is the stress contribution. All coefficients are nonnegative.

Differentiating gives

$$
\frac{\partial\log\operatorname{VIX}_T}{\partial z}
=\frac12\,
\underbrace{\frac{\int A_T(s)e^{\lambda_T(s)z}ds}{H_T(z)}}_{
\text{stress share of conditional variance}}
\underbrace{\frac{\int\lambda_T(s)A_T(s)e^{\lambda_T(s)z}ds}
{\int A_T(s)e^{\lambda_T(s)z}ds}}_{
\text{average stress sensitivity}}.
$$

The same formula extends to a positive integral over component positions by averaging over both $s$ and $u$, with sensitivities $bu e^{-a(s-T)}$.

To understand the shape, consider the one-sensitivity approximation

$$
H(z)=B+A e^{\lambda z},\qquad
p(z)=\frac{Ae^{\lambda z}}{B+Ae^{\lambda z}}.
$$

Then

$$
\frac{d\log\operatorname{VIX}}{dz}=\frac{\lambda}{2}p(z),
\qquad
\frac{d^2\log\operatorname{VIX}}{dz^2}
=\frac{\lambda^2}{2}p(z)(1-p(z)).
$$

The response bends most while the two contributions are comparable. Once the stress term dominates, $p$ approaches one and log VIX becomes almost linear in the Gaussian stress factor. A single exponential of a Gaussian produces a lognormal variable with a flat Black smile. The full model has a time integral and residual M2 randomness, so this is a mechanism explaining a tendency toward flattening, not a statement that its finite-strike smile must be exactly flat.

### What the payoff-weighted diagnostics show

For the matching two-component reference at the selected continuous parameters, weight each simulated state by its contribution to a call payoff, including its exact importance weight. This describes the states that matter to the option, rather than the most extreme unweighted simulation values.

| Strike | Payoff-weighted stress share of conditional variance | Payoff-weighted sensitivity of log VIX to a standard normal |
|---|---:|---:|
| 20 | 68.88% | 0.541 |
| 40 | 91.39% | 0.720 |
| 60 | 96.31% | 0.760 |
| 100 | 98.71% | 0.781 |

At strike 100, approximately 99.84% of the raw reference call payoff comes from states whose stress variance share exceeds 90%. The raw payoff effective sample size is about 1,478; the largest ten contributions account for about 1.13% of total payoff. These concentration measures are descriptive, not independent convergence certificates.

This is the strongest direct clue: **the far-right option is already sampling a region dominated by one stress response**. Adding intermediate components along the tested interpolation did not materially change that region in the selected fit. The mechanism has already used much of its transition from a mild response to a strong response.

The percentages above decompose conditional variance in payoff-contributing states of the matching reference. They are not a decomposition of the continuous model's option price into causal percentages. The continuous fit is close to that reference, but this distinction matters.

## 4 Shared parameters restrict the maturity shape

The stress strength $b$, persistence $a$ and measure $\mu$ are shared across the four expiries. The seven free variance levels change the variance profile, but they do not give each expiry an independent response shape.

In the main fit, $b\approx3.188$ is below its search ceiling of 5. Its wing error is therefore not simply a stress-strength parameter stuck at that ceiling. Increasing it may move other strikes or maturities adversely; that trade-off needs a constrained comparison. The retained broader-measure fit does reach the ceiling, which limits conclusions about that candidate.

Since $Z_0=0$, the conditional stress loading on a standard normal coordinate at expiry $T$, for a future variance date $T+\ell$, is

$$
L(T,\ell)=b\sqrt{1-e^{-2aT}}\,e^{-a\ell}.
$$

For the fitted $b\approx3.1877$ and 105.29-day half-life, at the illustrative middle of the VIX window, $\ell=15$ days:

| Expiry | $q(T)$ | Standard-normal stress loading $L(T,15)$ |
|---|---:|---:|
| 27 days | 0.2992 | 1.5796 |
| 55 days | 0.5153 | 2.0730 |
| 118 days | 0.7885 | 2.5645 |
| 153 days | 0.8666 | 2.6884 |

A smaller absolute loading at shorter expiry is **not sufficient** to explain an annualized IV slope. In the single-time pure-stress approximation, the Black volatility would be

$$
\sigma_{\mathrm{Black}}(T,\ell)
=\frac{L(T,\ell)}{2\sqrt T}.
$$

At the midpoint this proxy is about 290% for 27 days and 208% for 153 days. The factor $1/\sqrt T$ reverses the simple comparison of absolute loadings. These are illustrative proxies, not the actual fitted smile levels. In particular, the proxy exceeds the market's far-right 27-day IV: there is no demonstrated hard volatility ceiling below the market.

### A calculation illustrating the persistence trade-off

At fixed $b$ and one fixed future lag $\ell$, differentiating $\log L$ with respect to $a$ gives

$$
\frac{\partial\log L}{\partial a}
=\frac{T}{e^{2aT}-1}-\ell.
$$

The maximum occurs at

$$
a_* = \frac{\log(1+T/\ell)}{2T}.
$$

For $T=27$ and $\ell=15$ calendar days, the corresponding half-life is about 36.4 days. It raises this 27-day loading from approximately 1.580 to 1.92, while lowering the 153-day loading from approximately 2.688 to 2.39 at the same $b$. Thus changing one persistence parameter can help a front loading while moving later loadings the other way.

This is not a recalibration and does not prove that a 36-day half-life improves the full smile. The complete VIX uses all future dates, the floor, M2, the normalization and the forward constraints. It is a transparent illustration of why one shared persistence can create a compromise across maturities.

![The persistence trade-off in the single-time proxy](runs/m2_continuum_reduction_note_001/OU_midwindow_tradeoff.png)

[Bergomi's formulation](https://www.risk.net/media/download/920801/download), pp. 91–92, allows maturity-dependent mapping parameters. Our present response parameters are global, so the comparison does not establish that we have tested the same degree of maturity flexibility.

## 5 The first variance budget is also slightly too small

For any fit satisfying all archived bid–ask intervals and the exact observed 27-day forward $F=20.1$, the previously constructed static-payoff certificate implies

$$
\mathbb E[\operatorname{VIX}_T^2]\geq468.99933.
$$

The selected curve gives the model identity

$$
\mathbb E[\operatorname{VIX}_T^2]
=\frac{10000}{\Delta}\int_T^{T+\Delta}\xi(s)ds
=467.44959.
$$

The certificate was reverified in exact rational arithmetic, including its global payoff inequality and unchanged quote hash. The shortfall is approximately 1.54974, or 0.332% of the current second moment.

Conditional on the exact target forward, this corresponds to

$$
\operatorname{Var}(\operatorname{VIX}_T)=63.43959
\quad\text{versus a necessary value of at least }64.98933.
$$

Therefore this **fixed fitted curve** cannot place every option inside its archived spread while matching that exact population forward. The bound does not say that the entire model family is infeasible, and it does not locate the missing variance exclusively in the right wing.

The deficit is small relative to the visible shape error. Raising the average variance by at least $0.00015497$ would remove this particular moment obstruction, but would not establish a smile fit. Raising variance levels alone also moves the VIX forward; the response shape must adjust if the forward is to remain fixed.

Within days 27–57, the current curve has:

| Future interval | Annual variance level | Share of expected non-floor variance in the window |
|---|---:|---:|
| Days 27–29 | 0.0144, at the floor | Essentially zero |
| Days 29–55 | 0.047864 | 89.66% |
| Days 55–57 | 0.064544 | 10.34% |

The last two days also enter the 55-day option's future window. This creates some coupling between the first two curves. These levels were already free to fit; this is not a remaining three-level-versus-seven-level implementation mistake.

## 6 Which explanations are supported

| Explanation | Assessment |
|---|---|
| A forward-conversion plotting error creates the whole gap | Not supported: the measured wing-IV convention effect is at most 0.262 points, versus a 22.65-point average error; native far-call prices are also too low. |
| A handful of rare paths creates the observed tail fit | The weighted concentration checks do not indicate this in the current reference sample. They do not prove full convergence. |
| The 27-day slice receives negligible optimization weight | Not supported: it accounts for about 54% of the total option loss. |
| The current first-window second moment is insufficient for a complete within-spread fit | Established, conditional on the exact observed population forward, for this fixed curve by the payoff certificate. |
| Stress dominance and globally shared response parameters restrict the observed first wing | A well-supported mechanism and leading hypothesis. A controlled calibration comparison is still needed to establish how much of the mismatch it explains. |
| Every continuous integral model must fail | Not established. Only a particular interpolation path, Beta family and bounded search were tested. |

The old M2 cap is still part of the model. However, the current far-right reference states are overwhelmingly stress-driven. Simply identifying the old cap is not enough to explain this particular remaining gap.

## 7 The next discriminating experiments

The first useful experiment would be a **separate 27-day-only diagnostic calibration**, with its own exact numerical forward constraint, the same response family and bounds, and the three variance levels affecting days 27–57. The original four-maturity fit and all 217 quotes should remain intact. Reprice the other maturities afterward to expose the trade-off; do not treat that diagnostic as an adopted reduction of market scope.

If it can fit the first wing, that would support an across-maturity restriction. If it cannot, examine the marginal response family, floor/curve bounds and optimization more closely. Failure of a bounded optimizer alone would not prove impossibility.

If maturity flexibility is needed, two controlled extensions preserve a single stock-variance process:

1. **A deterministic variance-of-stress clock.** Let

   $$
   dY_t=-aY_tdt+\eta(t)dB_t,\qquad
   Q(t)=\int_0^t e^{-2a(t-r)}\eta(r)^2dr,
   $$

   and use $e^{uY_t-u^2Q(t)/2}$. The present stress factor corresponds to $Y_t=bZ_t$ and constant $\eta=b\sqrt{2a}$. A small number of time intervals for $\eta$ can change how quickly uncertainty accumulates. It is not unlimited independent expiry control: nonnegative innovation variance requires

   $$
   Q(T_2)\geq e^{-2a(T_2-T_1)}Q(T_1),\qquad T_2>T_1.
   $$

2. **A deterministic calendar-time measure $\mu_s$.** Define the same instantaneous variance $v_s$ using a positive normalized measure that changes with calendar date $s$. Each component remains mean one, so $\mathbb E[v_s]=\xi(s)$. Overlapping VIX windows still couple maturities. Assigning a different ad hoc measure to each option expiry without a common variance process would not preserve the present model.

A separate direction is to enlarge the component family itself. The current interpolation ties the M2 power $c$ and stress strength $\lambda$ through $c=1-u$, $\lambda=bu$. A genuinely larger family could integrate

$$
\frac{G(X_t)^c}{\mathbb E[G(X_t)^c]}
 e^{\lambda Z_t-\lambda^2q(t)/2}
$$

over a joint positive probability measure on $(c,\lambda)$, initially with $0\leq c\leq1$ and $0\leq\lambda\leq\Lambda<\infty$. That would release this tie. Its usefulness, integrability bounds, sampling and added calibration freedom would need their own tests. None of these extensions has been implemented or calibrated as part of this diagnosis.

## Evidence

- [Empirical reduction of the continuous model](CONTINUUM_EMPIRICAL_REDUCTION.md).
- [Calibration results and estimator limitations](CONTINUUM_RESPONSE_RESULTS.md).
- [Reproducible diagnostic script](diagnose_vix27_continuum.py).
- [Machine-readable diagnostics](runs/m2_continuum_reduction_note_001/diagnostics.json).
- [Payoff-weighted component contributions](runs/m2_continuum_reduction_note_001/payoff_weighted_reference_components.csv).
- [Matching-range maturity slopes](runs/m2_continuum_reduction_note_001/common_range_wing_slopes.csv), which avoid comparing different strike ranges across expiries.

Numerical audits, exact rational moment certificates, mathematical response identities and causal hypotheses have different scopes. No new Lean theorem is claimed here.
