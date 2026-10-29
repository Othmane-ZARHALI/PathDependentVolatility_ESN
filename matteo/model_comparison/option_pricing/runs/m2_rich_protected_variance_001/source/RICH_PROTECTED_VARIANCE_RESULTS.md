# Richer VIX-protected variance shapes — 29 October 2026

The richer shape improves aggregate SPX fit in this checking sample, but the observed downward skew remains far from market. No candidate is adopted. The production model and previous runs are unchanged. [Seven SPX smiles](visualizations/rich-protected-variance029/SPX-rich-protected.png), [time scaling](visualizations/rich-protected-variance029/local-scaling.png), full run `runs/m2_rich_protected_variance_001`.

## What this tests

The previous protected correction used a single global normalization and a product envelope. This experiment tests a more efficient local time normalization and two smooth state shapes. It keeps the same M2 readout, cap, Q normalizer, strictly positive spike shift, both stresses, mixture, stock theta=.8 and forward-variance curve. No new bank, curve refit or feedback of returns/variance into the Gaussian OU dynamics is added. Stock paths are recomputed from fresh primitive draws; the underlying exogenous factor construction is reusable.

Let the existing centered modes obey $dU=-\kappa Udt+\sqrt{2\kappa}dB$, $U_0=0$. Their half-lives are 3.329, 9.481 and 27.004 days. As before,

$$v_t=v^0_t+\frac f2\sum_j c_jh_j(t,U^j_t),\quad \sum_j|c_j|\le1,\quad f=0.0144,$$

$$q_j(t,u)=b_j(t)\phi_j(u),\qquad h_j=(\partial_t+\mathcal L_j)q_j=b'_j\phi_j+\kappa_jb_j(\phi_j''-u\phi_j').$$

Here $\phi$ is either $\tanh(u/s)$ or $\operatorname{erf}(u/(\sqrt2s))$, with widths $s=1,2$. The width changes how rapidly the response saturates as the Gaussian state moves away from zero; it adds no randomness.

## Local time scaling and positivity

Every quoted VIX-window endpoint is in $B=\{27,55,57,85,118,148,153,183\}/365$. On each successive interval $[l,r]$ take

$$a(t)=\sin^2(\pi(t-l)/(r-l)),\qquad A=\pi/(r-l).$$

Before the first endpoint, $a(t)=\cos^2(\pi t/(2T_1))$ with $A=\pi/(2T_1)$. After the last endpoint, use a 10-day $\sin^2$ ramp from 0 to 1, then keep it at 1, with $A=\pi/(2\tau)$ and $\tau=10/365$. These are boundary extensions of a single calendar construction; no expiry is given its own fitted factor or coefficient. It remains dependent on the finite quoted calendar, as the first test was.

For either shape $|\phi|\le1$ and $|\phi''-u\phi'|\le C_s$, with the analytic bounds

$$C_s^{\tanh}=\frac12+\frac{4}{3\sqrt3s^2},\qquad C_s^{\operatorname{erf}}=\sqrt{\frac2{\pi e}}(1+s^{-2}).$$

On each interval use the constant

$$K=\frac{\kappa C_s}{2}+\sqrt{A^2+\frac{\kappa^2C_s^2}{4}},\qquad b(t)=a(t)/K.$$

To see the bound, write $a=(1-\cos2\theta)/2$, $a'=A\sin2\theta$. For fixed $u$, the numerator $a'\phi+a\mathcal L\phi$ is a constant plus a sine and cosine. Its absolute maximum is at most $|\mathcal L\phi|/2+\sqrt{A^2\phi^2+(\mathcal L\phi)^2/4}$, which is at most $K$. The cosine boundary segment has the same bound, and the flat extension satisfies it too. Thus $|h_j|\le1$ globally and $v\ge f/2=0.0072$. This is a proved analytic bound under the stated OU generator, not a bound inferred from observed paths. No variance clipping is used. It preserves the same conservative positive floor as the first correction, rather than lowering it again.

Although $K$ changes between intervals, both $a$ and $a'$ vanish at their joins. Therefore $b$ is continuously differentiable, and the potential has no jump term. A large local scale cannot be inserted without this endpoint check. The two possible state shapes and their generators are implemented directly; there is no post-pricing smile smoothing.

## VIX protection and its limits

For every quoted expiry $T_i=27,55,118,153$days and $\Delta=30/365$, $q(T_i,u)=q(T_i+\Delta,u)=0$ for all states. Conditional Dynkin gives

$$E_{T_i}\int_{T_i}^{T_i+\Delta}h(s,U_s)ds=E_{T_i}q(T_i+\Delta,U_{T_i+\Delta})-q(T_i,U_{T_i})=0.$$

This protects model VIX values state by state, hence all quoted-date VIX options and forwards, for the unchanged full Gaussian filtration. It inherits the baseline VIX fit: loss approximately 0.05008 and maximum checking-forward gap 0.05517 points. No new claim of exact market/population forwards is made and no unchanged baseline exposure estimate is mislabelled as fresh simulation. Unquoted dates and the full VIX process law are not protected. Quoted-date joint VIX distributions are protected too, a stronger restriction than preserving only marginal distributions.

Both shapes are odd in the centered Gaussian state, so $E[h_t]=0$ and ideal $E[v_t]=\xi(t)$ is unchanged. Positivity and conditional cancellation are separate requirements; neither by itself guarantees a satisfactory SPX smile. No new continuous-stock martingale theorem or Lean proof is asserted.

## Selection and independent checking

The pilot uses 512 stock paths and 26 configurations: baseline, the previous selected correction, and 24 new choices from three OU modes, two state shapes, two widths and both signs, at the same maximum L1 budget. This is a discrete shape search, not full coefficient optimization. The best new choice by all-primary fine-grid SPX interval loss was frozen before 1,024 fresh stock paths: the error-function response on the existing approximately 9.5-day OU mode, with width 1. These target draws differ from the first experiment; all three checking configurations share draws for precise differences. The same existing 65,536-path auxiliary reference and market snapshot are reused, so this is not market out-of-sample validation.

All 1,751 eligible and 1,294 primary SPX quotes and seven expiries (29, 57, 119, 211, 302, 666 and 1,030 days) remain. All 217 VIX quotes at four dates are protected. The table is on the 256 steps per trading day grid. The 29-day wing gap is implied volatility near strike/forward1.1 minus that near0.9, in percentage points; market gap is -14.05pp.

| Model | SPX price loss | IV mean absolute error (pp) | 29-day wing gap (pp) | Signed negative primary prices |
|---|---:|---:|---:|---:|
| Baseline | 1945.52 | 3.209 | -5.36 | 0 |
| Previous protected correction | 1909.39 | 3.162 | -5.75 | 0 |
| Richer protected correction | 1884.31 | 3.137 | -5.85 | 0 |

The richer correction changes price loss by 3.15% relative to baseline and 1.31% relative to the previous correction. Approximate paired eight-batch bootstrap intervals for candidate-minus-reference loss are 

| Comparison | Loss change | Approximate 95% interval |
|---|---:|---:|
| Previous protected correction minus Baseline | -36.13 | [-51.24, -17.72] |
| Richer protected correction minus Baseline | -61.20 | [-92.38, -30.72] |
| Richer protected correction minus Previous protected correction | -25.07 | [-43.78, -8.51] |

Each uses 1,024 resamples of paired full-model batches with independent resampling of the common cheap-reference batches. They exclude deterministic grid error, fitting variability and market uncertainty. Signed negative prices are retained in the price score; unavailable IVs are reported and omitted only from the available-IV mean. IV means equally weight the seven maturities and are not a replacement for wing-by-wing inspection.

## Implementation and numerical evidence

Compiled/Python formula discrepancy 5.3e-15; conditional-window maximum residual 5.16e-10 across all 13 basis functions, four quoted dates, states -4, -1, +1, +4 and GH64/96; endpoint potential/derivative discrepancy 0; baseline and previous-correction replay discrepancy 0. Independent nonzero scalar stock recurrence error 0. An additional independent finite-difference generator check has maximum error 3.26e-08. These finite numerical checks support implementation correspondence; they do not establish a uniform error bound over all states.

Both 128- and 256-step stock grids are retained;

| Model | 128-step loss | 256-step loss | Median price change / quote scale | Maximum |
|---|---:|---:|---:|---:|
| Baseline | 1940.47 | 1945.52 | 0.089 | 1.565 |
| Previous protected correction | 1903.76 | 1909.39 | 0.086 | 1.471 |
| Richer protected correction | 1878.59 | 1884.31 | 0.082 | 1.594 |

The selected correction's GH96 conditional left-sum audit on actual grids and states -4, -1, 0, +1, +4 at all four dates has maximum correction-induced VIX-squared change 3.277e-05points squared. Exact cancellation belongs to the continuous-time model; the discretized correction has small nonzero error. Baseline VIX quadrature is a separate inherited error. Paired option-price errors, native-IV/stock-mean diagnostics by maturity, and raw tail payoff ESS/top-10 concentration are saved. No sample stock-mean renormalization was used.

## Research conclusion

The richer shape improves aggregate SPX fit in this checking sample, but the observed downward skew remains far from market. This test distinguishes richer state/time shape from adding another stochastic bank, with the same variance budget and positivity floor. The allowed corrections remain small compared with stressed baseline variance, and finite-date protection restricts their calendar shape. These are limitations of this tested family, not an impossibility theorem for M2 or all exogenous readouts. Do not repeat this exact shape grid as untested, and do not increase amplitudes past the established bound without a new positivity argument. Any next step should address the remaining wing-level mismatch explicitly rather than treating aggregate price improvement as a solved double calibration.

All seven SPX panels and the time-scaling figure were visually reviewed. The 29-day correction primarily lowers the right side; the large shortfall in downside implied volatility remains visible. The improvement does not supply the missing left-wing level.
