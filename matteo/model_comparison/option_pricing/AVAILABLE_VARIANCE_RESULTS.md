# Using existing stress variance as the protected-correction budget — 29 October 2026

## Question and construction

The earlier additive corrections were limited by a fixed variance budget f/2, regardless of how high current variance was. Here the budget grows with the existing stress contribution. This retains the full exogenous OU bank, old capped/Q-normalized M2 readout, positive spike shift, both stress factors, their mixture, stock driver and fitted forward-variance curve. It does not add actual-return or volatility feedback to the OU dynamics. There is no new bank or market-data refit.

Write the baseline as

$$v^0_t=f+c_M(t)R_t+c_s(t)A^s_t+c_f(t)A^f_t,$$

where $c_M=(\xi-f)(1-w)$, $c_s=(\xi-f)w(1-r)$, $c_f=(\xi-f)wr$, all nonnegative. Each stress follows $dZ=-aZdt+\sqrt{2a}dB$, $Z_0=0$, with $A(t,Z)=\exp(bZ-b^2Q(t)/2)$ and $Q(t)=1-e^{-2at}$. The selected pre-existing mode $U$ is driven independently of these two stress factors. Their stock correlations remain unchanged.

Simply setting h=current variance times the old correction would not generally retain cancellation. Instead scale the potential before applying the full generator. For one stress component define

$$y=Z-bQ(t),\quad g(t,Z)=e^{-\lambda y^2},\quad\phi(U)=\operatorname{erf}(U/\sqrt2),$$

$$P_j(t,U,Z)=c_j(t)\beta_j(t)\phi(U)A(t,Z)g(t,Z),\quad h_j=(\partial_t+\mathcal L)P_j.$$

The Gaussian factor g keeps the derivative terms bounded relative to available stress variance. Its center bQ(t) is the mean of Z under weighting by A, so it follows the states that contribute most to mean stress variance. It is a smooth readout weight, not a new state or clipping rule.

Put $D_\lambda(t,y)=a[(4\lambda^2+2\lambda)y^2-b(1+2\lambda Q(t))y-2\lambda]$. Direct differentiation gives

$$\frac{(\partial_t+\mathcal L_Z)(Ag)}{A}=gD_\lambda,$$

$$h_j=c_j A g\,[\beta'_j\phi+\beta_j\mathcal L_U\phi+\beta_j\phi D_\lambda],\qquad\mathcal L_U\phi=-2\kappa U\sqrt{2/\pi}e^{-U^2/2}.$$

There is no mixed U/Z derivative term because the existing OU mode and stress driver are independent. The date-dependent coefficients c_j are constant between fitted curve breakpoints. The potential and its first time derivative are set to zero at every curve breakpoint as well as every quoted VIX-window endpoint. Thus coefficient jumps introduce no omitted jump term. Where xi=f, c_s=c_f=0 and this correction is zero; this differs from the earlier f/2 correction, which was active even on floor intervals.

## Global positivity bound

We have $|\phi|\le1$, $|\mathcal L_U\phi|\le\kappa C$ with $C=2\sqrt{2/(\pi e)}$. Since $Q(t)\in[0,1]$ and b>0,

$$|gD_\lambda|\le D_*:=a\left[\frac{4\lambda+2}{e}+\frac{b(1+2\lambda)}{\sqrt{2\lambda e}}+2\lambda\right].$$

This uses $|y|e^{-\lambda y^2}\le1/\sqrt{2\lambda e}$ and $y^2e^{-\lambda y^2}\le1/(\lambda e)$. On each interval between endpoints use the same sine-squared clock as in the preceding experiment, with $a_0=\sin^2(\pi(t-l)/(r-l))$ and $A_0=\pi/(r-l)$. Before the first endpoint use a cosine-squared boundary segment. After the last use a 10-day sine-squared ramp then a flat extension. Set

$$B_j=\kappa C+D_{*,j},\quad K_j=B_j/2+\sqrt{A_0^2+B_j^2/4},\quad\beta_j=a_0/K_j.$$

The same sine/cosine amplitude bound used previously now gives $|h_j|\le c_jA^j$. Consequently, with

$$v_t=v^0_t+\eta(h_s+h_f),\qquad |\eta|\le1,$$

$$v_t\ge f+c_MR+(1-|\eta|)(c_sA^s+c_fA^f)\ge f.$$

This retains the ORIGINAL floor f=0.0144, stronger than the f/2 guarantee in the previous experiment. It can exceed the old absolute correction budget f/2 in high-variance states. The inequality is analytic, not inferred from simulated minima. It covers this stress-scaled construction, not arbitrary multiplication of corrections by variance or the entire M2 component. No variance clipping is used. Gaussian damping makes the relevant potentials and derivatives integrable on finite horizons.

## Conditional VIX protection and mean variance

All potentials vanish at both ends of all four quoted VIX windows. For the unchanged full Gaussian filtration,

$$E_{T_i}\int_{T_i}^{T_i+\Delta}(h_s+h_f)dt=0,\quad T_i=27,55,118,153\text{ days},\quad\Delta=30/365.$$

Conditional Dynkin's identity therefore protects model VIX values state by state, hence all 217 quoted VIX options and four forwards at these dates. Baseline VIX loss approximately 0.05008 and max checking-forward gap 0.05517 points are inherited, not newly corrected or certified zero. Unquoted dates and the full VIX process law are not protected. Quoted-date joint VIX distributions are unchanged too. Oddness in the centered Gaussian U and independence from stresses give E[h]=0, retaining ideal E[v_t]=xi(t). Finite numerical normalization errors of the original M2 remain.

## Numerical translation and evidence scope

The compiled formula is compared against independent Python algebra, finite differences of the full time/U/Z potential, and a Gaussian conditional-expectation oracle. For the latter, exponential tilting and Gaussian damping give E[Ag] and E[AgD] analytically; E[erf(U/sqrt2)] and E[L_U phi] are analytic too. The product of these conditional moments is checked against compiled two-dimensional GH64 quadrature. Time integration splits at every curve endpoint and tests all three existing modes, both damping values, both stress factors, four VIX dates and four pairs of initial U/Z states. The actual fitted curve coefficient is included. Conditional residuals are scaled by initial stress and maximum curve level, since raw errors at extreme stress states scale with A. The analytic identity supplies the general cancellation claim; finite probes do not certify a uniform numerical bound.

The new construction is compared with the previously selected rich correction. This is not an ablation changing only amplitude: the calendar now includes every curve breakpoint, the stress readout is damped, and the new family retains the stricter original floor. Keep these differences when interpreting the outcome. No new Lean theorem, continuous-stock martingale theorem or universal model-family impossibility is asserted.

## Experiment and outcome

512 pilot paths test 38 configurations: baseline, the previous rich correction, and 36 new choices from three existing modes, damping lambda 0.125/0.5, strengths 0.25/0.5/1 and both signs. The best new primary-SPX interval-loss choice was frozen before 1,024 fresh stock paths: half-life3.329days, lambda0.125, eta1.0. All 1,751 eligible and 1,294 primary SPX quotes and seven expiries are retained. Signed prices and invalid IVs stay in the artifacts. The auxiliary 65,536-path reference and market snapshot are reused, so this is not market out-of-sample validation. New target draws are independent of previous experiments; all checking candidates share these draws.

Results on 256 steps per trading day; left/right IVs are near strike/forward 0.9/1.1 at 29 days. The market levels are26.90% and12.85%, giving a gap-14.05pp.

| Model | SPX price loss | IV MAE (pp) | 29d left IV (%) | 29d right IV (%) | Gap (pp) |
|---|---:|---:|---:|---:|---:|
| Baseline | 1948.66 | 3.251 | 18.97 | 13.25 | -5.72 |
| Previous rich correction | 1888.60 | 3.177 | 19.13 | 12.44 | -6.69 |
| Stress-budget correction | 1925.31 | 3.211 | 19.26 | 13.12 | -6.14 |

Price-loss reduction versus baseline:1.20%; versus previous rich correction:-1.94%. Approximate 1,024-resample eight-batch paired bootstrap intervals are below; they resample full-model batches together across candidates and resample cheap-reference batches independently. They exclude deterministic error, fitting variability and market uncertainty.

| Comparison | Loss change | Approximate 95% interval |
|---|---:|---:|
| Previous rich correction minus Baseline | -60.06 | [-86.44, -35.81] |
| Stress-budget correction minus Baseline | -23.35 | [-120.12, 76.59] |
| Stress-budget correction minus Previous rich correction | 36.70 | [-46.34, 128.55] |

Compiled/Python relative formula error5.14e-15; full-potential finite-difference generator error scaled by 1+A:5.33e-08; conditional scaled-integral error3.41e-13; analytic conditional oracle versus compiled GH64 error4.76e-16; baseline/previous-rich replay error0; independent signed stock recurrence error0. The scalar recurrence uses a synthetic clock to cover all time features; it is a finite correspondence audit, not a separate pricing model.

Both 128/256 stock grids and all per-maturity native-IV/stock-mean diagnostics are saved in `visualizations/available-variance029/summary.json` and `maturity_metrics.csv`. The selected correction's conditional left-sum audit on actual grids is `runs/m2_available_variance_001/discrete_window_audit.json`; maximum error relative to initial stress and maximal curve budget is2.61e-08years. Exact cancellation belongs to the continuous-time model. Discrete extreme-state errors must be assessed relative to stress size rather than presented as a universal absolute bound. Direct-Q variance-budget probes and raw option-tail payoff concentration/ESS are saved; their finite samples do not certify rare-tail convergence. No sample stock-mean renormalization or smile smoothing is used.

## Conclusion

The main mathematical answer is yes: this explicit correction can use the current stress variance while preserving quoted-date VIX and the original positive floor. The numerical table assesses whether this particular construction supplies the missing SPX shape; do not identify feasibility with calibration success. Both the remaining downside IV shortfall and any wing deterioration must be read from the plots, not just the aggregate price score. No candidate is adopted. Preserve this search rather than rerunning it as untested. Future changes must retain the full-generator terms, the curve-jump endpoint conditions and the positivity argument.

[Seven SPX smiles](visualizations/available-variance029/SPX-available-variance.png). Full artifacts: `runs/m2_available_variance_001`.

## Decision from this test

The new correction has a 1.20% lower sampled price loss than baseline, but its paired loss-change interval includes zero: [-120.12, 76.59]. Its sampled loss is higher than the previous rich correction; that paired interval also includes zero: [-46.34, 128.55]. This sample does not establish an improvement over either comparator. At 29 days, left-wing IV rises from 18.97% to 19.26%, while the market is 26.90%. The much larger downside-volatility shortfall remains.

The construction actually exceeds the old absolute variance-correction limit in some states: the day7 direct-Q probe has a 99th percentile 0.02714 versus the old fixed limit 0.0072. Probes at curve endpoints are zero by construction and cannot describe typical interior amplitudes. These finite probes do not certify rare-tail convergence.

The feasibility question has a constructive positive answer, but this numerical specification is not a successful calibration. No candidate is adopted. All seven SPX panels were visually reviewed.
