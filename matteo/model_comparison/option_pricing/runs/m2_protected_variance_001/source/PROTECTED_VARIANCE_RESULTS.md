# VIX-protected additive M2 variance correction — 29 October 2026

The first test produces a modest SPX improvement while analytically protecting the four quoted VIX expiry values. It does not solve the missing SPX skew. No candidate is adopted; the production model and all previous runs are unchanged. Figures: [all seven SPX smiles](visualizations/protected-variance029/SPX-protected-variance.png), [correction envelope](visualizations/protected-variance029/correction-envelope.png). Full artifacts: `runs/m2_protected_variance_001`.

## What changed, in plain language

Keep the existing M2 component, both stress factors, their mixture, cap, normalizers, positive spike shift, forward-variance curve and stock driver. Add a small signed correction directly to variance. It can redistribute volatility along stock paths, but its conditional expected total over each quoted VIX window is zero. This resembles moving water within a container: its distribution can change while its expected total remains fixed. That cancellation is conditional on the entire current model state, stronger than cancelling only in the unconditional average.

Let the baseline variance be $v^0_t$. Existing centered OU modes satisfy

$$dU^j_t=-\kappa_j U^j_tdt+\sqrt{2\kappa_j}dB^j_t,\quad U^j_0=0.$$

Use the same geometric feedback bank, modes 7,6,5 with half-lives [3.3291252855771765, 9.48146822384038, 27.0035615869317]. There is no new bank or feedback of returns or variance into these dynamics. Stock paths are recomputed; Gaussian factors/random draws remain reusable.

$$v_t=v^0_t+\frac f2\sum_j c_jh_j(t,U^j_t),\qquad \sum_j|c_j|\le1,\quad f=0.0144.$$

$$a(t)=\prod_{b\in B}\tanh^2((t-b)/\tau),\quad B=\{27,55,57,85,118,148,153,183\}/365,\quad \tau=10/365.$$

$$q_j(t,u)=a(t)\tanh(u)/K_j,\quad K_j=1/\tau+\kappa_j[1/2+4/(3\sqrt3)],$$

$$h_j=(\partial_t+\mathcal L_j)q_j=\frac{a'\tanh u-\kappa_j a(u+2\tanh u)\operatorname{sech}^2u}{K_j}.$$

The common calendar envelope treats every quoted expiry by the same rule. This is deliberately an initial finite-expiry protection test, not a universal calendar-independent VIX model. Smoothness here does not mean a small effect is guaranteed or that smile fit follows automatically.

## Why VIX is protected

For each quoted $T_i=27,55,118,153$ days and $\Delta=30/365$, both endpoint functions vanish. Conditional Dynkin's identity gives

$$E_{T_i}\int_{T_i}^{T_i+\Delta}h_j(s,U^j_s)ds
=E_{T_i}q_j(T_i+\Delta,U^j_{T_i+\Delta})-q_j(T_i,U^j_{T_i})=0.$$

Thus $\mathrm{VIX}^2_{T_i}=10000\Delta^{-1}E_{T_i}\int_{T_i}^{T_i+\Delta}v_sds$ is unchanged state by state for the unchanged full Gaussian filtration. All VIX strikes and forwards at these four dates inherit their baseline values. Baseline VIX loss is approximately0.05008, with maximum checking-forward gap0.05517points: protection does not newly make those baseline market errors zero. Baseline VIX estimates are reused under this identity, not presented as a new independent calibration. Unquoted dates and the full VIX process law are not protected. At the protected dates even their joint law is unchanged, so this is more restrictive than the paper's approach of fixing only expiry marginal distributions.

Oddness of $h_j$ in centered Gaussian $U^j$ also preserves $E[v_t]=\xi(t)$ in the ideal normalized baseline. A global bound $|h_j|\le1$ and the amplitude budget imply $v_t\ge f/2=0.0072$: there is no clipping. This retains strict positivity but lowers the guaranteed floor from the original $f$; it does not retain the original exact hard floor. Since the correction is additive, it is not multiplied by $\xi-f$ and is active even on baseline floor intervals.

## Actual experiment and results

The pilot evaluated baseline plus twelve axis configurations: each of three modes, both signs, and strengths0.5/1 under the budget. There are512pilot stock paths. The minimum all-primary SPX interval-loss configuration was frozen before1024fresh stock paths: $c=(0,0,1)$, the27-day mode at its allowed amplitude. No curve refit or VIX-forward fit was performed. This is a small discrete parameter search, not optimization over all combinations or a proof of optimality.

Keep all1751eligible SPX quotes,1294primary quotes and seven expiries29/57/119/211/302/666/1030days. Four VIX expiries and217quotes remain protected. Signed prices and all native-IV failures remain in artifacts. The table uses fresh paths and256steps per trading day. The last column is IV near strike/forward1.1 minus IV near0.9 at29days, in percentage points: more negative means a stronger downward skew.

| Model | SPX price loss | IV mean absolute error (pp) | 29-day right minus left IV (pp) |
|---|---:|---:|---:|
| Baseline | 1841.07 | 3.169 | -6.00 |
| Protected correction | 1801.18 | 3.121 | -6.59 |

The market29-day gap is-14.05pp. Price loss improves2.17%. The approximate paired1024-resample eight-batch95% bootstrap interval for candidate-minus-baseline loss is[-51.17,-28.63]. It resamples shared full-model batches and independently resamples the shared65536-path auxiliary-reference batches. The new target simulations are fresh; the auxiliary reference and market snapshot are reused. These are sampling diagnostics, not a market out-of-sample certificate, and exclude deterministic error and fitting uncertainty.

The result separates two questions: this construction can change SPX prices while protecting quoted-date VIX; this first conservative correction is too weak to recover the observed SPX skew. A small aggregate loss improvement is not a satisfactory double calibration. The conservative bound constrains amplitude globally even where the baseline variance has more room. Also, the common envelope is suppressed near the densely spaced early VIX endpoints. Those are plausible limitations of this construction, not numerically established explanations or a no-go theorem for exogenous M2.

## Verification and limits

Compiled/Python formula maximum discrepancy5.52e-14. Conditional cancellation probes cover three modes, four dates, five conditioning states and GH64/96; maximum integral residual2.72e-10. Zero-amplitude compiled stock replay error0. Independent scalar recurrence on a97-step synthetic clock spanning envelope features, five configurations, both antithetic partners and both grids gives maximum stock discrepancy0. These are numerical correspondence checks, alongside the stated analytic identity, with no new Lean proof or continuous-stock martingale theorem.

Both128/256stock grids are retained. Candidate-minus-baseline loss changes are -39.96 on the coarse grid and -39.89 on the fine grid, so the observed benefit is stable in this paired refinement, despite nonzero absolute grid error. Grid diagnostics are in `visualizations/protected-variance029/summary.json`. A separate GH96 conditional left-sum audit on both actual grids, selected mode and five states at all four dates is in `discrete_window_audit.json`; its largest correction-induced VIX-squared change is3.17054e-05VIX-points squared. Exact cancellation is a continuous-time claim; the grid sum has small nonzero error. This audit concerns the new correction, not the baseline VIX quadrature. No finite probe establishes a uniform discretization bound for every state.

Per-maturity native-IV/skew/stock-mean diagnostics, paired price errors and raw tail payoff concentration are saved. Both fresh fine-grid candidates have five signed negative primary option estimates with unavailable implied volatilities; these remain in the data and price loss, and are excluded only from the available-IV average. All seven SPX panels and the correction-envelope figure were visually reviewed. No stock-mean renormalization or post-pricing smile smoothing was used. The old valid conditional-lognormal auxiliary reference and fixed deterministic delta controls remain unchanged.

## Research conclusion

Retain the protected-correction idea as a demonstrated way to separate part of SPX path dependence from four quoted-date VIX constraints, but do not adopt this first bounded-tanh specification as a solution. Do not repeat this exact mode/strength grid as untested. Before adding banks, a useful next design would investigate better use of the available positivity budget and less restrictive time/state bases while preserving the explicit conditional cancellation. Any enlarged correction requires a new positivity argument and implementation checks; simply increasing the coefficient beyond the tested bound is not justified by this experiment.
