# Protecting the successful VIX fit while improving SPX pricing

28 October 2026. This continues the user's authorized two-speed SPX/VIX investigation. The numerical experiments below keep the successful VIX-only model fixed. They do not adopt the unsatisfactory joint candidate, change correlations, add factors or refit option prices.

## What stays fixed

The variance remains

$$v_t=f+(\xi_t-f)[(1-w)R_t+w\{(1-r)A_t^s+rA_t^f\}].$$

The slow/fast stress factors, original M2 response and normalizer, positive spike shift, risk premia, initial state, floor and all16 variance levels retain the successful VIX-only values. In particular, stress weight is0.922838 and the volatility floor is12%. The model artifact in `runs/m2_two_speed_001/experimental_model_two_speed.json` is unchanged. Its existing VIX fit and forward accuracy are therefore unchanged; this does not upgrade sampled forward matching into exact population equality.

All seven SPX maturities29/57/119/211/302/666/1030days remain active. The full experiment includes all1751 eligible native SPX quotes, with1294 primary quotes under the agreed cutoff3. All217 VIX quotes and four forward targets remain the protected VIX scope. Signed SPX estimates and undefined implied volatilities are retained and flagged, with no clipping or smoothing.

## A simpler reference helps price the actual model

The new stress factors are independent of the stock's Brownian driver. An auxiliary model replaces the old M2 response by1:

$$\bar v_t=f+(\xi_t-f)[(1-w)+w\bar H_t].$$

The auxiliary stress is held on a daily-ish subset of the existing stock grid, split additionally at each represented variance-curve knot. The actual model is still simulated at64steps per trading day. A cheap reference can use its own coarser rule because we subtract its sampled payoff and add back its expectation. It does not replace the actual variance process.

Conditional on the auxiliary stress path, its variance curve is known and its option expectation is an ordinary Black price. Write its total variance as I, its normalized stock as Sbar, the actual normalized stock as S, and the native-call payoff as g. For the old-path likelihood L,

$$E_Q[g(S)]=E_{stress}[B(1,K,I)]+E_{proposal,stress}[L\{g(S)-g(\bar S)\}].$$

In plain language: accurately price most of the option using a simpler model, then simulate only the difference between the actual model and that reference. The stress-only Black expectation is computed separately on new cheap paths. It is not multiplied by the old-path importance weight. The expensive correction retains the full old-path likelihood.

The first pilot also tests terminal-delta and daily dynamic hedge controls. These have mean zero conditional on the entire independent stress path. Their expectation is zero even though they need not replicate the payoff exactly. See [TWO_SPEED_PRECISION_PLAN.md](TWO_SPEED_PRECISION_PLAN.md) and `runs/m2_two_speed_precision_002`.

## A stronger, one-step conditional control

We then integrate each stock step analytically. This is an expectation identity for the existing discrete stock law; it does not change its time grid or approximate the actual model by the reference.

Let J_j be remaining reference integrated variance at stock date j, and define

$$U_j(s)=B(s,K,J_j),\qquad J_j-J_{j+1}=\bar v_j\Delta t_j.$$

Given the current full M2 state and the independent stress path, the next actual stock step is lognormal. Therefore

$$E[U_{j+1}(S_{j+1})\mid\mathcal F_j,stress]=B(S_j,K,J_{j+1}+v_j\Delta t_j).$$

Define the correction

$$D_j=B(S_j,K,J_{j+1}+v_j\Delta t_j)-B(S_j,K,J_j).$$

Telescoping the Black values and taking expectations gives

$$E_Q[g(S_T)]=E_{stress}[B(1,K,J_0)]+E_{proposal,stress}\left[L\sum_jD_j\right].$$

Thus the noisy final payoff is replaced by a sum of small conditional-price differences. The calculation still accounts for the actual M2/equity dependence. When actual and auxiliary variances agree, every correction is zero. Conditional on the entire independent stress path, the old equity driver retains its Brownian law, so this conditioning is legitimate. The argument would require review if stress/equity correlations were introduced later.

The cheap reference additionally uses eight independent scrambled Sobol nets with an OU bridge. The bridge generates large-scale stress movements first, then fills the intermediate dates with their exact conditional Gaussian distributions. Each net has4096points, giving32768reference paths per expiry and role. Uncertainty uses independent net means, not individual Sobol points.

## Numerical results

The generated tables and plots below compare the same frozen model. Lower pricing noise does not mean the model fits the market better.

| Days | Previous median SE / scale | New median SE / scale | Reduction | Expensive component | Cheap component |
|---:|---:|---:|---:|---:|---:|
| 29 | 4.681 | 0.579 | 8.1x | 0.574 | 0.088 |
| 57 | 15.029 | 0.550 | 27.3x | 0.506 | 0.201 |
| 119 | 12.173 | 0.903 | 13.5x | 0.854 | 0.297 |
| 211 | 12.756 | 1.153 | 11.1x | 0.939 | 0.646 |
| 302 | 13.082 | 1.239 | 10.6x | 1.192 | 0.320 |
| 666 | 5.738 | 0.435 | 13.2x | 0.372 | 0.204 |
| 1030 | 5.005 | 0.174 | 28.8x | 0.162 | 0.059 |

All figures are checking-role results for the same frozen VIX-only model. The estimated median SE reduction ranges from8.1x to28.8x across maturities. That is a reduction in SE, not in variance or model error. Only1/1294primary quotes passes the0.125SE/scale screen. Negative native checking estimates:0/1751. No option optimizer was launched.

[Precision comparison](visualizations/two-speed-precision027/precision-comparison.png) and [all seven native SPX smiles](visualizations/two-speed-precision027/frozen-model-spx-smiles.png). The corrected IV plot uses percentage points directly from the existing inversion routine; it does not multiply percentages by100 again. Its CSV retains every primary quote and flags invalid model IVs.

The clearer SPX plots still show a substantially more U-shaped model smile than the observed downward skew. This is a fixed-parameter discrepancy; it is not a proof that refitting the nine allowed curve entries cannot help, or that this two-speed model family is infeasible. The small plotted sampling bars do not explain the broad current smile mismatch.

The expensive correction is the main remaining sampling bottleneck. Taking only the cheap reference to infinite precision would leave the median expensive-component SE above the screen at every maturity. A simple equal-variance-budget projection gives:

| Days | Projected median expensive paths | Projected worst-quote expensive paths |
|---:|---:|---:|
| 29 | 21,632 | 158,080 |
| 57 | 16,832 | 248,832 |
| 119 | 47,872 | 1,203,456 |
| 211 | 57,856 | 231,552 |
| 302 | 93,184 | 416,896 |
| 666 | 9,088 | 25,728 |
| 1030 | 1,792 | 10,496 |

These are pilot-based sqrt(N) projections using the same expensive estimator/proposal, not guaranteed requirements or achieved precision. Each component is assigned SE0.125/sqrt(2). For the cheap reference, the companion CSV projects more independently scrambled4096-point nets; it does not assume a convergence rate when enlarging one Sobol net. The worst-quote expensive projection is about1.2million paths, versus512 currently stored. This makes blind expansion of the full fine-path cache unattractive: use streaming, adaptive allocation and a coupled multilevel estimator before committing to that cost. The projections exclude timestep and normalizer bias and can themselves be unstable for concentrated tails.

Primary native interval loss (equal maturity weight) of the frozen model with the new estimator: 6304.131. This is a diagnostic of the unchanged model, not a new calibration result.



## What the tests establish

Independent checks cover replay of the original target stock, a scalar Python reconstruction of the auxiliary stock and its integrated variance, all16 curve derivative directions on a1030day probe that actually visits all16 intervals, exact covariance reconstruction of both OU bridges, an independent Python conditional-price correction, importance-weighted finite-law control-mean checks, and replay of one stored reference scramble at every expiry and role. Source/input and artifact hashes are checked separately.

These are implementation and finite numerical checks. No Lean theorem or uniform price-error certificate is supplied. The existing expensive search/check streams are reused, so this is not untouched simulation validation or market out-of-sample testing. Eight-block sampling errors omit normalizer error and stock timestep error. The earlier coupled64/128/512stock study left timestep precision unresolved; the present control identities do not remove that bias.

The preliminary run `m2_two_speed_precision_001` stopped in preflight: its identical-model fixture incorrectly combined fine target stress with held reference stress. The fixture was corrected to use the same stress grid in both models; no price experiment completed in that failed preflight. Run002 preserves the corrected source and passed checks. The35-quote conditional pilot is preserved as `m2_two_speed_precision_conditional_001`; its full-strike successor is002. A read-only first-order integrated-old-variance control probe did not consistently reduce residual SE and was not adopted.

## Next permitted calibration

The nine off-window curve indices are0,5,6,9,11,12,13,14,15. Varying only these entries leaves the exact VIX parameter mapping unchanged in the current exogenous model. All stress parameters, the floor and VIX-sensitive levels remain fixed. That preserves the current VIX law, not merely its forward means.

The primary sampling screen is $2SE\le price\_scale/4$, or $SE/price\_scale\le0.125$. A price optimizer remains disabled until the numerical gates are met. Further work should stream fresh full-M2 paths for the expensive correction, with independently checked batches and a coupled stock-grid refinement. More cheap stress paths alone cannot remove the expensive correction's remaining uncertainty. Then fit only the nine protected curve levels, check all seven SPX maturities and retain the VIX forward constraints and smile protection for any broader joint experiment. Stress/equity correlations remain a subsequent proposed experiment, not an implemented change.

For restarting, use the current top entries of [CONTINUE_FITTING_PROMPT.md](CONTINUE_FITTING_PROMPT.md), [NEXT_STEPS.md](NEXT_STEPS.md), and [RESEARCH_LEDGER.md](RESEARCH_LEDGER.md). A new chat is optional; the repository records preserve the scientific state.
