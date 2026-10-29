# SPX skew, VIX consistency and fast exogenous-model experiments

28 October 2026. Literature and source review, with a restricted new Lean check. This note incorporates the user's explicit constraint: **keep the OU/state dynamics exogenous; do not introduce dependence of their coefficients on volatility or simulated returns.** Preserve the existing Zumbach-motivated architecture and reuse simulated Gaussian banks. Othmane's endogenous bank is studied as an explanation of his results, not proposed for adoption.

No new stock-driver parameter, pricing engine, calibration or ablation sweep is implemented in this review. The new Lean results below check algebraic bridges for the proposed experiment, not a completed stochastic-model implementation or market fit.

## 1. The mechanism to target

A downward SPX skew means that downside options have higher implied volatility than comparable upside options. A standard mechanism links falling stock prices with increasing future variance. In a diffusion with stock loading sqrt(v), variance loading a(v) and Brownian correlation rho, the instantaneous stock-log/variance covariance is rho*sqrt(v)*a(v). With a nonnegative variance loading, negative rho contributes negatively.

The relation between this dependence and ATM implied-volatility skew is quantified by small-vol-of-vol expansions; it is not a theorem that every finite-maturity wing is monotone in rho. [Guyon, revisiting the Bergomi–Guyon expansion](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3956786) gives the general framework. [Abi Jaber–Illand–Li, section 2](https://arxiv.org/html/2212.08297v2) explicitly connects first-order ATM skew to integrated spot/forward-variance covariance in their Gaussian polynomial class. The same paper shows that kernel choices can obstruct joint fits even at very negative correlation. More negative correlation is therefore a targeted experiment, not a guaranteed solution.

Our successful VIX fit puts weight w=0.922838 on stress responses independent of the old equity shock. That is a mixture coefficient, not a statement that exactly 92% of variance comes from stress on every path. The old M2 component still supplies leverage. The new stress shocks have no direct stock-shock covariance. The tested curve-only fit lowers SPX levels but does not supply the required downward skew; see [the durable empirical conclusion](PROTECTED_CURVE_CONCLUSION.md). Capturing a Zumbach-type time asymmetry and matching an option-smile skew are different requirements.

## 2. Relevant literature and what it establishes

| Work | Relevant evidence | Scope/limit |
|---|---|---|
| [Abi Jaber, Illand and Li, Gaussian polynomial volatility (2025)](https://arxiv.org/html/2212.08297v2) | Joint SPX/VIX calibration with exogenous Gaussian factors; kernel/readout comparison; quantization and neural pricing. | The extensive daily-data experiment emphasizes short expiries. Their empirical success does not certify our capped readout or our longer maturities. |
| [Abi Jaber and Li, two-factor Quintic OU (2026 revision)](https://arxiv.org/html/2503.14158v2) | Two OU speeds driven by a common Brownian motion; reported joint SPX/VIX and smile-dynamics results beyond two years. | This is evidence for exogenous shared-driver constructions. Its polynomial readout is not introduced into our model. |
| [Fouque and Saporito, stochastic vol-of-vol](https://arxiv.org/html/1706.00873v1) | A stochastic-vol-of-vol extension of Heston with joint calibration examples and efficient perturbative pricing. | Prices are approximations with stated asymptotic scope; adding random vol-of-vol alone does not guarantee our skew. |
| [Gazzani and Guyon, four-factor PDV pricing/calibration (2025)](https://arxiv.org/abs/2406.02319) | An endogenous return-history model, with a learned conditional VIX map for joint calibration. | Useful computational precedent; its return/volatility feedback is outside our current model scope. |
| [Baschetti, Bormetti and Rossi, joint deep calibration (2025)](https://papers.ssrn.com/sol3/Delivery.cfm/5349810.pdf?abstractid=5349810&mirid=1) | A trained pricing surrogate makes repeated PDV calibration evaluations inexpensive. | Offline simulation/training and independent accuracy checks remain necessary. It does not make an incorrect model correct. |

The practical lesson is to separate factor-distribution/readout flexibility, memory speeds and spot/variance coupling. Changing one mechanism can be more informative than adding several factors simultaneously. Quintic readouts, jumps and endogenous feedback remain unadopted.

## 3. What Othmane actually did

The relevant option report is the 19 October 2026 `cs_doc`, not the historical/P calibration in `gs_doc` or `jc_doc`. Its equations and material numbers were rechecked against page images 3,7,8. The earlier 14 October email figures are a different evidence item; an exact matching implementation/configuration for those figures has not been established here. No independent repricing of Othmane's October implementation is claimed.

His volatility is

$$\sigma_t=\widetilde\sigma(t)e^{\psi(\zeta_t)}F(P_t).$$

The main SPX-skew mechanisms are:

- A rough-bank shock direction strongly correlated with the stock shock, fitted correlation 0.9264.
- A decreasing rough-factor readout psi: negative factor moves raise volatility. The positive quoted correlation is therefore compatible with negative stock/volatility leverage; sign conventions must be read together with the readout.
- A decreasing, asymmetric feedback readout F of filtered returns. The report gives F(-8)=12.2 and roughly 0.21 after large gains; it attributes much of SPX skew to this feedback and leverage.
- Fitted initial states and deterministic volatility term structure.

His feedback bank contains

$$ds_t^j=-\kappa_js_t^jdt+\sqrt{2\kappa_j}\frac{\sigma_t}{\widetilde\sigma(t)}dW_t^X.$$

This is precisely the type of endogenous coefficient excluded by the user's current instruction. Changing its readout generally changes later s states, so complete state paths cannot simply be held fixed during such comparisons.

For VIX flexibility he adds Gamma to the rough-bank diffusion:

$$dr_t^i=-\lambda_ir_t^idt+\sqrt{2\lambda_i\Gamma_t}\,dW_t^i.$$

Gamma is a positive mean-reverting jump modulation, not a final independent volatility multiplier. Its fitted initial level is below its long-run mean. The modulation dynamics were fixed from the historical architecture selection. This is a different mechanism from our mean-one stress readout.

His reported SPX RMSEs are 0.63/0.49/0.51/1.12 volatility points across four expiries. VIX near-ATM fits are also useful. However, the report explicitly retains too-flat VIX right wings at 2–5 months; all-strike VIX RMSEs are 7.3/10.1/4.7/4.1 points. The VIX estimator uses a regression approximation, with a nested Monte Carlo check left pending. He fits a narrower VIX strike range with strong near-ATM weights, and fewer SPX expiries than our present scope. Those differences prevent a direct numerical ranking against our full target. See [page7](../othmane_docs_20261026.prepared/cs_doc/pages/page-0007.md), [page8](../othmane_docs_20261026.prepared/cs_doc/pages/page-0008.md), and [limitations, page13](../othmane_docs_20261026.prepared/cs_doc/pages/page-0013.md).

## 4. A compatible stock-shock rotation

Keep every OU equation and the entire joint variance-factor law unchanged. Let W 0 be the old equity/M2 driver, and let Bs,Bf be the existing independent slow/fast stress drivers. With the existing r fixed, define

$$dB^*=\sqrt{1-r}\,dB^s+\sqrt r\,dB^f,$$

$$dW^{stock}_\theta=\sqrt{1-\theta^2}\,dW^0-\theta dB^*,\qquad0\le\theta<1.$$

Only the stock equation changes. The old M2 and stress paths remain the same. This introduces statistical dependence between stock and stress shocks without introducing volatility/return feedback into any OU equation. Old stock/M2 correlations weaken by sqrt(1-theta^2); this is an intended change, not a claim that every stock correlation is held fixed.

For the current stress exponentials, their diffusion loadings in variance are

$$k_s=(\xi_t-f)w(1-r)b_s\sqrt{2a_s}A_t^s,\qquad k_f=(\xi_t-f)wrb_f\sqrt{2a_f}A_t^f.$$

They are nonnegative under the existing parameter restrictions. The new stress contribution to instantaneous log-stock/variance covariance is

$$-\theta\sqrt{v_t}\{\sqrt{1-r}\,k_s+\sqrt r\,k_f\}\,dt\le0.$$

The old M2 covariance contribution is multiplied by sqrt(1-theta^2). This calculation identifies the proposed mechanism; it does not prove that the complete finite-maturity implied skew is monotone in theta.

Because the full joint Q law of the exogenous factors and the variance readout is unchanged, conditional future variance—and hence the VIX distribution and its derivatives—remains unchanged. This is a model-level conditional-law argument, not just an unconditional-moment argument. Integrability, the continuous stock martingale property, joint Gaussian implementation and numerical accuracy remain separate obligations. See [the detailed proposal](STOCK_DRIVER_ROTATION_PROPOSAL.md).

## 5. Fast response scans without resimulating OU paths

For fixed variance/readout parameters, retain at every expiry only

$$I_T=\sum_jv_j\Delta t_j,\quad U_T=\sum_j\sqrt{v_j}\Delta W_j^0,\quad V_T=\sum_j\sqrt{v_j}\Delta B_j^*.$$

The discrete forward-normalized log-stock recursion then gives exactly, over real arithmetic,

$$\log(S_T^\theta/S_0)=-\tfrac12I_T+\sqrt{1-\theta^2}U_T-\theta V_T.$$

Thus a whole theta grid reuses three stored numbers per path/expiry. For slow-only, fast-only and mixed leverage directions, store the Bs and Bf integrals separately: four numbers in total. This avoids both OU resimulation and repeated time-stepping for a fixed-variance theta/direction scan. It still requires payoff evaluation, uncertainty analysis and a valid joint Brownian/OU sample.

The existing normalized OU innovations are not themselves Brownian increments. Complete each stress increment with its correctly correlated Brownian component before forming V. Cached OU states/innovations can be retained; new independent conditional Gaussian residuals supply the missing Brownian information. The proposal records the exact covariance formula. Reusing a changed-law control estimator is not justified: current independent-stress conditional-delta arguments must be replaced or rederived for the rotation.

If xi, weights or readout coefficients change, rebuild v and these downstream sums from the same OU paths. If mean-reversion rates or OU diffusion dynamics change, the old state paths are generally not reusable as-is. To retain a fixed bank during optimization, keep its rates and state dynamics fixed and fit readouts/couplings. The current experiment keeps both stress speeds fixed.

## 6. A minimal ablation programme

First freeze all 16 variance levels at the protected candidate and all factor/readout parameters. Do not start with another large joint optimization.

| Variant | Change | Variance-factor/VIX law |
|---|---|---|
| Baseline | theta=0 | Unchanged |
| Mixed leverage scan | theta in 0,.2,.4,.6,.8,.95; direction from fixed r | Unchanged |
| Fast-only leverage | Same theta grid, couple stock only to Bf | Unchanged |
| Slow-only leverage | Same theta grid, couple stock only to Bs | Unchanged |

These ablations remove a leverage channel, not the variance factor itself. Removing a variance factor would generally change VIX and must be a separate experiment with explicit renormalization and full VIX checking.

For all seven retained SPX maturities report native bid–ask price loss, native-IV errors, predeclared left/right quoted-strike IV gaps, invalid IV flags, paired uncertainty, raw weighted payoff concentration and stock-grid sensitivity. Keep all native quotes, signs and existing scope. Do not select a model from price loss alone. Paired sweeps isolate responses at fixed parameters; a later refit measures how much remaining flexibility can compensate. These answer different questions.

Only if the scan shows a useful direction should theta and the nine protected xi entries be optimized, followed by fresh independent checking. If the skew term structure requires different speed contributions, evaluate a second direction parameter only after the fast/slow ablations—not another stress speed by default. No endogenous feedback is proposed.

## 7. Lean: what is now checked and what is not

New file [LeverageRotation.lean](formalization/LeverageRotation.lean) has three checked theorems, with warnings treated as errors and the pinned existing Mathlib environment:

1. The squared coefficients of the proposed rotation sum to 1 when theta and r are in[0,1]. This checks the covariance normalization algebra; Gaussian-process construction is a further correspondence step.
2. Nonnegative stress diffusion loadings give a nonpositive instantaneous stress-covariance contribution under the proposed negative rotation. Identifying those loadings with the SDE remains explicit; this is not a theorem of finite-wing implied-volatility monotonicity.
3. The declared finite-step log-stock recurrence is reconstructed exactly by the three accumulated sums. This is the reusable computational bridge, not certification of floating-point sampling or code.

Audit: `formalization/leverage_audits/attempt_001/audit.json`; no admitted proof or custom axiom. Existing Gaussian-mean-one/control-reuse and wing-moment proof sources also match their previously successful audit hashes. They do not already prove the new rotation's full VIX invariance or new control identities.

Lean can help rule out algebraic mistakes and certify exact shortcuts and necessary constraints. It does not generate empirical smile responses, search parameter space or remove rare-event sampling error. The efficient workflow is to derive mechanisms, prove the reusable structural claims, screen with paired fast numerics, and validate the few surviving candidates independently.
