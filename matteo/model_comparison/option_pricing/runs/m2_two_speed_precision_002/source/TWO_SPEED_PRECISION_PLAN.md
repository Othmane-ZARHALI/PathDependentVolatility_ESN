# Protected two-speed SPX pricing experiment

28 October 2026. The user authorized the proposed sequence: improve SPX pricing precision with the successful VIX model fixed; then fit only the nine off-window variance levels if prices become reliable; subsequently consider stress/equity correlation if needed. Continue in the existing chat. No model extension is introduced in this experiment.

## Frozen model and scope

Use the baseline vector in `runs/m2_two_speed_joint_001/design.json`, which embeds the selected VIX-only model from `runs/m2_two_speed_001`. Keep the original M2 dynamics, readout, normalizer, initial state, premia, positive spike shift, floor, both stress factors and all16 variance levels fixed for this precision pilot. Preserve all1751 eligible SPX quotes, with1294 primary quotes at29/57/119/211/302/666/1030 days. Keep all217 VIX quotes and four forward targets; they need no repricing because no model parameter changes.

Later, only xi indices0,5,6,9,11,12,13,14,15 may move in the protected fit. The exact VIX mapping ignores those indices. In this exogenous model, changing them cannot affect any conditional variance inside the four VIX windows. This preserves the existing VIX law and existing forward accuracy, not a new certificate of exact population forwards.

## Independent-stress pricing control

Write the actual variance as

$$v_t=f+(\xi_t-f)[(1-w)R_t+wH_t],\qquad H_t=(1-r)A_t^s+rA_t^f.$$

Construct an auxiliary reference with the original M2 response replaced by its mean1:

$$\bar v_t=f+(\xi_t-f)[(1-w)+w\bar H_t].$$

Here the reference stress is held at a daily-ish subset of the existing64-step/trading-day stock grid, additionally split at every represented curve-segment change. This is a control only: the target variance and stock discretization are unchanged. The reference is driven by the same equity increments as the target. Conditional on its stress path, its normalized terminal stock is lognormal with total variance $I=\sum\bar v_t\Delta t$. Its call expectation is the Black price $B(I,K)$.

For old-path importance weight L, payoff g, actual stock S and auxiliary stock Sbar,

$$E_Q[g(S)]=E_{\rm stress}[B(I,K)]+E_{\rm proposal,stress}[L\{g(S)-g(\bar S)\}].$$

The Black term is estimated on fresh independent cheap stress paths. It is not multiplied by the old-path likelihood. The expensive correction retains the original full-path likelihood; the two noises are independent. This identity concerns the existing discrete stock model, not a continuous-time approximation or proof of numerical accuracy.

Compare three predeclared corrections: unhedged; terminal Black-delta control; and a daily-ish dynamic difference hedge. Conditional on the entire independent stress path, the old equity driver remains Brownian. The target and reference discrete stocks have conditional-mean-zero increments. Thus terminal delta times their terminal difference and predictable Black-delta stochastic sums have zero mean. Dynamic correction is

$$L\left[g(S_T)-g(\bar S_T)-\sum_j\{\delta_j(S_j)(S_{j+1}-S_j)-\delta_j(\bar S_j)(\bar S_{j+1}-\bar S_j)\}\right].$$

The deltas use remaining reference integrated variance. Reference and target stocks are observed at the same subset grid; the target evolves on every fine step between hedge dates. This is a discrete zero-mean control, not an assertion of exact replication. The dynamic estimator is the prespecified main candidate; do not select a method using checking fit loss. Native puts are reconstructed using the model's mathematical unit forward, with signed estimates preserved.

## Checking and decision gate

Use both existing512-path search/check old-state caches and their separate stress streams. The caches have already been reused: checking is numerical replication, not untouched market or simulation validation. Eight expensive blocks preserve the earlier uncertainty convention. Use eight independently seeded cheap batches of4096 antithetic pairs each (65536 stress paths per expiry/role); a pair is the sampling unit. The combined SE adds independent expensive-correction and cheap-reference variances. Store both components separately, all prices, blocks, seeds, integrals, hashes and costs.

Verify target-stock replay against the original kernel, an independent Python reference path and integrated variance, all16 curve derivatives, the identical-model zero-residual limit, and invariance of the VIX parameter map under all nine permitted changes. Follow up with independent cheap-RNG replay and finite-law control-mean checks before interpreting a gain.

The existing primary sampling gate is $2SE\le price\_scale/4$. Every primary quote must pass before this experiment permits calibration. This is an empirical sampling screen only; stock timestep and normalizer errors require separate checks. If the gate fails, quantify whether the remaining uncertainty comes from expensive corrections or cheap stress integration and improve that component. Do not launch another unrestricted joint fit. No new factor, correlation, optimizer, production baseline replacement, clipping, smoothing or maturity-specific parameter is authorized by this pricing pilot.
