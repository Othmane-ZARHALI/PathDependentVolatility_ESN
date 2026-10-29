# A modest bank of fixed geometric stress speeds

Calibration: 27 October 2026. Audited and recorded: 28 October 2026. All 217 VIX quotes and all four maturities fitted jointly; SPX remains paused.

Checking option loss is **0.062160**, compared with **0.050081** for the previously selected fitted two-speed model. Quotes inside bid–ask: **150/217**, versus 153/217. Training loss is 0.032861. This compares two bounded fits on one stored market snapshot, not performance across market regimes.

The four-node bank is retained as an experimental candidate, not adopted as a reliable exact-forward calibration. Training forward means match, but their block SEs are much larger than the equality residuals, and checking errors reach 0.216 VIX points. The current numerical experiment does not establish an improvement over the two-speed model, or an inability of this bank family to fit the forwards. Improve forward estimation before increasing the bank further.

## Choosing the half-lives

We fixed the half-lives before calibration:

$$h_j=7\,4^{j-1}\text{ days},\qquad j=1,2,3,4,$$

| Half-life | Calendar interpretation |
|---:|---|
| 7 days | One week |
| 28 days | Four weeks, roughly one month |
| 112 days | Sixteen weeks, roughly four months |
| 448 days | Sixty-four weeks, roughly fifteen months |

Each step multiplies the persistence time by four. Equal spacing on a logarithmic axis covers short and persistent stress with only four nodes. These calendar scales were not selected to reproduce the earlier fitted 71/365-day pair. This is a transparent modelling choice, not a theorem that these endpoints are optimal.

The 7-day node can describe shocks that fade within the 30-day VIX averaging window. The 28/112-day nodes span intermediate persistence. The 448-day node allows stress to outlast the longest priced expiry and extends beyond the previous 365-day search ceiling. A half-life measures shock decay; it is not an option expiry. Every node affects all maturities.

Long half-lives can be difficult to distinguish over a short observation/pricing horizon, especially when sensitivities are also fitted. A fitted weight on a node therefore does not prove the market contains precisely that persistence. A later grid-sensitivity test would be separate; this experiment does not choose a new grid after viewing checking prices.

## The shared variance process

$$dZ_t^j=-a_jZ_t^jdt+\sqrt{2a_j}dB_t^j,\quad Z_0^j=0,\quad a_j=\log2/h_j,$$

$$A_t^j=\exp\left(b_jZ_t^j-\tfrac12b_j^2(1-e^{-2a_jt})\right),$$

$$v_t=f+(\xi(t)-f)\left[\alpha_0R_t+\sum_{j=1}^4\alpha_jA_t^j\right],\quad\alpha_j\ge0,\quad\sum_{j=0}^4\alpha_j=1.$$

Here $R_t$ is the original normalized M2 response. The four stress Brownian motions are independent of each other and the old M2 bank, matching the independence assumption used in the previous two-speed test. Each response has mean one, so $E[v_t]=\xi(t)$ and $v_t\ge f$ when $\xi(t)\ge f$. Mixture weights apply to the normalized non-floor variance response, not to option prices.

Only speeds are fixed. Four sensitivities, four independent mixture weights, seven active variance levels and the floor are fitted: 16 parameters, versus 14 previously. Nine inactive curve levels remain unchanged. The amplitude bound is 5, active variance upper bound 0.2, floor volatility 5–12%. Positive weights use stick-breaking internally so every trial remains a positive mixture. There is no sparsity penalty or smoothness prior on the weights.

The old M2 dynamics, readout, positive spike shift, premiums and smooth cap remain frozen. Stress exponentials have their Gaussian normalization and no hard state cap. The bank is a finite positive measure over speeds; a continuous speed measure is not simultaneously calibrated.

## Fitted bank

| Component | Sensitivity | Weight in normalized response |
|---|---:|---:|
| M2 | Frozen | 1.5257% |
| 7-day stress | 4.643263 | 7.4820% |
| 28-day stress | 4.495896 | 42.7015% |
| 112-day stress | 3.576861 | 24.0654% |
| 448-day stress | 2.686159 | 24.2255% |

Floor volatility: 12.0000%. Active variance levels: 0.01440000, 0.04093954, 0.18585615, 0.05728694, 0.20000000, 0.07021709, 0.08372825.

Weights near zero mean the corresponding node contributes little in this fitted snapshot; its fitted amplitude may then be weakly identified. Conversely, positive weights do not establish unique economic factors. Bound hits and optimizer records are retained in the run folder.

## All four maturities

| Days | Two-speed loss contribution | Bank loss contribution | Two-speed wing IV MAE | Bank wing IV MAE |
|---:|---:|---:|---:|---:|
| 27 | 0.014328 | 0.019718 | 4.446 pp | 6.356 pp |
| 55 | 0.011342 | 0.019723 | 5.247 pp | 7.114 pp |
| 118 | 0.020561 | 0.015517 | 2.463 pp | 2.158 pp |
| 153 | 0.003850 | 0.007203 | 2.160 pp | 2.051 pp |

Wing means use native quotes with strike/market forward >2. Strike ranges differ by maturity. Native IVs use each model forward; canonical-call plots use the market forward. Because checking forward errors are material, their left-wing IVs can be unavailable even for a positive call estimate below market-forward intrinsic value. The separate native-IV panels use each model forward and avoid that parity-display issue; they still inherit the forward-estimation uncertainty. Signed prices and unavailable IVs are retained without clipping or smoothing.

| Days | Target forward | Bank checking error ± block SE |
|---:|---:|---:|
| 27 | 20.1000 | -0.093494 ± 0.037842 |
| 55 | 20.7300 | +0.178419 ± 0.071646 |
| 118 | 21.1900 | -0.216210 ± 0.083138 |
| 153 | 21.7200 | +0.194153 ± 0.077991 |

Maximum training forward gap: 8.48e-12 VIX points. These are equalities on the calibration estimator. Checking residuals and SEs do not certify exact population forward matching.

Paired change in historical total loss (options plus checking-forward contribution): +0.131919; eight-block jackknife SE 0.081511. Shared old-state streams are reflected in the pairing. This uncertainty excludes fitting variability, grid uncertainty and model error.

Training forward block SEs range from about 0.049 to 0.154 VIX points; checking SEs range from 0.038 to 0.083. The checking gaps are roughly 2.5 checking SEs, and are smaller relative to combined training/checking uncertainty. This is a numerical reliability limitation rather than evidence that the model family cannot match the forwards.

Option-only paired checking loss change: +0.012079, jackknife SE 0.024117. The apparent small ranking difference should not be treated as established model superiority.

## Adding further half-lives

More finite nodes are mathematically admissible, but independent drivers can diversify away stress as their weights spread. A denser independent bank is not simply a finer quadrature of the same random continuous spectrum. For genuinely smooth refinement, correlated factors driven by a small number of common shocks are a possible separate comparison. See [the refinement note](STRESS_BANK_REFINEMENT.md). The present independent four-node model was kept unchanged during fitting.

## Cached bank and numerical checks

For European VIX prices we cache terminal standard-normal coordinates at each expiry and compute future conditional variance analytically. Complete joint time paths are not needed for these marginal prices. Speeds are fixed, while amplitudes and weights rebuild responses, variance integrals, VIX and payoffs. Reusing coordinates does not mean reusing completed prices.

Training reuses 32,768 old/node2/node3 coordinates per expiry from earlier work, and checking reuses 65,536. Nodes0/1 are newly drawn jointly from a defensive Gaussian mixture, independently permuted within each block. The target has four independent Gaussian stress coordinates. Product importance likelihoods are used without sample normalization. Checking is not wholly fresh and uses the same market snapshot.

The multilevel identity is $E[P_{bank}]=E[P_{reference}]+E[P_{bank}-P_{reference}]$. The reference replaces the first three stress responses by their unconditional mean one and retains the slowest node and M2. The large term uses original likelihoods; both paired terms use identical product likelihoods and controls. Both sides form conditional variance before the VIX square root. This changes the estimator, not the bank model.

All 16 derivative directions and the exact reference boundary replay pass. The independent audit reconstructs likelihoods, replays all 32 additional checking seeds, checks Gaussian mean-one normalization and conditional expectations with GH80, and reconstructs VIX on payoff-relevant states using adaptive scalar stress integrals. Existing M2 conditional tables and their earlier audits are reused. These numerical checks are not a new Lean theorem or uniform error certificate.

| Days | Weighted tail payoff ESS | Top 10 tail share | Likelihood mass ± SE |
|---:|---:|---:|---:|
| 27 | 1310.3 | 1.766% | 1.000651 ± 0.005734 |
| 55 | 1755.9 | 1.511% | 0.997385 ± 0.006492 |
| 118 | 2257.3 | 1.311% | 1.010445 ± 0.006326 |
| 153 | 2427.3 | 1.209% | 0.998797 ± 0.006158 |

Negative native checking price estimates: bank 1; previous two-speed 1. They are preserved.

The joint SLSQP stage is preserved in runs/m2_stress_bank_001. A finishing least-squares step adjusts only the seven already-free variance levels at fixed stress parameters, restoring the four forward equalities for each candidate. Original and restored candidates are retained and selection uses training loss only. New additional checking coordinates are generated after this final selection. This changes optimization completion, not the variance model, grid, quote set or forward targets.

The two-speed reference is the earlier selected fit and stored checking blocks; it is not refitted here. The new fixed-speed family does not contain those exact fitted speeds, so the comparison is not a nested-model test. The winning candidate was the original converged third start, already matching forwards within8.5e-12; the finishing step did not change the selected model. Three starts were attempted; all successes/failures and feasibility gaps are retained. Neither the production baseline nor the working VIX candidate is replaced.

## Artifacts

- [Native IV panels, using each model forward](visualizations/stress-bank025/vix-stress-bank-native-IV.png)
- [Four smile panels](visualizations/stress-bank025/vix-stress-bank-smiles.png)
- [Right wings](visualizations/stress-bank025/vix-stress-bank-right-wing.png)
- [Call prices](visualizations/stress-bank025/vix-stress-bank-prices.png)
- [Checking forwards](visualizations/stress-bank025/vix-stress-bank-forwards.png)
- [Persistence grid and fitted weights](visualizations/stress-bank025/stress-bank-persistence-weights.png)
- [Experiment plan](STRESS_BANK_PLAN.md)
- Run: `runs/m2_stress_bank_002`.
