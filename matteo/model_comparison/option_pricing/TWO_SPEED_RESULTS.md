# Shared fast and slow stress: joint VIX experiment

27 October 2026. All four maturities fitted together; no expiry-specific parameters.

The bounded pilot changes checking option loss from **1.725555 to 0.050081** (+97.10% improvement), with 83 → 153 of 217 native quotes inside bid–ask. This is a numerical result on the same market snapshot, not evidence across market regimes.

Training option loss: 1.733264 → 0.043630. Maximum training forward residual: 2.49e-10 VIX points. Checking forward residuals are shown below; training equality does not establish population equality.

## The model in plain language

The old M2 response describes the original variance mechanism. The stress part can now react on two different timescales: one forgets shocks quickly, while the other keeps them for longer. Every option expiry uses the same states, weights, amplitudes and half-lives. Variance is mixed first; the VIX square root and option payoff come afterwards.

$$v_t=f+(\xi(t)-f)\left[(1-w)R_t+w\{(1-r)A_t^s+rA_t^f\}\right],$$

$$A_t^i=\exp(b_i Z_t^i-\tfrac12b_i^2q_i(t)),\quad dZ_t^i=-a_iZ_t^i dt+\sqrt{2a_i}dB_t^i,\quad a_i=\log2/h_i.$$

The two stress Brownian innovations are independent of each other and the old M2 bank. Each normalized response has mean one, so $E_Q[v_t]=\xi(t)$. The seven active variance levels and the floor are still fitted. See [the full plan](TWO_SPEED_PLAN.md) for normalization, conditional expectations, bounds, and estimator details.

This tests the simpler endpoint representation motivated by the previous continuous-measure experiment, with an additional stress state. It does not refit the Beta measure at the same time. Three variance branches now replace the old M2-plus-single-stress pair; this is not an arbitrary continuous measure over speeds.

## Fitted shared parameters

| Parameter | Fitted value |
|---|---:|
| Slow half-life | 365.0000 days |
| Fast half-life | 71.3905 days |
| Slow amplitude | 2.285440 |
| Fast amplitude | 3.535495 |
| Total stress fraction $w$ | 0.922838 |
| Fast fraction of stress $r$ | 0.658093 |
| Floor volatility | 12.0000% |

The slow half-life and floor volatility reach their upper bounds, as do two active variance levels; the first active level reaches the floor. Three starts converge to essentially the same best fit. The fourth reaches its iteration limit and misses forward feasibility, so it is excluded and preserved in the optimizer ledger. These results do not establish that the bounds are economically optimal.

The old M2 smooth cap remains unchanged (level 5.5027800568, sharpness 2). The stress exponentials have no hard state cap; their Gaussian normalization makes their means finite. The amplitude bound of 5 is a calibration parameter bound, not a cap on realized variance.

The labels fast/slow enforce ordering, not separate rules for short/long options. Either component affects every maturity. Amplitudes are bounded by 5, half-lives by 1–365 days, floor volatility by 5–12%, and free variance levels by floor–0.2. Four starts and a refitted nested boundary were compared; this does not prove a global optimum.

## All four maturities

| Days | Single loss contribution | Two-speed loss contribution | Single wing IV MAE | Two-speed wing IV MAE |
|---:|---:|---:|---:|---:|
| 27 | 0.927941 | 0.014328 | 22.843 pp | 4.446 pp |
| 55 | 0.378273 | 0.011342 | 8.501 pp | 5.247 pp |
| 118 | 0.062100 | 0.020561 | 1.325 pp | 2.463 pp |
| 153 | 0.357241 | 0.003850 | 1.297 pp | 2.160 pp |

Wing means use native quotes with strike/market forward > 2. Their strike ranges differ by maturity, so compare models within each row, not slopes across different strike ranges. Native IVs use each model's own forward; plot IVs use canonical calls and market forward. Signed estimates are preserved and unavailable IVs are not replaced.

| Days | Target forward | Single error ± SE | Two-speed error ± SE |
|---:|---:|---:|---:|
| 27 | 20.1000 | +0.026777 ± 0.003513 | -0.026838 ± 0.013620 |
| 55 | 20.7300 | -0.011981 ± 0.005517 | +0.055170 ± 0.030998 |
| 118 | 21.1900 | -0.002539 ± 0.005007 | -0.011069 ± 0.029063 |
| 153 | 21.7200 | +0.008292 ± 0.007060 | +0.035924 ± 0.024641 |

Paired change in the historical total loss (option loss plus forward-error contribution): -1.671233, jackknife SE 0.028181. This conditional evaluation uncertainty excludes calibration variability and model error.

## Numerical checks and interpretation limits

All 14 pricing/forward derivative directions and the exact nested boundary replay passed. The additional audit independently reconstructs fast importance weights, checks conditional stress expectations with GH80, compares stress time integration with adaptive scalar quadrature on payoff-relevant states, and reports weighted tail ESS/concentration. M2 conditional tables and their earlier audits are reused.

| Days | Weighted tail ESS | Top 10 tail payoff share | Likelihood mass ± SE |
|---:|---:|---:|---:|
| 27 | 3844.9 | 0.501% | 1.000052 ± 0.003001 |
| 55 | 3876.0 | 0.503% | 0.997768 ± 0.003836 |
| 118 | 3939.0 | 0.492% | 1.002065 ± 0.003666 |
| 153 | 4070.6 | 0.488% | 1.001605 ± 0.005773 |

The estimator is $E[P_{slow}]+E[P_{two}-P_{slow}]$. Training reuses 32,768 old/slow states per expiry; checking reuses 65,536 old/slow states from the earlier study and draws new independent fast coordinates. Independent within-block permutations avoid treating matching one-dimensional net ranks as a joint QMC design. Both paired terms use the product likelihood; eight whole blocks determine SEs. Checking is therefore not wholly fresh. The independent permutation can increase variance; compare forward gaps with SEs before considering adoption.

Negative native checking estimates: single 1, two-speed 1; retained without clipping.

Extra timescales do not guarantee a steeper wing. If one exponential dominates payoff-contributing states, the response can remain approximately lognormal. Any improvement must be assessed across all four maturities and along with the checking forwards.

Original production and working VIX baselines remain unchanged. SPX stays paused. No new Lean theorem or population forward certificate is claimed.

## Artifacts

- [Four smile panels](visualizations/two-speed024/vix-two-speed-smiles.png)
- [Right wings](visualizations/two-speed024/vix-two-speed-right-wing.png)
- [Call prices](visualizations/two-speed024/vix-two-speed-prices.png)
- [Checking forwards](visualizations/two-speed024/vix-two-speed-forwards.png)
- Run: `runs/m2_two_speed_001` (source snapshots, inputs hashes, optimizer starts, training/checking tables, blocks, paths, audit).
