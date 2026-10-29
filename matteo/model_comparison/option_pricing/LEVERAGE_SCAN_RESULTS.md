# Stock/stress leverage scan — 28 October 2026

The experiment introduces negative stock/stress dependence while keeping the full variance-factor law fixed. It materially improves both price loss and native implied-volatility error. It does not yet give a satisfactory double calibration.

## What changed

For fixed stress direction (cs,cf), use

$$dW^{stock}_\theta=\sqrt{1-\theta^2}dW^0-\theta(c_s dB^s+c_f dB^f).$$

A positive stress shock raises variance and now tends to push the stock down. We tested theta = 0, .2, .4, .6, .8, .95 with mixed, slow-only and fast-only directions. Mixed uses (sqrt(1-r),sqrt(r)); the other directions isolate each factor. All OU dynamics, their joint law, M2 readout, floor and protected-candidate variance curve remain fixed. Original positive spike shift remains 0.369729973. There is no endogenous volatility/return feedback. The old stock/M2 correlations shrink by sqrt(1-theta²).

The variance-factor process and conditional future variance law are unchanged. Thus the VIX distribution, its options and forwards are protected algebraically. Existing finite-sample VIX forward error is inherited; this is not a new population-forward equality certificate. Baselines and selected models from earlier runs are not overwritten.

## Fresh numerical comparison

Selection was frozen before fresh expensive checking: **mixed, theta=0.8**, based on the replay sample's fine-grid equal-maturity price loss. All 1751 eligible SPX quotes / all seven maturities remain in CSV; the original 1294 primary quotes determine scores.

| Metric, 128 steps/trading day | Theta 0 | Selected theta 0.8 |
|---|---:|---:|
| Equal-maturity price interval loss | 3788.48 | 1861.38 |
| Equal-maturity IV MAE, percentage points | 4.493 | 3.199 |
| Primary prices inside bid–ask | 7 / 1294 | 18 / 1294 |
| Invalid primary IVs | 0 | 0 |
| Negative primary price estimates | 0 | 0 |

Price loss falls **50.9%**. Approximate paired 95% batch-bootstrap loss-change interval: [-2182.804800929469, -1622.6946831120715]. This interval covers sampling of expensive and cheap blocks; it excludes timestep, normalizer and market errors. Only eight independent batches underpin this approximation. IV MAE omits unavailable IVs, whose counts remain explicit. No signed estimate is clipped or smoothed.

Fresh full-grid scores are in [scores.csv](runs/m2_leverage_scan_001/fresh/scores.csv). Looking across that grid after checking is exploratory; it is not an independent validation of a newly selected winner. Some maturities remain visibly mismatched; inspect the seven maturity panels and maturity_metrics.csv rather than treating the aggregate gain as a good calibration.

## What still misses the market

The improvement occurs at every maturity, but the left wing remains too low. For 29 days, compare quotes nearest strike/forward .9 and 1.1: the model IV(right)-IV(left) changes from -1.52 to -5.58 percentage points, against market -14.05. Negative stress leverage tilts the smile in the correct direction but does not yet produce enough slope. At 57 days the corresponding model gap changes -0.63 to -6.29 pp versus -13.12 in the market. Increasing theta mainly lowers the right wing and modestly raises the left; the plots make the residual mismatch clear.

Only 18 of 1294 primary prices are inside bid-ask. Median selected Monte Carlo SE is 1.36 quote scales; maturity-wise median 64/128 grid changes range .057 to .379 scales, with a largest quote change 1.58 scales. Final precision gates are unresolved. The selected auxiliary stock-forward sample means range .999870 to 1.000028, with batch SE .000021 to .000249; these are diagnostics, not exact numerical forward matching. The largest single conditional-forward contribution is below .0041%.

Mixed theta=.95 has a fresh price loss of 1787.57 and IV MAE 3.120 pp, slightly better than .8. It was not selected before fresh checking, and the replay sample preferred .8, so no uniquely determined optimal theta is claimed.

## Computation and checks

Each factor stream is built once, then caches I=sum(v dt), U0=sum(sqrt(v)dW0), Us=sum(sqrt(v)dBs), Uf=sum(sqrt(v)dBf). For every theta/direction,

$$\log S_T=-I/2+\sqrt{1-\theta^2}U_0-\theta(c_sU_s+c_fU_f).$$

Replay uses 512 earlier paths; checking uses 1024 fresh full-model paths (512 antithetic pairs). Both reuse the same 65536-path independent auxiliary stress reference, with its separate batch uncertainty included. Thus the expensive checking streams are fresh, while the complete estimator is not wholly independent of screening. These are the same market quotes, not market out-of-sample data. Coupled 64/128 grids diagnose timestep sensitivity. Total scan wall time: 210.4 seconds, including compile, both factor samples, reference simulation, pricing all configurations, IV inversion and paired bootstrap.

Stress Brownian increments are sampled conditionally on OU innovations using their exact joint covariance, including the independent residual. Full replay at theta=0 matches stored stocks to 3.6e-13. Independent scalar target/reference sufficient-statistic oracles agree to 6.94e-18 and 2.44e-15; standardized stress covariance error is 1.86e-16.

The old future-stress-dependent delta control is not reused. The new auxiliary call is integrated conditionally with log mean -I/2-theta(cs Us+cf Uf) and log variance (1-theta²)I. Add the paired full/reference payoff difference and subtract a deterministic strike-dependent Black delta times their stock difference. Each finite predictable log-Euler stock has mathematical mean one; no sample-forward renormalization is applied. Continuous-time martingality and tail convergence are separate questions, not established here. Stock/reference forward diagnostics are saved, not used to repair prices.

Existing Lean results prove coefficient normalization, a covariance-sign identity and an abstract discrete reconstruction identity. These numerical implementation checks do not turn them into a full simulator, VIX-invariance or continuous-martingale theorem.

## What to do next

Keep the mixed-leverage construction as a promising experimental extension. Before fitting additional parameters, improve independent reference/tail precision for the promising theta range and inspect remaining maturity-specific shape errors. Then a joint SPX/VIX fit can vary the previously authorized curve/readout parameters within explicit VIX protection. Do not introduce endogenous OU feedback or adopt a new production model from this exploratory scan alone.

Figures: [theta response](visualizations/leverage028/SPX-response.png), [all seven mixed-direction smiles](visualizations/leverage028/SPX-theta-fan.png), [stress-direction comparison](visualizations/leverage028/SPX-stress-directions.png). Detailed maturity, grid and reference-forward diagnostics accompany them.
