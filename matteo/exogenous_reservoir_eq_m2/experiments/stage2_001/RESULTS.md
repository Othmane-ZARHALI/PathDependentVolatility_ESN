# Targeted Stage-2 M2/M3/M2+O/M4 validation

## Status and interpretation boundary

This is a physical-measure capacity study on synthetic model paths. The acceptance bands are declared proxies, not empirical SPX targets. Accordingly, acceptance rates compare the four cells under one declared sampling measure; they do not yet validate the model against market data.

## Experimental design

- Cache seeds: `[196613, 221251, 262147, 294001, 324503]`.
- Full readout draws: `32` using `sobol` with design seed `20261021`.
- Parameter space: `targeted_stage2`.
- Structural scenarios: `4`.
- Simulation: `{'years': 5.0, 'burn_years': 1.0, 'paths': 16, 'observations_per_year': 252, 'steps_per_observation': 2}`.
- Cells: M2 (neither optional loading), M3 (echo only), M2+O (orthogonal energy only), and M4 (both).
- The echo consumes only the common feedback, spike, and clustering banks, so the orthogonal switch cannot silently redefine M3.

## Latent architecture calibration

- Optimizer seeds: `[20260902, 20260903, 20260904]`; selected seed: `20260902`.
- Objective: `0.03039`.
- Objective range: `[0.03038740521353733, 0.030387405276293338]`.
- Optimizer status: `Optimization terminated successfully.`.
- Feedback rates: `[0.05000000000000001, 5000.000000000004]`.
- Clustering rates: `[0.5, 20.710215785426232]`.
- Search-bound hits: `['feedback_low:lower', 'feedback_high:upper', 'cluster_low:upper']`.

This fit targets only the finite-band rough variogram and latent clustering curve. It does not determine leverage, Zumbach, tails, or Taylor effects. Endpoint hits require a wider-box or regularised follow-up before treating the architecture as identified.

## Acceptance by independent cache seed

| seed | M2 | M3 | M2+O | M4 |
| --- | --- | --- | --- | --- |
| 196613 | 0.0% | 0.0% | 0.0% | 0.0% |
| 221251 | 0.0% | 0.0% | 0.0% | 0.0% |
| 262147 | 0.0% | 0.0% | 0.0% | 0.0% |
| 294001 | 0.0% | 0.0% | 0.0% | 0.0% |
| 324503 | 0.0% | 0.0% | 0.0% | 0.0% |

## Acceptance by structural scenario

| scenario | M2 | M3 | M2+O | M4 |
| --- | --- | --- | --- | --- |
| 0 | 0.0% ± 0.0pp | 0.0% ± 0.0pp | 0.0% ± 0.0pp | 0.0% ± 0.0pp |
| 1 | 0.0% ± 0.0pp | 0.0% ± 0.0pp | 0.0% ± 0.0pp | 0.0% ± 0.0pp |
| 2 | 0.0% ± 0.0pp | 0.0% ± 0.0pp | 0.0% ± 0.0pp | 0.0% ± 0.0pp |
| 3 | 0.0% ± 0.0pp | 0.0% ± 0.0pp | 0.0% ± 0.0pp | 0.0% ± 0.0pp |

## Stylised-fact comparison

Entries are fixed-design means across seeds; the value after `±` is the between-seed standard deviation.

| model | max \|return ACF\| | mean log-variance ACF | rough H | rank leverage | rank Zumbach | excess kurtosis | absolute-tail ratio | Taylor gap | P cap fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| M2 | 0.07845 ± 0.002537 | 0.2006 ± 0.02037 | 0.04081 ± 0.002726 | -0.005818 ± 0.003356 | 0.01524 ± 0.007758 | 8.499 ± 0.4154 | 3.653 ± 0.03248 | 0.0209 ± 0.002974 | 0.002919 ± 0.0001278 |
| M3 | 0.07847 ± 0.002542 | 0.2008 ± 0.02034 | 0.04106 ± 0.002722 | -0.005819 ± 0.003329 | 0.01483 ± 0.00772 | 8.51 ± 0.4117 | 3.654 ± 0.03171 | 0.02093 ± 0.002959 | 0.00294 ± 0.000131 |
| M2+O | 0.07986 ± 0.002133 | 0.2172 ± 0.01993 | 0.04823 ± 0.003101 | -0.005309 ± 0.003162 | 0.01489 ± 0.003312 | 8.751 ± 0.5254 | 3.672 ± 0.02525 | 0.02458 ± 0.003114 | 0.002986 ± 0.0001746 |
| M4 | 0.07989 ± 0.00214 | 0.2174 ± 0.01993 | 0.04842 ± 0.003097 | -0.005312 ± 0.00314 | 0.01447 ± 0.003248 | 8.762 ± 0.5228 | 3.672 ± 0.02574 | 0.0246 ± 0.003098 | 0.003011 ± 0.0001794 |

## Individual gate pass rates

These rates identify the bottlenecks hidden by joint acceptance.

| model | max_abs_return_acf | mean_log_variance_acf | rough_hurst | leverage_rank | zumbach_rank | excess_kurtosis | taylor_gap |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M2 | 100.0% ± 0.0pp | 53.1% ± 3.6pp | 8.1% ± 1.6pp | 34.2% ± 19.8pp | 50.3% ± 22.6pp | 90.6% ± 1.7pp | 100.0% ± 0.0pp |
| M3 | 100.0% ± 0.0pp | 53.9% ± 3.2pp | 8.1% ± 2.5pp | 35.5% ± 19.4pp | 48.4% ± 19.8pp | 90.2% ± 1.9pp | 100.0% ± 0.0pp |
| M2+O | 100.0% ± 0.0pp | 57.0% ± 6.1pp | 17.8% ± 1.0pp | 29.7% ± 17.0pp | 46.2% ± 9.7pp | 91.9% ± 2.0pp | 100.0% ± 0.0pp |
| M4 | 100.0% ± 0.0pp | 58.0% ± 6.3pp | 16.7% ± 2.2pp | 29.5% ± 17.4pp | 42.7% ± 10.0pp | 91.7% ± 2.1pp | 100.0% ± 0.0pp |

## Post-hoc threshold sensitivity

These deliberately stricter thresholds are diagnostics, not preregistered acceptance criteria or empirical SPX estimates. They test whether the declared proxy gates are masking weak magnitudes.

| model | rough H in [0.08, 0.12] | rank leverage < -0.01 | rank Zumbach > 0.02 | rank Zumbach > 0.05 | excess kurtosis > 3 | excess kurtosis > 5 | joint H/leverage/Z>0.02/kurtosis>3 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M2 | 4.2% ± 0.7pp | 24.8% ± 16.4pp | 40.2% ± 18.3pp | 1.7% ± 1.7pp | 90.6% ± 1.7pp | 81.6% ± 3.9pp | 0.0% ± 0.0pp |
| M3 | 3.8% ± 1.4pp | 24.7% ± 15.0pp | 36.9% ± 17.3pp | 3.4% ± 2.1pp | 90.2% ± 1.9pp | 81.6% ± 3.6pp | 0.0% ± 0.0pp |
| M2+O | 11.6% ± 0.9pp | 22.3% ± 15.1pp | 34.7% ± 9.0pp | 1.6% ± 1.5pp | 91.9% ± 2.0pp | 83.3% ± 2.0pp | 0.0% ± 0.0pp |
| M4 | 10.6% ± 1.6pp | 22.2% ± 13.5pp | 30.9% ± 10.3pp | 3.3% ± 2.4pp | 91.7% ± 2.1pp | 83.0% ± 2.3pp | 0.0% ± 0.0pp |

## Exploratory M2 feedback-curvature screen

This post-hoc quartile screen is descriptive, not causal. It is included to choose the next parameter box rather than to validate a fitted effect.

| feedback curvature | rough_hurst | excess_kurtosis | leverage_rank | zumbach_rank | accepted |
| --- | --- | --- | --- | --- | --- |
| [0.1012, 0.1428] | 0.07276 ± 0.003217 | 5.947 ± 0.1675 | -0.002025 ± 0.003686 | 0.02114 ± 0.004199 | 0 ± 0 |
| [0.1549, 0.2194] | 0.04063 ± 0.002865 | 7.018 ± 0.3537 | -0.005142 ± 0.003453 | 0.0177 ± 0.007689 | 0 ± 0 |
| [0.2325, 0.3329] | 0.03241 ± 0.00268 | 8.873 ± 0.5433 | -0.006829 ± 0.003422 | 0.01162 ± 0.009534 | 0 ± 0 |
| [0.3395, 0.4847] | 0.01743 ± 0.002273 | 12.16 ± 0.6406 | -0.009277 ± 0.003133 | 0.01051 ± 0.01067 | 0 ± 0 |

## Paired optional-component effects

Contrasts use identical primitive paths, structural scenarios, and core readout draws. Positive acceptance effects are improvements; for distance, negative effects are improvements. Other signs retain the metric's meaning.

| effect | acceptance_rate | distance | rough_hurst | leverage_rank | zumbach_rank | excess_kurtosis |
| --- | --- | --- | --- | --- | --- | --- |
| echo_at_O0 | 0 ± 0 | -0.001505 ± 0.000288 | 0.0002495 ± 4.627e-06 | -2.399e-07 ± 3.308e-05 | -0.0004154 ± 0.0001376 | 0.01115 ± 0.005017 |
| echo_at_O1 | 0 ± 0 | -0.001815 ± 0.000287 | 0.0001845 ± 4.976e-06 | -3.708e-06 ± 3.15e-05 | -0.0004261 ± 0.0001489 | 0.01123 ± 0.003426 |
| echo_orthogonal_interaction | 0 ± 0 | -0.0003104 ± 0.0001628 | -6.497e-05 ± 4.628e-06 | -3.468e-06 ± 7.841e-06 | -1.073e-05 ± 5.867e-05 | 7.944e-05 ± 0.001682 |
| orthogonal_at_echo0 | 0 ± 0 | -0.03924 ± 0.006905 | 0.007423 ± 0.001283 | 0.0005099 ± 0.0008207 | -0.0003495 ± 0.005286 | 0.2526 ± 0.1324 |
| orthogonal_at_echo1 | 0 ± 0 | -0.03955 ± 0.006854 | 0.007358 ± 0.001279 | 0.0005064 ± 0.0008141 | -0.0003602 ± 0.005231 | 0.2526 ± 0.1338 |

## Evidence summary

- Joint acceptance was M2 `0.0%`, M3 `0.0%`, M2+O `0.0%`, and M4 `0.0%` under the declared proxy bands.
- The echo-only contrast changed distance by `-0.001505 (95% seed CI [-0.001863, -0.001147])`. It changed rank leverage by `-2.399e-07 (95% seed CI [-4.131e-05, 4.083e-05])` and rank Zumbach by `-0.0004154 (95% seed CI [-0.0005862, -0.0002446])`; signs must be judged against the desired negative leverage and positive Zumbach directions.
- The orthogonal-only distance contrast was `-0.03924 (95% seed CI [-0.04782, -0.03067])`. An effect of this size should not justify the additional bank unless it becomes stable under a broader loading range or a targeted conditional test.
- At the row level, the echo changed `0` failures into passes and `0` passes into failures. The orthogonal block changed `0` failures into passes and `0` passes into failures. This distinguishes a small net rate change from uniformly better candidates.
- For M2, the Zumbach gate pass rate averaged `50.3%` and ranged from `27.3%` to `87.5%` across seeds. This is a warning against relying on one simulation seed.
- The tighter post-hoc checks are substantially less favourable: M2 reached H in [0.08, 0.12] in `4.2%` of rows, excess kurtosis above 5 in `81.6%`, and the stricter joint check in `0.0%`. Therefore the proxy joint rate must not be read as evidence that all target magnitudes have already been obtained.
- The lowest-distance M2 feedback-curvature bin is the clearest next search region: mean H was `0.07276`, excess kurtosis `5.947`, rank leverage `-0.002025`, and joint proxy acceptance `0.0%`. Rank Zumbach was `0.02114`, so this is a promising narrowing direction, not a complete solution.

## Files and continuation rule

Raw rows, per-seed checkpoints, seed/scenario/gate summaries, metric summaries, paired effects, and the complete manifest are stored beside this report. Keep an optional component only if its improvement exceeds between-seed uncertainty and does not materially degrade leverage, Zumbach, return ACF, or cap behaviour.

With only 5 seeds, uncertainty intervals are necessarily wide. A later market-data study must replace the proxy bands and retain an untouched validation period.

## Pilot-to-Stage-2 validation comparison

Pilot 1 was used only to choose the targeted box. Stage 2 uses a new scrambled design, new structural design, and disjoint cache seeds. Both studies below are re-scored under the same Stage-2 proxy gates. Because the sampling boxes differ, this compares targeted capacity, not global parameter-space volume.

| model | stage2_proxy_acceptance | rough_hurst | leverage_rank | zumbach_rank | excess_kurtosis |
| --- | --- | --- | --- | --- | --- |
| M2 | 0 → 0 | 0.2091 → 0.04081 | -0.007353 → -0.005818 | 0.0187 → 0.01524 | 2.372 → 8.499 |
| M3 | 0.006944 → 0 | 0.2098 → 0.04106 | -0.006775 → -0.005819 | 0.01538 → 0.01483 | 2.398 → 8.51 |
| M2+O | 0 → 0 | 0.2092 → 0.04823 | -0.007518 → -0.005309 | 0.01874 → 0.01489 | 2.373 → 8.751 |
| M4 | 0.006944 → 0 | 0.21 → 0.04842 | -0.006929 → -0.005312 | 0.0154 → 0.01447 | 2.399 → 8.762 |

## Configuration stability across Stage-2 seeds

A configuration is one structural scenario, readout candidate, and model cell. The rates below show how often the same configuration passes across seeds; this is stricter than pooling all rows.

| model | configs | pass ≥3/5 | pass ≥4/5 | pass 5/5 | best rate |
| --- | --- | --- | --- | --- | --- |
| M2 | 128 | 0.0% | 0.0% | 0.0% | 0.0% |
| M3 | 128 | 0.0% | 0.0% | 0.0% | 0.0% |
| M2+O | 128 | 0.0% | 0.0% | 0.0% | 0.0% |
| M4 | 128 | 0.0% | 0.0% | 0.0% | 0.0% |

## Stage-2 reading

- `M2` has the largest mean Stage-2 proxy pass rate at `0.0%`; the corresponding ≥4/5-seed stable-configuration rate is `0.0%`.
- Adding standardized orthogonal energy to M2 changes acceptance by `0` with a seed interval `[0, 0]`.
- Adding the echo to M2 changes acceptance by `0` with a seed interval `[0, 0]`.
- Retain an optional block only if its paired benefit is directionally stable and its leverage, Zumbach, roughness, tail, return-ACF, and cap diagnostics do not reveal a compensating deterioration.

## Selection boundary

The frozen latent architecture came from `experiments/pilot_001` and was not refitted. Stage-2 results may guide a later experiment, but any newly selected subregion or best candidate requires another untouched design before it can be called validated. Option-surface calibration remains untested because no SPX/VIX option data or pricing loss was supplied.
