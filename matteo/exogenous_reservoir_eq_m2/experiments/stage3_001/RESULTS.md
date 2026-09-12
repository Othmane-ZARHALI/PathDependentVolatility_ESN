# Stage-3 confirmation of the slower-spike/high-shift ridge

## Status and interpretation boundary

This is a physical-measure capacity study on synthetic model paths. The acceptance bands are declared proxies, not empirical SPX targets. Accordingly, acceptance rates compare the four cells under one declared sampling measure; they do not yet validate the model against market data.

## Experimental design

- Cache seeds: `[360007, 390001, 420001, 450007, 480013]`.
- Full readout draws: `32` using `sobol` with design seed `20261201`.
- Parameter space: `stage3_local_slow_spike_ridge`.
- Structural scenarios: `4`.
- Simulation: `{'years': 6.0, 'burn_years': 1.0, 'paths': 24, 'observations_per_year': 252, 'steps_per_observation': 2}`.
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
| 360007 | 82.8% | 82.8% | 87.5% | 87.5% |
| 390001 | 71.9% | 71.9% | 71.9% | 71.9% |
| 420001 | 3.9% | 3.9% | 5.5% | 5.5% |
| 450007 | 75.8% | 75.8% | 80.5% | 80.5% |
| 480013 | 75.8% | 75.8% | 78.1% | 78.1% |

## Acceptance by structural scenario

| scenario | M2 | M3 | M2+O | M4 |
| --- | --- | --- | --- | --- |
| 0 | 63.1% ± 34.0pp | 63.1% ± 34.0pp | 66.2% ± 34.2pp | 66.2% ± 34.2pp |
| 1 | 63.7% ± 34.3pp | 63.7% ± 34.3pp | 65.6% ± 35.6pp | 65.6% ± 35.6pp |
| 2 | 61.3% ± 32.6pp | 61.3% ± 32.6pp | 64.4% ± 32.8pp | 64.4% ± 32.8pp |
| 3 | 60.0% ± 30.2pp | 60.0% ± 30.2pp | 62.5% ± 31.9pp | 62.5% ± 31.9pp |

## Stylised-fact comparison

Entries are fixed-design means across seeds; the value after `±` is the between-seed standard deviation.

| model | max \|return ACF\| | mean log-variance ACF | rough H | rank leverage | rank Zumbach | excess kurtosis | absolute-tail ratio | Taylor gap | P cap fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| M2 | 0.07535 ± 0.004896 | 0.209 ± 0.02073 | 0.07543 ± 0.007044 | -0.01746 ± 0.006087 | 0.04158 ± 0.01294 | 9.837 ± 1.341 | 3.67 ± 0.0612 | 0.02899 ± 0.005383 | 9.413e-05 ± 3.558e-05 |
| M3 | 0.07535 ± 0.004896 | 0.209 ± 0.02073 | 0.07543 ± 0.007044 | -0.01746 ± 0.006087 | 0.04158 ± 0.01294 | 9.837 ± 1.341 | 3.67 ± 0.0612 | 0.02899 ± 0.005383 | 9.413e-05 ± 3.558e-05 |
| M2+O | 0.07572 ± 0.004265 | 0.2176 ± 0.02184 | 0.07772 ± 0.006648 | -0.01618 ± 0.005748 | 0.03953 ± 0.011 | 9.904 ± 1.404 | 3.674 ± 0.0538 | 0.03076 ± 0.005286 | 9.744e-05 ± 4.091e-05 |
| M4 | 0.07572 ± 0.004265 | 0.2176 ± 0.02184 | 0.07772 ± 0.006648 | -0.01618 ± 0.005748 | 0.03953 ± 0.011 | 9.904 ± 1.404 | 3.674 ± 0.0538 | 0.03076 ± 0.005286 | 9.744e-05 ± 4.091e-05 |

## Individual gate pass rates

These rates identify the bottlenecks hidden by joint acceptance.

| model | max_abs_return_acf | mean_log_variance_acf | rough_hurst | leverage_rank | zumbach_rank | excess_kurtosis | taylor_gap |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M2 | 100.0% ± 0.0pp | 98.8% ± 1.7pp | 66.7% ± 25.3pp | 95.8% ± 7.8pp | 97.0% ± 5.4pp | 100.0% ± 0.0pp | 100.0% ± 0.0pp |
| M3 | 100.0% ± 0.0pp | 98.8% ± 1.7pp | 66.7% ± 25.3pp | 95.8% ± 7.8pp | 97.0% ± 5.4pp | 100.0% ± 0.0pp | 100.0% ± 0.0pp |
| M2+O | 100.0% ± 0.0pp | 98.8% ± 1.7pp | 72.5% ± 21.6pp | 92.5% ± 12.7pp | 97.5% ± 4.4pp | 100.0% ± 0.0pp | 100.0% ± 0.0pp |
| M4 | 100.0% ± 0.0pp | 98.8% ± 1.7pp | 72.5% ± 21.6pp | 92.5% ± 12.7pp | 97.5% ± 4.4pp | 100.0% ± 0.0pp | 100.0% ± 0.0pp |

## Post-hoc threshold sensitivity

These deliberately stricter thresholds are diagnostics, not preregistered acceptance criteria or empirical SPX estimates. They test whether the declared proxy gates are masking weak magnitudes.

| model | rough H in [0.08, 0.12] | rank leverage < -0.01 | rank Zumbach > 0.02 | rank Zumbach > 0.05 | excess kurtosis > 3 | excess kurtosis > 5 | joint H/leverage/Z>0.02/kurtosis>3 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M2 | 30.8% ± 21.3pp | 88.3% ± 14.4pp | 92.8% ± 11.9pp | 28.9% ± 32.0pp | 100.0% ± 0.0pp | 100.0% ± 0.0pp | 24.4% ± 23.1pp |
| M3 | 30.8% ± 21.3pp | 88.3% ± 14.4pp | 92.8% ± 11.9pp | 28.9% ± 32.0pp | 100.0% ± 0.0pp | 100.0% ± 0.0pp | 24.4% ± 23.1pp |
| M2+O | 40.2% ± 19.3pp | 83.0% ± 23.3pp | 92.2% ± 10.1pp | 23.8% ± 25.9pp | 100.0% ± 0.0pp | 100.0% ± 0.0pp | 32.0% ± 23.0pp |
| M4 | 40.2% ± 19.3pp | 83.0% ± 23.3pp | 92.2% ± 10.1pp | 23.8% ± 25.9pp | 100.0% ± 0.0pp | 100.0% ± 0.0pp | 32.0% ± 23.0pp |

## Exploratory M2 feedback-curvature screen

This post-hoc quartile screen is descriptive, not causal. It is included to choose the next parameter box rather than to validate a fitted effect.

| feedback curvature | rough_hurst | excess_kurtosis | leverage_rank | zumbach_rank | accepted |
| --- | --- | --- | --- | --- | --- |
| [0.06014, 0.06624] | 0.08508 ± 0.007636 | 8.803 ± 1.245 | -0.01564 ± 0.006505 | 0.05176 ± 0.01328 | 0.7812 ± 0.4227 |
| [0.06692, 0.07175] | 0.07857 ± 0.007252 | 9.419 ± 1.317 | -0.01672 ± 0.006248 | 0.04536 ± 0.01308 | 0.825 ± 0.3913 |
| [0.07252, 0.07862] | 0.07138 ± 0.006818 | 10.22 ± 1.388 | -0.01828 ± 0.005915 | 0.03743 ± 0.01288 | 0.475 ± 0.2825 |
| [0.07949, 0.08432] | 0.06668 ± 0.006481 | 10.9 ± 1.422 | -0.01919 ± 0.005689 | 0.0318 ± 0.01275 | 0.4 ± 0.2279 |

## Paired optional-component effects

Contrasts use identical primitive paths, structural scenarios, and core readout draws. Positive acceptance effects are improvements; for distance, negative effects are improvements. Other signs retain the metric's meaning.

| effect | acceptance_rate | distance | rough_hurst | leverage_rank | zumbach_rank | excess_kurtosis |
| --- | --- | --- | --- | --- | --- | --- |
| echo_at_O0 | 0 ± 0 | 0 ± 0 | 0 ± 0 | 0 ± 0 | 0 ± 0 | 0 ± 0 |
| echo_at_O1 | 0 ± 0 | 0 ± 0 | 0 ± 0 | 0 ± 0 | 0 ± 0 | 0 ± 0 |
| echo_orthogonal_interaction | 0 ± 0 | 9.675e-24 ± 1.913e-23 | 0 ± 0 | 1.355e-21 ± 3.03e-21 | 0 ± 0 | 0 ± 0 |
| orthogonal_at_echo0 | 0.02656 ± 0.02037 | -0.001224 ± 0.001547 | 0.002286 ± 0.0007494 | 0.00128 ± 0.0008811 | -0.002059 ± 0.002262 | 0.06707 ± 0.1938 |
| orthogonal_at_echo1 | 0.02656 ± 0.02037 | -0.001224 ± 0.001547 | 0.002286 ± 0.0007494 | 0.00128 ± 0.0008811 | -0.002059 ± 0.002262 | 0.06707 ± 0.1938 |

## Evidence summary

- Joint acceptance was M2 `62.0%`, M3 `62.0%`, M2+O `64.7%`, and M4 `64.7%` under the declared proxy bands.
- The echo-only contrast changed distance by `0 (95% seed CI [0, 0])`. It changed rank leverage by `0 (95% seed CI [0, 0])` and rank Zumbach by `0 (95% seed CI [0, 0])`; signs must be judged against the desired negative leverage and positive Zumbach directions.
- The orthogonal-only distance contrast was `-0.001224 (95% seed CI [-0.003145, 0.0006966])`. An effect of this size should not justify the additional bank unless it becomes stable under a broader loading range or a targeted conditional test.
- At the row level, the echo changed `0` failures into passes and `0` passes into failures. The orthogonal block changed `25` failures into passes and `8` passes into failures. This distinguishes a small net rate change from uniformly better candidates.
- For M2, the Zumbach gate pass rate averaged `97.0%` and ranged from `87.5%` to `100.0%` across seeds. The sign gate is comparatively stable across the declared seeds.
- The tighter post-hoc checks are substantially less favourable: M2 reached H in [0.08, 0.12] in `30.8%` of rows, excess kurtosis above 5 in `100.0%`, and the stricter joint check in `24.4%`. Therefore the proxy joint rate must not be read as evidence that all target magnitudes have already been obtained.
- The lowest-distance M2 feedback-curvature bin is the clearest next search region: mean H was `0.08508`, excess kurtosis `8.803`, rank leverage `-0.01564`, and joint proxy acceptance `78.1%`. Rank Zumbach was `0.05176`, so this is a promising narrowing direction, not a complete solution.

## Files and continuation rule

Raw rows, per-seed checkpoints, seed/scenario/gate summaries, metric summaries, paired effects, and the complete manifest are stored beside this report. Keep an optional component only if its improvement exceeds between-seed uncertainty and does not materially degrade leverage, Zumbach, return ACF, or cap behaviour.

With only 5 seeds, uncertainty intervals are necessarily wide. A later market-data study must replace the proxy bands and retain an untouched validation period.

## What Stage 3 changed

Stage 2 selected a slower spike grid (10–500 inverse years), spike projection power -0.5, and a high-shift/low-feedback-curvature neighborhood. Stage 3 freezes those choices, uses a new Sobol design and disjoint cache seeds, and raises the Q-normalisation ensemble to 24 paths. The echo loading is fixed at zero, so M3 duplicates M2 and M4 duplicates M2+O by construction.

## Same-gate Stage-2/Stage-3 comparison

| model | stage2_proxy_acceptance | rough_hurst | leverage_rank | zumbach_rank | excess_kurtosis |
| --- | --- | --- | --- | --- | --- |
| M2 | 0 → 0.6203 | 0.04081 → 0.07543 | -0.005818 → -0.01746 | 0.01524 → 0.04158 | 8.499 → 9.837 |
| M2+O | 0 → 0.6469 | 0.04823 → 0.07772 | -0.005309 → -0.01618 | 0.01489 → 0.03953 | 8.751 → 9.904 |

## Stability of the same configuration across seeds

| model | configs | pass ≥3/5 | pass ≥4/5 | pass 5/5 | best rate |
| --- | --- | --- | --- | --- | --- |
| M2 | 128 | 75.8% | 69.5% | 3.9% | 100.0% |
| M2+O | 128 | 78.9% | 70.3% | 5.5% | 100.0% |

## Strong-Zumbach sensitivity

| model | Zumbach > 0.03 | Zumbach > 0.05 | joint with Zumbach > 0.03 | joint with Zumbach > 0.05 |
| --- | --- | --- | --- | --- |
| M2 | 75.2% | 28.9% | 54.1% | 26.2% |
| M2+O | 73.3% | 23.8% | 55.5% | 21.2% |

## Best untouched configurations

| model | candidate | scenario | seed pass | H | leverage | Zumbach | kurtosis |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M2 | 18 | 0 | 100.0% | 0.08445 | -0.01513 | 0.04805 | 10.28 |
| M2+O | 18 | 0 | 100.0% | 0.08484 | -0.01473 | 0.04718 | 10.26 |

## Confirmation reading

- M2: row acceptance `62.0%`, neighborhood ≥4/5-seed stability `69.5%`, best configuration `100.0%`.
- M2+O: row acceptance `64.7%`, neighborhood ≥4/5-seed stability `70.3%`, best configuration `100.0%`.
- The paired standardized-O acceptance effect is `0.02656` with seed interval `[0.001267, 0.05186]`.
- Recommend the simpler M2 unless M2+O shows a stable joint-rate gain rather than only a marginal roughness increase. M3 remains rejected for this stage.

The gates remain synthetic capacity proxies. They do not replace empirical estimation, option-surface calibration, conditional VIX pricing validation, or an out-of-time market-data test.
