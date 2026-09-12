# Exogenous-reservoir exponential–quadratic volatility model

This package is the executable research version of the current model descended from ESN-A2. It combines exact multirate OU state transitions, an optional fixed contractive echo-state layer, and a smoothly tempered exponential–quadratic variance readout.

The proposed experiment is rational, with one necessary qualification: full return and volatility stylised facts are generally **not** functions of hyperparameters alone. Rate grids and fixed projection profiles can be selected from latent roughness and clustering diagnostics, but tails, leverage, Zumbach, Taylor, and even measured post-readout roughness depend on the readout and its interactions. The package therefore uses a hierarchical design.

1. Calibrate architecture-level latent feature shapes.
2. Draw one cache of primitive Gaussian randomness.
3. Optionally sample state-affecting correlations and deterministic risk premia in an outer design.
4. Sample cheap exponential–quadratic readout parameters in an inner design.
5. Report the **joint acceptance rate and accepted region**, not only the best parameter vector.
6. Repeat the entire design over independent cache seeds before drawing a model-level conclusion.

The default numerical bands and parameter ranges are illustrative research placeholders. They are not empirical SPX targets.

## Implemented model

For feedback, spike, clustering, and optional orthogonal banks, each standardised OU coordinate satisfies under the physical measure

\[
dX_{j,t}^{b,P}=-\lambda_j^b X_{j,t}^{b,P}\,dt+\sqrt{2\lambda_j^b}\,dW_{j,t}^{b,P}.
\]

The feedback and spike banks use drivers partly correlated with the equity Brownian motion. Slow clustering modes and optional orthogonal modes use independent primitives by default. The code samples the joint equity/OU transition exactly rather than with Euler discretisation. Variance at a fine step uses the state at the **start** of the step, so the stochastic integrand is adapted.

Fixed projections produce \(Y_t\), \(F_t\), and \(C_t\), each normalised at the stationary P reference point. For the optional orthogonal bank, the recommended energy feature is centered and stationary-standardised,

\[
E_t^O=\frac{N_O^{-1}\sum_{j=1}^{N_O}(X_{j,t}^O)^2-1}{\sqrt{2/N_O}}.
\]

The score is

\[
\eta_t=c_C C_t+q_Z(Y_t-\beta)^2+q_JF_t^2+q_OE_t^O+\delta Z_t.
\]

The recommended M2 model uses \(\delta=q_O=0\). M2+O activates only \(q_O\), M3 activates only \(\delta\), and M4 activates both. If enabled, \(Z_t\) is produced by a fixed bounded echo-state layer satisfying the sufficient contraction bound

\[
(1-\alpha)+\alpha\lVert A\rVert_2<1.
\]

The production multiplier and variance are

\[
h_{L,\kappa}(x)=L-\kappa^{-1}\log(1+e^{\kappa(L-x)}),\qquad
G(x)=e^{2h_{L,\kappa}(x)},
\]

\[
V_t=\underline v+(\xi_0(t)-\underline v)\frac{G(\eta_t)}{M_Q(t)},
\qquad M_Q(t)=E_0^Q[G(\eta_t)].
\]

`M_Q` is estimated cross-sectionally from the matched Q path ensemble. The simulator records the maximum forward-mean normalisation error and the P/Q cap exceedance fractions for every candidate.

A constant deterministic Girsanov kernel changes only OU conditional means. The code applies the exact mean shift from the common initial state and reuses the same centred Gaussian innovations under P and Q.

## Parameter hierarchy and path reuse

| Class | Examples | What is reused? |
|---|---|---|
| Architecture hyperparameters | dimensions, OU rates, projection powers, fixed ESN matrices | Changing these rebuilds features from the same primitive random numbers. |
| Outer structural parameters | feedback/spike correlations | The equity primitive path remains bitwise identical; affected OU states are transformed again. |
| Measure-change parameters | deterministic equity and unspanned risk prices | No new randomness; Q means and, if active, Q echo features are recomputed. |
| Inner readout parameters | \(c_C,q_Z,\beta,q_J,q_O,\delta,L,\kappa\) | P/Q OU and echo features are bitwise unchanged; only score, variance, returns, and diagnostics are recomputed. |
| Daily inputs | forward variance, rates, dividends | No new randomness; pricing outputs are recomputed. |

Path reuse never means that final stock, variance, VIX, or payoff paths remain unchanged after a volatility parameter changes. It means that the expensive exogenous randomness and permitted state features are reused.

## Installation and tests

From this directory:

```bash
python -m pip install -e .
python -m unittest discover -s tests -p 'test_*.py' -v
```

The suite contains analytical unit tests for exact covariance, deterministic P/Q shifts, cap stability, echo-state contraction, sampling reproducibility, and statistic orientation, plus end-to-end integration tests for path reuse and hierarchical exploration.

## Run a reproducible exploration

After installation:

```bash
esn-eq-explore \
  --samples 64 \
  --structural-samples 8 \
  --years 8 \
  --burn-years 2 \
  --paths 8 \
  --sampler latin_hypercube \
  --seed 2026 \
  --output exploration_results
```

This evaluates \(8\times64=512\) joint candidates. The command writes:

- `candidates.csv`: every candidate, its parameters, diagnostic means, Monte Carlo standard errors, audit metadata, joint acceptance flag, and continuous distance;
- `rank_sensitivities.csv`: marginal Spearman screens for the outer and inner parameters;
- `summary.json`: candidate count and joint acceptance rate.

Latin hypercube is the default. Scrambled Sobol sampling is also available. Both give better domain coverage than naive independent uniform draws at the same budget.

`ParameterSpace.plausible_m2()` fixes the inactive ESN and orthogonal loadings at zero. Use `plausible_m3()` only with an enabled echo layer, and `plausible_m4()` only after adding a non-empty orthogonal bank. This avoids wasting sampling dimensions on parameters that cannot affect the selected architecture.

`ParameterSpace.targeted_stage2()` reproduces the deliberately aggressive curvature study. `ParameterSpace.stage3_local()` is the independently validated local region around the slower-spike/high-shift solution; it fixes the rejected echo loading at zero.

## Run the paired M2/M3/M2+O/M4 experiment

The factorial runner is the preferred way to decide whether the optional mechanisms earn their complexity:

```bash
esn-eq-factorial \
  --candidate-samples 16 \
  --structural-samples 3 \
  --seeds 104729 130363 155921 \
  --years 5 \
  --burn-years 1 \
  --paths 8 \
  --steps-per-observation 2 \
  --output factorial_pilot
```

One scrambled full-M4 parameter design is generated and held fixed across cache seeds. For every full draw, the runner evaluates exact-zero controls:

| Cell | Echo loading | Orthogonal curvature |
|---|---:|---:|
| M2 | (0) | (0) |
| M3 | sampled (delta) | (0) |
| M2+O | (0) | sampled (q_O) |
| M4 | sampled (delta) | sampled (q_O) |

All four cells use the same primitive paths, structural scenario, and core readout parameters. The echo feature consumes only the common feedback, spike, and clustering banks by default. Consequently, enabling the orthogonal bank does not redefine the echo map, and the two main effects and their interaction have a clean factorial interpretation.

The output directory is a cumulative experiment ledger:

- `checkpoints/seed_<seed>.csv` is written after each completed seed;
- `runs.csv` retains every seed, scenario, candidate, and model cell;
- `seed_summary.csv` exposes simulation-seed variation;
- `metric_summary.csv` reports seed-level means, standard deviations, Student intervals, and pooled parameter-design quantiles;
- `paired_effects.csv` and `effect_summary.csv` report conditional echo effects, conditional orthogonal effects, and their interaction;
- `criterion_summary.csv`, `scenario_summary.csv`, and `posthoc_threshold_sensitivity.csv` expose gate bottlenecks, structural sensitivity, and explicitly labelled stricter diagnostic thresholds;
- `feedback_curvature_screen.csv` provides a non-causal M2 quartile screen for narrowing the next readout parameter box;
- `manifest.json` records the complete architecture, simulation, parameter design, structural scenarios, and acceptance bands;
- `RESULTS.md` gives a compact comparison while preserving the interpretation boundary.

Three seeds are a sensible pilot minimum, not a final uncertainty analysis. Keep the parameter design fixed across these seeds to isolate Monte Carlo variation. After narrowing the parameter box, a confirmatory study should also change the scrambled design seed and increase both the cache-seed count and Q-normalisation path count.

## Current multi-seed evidence

The broad Pilot 1 and aggressive Stage 2 boxes did not jointly satisfy the more demanding proxy gates. Stage 2 generated strong tails but made volatility too rough at the tested daily band. Training-only mechanism probes then identified two complementary changes:

- slow the spike lift from rates 20–2000 to 10–500 inverse years and weight its slower coordinates with power -0.5;
- keep feedback curvature low while increasing the shift, which preserves the rough component while supplying a sufficiently negative linear leverage coefficient.

Stage 3 froze those choices before using a new Sobol design, four new structural scenarios, five disjoint cache seeds, 24 paths, and six simulated years. Under the unchanged synthetic proxy gates, M2 achieved 62.0% row acceptance and M2+O 64.7%. For M2, 75.8% of configurations passed at least three of five seeds, 69.5% passed at least four, and five configurations passed all five. Fixed-design means were:

| Model | Rough H | Rank leverage | Rank Zumbach | Excess kurtosis | Max return ACF |
|---|---:|---:|---:|---:|---:|
| M2 | 0.0754 | -0.0175 | 0.0416 | 9.84 | 0.0754 |
| M2+O | 0.0777 | -0.0162 | 0.0395 | 9.90 | 0.0757 |

M2 remains the recommended baseline. The standardized orthogonal block is a defensible optional roughness stabilizer: it raised acceptance by 2.66 percentage points in the paired study, but slightly weakened leverage and Zumbach and only modestly increased stable-region volume. The sampled echo block had negligible benefit in Stage 2 and is not recommended.

Reproduce the two latest studies with:

```bash
esn-eq-stage2 --output experiments/stage2_001
esn-eq-stage3 --output experiments/stage3_001
```

The experiment ledgers contain every raw row, seed checkpoint, sampling bound, paired contrast, stronger-Zumbach sensitivity, and same-configuration stability result.

The differential-evolution latent calibration is also replicated over three optimizer seeds by default. The selected architecture is the smallest-objective replication, while `latent_calibration.json` and `latent_calibration_replicates.csv` retain every fit, convergence message, fitted coordinate, and search-bound hit. This prevents a stable simulation study from being built on an unreported optimizer-seed accident.

## Python workflow

```python
from esn_eq import (
    AcceptanceCriteria,
    ArchitectureConfig,
    HierarchicalExplorer,
    LatentArchitectureCalibrator,
    LatentTargets,
    MarketEnvironment,
    ParameterSpace,
    QmcParameterSampler,
    QmcStructuralSampler,
    SimulationConfig,
)

# 1. Select latent feature-engine hyperparameters.
base = ArchitectureConfig()
latent_fit = LatentArchitectureCalibrator().fit(
    base,
    LatentTargets(rough_hurst=0.10),
    seed=100,
)
architecture = latent_fit.architecture

# 2. Build space-filling outer and inner designs.
structures = QmcStructuralSampler().sample(8, seed=101)
readouts = QmcParameterSampler(ParameterSpace.plausible_m2()).sample(
    64,
    seed=102,
    method="latin_hypercube",
)

# 3. Replace these illustrative bands with market-estimated intervals.
criteria = AcceptanceCriteria.illustrative()
experiment = HierarchicalExplorer(
    architecture,
    SimulationConfig(years=8, burn_years=2, paths=8),
    MarketEnvironment(forward_variance=0.20**2),
    criteria,
)
result = experiment.run(structures, readouts, seed=103)
result.save("exploration_results")
```

For daily output, leave `observations_per_year=252`. For intraday output, set it to the desired number of observations per year and express all OU rates consistently in inverse years. The empirical Zumbach target should use intraday realised variance; daily squared returns are not an adequate substitute.

Frequency-dependent lags are explicit in `DiagnosticConfig`. `aggregate_result` can reuse one fine-frequency simulation to examine tails and dependence after coarser aggregation, provided the aggregation factor divides `observations_per_year`.

## Diagnostics

Each path produces:

- maximum absolute return autocorrelation;
- mean log-variance autocorrelation;
- finite-band roughness estimate;
- Pearson and rank leverage;
- Pearson and rank Zumbach asymmetry using matched realised-variance windows;
- excess kurtosis, a robust absolute-return tail ratio, and multi-threshold Hill summary;
- Taylor gap and fraction of positive Taylor lags;
- mean and 99.5% volatility;
- Monte Carlo standard errors across independent simulated paths.

The Hill estimate is included as a finite-range diagnostic only. Tempering ultimately truncates the variance mixture, so it would be wrong to interpret it as proof of an asymptotic power law.

## How to judge whether the model is “generally good”

A large parameter box can make any pass rate look poor, while a hand-picked narrow box can make it look excellent. The research conclusion should therefore report all of the following:

- parameter bounds and sampling measure;
- joint pass rate with confidence intervals across cache seeds;
- connectedness and volume of accepted regions;
- metric trade-offs, especially tails/Taylor versus leverage/Zumbach;
- cap-binding rate and Q-forward normalisation error;
- sensitivity to observation frequency and aggregation horizon;
- out-of-sample stability after hyperparameter selection;
- ablations M0, M1, M2, then M3/M4 only if justified.

Do not calibrate hyperparameters and assess capacity on the same simulated or empirical window. Use nested or rolling validation. A useful next stage is to compare M0–M2 on synthetic surfaces generated by a quadratic-rough-Heston and a quintic-Gaussian benchmark before real SPX/VIX option data are available.

## Current limitations

- The package is a physical-measure capacity and path-reuse research engine, not yet a production SPX/VIX option pricer.
- `M_Q` is a matched-path Monte Carlo normaliser. Conditional VIX valuation still needs quadrature, regression, quantisation, or a controlled surrogate.
- The default target bands and parameter ranges are placeholders, not estimates from market data.
- The successful Stage-3 evidence is daily-frequency synthetic capacity evidence; intraday robustness has not yet been tested.
- A finite OU lift is rough only on its designed scale band; it has Brownian asymptotics at sufficiently small lags.
- A positive Zumbach statistic is a numerical mechanism result, not a theorem guaranteeing the empirical magnitude.
- The reported Q forward-mean error is an in-sample algebraic check because the same paths estimate and verify the cross-sectional normaliser; it is not an out-of-sample integration-error estimate.
- The nonlinear ESN layer and orthogonal dispersion block are disabled by default and should be activated only after ablation evidence.
