# Paired Pearson-Zumbach mechanism ablation

This experiment supplies the reproducible results reported in Section 4.5 of the model
specification. It uses ten declared seeds, 32 paths per seed, ten simulated years with a
two-year burn-in, 252 daily observations per year, and two integration steps per day.

For every model path, the ordinary Pearson Zumbach statistic is evaluated at 5-, 10-, and
20-day windows using model-integrated daily average variance as the noise-free realized-
variance limit. The three values are averaged within path, the 32 paths are averaged within
seed, and uncertainty is computed across the ten seed means. Model ablations reuse identical
primitive Gaussian draws.

The table below reports the across-seed mean, sample standard deviation, and two-sided 95%
Student-t confidence interval.

| Configuration | Mean | Between-seed SD | 95% CI |
| --- | ---: | ---: | ---: |
| Recommended model | 0.180227 | 0.013944 | [0.170252, 0.190202] |
| Feedback quadratic removed | 0.254766 | 0.009016 | [0.248316, 0.261215] |
| Fast quadratic removed | 0.028513 | 0.011028 | [0.020624, 0.036401] |
| Both correlated quadratics removed | -0.002658 | 0.006025 | [-0.006968, 0.001652] |
| Return-state correlations removed | -0.001116 | 0.015253 | [-0.012027, 0.009796] |
| GARCH(1,1) positive control | 0.156712 | 0.008153 | [0.150880, 0.162544] |

Run from the package root with a Python 3.10+ environment containing the declared project
dependencies:

```bash
PYTHONPATH=src python experiments/pearson_zumbach_001/run.py
```

The runner rewrites `per_seed.csv`, `summary.csv`, and `manifest.json` deterministically.
