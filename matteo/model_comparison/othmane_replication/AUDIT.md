# Initial source audit and SP500 baseline

Date: 2026-09-30. This records source observations and a small reproducible
baseline. It is not a claim about which model fits historical markets better.

## Inputs and reproducibility

All six ignored `data/processed/bench9_real_*.json` files were rebuilt in
memory from the local raw OHLC files with `matteo.preprocessing` and compared
with the saved values. The JSON objects matched exactly. The six stock series
have 8,054 OHLC rows each; the three index series have 7,521 each. No target
file was overwritten. The input hashes for the SP500 run are in
`sp500_exploratory_001.json`.

The nine-asset study is an actual historical-target fit. Its source is
`SimultanuousAllFeatures/joint_calibration.py`, the saved fits are in
`joint_calibration_final.json`, and the published write-up is
`joint_calibration_protocol.tex`. The separate `CrossSectionalData/` scripts
use a 2026-02-19 SPX/VIX option snapshot. M2's checked-in Stage-3 work uses
synthetic data, so it is a model-capacity result rather than an empirical fit.

## Confirmed source defects

| Finding | Evidence | Effect |
| --- | --- | --- |
| Five values are passed to a six-dimensional optimizer | `starting_point()` returns `(H0, 0.4, 0.03, lh0, 0.0)`; `calibrate_joint()` unpacks six and defines six lower/upper bounds. A direct call with five values raises `ValueError: Inconsistent shapes between bounds and x0`. | The checked-in `__main__` cannot regenerate the final JSON. |
| Residual dimension is mislabeled | The blocks contain `1+40+6+6+5+1=59` values. Direct evaluation returns 59. The invalid-`kappa0` branch constructs 53 values, and the report calls the residual 53-dimensional. | The branch would give an inconsistent shape if visited. Under the current negative rough orientation and positive rough scale, `kappa0` appears strictly negative, so this branch is currently unreachable. The reported dimension is still wrong. |
| Stage-2 kurtosis step uses the wrong feedback setting | `bisect_m1_for_kurtosis()` has `alpha=1.0` as its default. The call in `__main__` does not pass fitted `c["alpha_fit"]`, although saved final alpha values are below 1. | If run as written after fixing the starting point, Stage 2 would optimize `m1` for a different simulator from the fitted `(m1, alpha)` pair. The saved JSON has uncertain code provenance, so its `m1` cannot be attributed to this call without a rerun. |
| Final serializer omits fitted alpha | The current `final[name]` dictionary has only five fitted fields. The archived `joint_calibration_final.json` has a sixth, `alpha_fit`. | The archived final JSON came from a different script state or manual edit; the checked-in entry point cannot produce it verbatim. |
| Roughness lag bands differ | `matteo.preprocessing.hurst.REAL_LAGS` is 2–40; `HURST_LAGS_ESN` is 1–40. | Historical and model H estimates use different finite bands. This may be deliberate proxy-noise handling, but it needs sensitivity analysis before fair fitting. |

## Exploratory paired evaluation

Command: `compare.py --asset SP500 --models both --days 1000 --burn-days 252
--paths 8 --seed 93271 --output sp500_exploratory_001.json`.

| Saved parameter set | Weighted squared error | H estimate | Excess kurtosis |
| --- | ---: | ---: | ---: |
| Historical target | — | 0.08576 | 10.31395 |
| Othmane archived SP500 fit | 12.5758 | 0.07912 | 4.08142 |
| M2 Stage-3 synthetic anchor, candidate 18/scenario 0 | 21.0944 | 0.07649 | 3.02169 |

The values above come from eight simulated paths of 1,000 retained days,
using Othmane's statistic extractor and objective for both models. The M2
anchor has not been fitted to historical SP500. Its Q variance normalizer is
estimated on the same eight paths, a small budget. Othmane's published
uncertainty bands used longer paths and more replicates. Neither model's
current variance output is the OHLC Garman–Klass proxy used to build the
historical H, leverage and Zumbach targets. These scores serve only as an
execution and convention check. The individual block errors are in the JSON.

The baseline was run with a Python 3.12 environment containing NumPy 2.0.2,
SciPy 1.18.1 and pandas 3.0.5. The in-memory target rebuild additionally used
openpyxl 3.1.5 for the index workbook. The source target files and archived
model outputs were not changed.
