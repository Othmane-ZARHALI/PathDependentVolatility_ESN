# Historical physical-statistic stability

`matteo.preprocessing.physical_stability` analyzes any of Othmane's six selected
stocks and three indices. It reads the ignored raw OHLC files, leaves them
unchanged, and writes one directory per asset here. The default comparison is
six non-overlapping five-calendar-year periods from 1994 through 2023. The
2019-2023 index period ends on the last available date, 2023-11-14.

From the repository root, with NumPy, pandas and openpyxl installed:

```sh
python -m matteo.preprocessing.physical_stability --assets SP500 \
  --replicates 1000 --block-lengths 126 63 252 \
  --rolling --rolling-window-years 5 10 --rolling-step-years 1 \
  --rolling-replicates 500 --rolling-block-length 126
```

Use `--assets all` for all nine assets or list any subset, for example
`--assets SP500 RUSSELL2000 NASDAQ100`. Runs for one asset use the same random
stream whether that asset is selected alone or with others. The seed, source
file hash, period definitions and bootstrap settings are recorded in each
`run.json`. `--first-year`, `--last-year`, `--period-years` and
`--include-partial` support other non-overlapping partitions. Compare assets
on a shared date range before interpreting differences. To examine 2024-2025
for stocks, run them separately with `--first-year 2024 --last-year 2025
--period-years 2`; this is too short for a historical stability test.

## Measurements

Close-to-close log returns are paired with same-day Garman-Klass variance from
open, high, low and close. Hurst uses the established four moment orders and
lags 2-40, applied to log Garman-Klass volatility. Excess kurtosis is the
existing finite-sample corrected daily-return statistic. Leverage is the mean
return-to-future-variance correlation over lags 1-5. Zumbach is the difference
of the two Pearson correlations between squared cumulative returns and mean
daily Garman-Klass variance at horizons 5, 10 and 20 trading days. Taylor is
the mean gap between absolute-return and squared-return ACFs at lags 1-20.
Clustering is the mean log-volatility ACF at lags 1-20. Return ACF is lag 1.

`proxy_sensitivity.csv` recomputes the same statistics with (a) Garman-Klass
plus the squared overnight open-to-previous-close log gap and (b) Parkinson's
high/low-only variance. These are measurement sensitivities, not intraday
realised variance. Return-only statistics are identical under all proxies.
`source_years.csv` tracks the fraction of days with open equal to the previous
close, which is unusually high in early index observations; this matters when
comparing open-dependent estimators across years.

## Uncertainty

For each period, the script resamples fixed-length contiguous blocks of
aligned daily returns and OHLC log ratios. It recomputes variance from the
resampled OHLC ratios. Lagged calculations exclude every pair or Zumbach
window that crosses an artificial block join. The 95% intervals in
`period_estimates.csv` are percentile intervals from the resampled estimates.
`bootstrap_covariance.csv` contains within-period sampling covariance of the
statistics for each block length. `source_checks.csv` gives actual dates,
observation counts, variance floors, recorded-open checks and the overnight
contribution.

With `--rolling`, `rolling_estimates.csv` contains every metric in annually
stepped 5- and 10-year windows with 95% **pointwise** bootstrap intervals.
`rolling_proxy_sensitivity.csv` repeats the point estimates under the two
other variance proxies, and `rolling_source_checks.csv` tracks recorded-open
and overnight-gap behavior in each window. Each rolling window is resampled
separately, so the output does not contain covariance between overlapping
windows or a global test of a flat rolling curve. Adjacent windows share most
of their observations; their intervals must not be treated as independent or
as additional calibration observations. Specify other window lengths, step or
replicate count with the corresponding `--rolling-*` options.

With Plotly installed, generate the standalone interactive chart after a
rolling run:

```sh
python -m matteo.preprocessing.plot_physical_stability --asset SP500
```

This writes `SP500/rolling_estimates.html`. The chart displays seven
statistics with pointwise interval bars and the recorded-open source check;
the CSV retains all nine statistics. Use `--input` to read another experiment
directory or `--output` to choose the HTML path.

`heterogeneity.csv` reports an inverse-sampling-variance Q statistic and a
centered bootstrap p-value for the hypothesis that all periods have the same
statistic. The test treats the disjoint periods as independent and assumes each
period is approximately stationary. The p-value is exploratory, particularly
for heavy tails and long volatility episodes. Failure to reject equal values
does **not** establish practical stability. Calibration needs a declared
tolerance for each statistic, proxy agreement, and a target covariance for the
selected stable vector. The physical-measure series ends before the 2026-02-19
option snapshot, so the study does not supply that date's latent state.

The daily OHLC Zumbach measure is a proxy for the PDF's intraday realised-
variance statistic. Daily OHLC roughness can also be affected by volatility
measurement error. These require intraday data for stronger conclusions.

## Three-year rolling sensitivity

The [combined three-, five- and ten-year chart](rolling_3y/SP500/rolling_estimates.html)
overlays all three window lengths with shaded pointwise 95% intervals. Its
[three-year estimates](rolling_3y/SP500/rolling_estimates.csv) recompute each
statistic from the daily data in 28 annually stepped three-year windows,
rather than averaging estimates from longer windows. The five-/ten-year series
come from the main output above. Reproduce the combined chart without replacing
that output:

```sh
python -m matteo.preprocessing.physical_stability --assets SP500 \
  --replicates 1000 --block-lengths 126 63 252 \
  --rolling --rolling-window-years 3 --rolling-step-years 1 \
  --rolling-replicates 500 --rolling-block-length 126 \
  --output matteo/exogenous_reservoir_eq_m2/experiments/physical_stability_001/rolling_3y
python -m matteo.preprocessing.plot_physical_stability --asset SP500 \
  --input matteo/exogenous_reservoir_eq_m2/experiments/physical_stability_001/rolling_3y \
  --additional-input matteo/exogenous_reservoir_eq_m2/experiments/physical_stability_001 \
  --interval-style band
```

The shaded intervals are pointwise. Adjacent three-year windows share roughly
two years of observations, and 126-day bootstrap blocks give only about six
blocks per window. Tail and time-asymmetry estimates can therefore be
especially imprecise; use the longer windows and block-length sensitivity for
calibration decisions.
