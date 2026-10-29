# SPX physical-statistic uncertainty, 1994-2023

Run: 1,000 moving-block bootstrap replicates per disjoint five-calendar-year
period, with 126-day blocks as the displayed interval and 63/252-day block
sensitivities. The paired daily close-to-close return and OHLC log ratios were
resampled together. Lagged pairs and Zumbach windows crossing resampled-block
joins were excluded. The primary volatility proxy is Othmane's daily
Garman-Klass variance. See `SP500/run.json` for the source hash and full run
configuration and `README.md` for definitions. The last period ends on
2023-11-14, the final index date.

## SPX estimates and 95% bootstrap intervals

| Period | H from log GK volatility | Daily excess kurtosis | Leverage mean, lags 1-5 | 10-day Zumbach proxy |
| --- | ---: | ---: | ---: | ---: |
| 1994-1998 | 0.045 [0.017, 0.076] | 8.49 [1.33, 12.86] | -0.141 [-0.186, -0.063] | -0.050 [-0.116, 0.144] |
| 1999-2003 | 0.067 [0.044, 0.097] | 1.22 [0.68, 1.73] | -0.110 [-0.152, -0.067] | 0.187 [0.046, 0.270] |
| 2004-2008 | 0.077 [0.031, 0.121] | 15.72 [1.31, 23.64] | -0.171 [-0.210, -0.102] | 0.354 [-0.112, 0.464] |
| 2009-2013 | 0.080 [0.055, 0.116] | 4.16 [1.56, 5.52] | -0.103 [-0.163, -0.044] | 0.125 [-0.202, 0.342] |
| 2014-2018 | 0.097 [0.060, 0.128] | 3.79 [1.67, 5.81] | -0.203 [-0.227, -0.136] | 0.129 [-0.179, 0.325] |
| 2019-2023 | 0.128 [0.062, 0.213] | 13.99 [1.28, 19.41] | -0.187 [-0.234, -0.097] | 0.278 [-0.037, 0.321] |

The intervals are percentile bootstrap intervals. Skewed bootstrap
distributions, especially for kurtosis and Zumbach, are visible in the
`bootstrap_mean` column of `period_estimates.csv`; the interval should not be
interpreted as symmetric around the point estimate. The 2004-2008 kurtosis
estimate, for example, is 15.72 while its bootstrap mean is 8.32 because
resampling sometimes omits influential extreme episodes.

## Between-period evidence

The table reports centered-bootstrap p-values for a common statistic across
the six periods. It is an exploratory heterogeneity test, not a test of
practical equivalence. Its assumptions include approximate stationarity
within each five-year period and weak enough dependence for block resampling.

| Statistic | 63-day blocks | 126-day blocks | 252-day blocks | Reading |
| --- | ---: | ---: | ---: | --- |
| Excess kurtosis | 0.001 | 0.002 | 0.001 | Strong evidence against one constant magnitude. |
| Leverage, lags 1-5 | 0.058 | 0.028 | 0.016 | Negative in every block; magnitude likely varies. |
| Lag-1 return ACF | 0.005 | 0.013 | 0.008 | The near-zero approximation is period dependent. |
| Roughness H | 0.305 | 0.213 | 0.044 | Apparent rise; conclusion depends on block length. |
| 10-day Zumbach proxy | 0.317 | 0.080 | 0.037 | Sign and magnitude are uncertain; block sensitive. |
| Taylor gap, lags 1-20 | 0.134 | 0.062 | 0.006 | Block sensitive; five-year estimates do not support a fixed gap. |
| Log-volatility clustering, lags 1-20 | 0.315 | 0.377 | 0.240 | Positive throughout; equal magnitude is not established. |

There are nine related statistics and multiple block lengths in the full
output, so the p-values are descriptive rather than multiplicity-adjusted
discoveries. Failure to reject a common value does not show that periods are
equivalent within a calibration tolerance.

## Source-data and volatility-proxy checks

The recorded SPX open equals the previous close on 96% of days in 1999-2003,
51% in 2004-2008 and less than 1% in 2014-2023. The transition occurs mainly
around 2006-2008; `SP500/source_years.csv` shows every year. The squared
overnight gap contributes about 0.1% of GK-plus-gap variance in 1999-2003,
12% in 2014-2018 and 33% in 2019-2023. This pattern suggests a change in
recorded-open conventions or construction, but source provenance is needed to
establish why. It makes open-dependent volatility estimates less comparable
over the full 30 years. No SPX Garman-Klass value was floored in these blocks.

A high/low-only Parkinson variance gives H estimates within 0.006 of GK in
each block. Its leverage estimates differ by at most 0.013 and its 10-day
Zumbach values by at most 0.047. Adding the recorded overnight-gap square
changes 2019-2023 leverage from -0.187 to -0.155 and H from 0.128 to 0.133;
earlier leverage estimates hardly move. These sensitivities do not validate
the early index bars: both range-based proxies use the same high/low data.
They also do not reproduce an intraday realised-variance estimate.

## Implication for joint calibration

The negative leverage *sign* and positive volatility clustering are plausible
shared physical constraints. Their exact magnitudes need uncertainty bands;
leverage already shows period heterogeneity. A single SPX excess-kurtosis
target over 1994-2023 is poorly supported. The current roughness and Zumbach
proxy results are insufficient for narrow, fixed targets, especially given
measurement risk and block-length sensitivity. The full statistic vector is
correlated, so any physical calibration loss should use a regularised joint
covariance of its selected components, rather than nine independent weights.
`bootstrap_covariance.csv` records within-period covariance for this purpose;
it does not itself justify pooling across periods.

The historical index series ends in 2023 while the available SPX/VIX option
quotes are a single 2026-02-19 snapshot. Historical physical statistics can
constrain shared model parameters but do not identify the latent reservoir
state on that option date. A common-target calibration also requires explicit
practical tolerances, additional source-data checks and, for definitive
roughness and Zumbach targets, a consistent intraday realised-variance series.

Method references: [Garman and Klass (1980)](https://www-2.rotman.utoronto.ca/~kan/3032/pdf/FinancialAssetReturns/Garman_Klass_JB_1980.pdf),
[Parkinson (1980)](https://www.cmegroup.com/trading/fx/files/michael_parkinson.pdf),
and the [model specification](../../documentation/esn_exponential_quadratic_model.pdf).
