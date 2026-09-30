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

## Rolling estimates: annual steps, five- and ten-year windows

The [interactive rolling plot](SP500/rolling_estimates.html) shows all nine
statistics and their 95% pointwise intervals, together with the share of days
whose recorded open equals the previous close. The complete numbers are in
[rolling_estimates.csv](SP500/rolling_estimates.csv); corresponding
[proxy sensitivities](SP500/rolling_proxy_sensitivity.csv) and
[source checks](SP500/rolling_source_checks.csv) permit the same measurement
checks by window. The five-year series has 26 windows ending in 1998-2023;
the ten-year series has 21 ending in 2003-2023. Each moves by one calendar
year. The final windows end on 2023-11-14. Within each window, 500
moving-block bootstrap replicates of length 126 trading days give the 95%
percentile intervals. The paired returns and OHLC ratios and the exclusion of
artificial block joins are the same as in the disjoint-period analysis above.

### Five-year point estimates

The recorded-open column is a source-data diagnostic, not a stylised fact.
Values here are rounded; use the CSV for intervals and unrounded values.

| Window | H from log GK volatility | Leverage, lags 1-5 | 10-day Zumbach proxy | Excess kurtosis | Open = previous close |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1994-1998 | 0.045 | -0.141 | -0.050 | 8.49 | 82.6% |
| 1995-1999 | 0.044 | -0.138 | -0.056 | 5.60 | 93.0% |
| 1996-2000 | 0.052 | -0.126 | 0.015 | 3.57 | 97.1% |
| 1997-2001 | 0.063 | -0.128 | 0.064 | 2.66 | 96.6% |
| 1998-2002 | 0.080 | -0.132 | 0.126 | 1.66 | 96.1% |
| 1999-2003 | 0.067 | -0.110 | 0.187 | 1.22 | 96.1% |
| 2000-2004 | 0.066 | -0.111 | 0.181 | 1.79 | 96.3% |
| 2001-2005 | 0.062 | -0.119 | 0.159 | 2.43 | 96.3% |
| 2002-2006 | 0.045 | -0.115 | 0.161 | 3.11 | 85.9% |
| 2003-2007 | 0.046 | -0.082 | 0.168 | 1.70 | 70.0% |
| 2004-2008 | 0.077 | -0.171 | 0.354 | 15.72 | 50.6% |
| 2005-2009 | 0.088 | -0.156 | 0.327 | 10.14 | 31.3% |
| 2006-2010 | 0.106 | -0.153 | 0.322 | 8.62 | 12.3% |
| 2007-2011 | 0.114 | -0.151 | 0.330 | 6.56 | 5.7% |
| 2008-2012 | 0.107 | -0.152 | 0.326 | 7.01 | 7.6% |
| 2009-2013 | 0.080 | -0.103 | 0.125 | 4.16 | 12.6% |
| 2010-2014 | 0.078 | -0.112 | 0.107 | 4.67 | 12.4% |
| 2011-2015 | 0.078 | -0.153 | 0.223 | 4.93 | 12.3% |
| 2012-2016 | 0.075 | -0.154 | 0.099 | 2.01 | 10.4% |
| 2013-2017 | 0.069 | -0.159 | 0.092 | 2.91 | 5.3% |
| 2014-2018 | 0.097 | -0.203 | 0.129 | 3.79 | 0.2% |
| 2015-2019 | 0.104 | -0.191 | 0.180 | 3.85 | 0.0% |
| 2016-2020 | 0.143 | -0.203 | 0.284 | 22.10 | 0.0% |
| 2017-2021 | 0.140 | -0.205 | 0.282 | 22.01 | 0.1% |
| 2018-2022 | 0.137 | -0.198 | 0.278 | 12.95 | 0.1% |
| 2019-2023 | 0.128 | -0.187 | 0.278 | 13.99 | 0.1% |

### What the rolling curves show

Leverage remains negative and log-volatility clustering remains positive in
every five- and ten-year window, supporting their **signs** as robust
constraints. Their magnitudes still move. For example, five-year leverage
ranges from -0.082 to -0.205, while the corresponding clustering mean ranges
from 0.224 to 0.625. The recent 2019-2023 leverage estimate is -0.187
[-0.234, -0.105] and clustering is 0.436 [0.274, 0.538]. These are pointwise
intervals for those particular windows, not bounds on a common level.

Excess kurtosis has large episode effects. The 2003-2007 five-year estimate
is 1.70, the 2004-2008 estimate is 15.72, and the 2015-2019 and 2016-2020
estimates are 3.85 and 22.10. Wider windows smooth some variation but do not
remove it: ten-year kurtosis changes from 11.13 in 2008-2017 to 5.02 in
2009-2018, and from 4.61 in 2010-2019 to 17.97 in 2011-2020. The adjacent
windows both add and remove a year; these shifts show the influence of
extreme episodes on a rolling target, not a smooth secular trend. The latest
five-year kurtosis interval is particularly wide: 13.99 [1.13, 20.97].

The roughness estimate H generally rises from early to recent windows:
0.045 in 1994-1998 versus 0.128 in 2019-2023, or 0.054 in 1994-2003
versus 0.115 in 2014-2023. Measurement comparability and the block-length
sensitivity noted above limit a fixed historical roughness target. The
ten-year Zumbach proxy falls from 0.301 in 2008-2017 to 0.120 in 2009-2018
and rises from 0.131 in 2010-2019 to 0.276 in 2011-2020. Its recent
five-year interval is 0.278 [-0.037, 0.323]; even a positive point estimate
does not give a precise five-year target. The Taylor gap and return ACF also
vary in the plot and CSV; the latter is near zero in some windows and notably
negative in windows containing 2020.

The recorded-open diagnostic changes sharply across the rolling history:
about 96% of observations in 1999-2003 have an open equal to the previous
close, compared with about 0.1% in 2019-2023. This reinforces the need to
check the [proxy sensitivity results](SP500/rolling_proxy_sensitivity.csv)
before interpreting open-dependent H, leverage or Zumbach differences as
economic change. The high/low-only proxy is a sensitivity check; it does not
resolve the source-convention question or turn daily OHLC into intraday
realised variance.

Annual overlap is useful for seeing *when* estimates change, but neighboring
five-year windows share roughly four years of data and neighboring ten-year
windows share roughly nine. The 47 windows are therefore **not 47 independent
replications**. Their intervals are pointwise, not a simultaneous band for a
flat curve; visual overlap or separation of two bars is not a valid test of
equality. The disjoint-period heterogeneity analysis above remains the
primary statistical check for constancy. For calibration, the rolling curves
help identify plausible recent ranges and influential episodes, while any
pooled target still needs an explicit tolerance and a covariance estimate
from a defensible time span.
