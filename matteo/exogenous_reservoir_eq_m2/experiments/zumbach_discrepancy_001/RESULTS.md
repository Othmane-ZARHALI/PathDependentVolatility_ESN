# Pearson/rank Strong-Zumbach discrepancy

## Conclusion

The discrepancy is genuine and interpretable. It is not caused by a coding error, the use of
log returns, or the fine time grid.

Pearson correlation detects a pronounced asymmetry in the *magnitude* of extreme episodes.
Rank correlation asks whether the broad ordering of observations is more predictive forward
than backward. For the Stage-4 primary, that bulk ordering is almost time-symmetric even though
the extreme forward response is much larger. The two statistics therefore answer different
questions.

## Directional decomposition

The table averages 5-, 10-, and 20-day windows, paths, and the five untouched validation seeds.

| model | estimator | forward leg | backward leg | difference |
| --- | --- | ---: | ---: | ---: |
| primary | Pearson | 0.2705 | 0.1369 | 0.1336 |
| primary | Spearman rank | 0.1708 | 0.1730 | -0.0022 |
| primary | Gaussian rank | 0.1931 | 0.1740 | 0.0191 |
| primary | Kendall | 0.1153 | 0.1162 | -0.0010 |
| anchor | Pearson | 0.3030 | 0.1398 | 0.1632 |
| anchor | Spearman rank | 0.2193 | 0.1983 | 0.0210 |
| anchor | Gaussian rank | 0.2385 | 0.1961 | 0.0424 |
| anchor | Kendall | 0.1485 | 0.1335 | 0.0151 |

For the primary, the Pearson forward correlation is nearly twice the backward correlation.
After replacing observations by uniform ranks, however, the two legs are essentially equal.
Normal-score ranks, which give more weight to the tails than uniform ranks, recover a small
positive difference, although its 95% interval `[-0.0192, 0.0575]` still includes zero.

## Tail sensitivity

For the primary, the fraction of positive covariance associated with an observation lying in
the upper 5% of either variable is `65.8%` in the forward leg and `51.5%` in the backward leg.
At the upper 1%, the corresponding shares are `42.2%` and `28.1%`. Extreme observations are
important to both correlations, but disproportionately important in the causal forward leg.

| treatment | primary Zumbach | 95% interval |
| --- | ---: | ---: |
| ordinary Pearson | 0.1336 | [0.1078, 0.1595] |
| winsorise above 99% | 0.1239 | [0.1001, 0.1477] |
| winsorise above 97.5% | 0.0965 | [0.0668, 0.1262] |
| winsorise above 95% | 0.0663 | [0.0270, 0.1056] |
| remove upper 1% | 0.0694 | [0.0379, 0.1010] |
| remove upper 2.5% | 0.0267 | [-0.0111, 0.0646] |
| remove upper 5% | -0.0083 | [-0.0468, 0.0302] |

Winsorisation preserves the information that a tail event occurred while limiting its
magnitude; the effect remains positive. Removing the observations eliminates the asymmetry.
Thus both tail-event incidence and tail magnitude carry the Pearson signal.

The conditional decile profiles tell the same story. In the top conditioning decile, the
primary's mean standardised response is `0.693` in the forward direction versus `0.316` in the
backward direction. Across the lower nine deciles the two profiles are much closer.

## Horizon decomposition

| window | Pearson Z | Spearman Z | Gaussian-rank Z |
| ---: | ---: | ---: | ---: |
| 5 days | 0.1606 | 0.0012 | 0.0218 |
| 10 days | 0.1513 | 0.0100 | 0.0332 |
| 20 days | 0.0890 | -0.0178 | 0.0023 |

Pearson asymmetry is positive at every horizon and becomes smaller as the window length grows.
The rank statistic is near zero at every horizon rather than being cancelled only by averaging
the three windows.

## Return convention

Using the matched Ito relative return changes the primary Pearson result from `0.1336` to
`0.1280` and the Spearman result from `-0.0022` to `-0.0075`. The discrepancy therefore does
not come from the second-order correction between `d log S` and `dS/S`.

## Implication for the experimental protocol

The model document already defines ordinary Pearson correlation as the primary Zumbach statistic
and describes rank correlation as a distinct robustness diagnostic. The Stage-2
acceptance code nevertheless promotes `zumbach_rank` to a mandatory gate. These two conventions
are internally inconsistent.

This study does not retroactively change the Stage-4 decision: under its declared rank gate,
the primary failed. For future experiments, however, the scientifically coherent protocol is:

1. use Pearson Zumbach as the primary stylised-fact gate because it is the statistic defined in
   the document and it captures the extreme-event magnitude asymmetry of interest;
2. retain Spearman, Gaussian-rank, and tail-trimmed results as robustness and mechanism
   diagnostics rather than requiring them to have the same magnitude or sign;
3. predeclare a tail-robust Pearson check, such as 99% winsorisation, to ensure that the result
   is not determined by a single observation;
4. report each horizon separately as well as their average.

Under such a protocol, both the primary and the anchor exhibit a robust positive Strong
Zumbach effect on the 64-step grid. The Stage-4 primary still has a weaker effect than the
anchor, but it combines that positive Pearson asymmetry with materially stronger leverage and
better-centred roughness.
