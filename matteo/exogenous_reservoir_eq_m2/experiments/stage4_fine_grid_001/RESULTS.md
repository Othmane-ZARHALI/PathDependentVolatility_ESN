# Stage-4 fine-grid recalibration results

## Decision

The targeted recalibration did not produce a candidate that can yet be advanced to the
64-versus-128 confirmation. The frozen primary passed all four search seeds but only two of five
untouched validation seeds. Under the unchanged seven-gate protocol, this is a failed
validation, not a successful recalibration.

The failure is specific: return autocorrelation, volatility persistence, roughness, leverage,
tails, and the Taylor effect passed on every validation seed. The rank-based Zumbach gate passed
on only two seeds.

## Frozen primary

The search-only rule selected candidate 57 under structural scenario 2:

| quantity | value |
| --- | ---: |
| cluster loading | 0.184230 |
| feedback curvature | 0.084458 |
| feedback shift | 3.096703 |
| spike curvature | 0.207564 |
| spike shift | 0.230388 |
| cap level | 5.565315 |
| cap sharpness | 9.854119 |
| feedback/equity correlation | 0.908481 |
| spike/equity correlation | 0.964497 |

It passed 4/4 search seeds. Its untouched validation pass rate was 2/5.

## Untouched validation

Values are means across the five validation seeds. The intervals are paired-seed or ordinary
95% Student intervals, as appropriate; five seeds still imply material uncertainty.

| diagnostic | primary | preceding anchor | reading |
| --- | ---: | ---: | --- |
| max absolute return ACF | 0.0690 | 0.0697 | both pass |
| mean log-variance ACF | 0.2295 | 0.2729 | both pass |
| rough H | 0.1174 | 0.1437 | primary is nearer the centre of the proxy band |
| rank leverage | -0.0236 | -0.0107 | primary is materially stronger |
| Pearson Zumbach | 0.1336 | 0.1632 | positive for every seed in both cases |
| rank Zumbach | -0.0022 | 0.0210 | primary fails the declared gate on 3/5 seeds |
| Ito-return Pearson Zumbach | 0.1280 | 0.1593 | close to the log-return statistic |
| Ito-return rank Zumbach | -0.0075 | 0.0169 | same instability as the log-return rank statistic |
| excess kurtosis | 7.3785 | 7.4490 | both pass |
| Taylor gap | 0.0264 | 0.0305 | both pass |

For the primary, the 95% interval for rank leverage is `[-0.0292, -0.0180]`, while the interval
for rank Zumbach is `[-0.0444, 0.0401]`. In contrast, Pearson Zumbach has interval
`[0.1078, 0.1595]`, and its Ito-return version has interval `[0.1022, 0.1538]`.

Relative to the anchor on matched validation seeds, the primary changes rank leverage by
`-0.0129` with paired interval `[-0.0148, -0.0110]`, but changes rank Zumbach by `-0.0232` with
paired interval `[-0.0383, -0.0081]`. The search therefore found a real leverage improvement,
but partly paid for it by weakening rank-based time-reversal asymmetry.

## Numerical stability

No volatility explosion or Q-normalisation failure occurred. The primary's mean P-score cap
exceedance fraction was `0.0182%`, with a maximum of `0.0394%` across validation seeds. The
failure is statistical/mechanistic rather than a fine-grid stability failure.

## Interpretation

The standard Pearson Strong Zumbach effect remains large, positive, and insensitive to whether
the return is represented as a log increment or as the matched Ito integral `dS/S`. The
second-order Ito correction is therefore not driving the conclusion at daily frequency.

The rank-based version behaves very differently and is highly seed-sensitive. Because the rank
gate was declared before the experiment, it cannot be discarded after seeing this result. The
current evidence instead says that M2 can robustly combine the other six proxy facts on the fine
grid, but this targeted box has not demonstrated robust rank-based Strong Zumbach asymmetry.

## Continuation rule

Do not run the 128-step confirmation for candidate 57. First decide, independently of these
results, whether the primary scientific object is the conventional Pearson Strong Zumbach
statistic or the rank-robustified statistic. If the rank version remains a required gate, the
next search should explicitly protect it while restoring leverage, and should use longer search
paths because the four-year search ranking did not generalise to the six-year validation paths.
