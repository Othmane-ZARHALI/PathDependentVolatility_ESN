# Fast-factor downside-asymmetry screen

## Interpretation boundary

This is a paired fresh-seed screen of shortlisted Stage-3 M2 configurations. The original seven stylised-fact gates are unchanged. The tail-event metrics measure direction, but are not empirical SPX acceptance targets.

## Results by shift

| spike shift | joint pass | H | leverage | Zumbach | kurtosis | negative day 0 | positive day 0 | day-0 gap | days 1-5 gap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.00 | 1 | 0.096416 | -0.017045 | 0.044228 | 9.8897 | 0.83336 | 0.48338 | 0.34998 | 0.29441 |
| 0.10 | 0.8375 | 0.098479 | -0.025098 | 0.026548 | 10.616 | 0.88156 | 0.49453 | 0.38703 | 0.37302 |
| 0.20 | 0.3375 | 0.10132 | -0.032211 | 0.0086878 | 11.469 | 0.92969 | 0.51117 | 0.41852 | 0.44412 |
| 0.30 | 0 | 0.10467 | -0.0384 | -0.0085307 | 12.431 | 0.9767 | 0.53495 | 0.44175 | 0.50819 |
| 0.40 | 0 | 0.10832 | -0.043707 | -0.024544 | 13.485 | 1.0223 | 0.55564 | 0.46661 | 0.56623 |
| 0.50 | 0 | 0.11212 | -0.0482 | -0.039086 | 14.616 | 1.0692 | 0.57971 | 0.48952 | 0.62053 |

Event responses are relative changes in volatility versus the preceding ten-day mean. Joint pass is the fraction satisfying all original Stage-3 proxy gates.

## Screen reading

- Highest joint pass rate: spike shift `0.00` with `100.0%`.
- Best non-zero shift on the old gates: `0.10` with `83.8%` joint acceptance and `6` of `16` configurations passing on every seed.
- At that non-zero shift, the only individual gate with failures was: `zumbach_rank`.
- Largest day-zero downside-minus-upside gap: spike shift `0.50` with `49.0%`.
- Positive-tail day-zero response fell for 0 of 5 tested non-zero shifts.
- Select a shift only after considering joint acceptance, leverage, Zumbach, return tails, cap use, and both asymmetry horizons together.

## Next evidential step

Estimate the same event-response statistics and uncertainty from matched SPX and intraday realised-variance data. Then run an untouched joint search around the best non-zero shift rather than treating this shortlist screen as final calibration.

Package version: `0.4.0`.
