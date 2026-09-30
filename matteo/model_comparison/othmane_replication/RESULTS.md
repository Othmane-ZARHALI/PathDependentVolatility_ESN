# M2 and ESN A2: first historical diagnostics

Run date: 2026-09-30. These are reproducible exploratory results, not a
publishable model-family ranking. [PLAN.md](PLAN.md) gives the remaining gates.

## Interactive Plotly figures

- [Nine-asset overview](plots/overview.html): total and return-only errors, H,
  and excess kurtosis.
- [Error contribution heatmap](plots/block_differences.html): which fact moves
  the full score for each asset.
- [Early readout screen](plots/readout_screen.html): pre-2015 selection and
  later-period return-only results.
- [Roughness lag sensitivity](plots/roughness_lag.html): common lag-band H
  estimates and their effect on the full score.
- Curves and 5th–95th simulated-path bands:
  [LHX](plots/curves_lhx.html), [BK](plots/curves_bk.html),
  [HUM](plots/curves_hum.html), [TSN](plots/curves_tsn.html),
  [IEX](plots/curves_iex.html), [KMB](plots/curves_kmb.html),
  [SP500](plots/curves_sp500.html),
  [RUSSELL2000](plots/curves_russell2000.html), and
  [NASDAQ100](plots/curves_nasdaq100.html).

The HTML figures use the bundled `plots/plotly.min.js` and work offline.

## Full-sample facts on fresh simulation seeds

The saved six-parameter ESN A2 fit was calibrated to the full historical
targets. M2 is one synthetic Stage-3 configuration, candidate 18/scenario 0;
it has not been fitted to these assets. Both models were simulated for 20
paths of 4,000 retained trading days after 504 burn-in days, starting at seed
720001. The score is Othmane's 59-term weighted squared error of **averaged
path statistics**. The return-only columns contain the two six-lag ACFs and
excess kurtosis; those facts avoid the direct historical-versus-model variance
proxy mismatch.

| Asset | Full ESN | Full M2 anchor | Return-only ESN | Return-only M2 anchor |
| --- | ---: | ---: | ---: | ---: |
| LHX | 2.81 | 30.39 | 1.37 | **0.70** |
| BK | 3.00 | 17.18 | **1.16** | 10.93 |
| HUM | 5.31 | 38.62 | 3.94 | **2.29** |
| TSN | 3.05 | 50.14 | 1.60 | **1.26** |
| IEX | 3.47 | 57.51 | **1.61** | 2.29 |
| KMB | 3.41 | 25.84 | 0.91 | **0.47** |
| SP500 | 2.96 | 14.23 | **1.13** | 8.91 |
| RUSSELL2000 | 9.37 | 17.00 | **6.24** | 11.67 |
| NASDAQ100 | 1.19 | 10.01 | **1.00** | 9.03 |

The full score favors the archived ESN on all nine assets. This is mainly a
roughness gap for stocks, amplified by the objective's 3× H weight: TSN's H
block contributes 0.51 for ESN and 46.63 for M2. The single M2 anchor gives
H≈0.075 on every asset, while the historical targets range from 0.035 to
0.086. A fitted per-asset M2 architecture could behave differently.

M2's leverage block is smaller on all nine assets under the *specified*
`(lag+1)^-2` weighting. For SP500, the target at lag 1 is −0.173, versus
−0.368 for ESN and −0.205 for M2. M2's index return persistence is too weak:
SP500 absolute-return ACF at lag 20 is 0.234 in the data, 0.162 in ESN,
and 0.054 in M2. This creates most of the return-only index gap.

The scoring convention matters. On SP500, squared error of the mean statistics
is 2.96 for ESN and 14.23 for M2; mean per-path squared error is 7.88 and
19.79. The bootstrap intervals in the overview resample the 20 simulated
paths, recompute the mean statistics, and then score them. They do **not**
express uncertainty in the historical data.

## Early-selection M2 screen

We screened all 32 existing M2 Stage-3 readouts with structural scenario 0,
using each asset's return ACFs and kurtosis through 2014-12-31. The screen
used 8 paths × 2,000 retained days and seed 850001. Each asset's lowest early
score selected candidate 3 or 12. We then evaluated the chosen candidate on
20 fresh paths × 4,000 days with seed 720001 against the post-2014 historical
target. The anchor was evaluated with the same budget. This is a coarse
readout choice; the two factor correlations and reservoir architecture were
not fitted to historical assets.

| Asset | Selected M2 | Anchor M2 | Archived ESN |
| --- | ---: | ---: | ---: |
| LHX | **0.93** | 1.04 | 1.39 |
| BK | **1.59** | 1.70 | 2.70 |
| HUM | **1.52** | 2.85 | 4.83 |
| TSN | **5.52** | 7.91 | 7.25 |
| IEX | **1.76** | 2.11 | 2.16 |
| KMB | **0.79** | 1.78 | 2.86 |
| SP500 | 13.43 | 15.16 | **3.89** |
| RUSSELL2000 | 8.98 | 9.57 | **7.97** |
| NASDAQ100 | 6.34 | 7.24 | **0.82** |

The early-selected M2 readout improves on the anchor for every late target.
It has a lower point estimate than archived ESN on all six stocks, especially
HUM and KMB, and a higher one on all three indices. The simulation-only 90%
bootstrap intervals overlap for several stock comparisons: for LHX they are
0.75–1.49 (selected M2) and 0.91–2.16 (ESN). HUM is 0.67–2.94 versus
3.04–6.62; KMB is 0.14–1.84 versus 1.91–3.98. See
[`late_return_scores.csv`](plots/late_return_scores.csv) for every interval.

The archived ESN fit used the **whole** historical sample, including these
late years. Its column therefore cannot serve as a clean out-of-sample
benchmark. The stocks run through 2025-12-31; the index workbook ends on
2023-11-14. The first return after each date cut is omitted to avoid a return
crossing the split.

## Source and estimator checks

The six saved historical target files exactly match an in-memory rebuild from
raw OHLC data; see [AUDIT.md](AUDIT.md) and
[`target_audit_001.json`](target_audit_001.json). The archived calibrator's
entry point cannot run as checked in. The
[`joint_calibration_reviewed.py`](joint_calibration_reviewed.py) copy adds the
sixth starting coordinate `alpha=0.5`, gives the invalid-region residual its
correct length of 59, passes fitted alpha to Stage 2, includes it in final
output, and writes to separate reviewed filenames. A short controlled run
confirmed seven finite evaluations from a valid six-coordinate start. The
invalid-region branch yielded the expected 59×25/2 = 737.5 cost.

The alpha omission is material on SP500. Holding the archived H, rough scale,
and decay rates fixed, a 12-path/2,500-day Stage-2 rerun gives `m1=0.196`
when `alpha=0.2204` is passed, but `m1=0.059` when the default `alpha=1`
is used. The latter m1, evaluated back at the fitted alpha, produces excess
kurtosis 5.13 instead of the target 10.31. This isolates the bug; it does not
identify how the archived final m1 was produced. The full nine-asset
multi-start recalibration has not been run.

The historical H target uses lags 2–40; the model estimator in Othmane's
objective uses 1–40. Recomputing both sides on common 2–40 lags changes
SP500 full scores from 2.96/14.23 (ESN/M2, original mixed bands) to
3.14/22.86. ESN remains lower on all nine assets under either common lag
band. The [lag audit](roughness_lag_sensitivity_001.json) records the values.

## Interpretation and next work

The return-only screen suggests that M2 already captures some stock return
persistence and tails at least as well as this archived ESN fit, while the
current M2 configurations underproduce persistence for broad indices. Its
leverage curve is closer at the first few lags, but the full objective gives
leverage much less numerical influence than H and the ACFs. None of these
results establishes a family-level advantage: M2 had only a 32-readout
screen, ESN used full-sample historical fitting, and variance observations
still differ.

The next decisive work is to construct the same simulated OHLC observation
for both models, fit both families on identical early years with several
starts, and evaluate later years with data-block uncertainty. Conditional
SPX/VIX pricing and the synthetic mechanism protocols from the plan also
remain open.

## Reproduce

From the repository root, install NumPy, SciPy, pandas, openpyxl and Plotly
into a Python 3.12 environment, then run:

```bash
python matteo/model_comparison/othmane_replication/audit_targets.py --output matteo/model_comparison/othmane_replication/target_audit_001.json
python matteo/model_comparison/othmane_replication/run_nine_assets.py
python matteo/model_comparison/othmane_replication/screen_m2.py
python matteo/model_comparison/othmane_replication/roughness_lag_audit.py
python matteo/model_comparison/othmane_replication/stage2_alpha_audit.py
python matteo/model_comparison/othmane_replication/make_plots.py
```

Run manifests, per-path statistics and source hashes are in
[`nine_assets_fresh_001.json`](nine_assets_fresh_001.json) and
[`m2_readout_screen_001.json`](m2_readout_screen_001.json). Scores are also
in [`nine_asset_scores.csv`](plots/nine_asset_scores.csv) and
[`late_return_scores.csv`](plots/late_return_scores.csv).
