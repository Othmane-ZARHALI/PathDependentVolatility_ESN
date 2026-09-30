# Othmane ESN A2 versus M2: empirical replication plan

## Assumptions and scope

The primary empirical study is Othmane's nine-asset historical calibration in
`SimultanuousAllFeatures/`. The six stocks are LHX, BK, HUM, TSN, IEX and KMB;
the three indices are SP500, RUSSELL2000 and NASDAQ100. The first comparison
uses the same five facts, horizons, smoothing and target files for both models.
The `CrossSectionalData/` SPX/VIX snapshot is a separate, second study. The
earlier `Roughness/`, `Taylor_Effect/`, `Leverage_Effect/`, and `Zumbach_Effect/`
folders use synthetic data-generating processes and are mechanism checks,
not historical fits.

The checked-in M2 Stage-3 results are synthetic capacity tests. They are not
historical calibrations. The exploratory runner in this folder evaluates an
existing Stage-3 M2 configuration against the historical targets; it does not
call that configuration fitted to SPX or to the other eight assets.

## What Othmane calibrated

The historical joint objective has six blocks: moment-based log-volatility
roughness (1 term), return-to-variance leverage at lags 1–40 (40), absolute
return ACF (6), squared return ACF (6), Pearson Zumbach at windows
2/10/20/30/40 (5), and excess return kurtosis (1). Lag curves are smoothed
over three adjacent *selected* lag entries. Each block is divided by a
cross-asset scale and the square root of its term count; roughness receives
an additional factor of 3. Leverage has `(L+1)^-2` weights. The model fits
six parameters per asset: `H`, `rough_scale`, `lam_lo`, `lam_hi`, `m1`, and
EWMA `alpha`; the reservoir architecture and other readout parameters are
fixed. Simulation uses fixed seeds within each least-squares run. A later
one-dimensional `m1` step targets kurtosis. Final fitted values and
simulation bands are in `joint_calibration_final.json` and `joint_ci.json`.

The real variance target is same-day Garman–Klass OHLC variance, paired with
close-to-close returns. Othmane's simulation evaluates roughness, leverage
and Zumbach on its *raw* variance, while its returns use EWMA variance.
This observation choice must be reported beside any fitted score.

## Bug and reproducibility audit, before interpreting fits

1. Recompute all six processed historical targets from the raw files and
   compare their arrays and source hashes. Preserve sample dates and the
   exact close-to-close/Garman–Klass alignment. Do this without overwriting
   the ignored `data/processed/` files.
2. Repair the historical calibrator in a reviewed copy, then verify a
   controlled rerun. Its current `starting_point()` returns five values but
   `calibrate_joint()` unpacks six. The stated/penalty residual length is 53,
   while its six blocks contain **59** terms. The checked-in final JSON is
   therefore not reproducible by invoking the current `__main__` unchanged.
   The fixed source must record a sixth starting value, correct the penalty
   length, and freeze the full parameter bounds and seed policy. Its Stage-2
   call also omits `alpha=c["alpha_fit"]`, so the `m1` bisection would default
   to `alpha=1` while a final fit retains a lower fitted `alpha`.
   Rerun Stage 2 with the fitted alpha and quantify the change. Its final
   serializer omits `alpha_fit` even though that field exists in the saved
   final JSON; restore it in the reviewed copy.
3. Check the calibrated objective against the report's equations, including
   lag indices, curve smoothing, roughness estimator, Zumbach direction,
   and scalar versus pathwise averaging. Test invalid-region handling and
   finite residuals before long optimizations. Historical roughness starts at
   lag 2 because of proxy noise, while model roughness starts at lag 1;
   report sensitivity to a common lag band.
4. Run Othmane's fitted nine-asset values on fresh held-out simulation seeds
   with a production budget. Compare exact termwise values and coverage with
   `joint_ci.json`; do not treat the archived bands as data uncertainty.
   Repeat fitting from multiple starts; record boundary hits and all losses.
5. Give M2 and ESN A2 the same historical estimator and observed targets.
   First validate the measurement layer: M2 currently outputs daily average
   instantaneous variance, whereas the data use Garman–Klass OHLC. Add
   simulated OHLC and apply the same proxy to both models, or label direct
   variance scores exploratory. Verify temporal alignment and estimator
   sensitivity with synthetic paths.
6. Estimate uncertainty from time blocks of the actual historical series,
   with at least 126/252/504-day block sensitivity. Fit a compact covariance
   or predeclared lag subset; the 40 leverage lags are correlated. Keep
   training and held-out years separate. Use identical data splits and
   reported simulation budgets for the two model families.

## Replication sequence

1. **Historical facts:** run the exploratory comparison CLI in this folder
   to check loading, conventions and a paired score table. Then fit M2
   readout parameters and the two equity/factor correlations on training
   windows, reusing exogenous path caches for cheap readout changes.
   Search several architectural grids because M2's synthetic Stage-3 grid
   was never optimized against these assets. Evaluate both models on untouched
   years and independent simulation seeds. Report every statistic, its curve,
   uncertainty, weighted objective contribution and parameter boundary.
2. **Synthetic mechanism protocols:** repeat the four DGP experiments with
   their original seeds, horizons and reference statistics. Compare M2, ESN A2
   and the DGP under one estimator and calibration budget; separately retain
   Othmane's original scoring convention for exact historical replication.
3. **SPX/VIX cross-section:** audit the 2026-02-19 quote preprocessing and
   reproduce Othmane's option plots and weighted price/IV errors. Implement
   conditional M2 SPX and VIX option valuation and an out-of-sample Q
   normalizer before attempting its analogous fit. The existing M2 package
   does not yet have those pricers. Freeze maturity/strike holdouts before
   fitting and report price errors against bid/ask intervals, rather than
   comparing unlike weighting schemes.

## Deliverables and gate

Keep code and run manifests in this folder. A publishable comparison requires
the historical source audit, matched variance observation, held-out historical
results, multi-seed simulation uncertainty and the separate option-pricing
validation. The current CLI result is an exploratory baseline only.
