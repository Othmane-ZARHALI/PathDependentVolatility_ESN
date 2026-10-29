# Independent amplitude after the original cap

The original capped normalized variance component is multiplied by an independent
mean-one OU exponential. The original curve, floor, reservoir and cap stay fixed.
Gamma is log-variance sensitivity; it is not a bounded multiplier range.

`vix-fresh-log` and `vix-fresh-linear` show new joint numerical draws, on both
strike axes. `vix-wing-shape` compares slope and curvature over K/marketF>2.
`spx-all-maturities-*` shows actual stock repricing. `vix-153-grid` retains all
nine training variants at the longest VIX expiry. `vix-forwards` checks forwards.

Lines connect raw price-derived Black call IVs using market forwards, without
smoothing/clipping. Missing IVs are marked at the foot of panels. Pointwise bands
are sampling SE only, excluding deterministic, pilot and selection uncertainty.
New numerical draws use the same market data; there is no market OOS claim.
See POSTCAP_VOV_RESULTS.md and the run's audit for results and limitations.
