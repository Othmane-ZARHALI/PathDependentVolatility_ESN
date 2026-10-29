# Point1 acceleration pilot with 4096 paths

Joint eleven-parameter calibration; full-price-guarded local search on the calibration simulations. The original OTM
observations are converted to calls by parity. Blue uses the new selected
parameters; orange uses the previous fitted vector on the SAME fixed training samples.

- [All seven SPX maturities](spx-all-maturity-fit.png)
- [All four VIX maturities](vix-all-maturity-fit.png)
- [Actual VIX forwards and path allocation](forwards-and-path-budgets.png)
- [SPX signed price errors](spx-price-errors.png)
- [VIX signed price errors](vix-price-errors.png)
- [Sampling precision by path budget](sampling-precision.png)

All figures also have SVG exports. Bands show ±2 pointwise sampling SE, not
parameter uncertainty or a simultaneous confidence band. Invalid IVs are marked;
signed prices remain in the tables and loss. Both market and model forward
conventions are preserved in the run; a common call plot does not establish
that the model VIX forward is correct. Original-OTM fit errors are also retained.
