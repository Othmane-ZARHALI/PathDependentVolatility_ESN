# Cap, weighted tails and continuous response experiments

Selected baseline remains the maturity-curve fit; none of these experiments is an adopted joint calibration. All original OTM inputs, fixed OU architecture, seven SPX/four VIX dates, four forwards and4096 paths per strike are preserved. No quintic/jump proposal.

Read [the numerical report](../../CAP_TAIL_RESULTS.md), [mathematics](../../CAP_TAIL_MATHEMATICS.md), and [stochastic vol-of-vol mechanisms and literature](../../STOCHASTIC_VOL_OF_VOL_RESEARCH.md).

- [Cap location and sharpness: all four VIX dates](vix-cap-location.png)
- [Alternative cap shapes: all four VIX dates](vix-cap-shapes.png)
- [Likelihood-weighted tail calls and pointwise SE](weighted-tail-calls.png)
- [Actual VIX forwards and pointwise SE](vix-forwards.png)
- [Tempering and candidate-specific normalization](cap-and-normalization.png)
- [Signed VIX price errors](vix-signed-price-errors.png)
- [Two VIX-focused fits](vix-focused-fits.png)
- [All seven SPX maturity tradeoffs](spx-all-maturities-iv.png)
- [Signed SPX price errors](spx-all-maturities-price-errors.png)
- [Separate previously used VIX checking streams](vix-separate-checking.png)
- [Two existing-cluster response-mixture profiles](vix-cluster-vol-of-vol.png)

Every PNG has an SVG export. Market/model IVs use observed market forwards with one call payoff. Original native model puts use model-forward parity in the saved scores. Curves are unsmoothed; unavailable IVs are marked red and signed prices remain in the CSVs. Colors are consistent across panels. Profile charts omit sampling bands for legibility; all SEs remain saved. Raw weighted-tail plots exclude signed control corrections; actual option-price plots include them.

VIX-focused fitting frees curve/floor without joint guards and is weakly identified outside VIX windows. Its SPX costs mix several changes and do not measure cap cost alone. All data remain in sample. Separate checking streams were used before and are not untouched holdout or market OOS. Bands show pointwise sampling uncertainty, excluding deterministic bias, frozen-pilot and parameter uncertainty. Large gaps and numerical second-moment issues remain.

`smile-shape.csv` reports ATM/right-wing IV levels, MAE, descriptive slope and quadratic curvature of RAW IV observations. These summaries do not smooth or alter the priced curves. `manifest.json` records source and plot hashes; run finalization records visual inspection.
