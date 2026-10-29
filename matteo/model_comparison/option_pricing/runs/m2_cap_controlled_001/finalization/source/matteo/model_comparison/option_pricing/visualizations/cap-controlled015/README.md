# A cleaner cap comparison

The selected baseline is unchanged. Gentler volatility brakes improve high-strike VIX curves, but both tested representatives worsen actual SPX fitting errors. Overall VIX price gains are smaller than the reported uncertainty on new draws.

Floor and all expected-variance levels outside VIX windows are exactly fixed. Only four small window adjustments are fitted. The OU reservoir, rates, grids, correlations, initial state and Q scenario remain unchanged. All 1511 original OTM quotes, seven SPX/four VIX maturities, four forwards and 4096 samples per strike remain active.

Start with:

- [VIX smiles, all four maturities](vix-all-maturities-iv.png)
- [Actual SPX tradeoffs, all seven maturities](spx-all-maturities-iv.png)
- [VIX on genuinely new checking draws](vix-fresh-checking.png)
- [Small, confined variance-curve changes](bounded-variance-curve.png)

Additional views:

- [Actual VIX forwards](vix-forwards.png)
- [Nine-shape tradeoffs](nine-shape-tradeoffs.png)
- [Signed VIX price errors](vix-all-maturities-signed-errors.png)
- [Signed SPX price errors](spx-all-maturities-signed-errors.png)
- [New-draw weighted tail calls](fresh-weighted-tail-calls.png)

Each PNG has an SVG export. Blue is the unchanged baseline; green raises the cap by .375; orange permits slow growth beyond the original transition. Gray dashed VIX curves use the original cap with the four-window fitting design. Market is black, with gray bid/ask bands where shown.

All curves use actual prices, with no smoothing or clipping. Market puts use market-forward parity; native model puts use model-forward parity. Bands show ±2 pointwise scramble sampling SE, excluding deterministic, pilot and parameter errors. Raw tail plots exclude signed control corrections; saved controlled prices include them. New draws check the SAME in-sample market observations, not market out-of-sample data. SPX has actual training-cache repricing, without a fresh SPX ensemble.

Small IV changes can create substantial fitting-error changes because quote price scales are narrow. Doubled loss does not mean doubled option prices. All nine local fits reach an active common-call safeguard; the experiment does not exhaust the model family.

[Full findings](../../CAP_CONTROLLED_RESULTS.md) · [Plan](../../CAP_CONTROLLED_PLAN.md) · [Run evidence](../../runs/m2_cap_controlled_001/run.json). The fresh-cache setup failure and correction are preserved; no fitting, SPX run or new random stream was repeated during recovery.
