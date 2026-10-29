# Research decision: variance levels versus SPX skew

Saved 28 October 2026 at the user's explicit request. Evidence: [full result](PROTECTED_CURVE_RESULTS.md), [plan](PROTECTED_CURVE_PLAN.md), run `runs/m2_protected_curve_001`, [SPX plot](visualizations/protected-curve028/SPX-protected-smiles.png).

**Remember both findings:** the protected nine-level fit lowers fresh SPX native price loss by 34.9% while leaving the VIX distribution unchanged, but equal-maturity IV MAE worsens 4.345→4.710 percentage points and the SPX downward skew remains missing. It is not a satisfactory joint calibration.

Only xi indices 0,5,6,9,11,12,13,14,15 moved. Both stress factors, their weights/speeds/sensitivities, the floor, old M2 and seven VIX-sensitive levels stayed fixed. All seven SPX maturities/1294 primary quotes were fitted; all 1751 eligible quotes checked. All 217 VIX quotes/four forwards remain protected. Two bounded starts; training selection before genuinely fresh 512-path paired Q validation. The 128-step paired loss-improvement interval excludes zero, but final sampling/grid/normalizer accuracy remains unresolved.

Do not spend another iteration rediscovering that lowering the variance curve can improve the price objective while worsening left-wing IVs. The next proposed mechanism is stock/stress leverage, not more stress speeds or another blind curve fit. [A driver-rotation construction](STOCK_DRIVER_ROTATION_PROPOSAL.md) can preserve the variance-factor joint law and VIX. It remains a proposal: no correlation parameter has been implemented or fitted.

This is a bounded empirical finding at one market snapshot, not proof that every possible curve, continuous stress measure or model in the family is incapable of joint calibration. New data or a clearly different hypothesis can justify revisiting it. Preserve price and IV metrics, numerical uncertainty and the successful VIX baseline when continuing.
