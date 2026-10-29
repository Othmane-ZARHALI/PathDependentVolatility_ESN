# VIX-only curve and post-cap factor recalibration, 26 October 2026

User authorized point 1: reassess the existing normalized post-cap OU multiplier
while allowing the variance curve to move. SPX stays paused; no new readout/bank.

Keep the existing 16-piece xi curve, including all original knots. Seven pieces
intersect [27,57], [55,85], [118,148], [153,183] calendar days: indices
1,2,3,4,7,8,10 (zero based). Optimize those seven independently. The other nine
have exactly zero sensitivity to these VIX contracts and retain their old values.
This is equivalent to freeing all 16 with the nine zero-gradient entries at their
initial values; it is not a three-piece curve or four broad window adjustments.
No monotonicity or smoothing penalty ties the seven levels together.

Fixed: old readout/cap, positive beta_J, OU banks, premiums, initial state, floor
0.00264. Each active xi is bounded between the floor and 0.20, the existing curve
upper bound. Compare old baseline; curve-only; gamma=1/half-life=7d with free
curve; and free curve plus gamma in [0,3], half-life in [1,365] calendar days.
Finite multistart least squares, at most 120 evaluations per start. This is an
exploratory bounded fit, not a global optimum or production-model adoption.

Use all 217 original OTM quotes, native put/call bid-ask residuals with existing
quote scales, equal weight per expiry, plus four forward residuals /0.5 point.
No SPX term, no artificial wing weighting. Report whole-wing shape separately.
Keep signed prices and unavailable IVs; no clipping/smoothing of prices or IVs.

Training reuses the previous post-cap study's *checking* draws and stored old-G
conditional ratios (4096 per expiry). These are now training data. Frozen pilot
controls remain unchanged. Freeze candidate parameters before generating fresh
37-dimensional exact importance-mixture draws: seed base 2026102691, eight
1024-point scrambles per expiry. No self-normalization. This numerical check is
on the same market quotes, not out-of-sample market validation.

Verify exact replay of the old baseline on reused draws; finite-difference price
Jacobian; exact interval-length second moments; fresh paired losses and SEs;
weighted moment/tail diagnostics; sampled GH32/time8 to GH48/time16 refinement.
Preserve source snapshots and hashes, all optimizer outcomes, prices and plots.
Preserve the previous selected model and all historic studies.
