# Joint parameter adjustment after uniform-path repair — 24 October 2026

User explicitly requests parameter adjustment, emphasizing remaining VIX misfit.
Keep all eleven agreed parameters, beta_J>0, fixed OU/Q/initial state/geometric
settings, all1294 SPX cutoff3 OTM inputs/seven dates,217 VIX inputs/four dates and
four forwards. One call payoff with observed market parity conversion, native
model OTM prices reconstructed with MODEL forward,4096 paths at every strike.
Preserve original family weights, quote scales, penalties, GH32/GH48/time8 and
stock512. Bounds match the last broader search; no automatic bound expansion.
Start at the last selected positive-shift vector, which is a retained incumbent.

Training: reuse uniform4096 cache and independently fitted frozen controls.
Independent numerical check: NEW4096 SPX paths and8x512 VIX points with new seed
2026102411. Regenerate the same known-zero control functions on those primitives;
reuse coefficients/pilot model, do not fit controls on checking draws. No checking
prices are read until final selection. This is not market out-of-sample testing.

Algorithm: bounded sequential local approximation of SPX, exact VIX evaluations.
At each anchor calculate the full joint4096-path price/Jacobian. Inside a box in
parameter units, use its first-order SPX prices and actual nonlinear4096-point VIX
prices/forwards, with the same bid/ask interval objective. Solve with bounded TRF
least squares and analytic derivatives, at most10 inner evaluations per cycle.
Every proposed outer step receives a FULL true4096-path joint repricing/Jacobian.
Only true prices can update the incumbent. Accept a step only if true joint loss
and VIX-option-plus-forward loss both improve, in BOTH canonical-call and
native-OTM conventions. This prevents a forward convention change from being
reported as a VIX improvement. Rejected steps shrink the box;
no surrogate score is reported as an actual price/fit. Keep all rejected proposals.
Maximum5 outer proposals, initial radius0.1*SCALE; shrink by0.5 on rejection,
grow at most1.5 up to0.3 on accepted agreement. Stop after90 optimization minutes
between evaluations; preparation/checking separate. This bounded attempt is not
a convergence/global optimum or identification certificate.

Before optimization: exact baseline replay against saved uniform repricing;
production VIX directional-Jacobian check; one-sided full-production directional
price check at the initial vector. Independent implementation tests for the new
local-approximation/acceptance logic. Retain signed prices and shape failures;
no clipping or IV smoothing. Existing original-OTM state-component accounting is
completed in the new fit wrapper; scalar prices/loss formulas are unchanged.

After selection: reprice incumbent and final on independent4096 checking streams,
report both common-call and native-OTM objectives, SPX/VIX/forward components,
every maturity, SEs, forward errors, shape/sign and raw-tail diagnostics. Paired
block sensitivity is descriptive; independent-pilot/deterministic error excluded.
Produce new all-maturity plots. Capture sources/config/seeds/evaluations/proposals,
full model vector, control coefficients, path arrays and source/input integrity.
Do not silently promote a failed check or claim VIX optimality after the budget.
