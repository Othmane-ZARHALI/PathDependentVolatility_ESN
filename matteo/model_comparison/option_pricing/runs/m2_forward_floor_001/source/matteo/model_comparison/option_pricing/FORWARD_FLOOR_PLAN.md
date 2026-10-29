# Tight VIX forwards, followed by a fitted minimum variance

Authorized: user requests points1and2 of the left-end proposals,26October2026.
Keep the successful additive M2/stress response; no new component or timescale.
Keep all217 native OTM quotes, four forwards,16 curve intervals with seven
independently fitted VIX-sensitive pieces, old readout/cap/banks and beta_J>0.
SPX stays paused. Preserve the prior selected model and every historical run.

v=f+(xi-f)[(1-w)R_old+w L_b], with the same independent normalized OU exponential
L_b, b in[0,5], w in[0,1], half-life in[1,365] calendar days.

Compare on common simulations:
- Prior additive price winner, unchanged, including its forward mismatches.
- Fit1: same fixed floor0.00264, all four VIX forwards as equality constraints.
- Fit2: same equality constraints, floor volatility free in[5%,12%], f=sigma_f².

xi_j=f+(0.20-f)y_j, y_j in[0,1], independently for the seven active intervals.
The other nine xi entries retain old values, all above the largest allowed floor.
This guarantees xi>=f for every curve piece without tying active levels together.
The lower5% bound preserves the agreed minimum-volatility requirement; the12%
upper bound is an exploratory bound below the lowest retained VIX strike12.5.
Both fits keep curve values and b/w/half-life free; fit2 does not force f upward.

Minimize native OTM bid–ask distance /existing quote scales with equal expiry
weights, subject to E[VIX_T]=the four archived forward targets. Six fixed-floor
starts, then seven free-floor starts; at most120SLSQP iterations/start. A bounded
curve-only least-squares prefit makes starts forward-feasible. Only candidates
with maximum training forward error<=1e-5point are eligible for selection. Save
all failures/outcomes. Reference native objective with old0.5point forward scale
is reported separately for comparisons; do not mix objective definitions.

The equality is imposed on the training estimator, not asserted as exact in the
population. Training reuses the prior mixture run's16384-point fresh draws and
stored Gaussian conditional ratios, now explicitly training inputs. Existing
independent pilot controls and exact importance weights remain frozen. Candidate
selection is frozen before new32768-point/expiry checking: eight4096-point nets,
the same68-component joint37D proposal with Z shifts0/2/4/6, seed base20261026083.
No self-normalization. Market targets stay in sample.

Check analytic native-price/forward Jacobians in all11unit-scaled directions,
including the moving floor and changing xi; replay the old additive model on the
training draws; independently replay saved path integrals, prices, parity,
forwards and exact second moments. Report signed prices, unavailable IVs, low-put
errors, raw-weighted tail concentration/ESS, right-wing changes and paired
eight-net uncertainty. Refine sampled GH32/time8 to GH48/time16, including large
weighted payoffs. Keep original market-forward call-IV plots and additionally
show native-put prices/IVs using each model's actual forward, without modifying
prices to hide residual forward error. No new Lean or convergence claim.
