# Controlled cap comparison: a cleaner experiment

25 October 2026. Completed bounded investigation. The selected baseline and every historical run remain unchanged. New evidence is in `runs/m2_cap_controlled_001`; all agreed maturities and 4096 paths per strike are retained. These are experimental comparisons, not adopted joint calibrations.

## In everyday terms

We eased the brake on extreme volatility while stopping the fit from changing unrelated parts of the expected-variance curve. This makes the comparison much easier to interpret than the previous VIX-only fit. The experiment separates whether the brake can help from the damage caused by moving variance at times that VIX options do not observe.

The model still rescales volatility by its own average multiplier. Consequently a gentler brake redistributes variance between ordinary and extreme states; it does not simply increase every option price. We kept that rescaling in every candidate.

## Exactly what changed

The variance floor stays at .00264 (5.138% instantaneous volatility), with positive beta_J unchanged. All original score coefficients, OU rates/geometric grids, initial state, correlations and Q risk-price scenario are fixed. SPX 29/57/119/211/302/666/1030 days, VIX 27/55/118/153 days, all four forwards and 1511 original OTM targets remain. Market puts use market-forward parity; native model puts use model forwards. No jumps, quintic proposal, smoothing or price clipping.

Only four broad adjustments can move expected variance: one over each VIX window [27,57], [55,85], [118,148], [153,183] days. Where the first two windows overlap, their adjustments are averaged. With B the normalized window-indicator matrix and s the four adjustments:

    xi_new = xi_baseline + (xi_baseline - floor) * (B @ s), |s_j| <= 0.15.

This keeps every off-window curve piece exactly fixed and changes excess variance by at most 15% per piece. xi remains above the floor, with no monotonicity condition. All boundaries already exist in the 16-piece curve. The 15% range is a declared experimental choice, not a user requirement.

Nine shapes cross cap-level increases 0/.375/.75 with continuation strengths a=0/.25/.5. Sharpness remains 2. At a=0 the model retains the logistic cap. At a>0 it adds a*log(1+softplus(k*(eta-L))/k) to the tempered score, permitting slow growth beyond the transition. Since a<=.5, the controlled-growth moment and stock-martingality arguments in [CAP_TAIL_MATHEMATICS.md](CAP_TAIL_MATHEMATICS.md) continue to apply under their stated assumptions.

Each shape fits only the four window adjustments using native VIX price-plus-forward loss. Diagnostic screens require forwards within .5 point of market and common-call VIX-plus-forward loss no worse than the unchanged baseline. Those screens are implementation choices. All nine local fits terminate successfully, but this finite grid and local optimization are not a global convergence result.

## All nine training outcomes

|Shape|Cap-level increase|Continuation a|Native VIX + forward loss|Right-wing IV MAE|Near-ATM IV MAE|Preferred screen|
|---|---:|---:|---:|---:|---:|:---:|
|Unchanged baseline|0|0|27.166865|33.646|12.389|reference|
|shape00|0|0|27.347687|33.606|13.127|no|
|shape01|0|0.25|26.438674|30.297|14.071|yes|
|shape02|0|0.5|26.910390|26.664|15.180|yes|
|shape03|0.375|0|26.865374|27.717|14.921|yes|
|shape04|0.375|0.25|28.099833|24.972|15.789|no|
|shape05|0.375|0.5|31.030736|21.889|16.913|no|
|shape06|0.75|0|30.490362|22.716|16.674|no|
|shape07|0.75|0.25|33.521607|20.384|17.602|no|
|shape08|0.75|0.5|38.442693|17.735|18.642|no|

IV MAE means average absolute implied-volatility error in percentage points, giving each maturity equal weight. Right wing uses K/F>2; near ATM uses |log(K/F)|<=.2. These summarize raw implied volatilities; they do not smooth prices. A preferred candidate must improve native loss, keep all wing IVs available, and worsen near-ATM IV MAE by no more than 5 points. Feasible candidates are ranked by wing MAE within the capped and continued families separately. This selects a useful shape comparison rather than merely the smallest aggregate loss.

Both representatives were frozen using TRAINING results before any new checking states were generated:

- **shape03 (capped)**: L increases by 0.375, a=0; four excess-variance adjustments are +1.950%, +2.242%, +1.670%, +1.324%. Floor/off-window levels are exactly unchanged.
- **shape02 (continued)**: L increases by 0, a=0.5; four excess-variance adjustments are +2.139%, +2.525%, +2.097%, +1.745%. Floor/off-window levels are exactly unchanged.

All nine fits reach the common-call loss safeguard, approximately 45.509602. Several also reach the +.5-point forward-error boundary at 118/153 days. Therefore this experiment does not show that the cap families have exhausted their flexibility: the chosen fitting safeguards are active. Larger cap/continuation changes improve wing errors further, but their fitted native price losses and near-ATM errors worsen in this comparison.

## Actual SPX tradeoffs and genuinely fresh VIX checking

|Candidate|Training SPX loss|Training native VIX + forwards|Training native joint loss|NEW-stream VIX + forwards|
|---|---:|---:|---:|---:|
|Unchanged baseline|43.777476|27.166865|70.977709|26.986598|
|shape03|87.400641|26.865374|114.299387|26.845187|
|shape02|98.795211|26.910390|125.738973|26.914668|

Every SPX value above comes from full nonlinear repricing at all seven dates, recomputing volatility and stock along all 4096 cached OU/equity paths per expiry. Only primitives and valid frozen controls are reused. No terminal stock path was treated as unchanged after a model change.

- shape03: training SPX loss change +43.623165, descriptive paired eight-block jackknife SE 5.958641; native joint change +43.321678, SE 5.892086. New-stream native VIX-plus-forward change -0.141411, paired eight-net SE 0.203733. Strict sampled decreases in all four prior joint/VIX conventions: False. The common-call VIX guard is effectively at equality, so tiny roundoff differences are not meaningful improvement.
- shape02: training SPX loss change +55.017735, descriptive paired eight-block jackknife SE 6.734290; native joint change +54.761264, SE 6.669611. New-stream native VIX-plus-forward change -0.071930, paired eight-net SE 0.228777. Strict sampled decreases in all four prior joint/VIX conventions: False. The common-call VIX guard is effectively at equality, so tiny roundoff differences are not meaningful improvement.

|Candidate|Native prices inside bid/ask /1511|Negative native prices|Unavailable common-call IVs|SPX monotonicity / convexity flags|
|---|---:|---:|---:|---:|
|baseline|139|0|0|1 / 66|
|shape03|94|0|0|1 / 68|
|shape02|87|0|0|1 / 75|

Shape flags concern finite controlled price estimates, not proven model-level arbitrage. The same fixed control functions can introduce finite-sample convexity artifacts; no plotted smoothing or replacement hides them.

The new-stream aggregate VIX gains are smaller than one reported paired jackknife SE. They do not establish a precise improvement in expected aggregate price loss, even though the wing curves visibly move closer to market. SPX loss deterioration persists across all eight leave-one-block comparisons for both candidates. These losses are fitting-error scores; a doubled score does not mean option prices doubled.

This supplies the actual tradeoff after removing the previous off-window/floor confounding. It remains a controlled combination of cap changes and the four small window adjustments, not a pure cap-only intervention. The earlier fixed-curve cap profiles provide that complementary comparison.

The checking experiment uses eight newly generated 512-point scrambled nets per VIX maturity, the unchanged defensive-mixture proposal and exact p/q weights, with no weight self-normalization. Training-pilot control coefficients remain frozen. Selection precedes new draw generation; checking results do not trigger reselection. Seeds are regenerated and checked against recorded historical production/pilot/diagnostic seeds. New random streams check numerical stability on the SAME in-sample market quotes; they are not market out-of-sample observations. SPX comparisons use the original training cache; no fresh SPX ensemble was run.

## Smile shape by maturity

|Candidate|Days|Wing IV MAE|Wing slope|Wing curvature|Near-ATM IV MAE|Unavailable IVs|
|---|---:|---:|---:|---:|---:|---:|
|baseline|27|28.430|35.398|-28.622|27.011|0|
|baseline|55|37.689|20.294|-16.928|13.867|0|
|baseline|118|35.591|12.325|-15.036|4.417|0|
|baseline|153|32.875|10.536|-11.238|4.260|0|
|shape03|27|20.358|39.068|-23.145|33.018|0|
|shape03|55|30.994|23.370|-13.010|18.270|0|
|shape03|118|30.803|14.892|-12.067|4.688|0|
|shape03|153|28.712|12.853|-8.203|3.710|0|
|shape02|27|19.182|39.915|-21.938|33.493|0|
|shape02|55|29.803|24.343|-11.622|18.650|0|
|shape02|118|29.832|15.679|-10.850|4.879|0|
|shape02|153|27.837|13.572|-7.026|3.699|0|

Slope and curvature are descriptive linear/quadratic regressions of raw IV against log strike within the right-wing region; they are not an asymptotic theorem or a new polynomial volatility model. Complete quote prices and pointwise SE remain saved. A steeper wing and a lower average price loss do not imply that every maturity or strike improved.

Fresh-stream equal-maturity IV error summaries:

|Candidate|New-stream wing IV MAE|New-stream near-ATM IV MAE|
|---|---:|---:|
|baseline|33.866|12.590|
|shape03|27.963|15.183|
|shape02|26.919|15.435|

## Weighted tails and second moments on new streams

|Candidate|Days|E[VIX²] estimate|Curve-implied target|Eight-net SE|Difference / SE|
|---|---:|---:|---:|---:|---:|
|baseline|27|477.5413|477.2295|2.1269|+0.15|
|baseline|55|536.2102|536.0269|1.6641|+0.11|
|baseline|118|581.2752|580.8043|1.5921|+0.30|
|baseline|153|608.6317|610.9326|2.2943|-1.00|
|shape03|27|486.3937|486.0736|2.1968|+0.15|
|shape03|55|547.6692|547.4027|1.7559|+0.15|
|shape03|118|590.5370|590.0612|1.6594|+0.29|
|shape03|153|616.2422|618.6692|2.4296|-1.00|
|shape02|27|487.2617|486.9397|2.2022|+0.15|
|shape02|55|549.0945|548.8296|1.7702|+0.15|
|shape02|118|592.9089|592.4311|1.6630|+0.29|
|shape02|153|618.6633|621.1318|2.4474|-1.01|

On the new baseline draws, day 153 gives 608.632 versus 610.933, SE 2.294, about 1.00 SE below target. That is closer than the historical reused-check estimate 606.208 versus 610.933 (SE 1.635), consistent with outer-draw sensitivity. One new eight-net ensemble is not a convergence pass.

The second-moment target changes only with each candidate’s four allowed curve adjustments. A new estimate closer to its target would not erase the historical 153-day discrepancy or certify convergence. Eight nets give limited information about rare tails and normalizer/conditional/time/stock errors remain separate.

|Candidate|Days|Highest strike|Proposal hits|Weighted probability|Raw call ± SE|Payoff ESS|Top10 share|
|---|---:|---:|---:|---:|---:|---:|---:|
|baseline|27|100|253|0.00048135|0.0106904 ± 0.00042|167.4|11.8%|
|baseline|55|170|141|6.7451e-05|0.0018068 ± 0.00016|85.5|20.5%|
|baseline|118|200|123|3.5883e-05|0.00109296 ± 6.2e-05|77.7|21.8%|
|baseline|153|200|100|3.8805e-05|0.00114614 ± 0.00011|53.7|30.6%|
|shape03|27|100|284|0.00067238|0.0177151 ± 0.00053|189.7|10.5%|
|shape03|55|170|173|0.00012184|0.00439378 ± 0.00022|108.2|16.5%|
|shape03|118|200|161|8.2593e-05|0.00311174 ± 0.00015|102.9|17.6%|
|shape03|153|200|130|7.7618e-05|0.00310301 ± 0.00019|70.0|25.4%|
|shape02|27|100|287|0.00070233|0.0191283 ± 0.00057|193.0|10.4%|
|shape02|55|170|182|0.00014131|0.00525775 ± 0.00023|114.0|15.7%|
|shape02|118|200|171|9.6081e-05|0.00392435 ± 0.0002|108.3|17.0%|
|shape02|153|200|138|9.0519e-05|0.0038838 ± 0.00021|74.6|24.0%|

These raw prices use likelihood-weighted payoff contributions. Actual controlled option prices also contain signed zero-mean control corrections. Raw hit counts alone do not establish sufficient risk-neutral tail mass. The table and plots expose payoff concentration and uncertainty without clipping or self-normalizing weights.

## Verification and limits

- Three standard-library tests check exact frozen coordinates, all adjustment-domain corners, overlap behavior, baseline nesting and basis derivatives. Production VIX objective sensitivities are checked for all nine shapes.
- Independent audit passes 15 saved score records, 139 fit-trial objective/forward reconstructions from archived exposures, 270 source snapshot hashes and all 32 complete new scramble/state/likelihood regenerations. It confirms the parent selection is unchanged.
- Representative direct same-GH32 conditional checks on new extreme states have maximum interpolation relative error 2.39e-10. This checks the fast integration table against that finite quadrature rule, not continuous-price convergence.
- GH32 conditional, GH48 normalization, time order 8 and the cached stock512 mesh remain approximations. Earlier representative order refinements are historical evidence, not a fresh whole-surface certificate for every candidate. Finite SE omits unsampled tails, deterministic bias, pilot and parameter uncertainty.
- Finite Gaussian moments of the multiplier and a true stock martingale establish different properties. They do not by themselves prove a finite second moment for long-maturity stock option payoffs or for importance-weighted estimators. The log-continuation candidate has no new population payoff-variance certificate; reported block/net SEs are descriptive numerical evidence.
- Signed prices, unavailable IVs, raw forward-ratio and shape flags remain in saved CSV/JSON/NPZ files. No new Lean/P-compatibility, global optimum or accepted production-fit claim is made.

The first attempt stopped before fresh model pricing because the engine constructor expected stock-cache metadata in the VIX-only directory. The setup was corrected; both fits, both completed SPX evaluations and all 32 saved fresh scrambles were retained. Only fresh VIX pricing resumed. The original failure log/source snapshots and separately captured correction/continuation are preserved; no random streams or candidate selections were replaced.

## Decision and how to read the result

Both tested representatives improve right-wing VIX shape but worsen the actual training SPX loss. Preserve the selected baseline. These candidates do not supply the desired joint improvement, even after unrelated curve changes are removed. The earlier catastrophic SPX losses were confounded by excessive curve freedom; the remaining tradeoff here is a more interpretable cap-plus-small-window-adjustment effect.

The result is not a proof that every cap modification fails. All nine fits hit an active common-call safeguard; stronger tail freedom reduces wing IV errors further but worsens the tested price/ATM compromise. Separate changing the model from changing the fitting requirements. A next mechanism experiment can target stochastic response strength, with persistence independently controlled by a dedicated continuous state, while retaining this clean four-window/fixed-floor design. That additional layer remains a separately specified proposal; it was not implemented or calibrated in this run.

Plots and reading guide: [cap-controlled015](visualizations/cap-controlled015/README.md). Plan: [CAP_CONTROLLED_PLAN.md](CAP_CONTROLLED_PLAN.md). Run: `runs/m2_cap_controlled_001`. Reproduce with modules `cap_controlled_study --output NEW_DIRECTORY` (current driver includes the VIX-only setup correction), `audit_cap_controlled --run RUN`, `plot_cap_controlled --run RUN --output PLOTS`, and `summarize_cap_controlled --run RUN` using the existing pdv-options-venv environment and package PYTHONPATH.

