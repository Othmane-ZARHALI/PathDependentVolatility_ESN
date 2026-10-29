# Dedicated stochastic volatility response: results

25 October 2026. Completed bounded fixed-curve experiment. The selected baseline remains `runs/m2_maturity_curve_fit_001/selected_model.json`. New run: `runs/m2_stochastic_vov_001`. No new model is adopted.

**Main finding:** the selected extra state brings right-wing IV levels closer to market, but makes the measured right-wing slope flatter at all four maturities under BOTH caps on fresh simulations. Near-ATM and joint price fit worsen. This specified spike-only layer has not solved the steepness problem.

## What we tried, in everyday language

The new state is a sensitivity dial. A high setting makes the existing spike factor move volatility more strongly; its half-life determines how long that episode tends to persist. This changes the distribution of future volatility while keeping its expected variance curve fixed through a new model-specific normalizer. It is an extra continuous independent OU state, with no jumps or changes to the existing reservoir.

We crossed 1.5×/2× sensitivity ranges with 7/30/90-day half-lives and two outer transformations: the original cap and the previously tested logarithmic continuation of strength .5. Both transformations also have a zero-strength reference. All16 variance levels, floor .00264, original readout coefficients, OU rates/geometric grids/state/correlations and risk-price scenario remain EXACTLY fixed. There is no curve or floor refit in this experiment.

## Equations and mathematical scope

Under Q:

```text
dZ = -aZ Z dt + sqrt(2aZ) dBZ, Z0=0, BZ independent;
A(Z)=exp(s tanh Z), aZ=log(2)*365/half_life_days;
eta=cC C+qY(Y-betaY)^2+qJ A(Z)(F-betaJ)^2;
G=exp(2 h(eta)), M(t)=E_Q[G(t)];
v=f+(xi(t)-f) G/M(t).
```

The original cap has h=L-softplus(k(L-eta))/k, with L=5.502780056785304 and k=2. The continuation adds .5 log(1+softplus(k(eta-L))/k). A stays between exp(-s) and exp(s). At s=0 the layer disappears exactly. The extra Brownian variance response is sqrt(2aZ)*partial_Z(v), where partial_Z(eta)=qJ A s sech²(Z)(F-betaJ)². Thus the new source directly changes spike sensitivity; cap saturation can still damp it.

A is centered at one in log space; its arithmetic mean need not equal one. The G/M normalization preserves expected variance xi, not the average spike coefficient. The comparison tests this specified modulation as a whole, including both its random sensitivity and its nonlinear average-response effects.

The capped multiplier remains bounded. With the continuation, bounded A gives G at most quadratic growth in the enlarged Gaussian state. Gaussian multiplier moments remain finite, and the finite-horizon linear-growth/nonexplosion stock-martingale argument in CAP_TAIL_MATHEMATICS.md extends to this layer. These mathematical arguments are distinct from quadrature/sampling checks and are not new Lean proofs. They do not certify finite long-maturity stock-payoff or importance-estimator variance for the continued model. The new Q mean is zero; no P calibration or added risk-premium identification is claimed.

## Complete VIX grid, entire curve fixed

|Case|Cap continuation|Sensitivity range|Half-life days|Native VIX + forwards|Wing IV MAE|Near-ATM IV MAE|Unavailable wing IVs|
|---|---:|---:|---:|---:|---:|---:|---:|
|cap0_ref|0|1×|—|27.166865|33.646|12.389|0|
|cap0_z1.5_h7|0|1.5×|7|41.784222|27.158|14.577|0|
|cap0_z1.5_h30|0|1.5×|30|67.626351|22.535|15.252|0|
|cap0_z1.5_h90|0|1.5×|90|72.314875|22.201|15.447|0|
|cap0_z2_h7|0|2×|7|82.078282|24.340|17.195|0|
|cap0_z2_h30|0|2×|30|228.413554|16.884|20.533|0|
|cap0_z2_h90|0|2×|90|315.722845|15.133|22.278|0|
|cap05_ref|0.5|1×|—|28.941170|27.509|12.467|0|
|cap05_z1.5_h7|0.5|1.5×|7|70.344886|19.342|14.833|0|
|cap05_z1.5_h30|0.5|1.5×|30|122.754792|14.928|15.612|0|
|cap05_z1.5_h90|0.5|1.5×|90|135.612072|14.645|16.044|0|
|cap05_z2_h7|0.5|2×|7|150.900985|15.390|17.680|0|
|cap05_z2_h30|0.5|2×|30|390.488913|13.088|20.686|0|
|cap05_z2_h90|0.5|2×|90|544.608862|11.986|22.418|0|

IV MAE is in percentage points, averaging maturities equally. Wing uses K/marketF>2; near ATM uses |log(K/marketF)|<=.2. Slopes and curvatures are descriptive regressions of raw implied volatilities, not a polynomial volatility model. Complete per-maturity values are saved in `grid_smile_metrics.csv` and `profiles.json`.

## Wing level is different from wing steepness

A lower wing IV error can come from lifting the whole smile, even when the right-wing slope becomes flatter. Selection by wing IV MAE does not select the best slope. Compare the actual slope and curvature below; they are measured against log strike over K/marketF>2.

|Case|Days|Mean wing IV|Wing slope|Wing curvature|
|---|---:|---:|---:|---:|
|Market|27|214.197|52.379|-35.649|
|Market|55|181.151|38.312|-17.605|
|Market|118|136.570|31.230|-3.911|
|Market|153|121.795|27.319|0.303|
|cap0_ref|27|185.767|35.398|-28.622|
|cap0_ref|55|143.462|20.294|-16.928|
|cap0_ref|118|100.979|12.325|-15.036|
|cap0_ref|153|88.921|10.536|-11.238|
|cap05_ref|27|193.743|40.538|-22.881|
|cap05_ref|55|150.307|24.846|-12.349|
|cap05_ref|118|106.121|15.952|-11.334|
|cap05_ref|153|93.508|13.769|-7.427|
|cap0_z2_h90|27|211.426|33.976|-23.423|
|cap0_z2_h90|55|175.145|11.433|-23.770|
|cap0_z2_h90|118|108.349|6.588|-19.509|
|cap0_z2_h90|153|105.330|4.491|-19.473|
|cap05_z2_h90|27|222.640|40.132|-11.252|
|cap05_z2_h90|55|186.321|19.012|-13.592|
|cap05_z2_h90|118|114.947|11.917|-12.914|
|cap05_z2_h90|153|112.380|9.653|-14.428|

The term “wing improvement” in MAE tables refers to closeness in IV level, not a claim of greater or better-matched steepness.

The native option-loss decomposition shows where the price-fit tradeoff occurs. These contributions sum to the saved VIX option loss; forward penalties are separate.

|Case|Near-ATM price loss contribution|Other strikes|Right-wing price loss contribution|
|---|---:|---:|---:|
|cap0_ref|6.3826|9.6002|10.2934|
|cap05_ref|11.3951|9.2636|8.0801|
|cap0_z2_h90|182.5170|119.1927|2.7971|
|cap05_z2_h90|298.2177|219.6744|2.8215|

The normalizer matters here: a larger raw sensitivity does not increase every normalized volatility outcome. The new model divides by its own average multiplier. For the capped case, increasing M lowers the normalized height of already saturated states. This can help explain why smile levels and wing slopes respond differently; it is not a proof of a unique causal decomposition.

|Cap|Selected sensitivity/half-life|M_new / M_no_layer at153 days|
|---|---|---:|
|cap0|z2_h90|2.4082|
|cap05|z2_h90|3.0743|

All six candidates and all four dates are saved in `normalizer_comparison.csv`. Expected variance still equals the same xi curve in the mathematical model; these ratios describe redistribution across states.

Training-only wing winners: `cap0_z2_h90`, `cap05_z2_h90`. Selection froze before SPX prices or fresh checking. The union of their sensitivity/half-life settings was crossed with BOTH caps, plus the two no-layer references. No aggregate, forward or near-ATM acceptance was imposed by this descriptive selection. The full-SPX subset does not establish the best joint point over the entire grid.

## Actual full SPX repricing

|Case|SPX loss|Native VIX + forwards|Native joint loss|Native prices inside bid/ask /1511|Negative native prices|Unavailable call IVs|
|---|---:|---:|---:|---:|---:|---:|
|cap0_ref|43.777476|27.166865|70.977709|139|0|0|
|cap05_ref|112.991814|28.941170|141.966352|84|0|0|
|cap0_z2_h90|561.510528|315.722845|877.266740|63|9|19|
|cap05_z2_h90|1376.764349|544.608862|1921.406579|60|9|22|

Every SPX entry above recomputes volatility and stock at all seven agreed maturities, using all4096 paths at every strike. Old exogenous OU states, equity increments, likelihood weights and frozen reference/hedge controls remain reusable; the terminal stock paths do not. The independent Z process is streamed at the SAME existing left-point mesh using exact OU transitions and parameter-independent counter-based Gaussian innovations. There is no additional multi-gigabyte all-time Z cache.

|Layer case|Comparison reference|SPX loss change ± paired block SE|Native joint change ± paired block SE|
|---|---|---:|---:|
|cap0_z2_h90|cap0_ref|+517.733052 ± 78.768516|+806.289032 ± 103.702773|
|cap0_z2_h90|cap0_ref|+517.733052 ± 78.768516|+806.289032 ± 103.702773|
|cap05_z2_h90|cap0_ref|+1332.986873 ± 132.995730|+1850.428870 ± 168.047704|
|cap05_z2_h90|cap05_ref|+1263.772535 ± 131.071933|+1779.440226 ± 165.845135|

Paired eight-block sensitivities are descriptive numerical diagnostics, not a population confidence or convergence certificate. Loss scores measure fitting error; changes in scores are not proportional changes in option prices. Native/common prices and finite monotonicity/convexity flags are retained in score/CSV artifacts, with no clipping or smoothing.

Negative native estimates here occur in SPX. VIX can separately have unavailable call IVs because the graph deliberately uses the observed market forward while the model forward differs: a positive model call can fall below the intrinsic bound associated with that market-forward Black inversion. These two kinds of failure are recorded separately; neither is hidden by changing the plotting forward.

The negative SPX estimates have magnitudes at most 0.816 of their reported sampling SE. They are compatible with finite-estimator noise around small positive prices; they are not evidence of negative theoretical option prices.

|Case|1030-day raw weighted stock mean (target1)|Raw iid SE|
|---|---:|---:|
|cap0_ref|0.962443|0.009266|
|cap05_ref|0.962786|0.009199|
|cap0_z2_h90|0.960550|0.009158|
|cap05_z2_h90|0.962104|0.009012|

The long-maturity raw stock-mean discrepancy already exists in the baseline. Exact mean-one controls remain mathematically valid, but do not certify sufficient finite-sample accuracy. New negative signed estimates or missing IVs must be read with the reported sampling errors and loss sensitivities, not as model-level arbitrage or a stock-martingale counterexample.

## Genuinely new joint VIX checking draws

|Case|New native VIX + forwards|Wing IV MAE|Near-ATM IV MAE|Unavailable wing IVs|Native VIX loss change vs baseline ± paired net SE|
|---|---:|---:|---:|---:|---:|
|cap0_ref|27.597211|33.841|12.670|0|-0.000000 ± 0.000000|
|cap05_ref|29.407494|27.798|12.792|0|+1.810283 ± 0.156578|
|cap0_z2_h90|272.754024|15.072|18.994|0|+245.156813 ± 6.418237|
|cap05_z2_h90|482.047297|11.733|18.555|0|+454.450086 ± 9.575264|

After selection,32 new complete scrambled nets were generated, eight512-point nets per maturity, in the enlarged37-dimensional Gaussian space. Old defensive-mixture shifts were extended by a zero Z component, retaining exact likelihood weights without self-normalization. New joint draws improve the dimension pairing over the training independent-1D augmentation; neither construction guarantees adequate sampling of joint high-response tails. Old independently fitted control coefficients remain frozen. Checking cannot trigger reselection. The market quotes are the SAME in-sample observations; no market out-of-sample claim is made. No fresh SPX ensemble was generated.

## Weighted tails and variance targets

|Case|Days|E[VIX²] estimate|Fixed curve target|Eight-net SE|Difference / SE|
|---|---:|---:|---:|---:|---:|
|cap0_ref|27|479.2426|477.2295|0.8269|+2.43|
|cap0_ref|55|536.0766|536.0269|1.3177|+0.04|
|cap0_ref|118|582.6876|580.8043|3.6484|+0.52|
|cap0_ref|153|606.4096|610.9326|2.4656|-1.83|
|cap0_z2_h90|27|478.2398|477.2295|3.9645|+0.25|
|cap0_z2_h90|55|522.7882|536.0269|10.0420|-1.32|
|cap0_z2_h90|118|566.9133|580.8043|13.8951|-1.00|
|cap0_z2_h90|153|608.1350|610.9326|12.2573|-0.23|
|cap05_ref|27|479.1784|477.2295|0.8880|+2.19|
|cap05_ref|55|536.0301|536.0269|1.4011|+0.00|
|cap05_ref|118|582.7892|580.8043|3.9063|+0.51|
|cap05_ref|153|606.3876|610.9326|2.6806|-1.70|
|cap05_z2_h90|27|477.9675|477.2295|4.0877|+0.18|
|cap05_z2_h90|55|520.4456|536.0269|10.6073|-1.47|
|cap05_z2_h90|118|564.9762|580.8043|14.9676|-1.06|
|cap05_z2_h90|153|608.2600|610.9326|13.3493|-0.20|

Because the entire xi curve stays fixed, the target is identical for every case at each maturity. Historical 153-day second-moment discrepancies are preserved; a new estimate closer to target is not a convergence pass.

|Case|Days|Highest strike|Proposal hits|Weighted tail probability|Raw call ± SE|Payoff ESS|Top10 payoff share|
|---|---:|---:|---:|---:|---:|---:|---:|
|cap0_ref|27|100|257|0.00049549|0.0105348 ± 0.00035|176.6|9.9%|
|cap0_ref|55|170|146|7.4576e-05|0.00222868 ± 9.1e-05|99.2|17.6%|
|cap0_ref|118|200|120|3.6452e-05|0.00103828 ± 9.7e-05|78.4|21.0%|
|cap0_ref|153|200|104|4.343e-05|0.00127605 ± 0.00014|56.3|30.0%|
|cap0_z2_h90|27|100|231|0.0015376|0.0415259 ± 0.0051|82.3|24.2%|
|cap0_z2_h90|55|170|119|0.00025027|0.0090017 ± 0.0012|44.2|37.7%|
|cap0_z2_h90|118|200|112|0.0001593|0.00530376 ± 0.00096|38.3|42.0%|
|cap0_z2_h90|153|200|79|0.00013586|0.00492899 ± 0.0011|28.1|50.5%|
|cap05_ref|27|100|279|0.00063961|0.0179716 ± 0.00057|200.6|9.0%|
|cap05_ref|55|170|181|0.00013155|0.00579324 ± 0.00024|122.2|14.8%|
|cap05_ref|118|200|157|8.0731e-05|0.00341842 ± 0.00026|102.1|16.5%|
|cap05_ref|153|200|142|9.7677e-05|0.00431284 ± 0.00039|82.5|22.8%|
|cap05_z2_h90|27|100|256|0.0022129|0.0725544 ± 0.0082|85.5|23.4%|
|cap05_z2_h90|55|170|150|0.0003918|0.0232637 ± 0.0034|49.6|34.8%|
|cap05_z2_h90|118|200|144|0.00041521|0.0201986 ± 0.0029|53.5|33.7%|
|cap05_z2_h90|153|200|114|0.00034707|0.0173605 ± 0.003|40.2|40.5%|

Raw hit counts are proposal events, not Q probabilities. Tail probabilities and calls use exact likelihood weights; raw calls exclude the signed control corrections used in actual option estimates. Low payoff ESS or concentrated contributions remain visible and limit interpretation.

## Verification and remaining numerical limits

- Three independent tests check the Gaussian-response integral against elementary and dense-grid oracles, counter-normal replay/OU variance, and streamed stock paths against a separate Python implementation. Zero-strength stock prices nest the old kernel for both cap choices.
- Original baseline full price replay maximum absolute difference: 1.36e-09. Independent saved-draw audit replays 22 price/SE/parity/loss records, 270 archived sources and all32 complete new joint scrambles. The selected baseline file is unchanged; selection timing is verified.
- Checked response interpolation uses degrees [32]; maximum sampled off-node relative error 7.94e-07. Six off-node factors per table include selected extreme old-state means. This checks finite quadrature/interpolation, not a uniform error bound over all states.
- A separate numpy four-dimensional Gaussian implementation checks 8 conditional/unconditional expectations, maximum relative error 5.63e-10.
- Representative time order8-to16 checks at ordinary/extreme153-day checking states change VIX by at most 1.29e-08 points.
- Representative original-state normalizer GH48-to96 (also refining new-state192-to256) changes M by at most 5.24e-08 relatively; all six added-state cases and both caps are checked at27/153/1030 days.
- cap0, new_Z_conditional_128_to_256: maximum representative relative change 8.91e-08.
- cap0, new_Z_normalizer_192_to_256: maximum representative relative change 1.48e-09.
- cap0, old_Gaussian_conditional_32_to_48: maximum representative relative change 3.82e-06.
- cap05, new_Z_conditional_128_to_256: maximum representative relative change 3.23e-08.
- cap05, new_Z_normalizer_192_to_256: maximum representative relative change 6.01e-10.
- cap05, old_Gaussian_conditional_32_to_48: maximum representative relative change 1.14e-06.
- Conditional original-state GH32, normalizer original-state GH48, additional-state GH128/GH192, time order8 and cached stock512 mesh are approximations. Sampled refinements are not whole-surface price-convergence proofs. Standard errors omit integration/stock biases, pilot/parameter uncertainty and unsampled tails. New-state terminal moment diagnostics are saved in the audit.

## Decision

Fresh simulations show flatter right-wing slopes at all four maturities for the extra state relative to its own no-layer cap reference. Improved wing IV levels therefore do not solve the original steepness question.

None of the matched layer representatives improves both sampled native SPX and VIX-plus-forward losses versus the selected baseline. Preserve the baseline. This bounded experiment does not establish that every continuous stochastic-response model fails.
For cap0, its wing-selected extra state improves training wing IV MAE from 33.646 to 15.133; near-ATM MAE changes from 12.389 to 22.278. Read this alongside actual forwards, price losses, SPX and new-stream evidence above.
For cap05, its wing-selected extra state improves training wing IV MAE from 27.509 to 11.986; near-ATM MAE changes from 12.467 to 22.418. Read this alongside actual forwards, price losses, SPX and new-stream evidence above.

This tests one spike-only bounded modulation and the predeclared two-cap grid, without fitting xi or the old coefficients. The earlier cluster experiment, a different modulation, and a possible future feedback modulation are distinct model questions. No full calibration convergence or structural impossibility follows.

Plan: [STOCHASTIC_VOV_PLAN.md](STOCHASTIC_VOV_PLAN.md). Figures and reading guide: [stochastic-vov016](visualizations/stochastic-vov016/README.md). Reproduce with the pinned pdv-options-venv and package PYTHONPATH using modules `stochastic_vov_study --output NEW_DIRECTORY`, `audit_stochastic_vov --run RUN`, `plot_stochastic_vov --run RUN --output PLOTS`, `summarize_stochastic_vov --run RUN`. Source snapshots and raw draws are in the run.

