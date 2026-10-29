# Independent amplitude after the cap: completed comparison

Completed 2026-10-26. Experimental fixed-curve comparison; original selected model unchanged.

**Main finding:** the strong 7-day multiplier improves the descriptive fresh VIX slope/curvature score by about31%, but the smile remains too flat and full SPX loss rises43.777 to121.147. The mild multiplier gives a small training joint-loss reduction70.978 to69.876, smaller than its paired sampling sensitivity (change -1.102, jackknife SE3.933), and does not improve fresh VIX price loss. Neither is an accepted joint fit.

## What changed

The old spike-only bounded sensitivity layer is absent. Starting from the original capped baseline, an independent Q OU state multiplies the normalized variance component AFTER the cap:

```text
dZ = -a Z dt + sqrt(2a) dB, Z0=0, q(t)=1-exp(-2at)
L(t)=exp(gamma Z(t)-gamma^2 q(t)/2)
v(t)=f+(xi(t)-f)*(G_base(X(t))/M_base(t))*L(t)
```

All original readout coefficients, OU banks, correlations, premiums, 16 variance levels, floor and cap stay fixed. gamma is a log-VARIANCE sensitivity, not a bounded multiplier range. The new factor is independent of every original driver. E[L]=1, so the original normalizer M_base is unchanged; there is no candidate-specific renormalization. Relative to the previous spike-only experiment, placement, factor shape and centering all change; this is not a placement-only causal test.

The analytic conditional factor is exp(gamma exp(-a(u-T)) Z_T - gamma^2 exp(-2a(u-T)) q(T)/2). Independence lets it multiply the original conditional expectation of G. No new conditional Gaussian quadrature is introduced.

## Full VIX grid on the original training simulations

|Case|Native VIX + forwards|Wing IV MAE|Near-ATM IV MAE|Wing slope/curvature score|Missing call IVs|
|---|---:|---:|---:|---:|---:|
|baseline|27.166865|33.646|12.389|0.460031|0|
|g0.25_h7|27.415521|33.124|12.519|0.452340|0|
|g0.25_h30|28.779218|32.607|12.866|0.425743|0|
|g0.25_h90|28.174556|32.730|12.838|0.420540|0|
|g0.5_h7|28.046353|32.068|12.646|0.342401|0|
|g0.5_h30|37.465191|30.314|14.067|0.283466|0|
|g0.5_h90|35.714005|30.606|14.075|0.275651|0|
|g1_h7|31.852273|28.975|12.839|0.183825|0|
|g1_h30|122.684181|22.920|17.844|0.223776|4|
|g1_h90|132.654771|23.352|18.910|0.244204|5|

MAE is in IV percentage points, averaging the four maturities equally. Wing: K/marketF>2; near ATM: abs(log(K/marketF))<=0.2. The descriptive shape score averages squared slope and curvature discrepancies, each divided by max(10,abs(market coefficient)) squared. Lower is closer; it is not an arbitrage condition or a fitted polynomial volatility model.

Price selection: `g0.25_h7`. Shape selection: `g1_h7`. These selections were frozen before SPX prices or fresh VIX draws. Neither selection implies joint acceptance. All nine nonzero variants and the zero-factor reference are retained.

## Actual full SPX repricing

|Case|SPX loss|Native VIX + forwards|Native joint loss|Inside original bid-ask /1511|Negative native estimates|Missing call IVs|
|---|---:|---:|---:|---:|---:|---:|
|baseline|43.777476|27.166865|70.977709|139|0|0|
|g0.25_h7|42.426757|27.415521|69.875646|155|2|2|
|g1_h7|121.147012|31.852273|153.032653|79|9|9|

All seven agreed SPX dates are repriced with 4096 samples at every strike. Original exogenous paths and equity shocks are reused; variance and stock paths are recomputed. Independent OU shocks use the same parameter-independent counter streams across candidates. No large additional factor-path cache is created. Original independent pilot controls stay frozen.

- `g0.25_h7`: paired SPX loss change -1.350719, descriptive eight-block jackknife SE 3.694715; native joint change -1.102062, SE 3.933472.
- `g1_h7`: paired SPX loss change +77.369536, descriptive eight-block jackknife SE 21.704057; native joint change +82.054944, SE 23.433859.

## Fresh joint VIX simulations after selection

|Case|Native VIX + forwards|Wing IV MAE|Near-ATM IV MAE|Shape score|
|---|---:|---:|---:|---:|
|baseline|26.852498|33.997|12.718|0.444949|
|g0.25_h7|26.932168|33.534|12.709|0.491890|
|g1_h7|29.700033|28.615|12.729|0.308749|

|Case|Days|Wing slope|Wing curvature|
|---|---:|---:|---:|
|Market|27|52.379|-35.649|
|Market|55|38.312|-17.605|
|Market|118|31.230|-3.911|
|Market|153|27.319|0.303|
|baseline|27|35.836|-8.272|
|baseline|55|19.928|-23.542|
|baseline|118|11.859|-12.834|
|baseline|153|9.709|-9.339|
|g0.25_h7|27|35.503|-1.392|
|g0.25_h7|55|20.701|-23.130|
|g0.25_h7|118|12.036|-13.178|
|g0.25_h7|153|9.908|-9.471|
|g1_h7|27|36.907|-7.909|
|g1_h7|55|28.538|-11.833|
|g1_h7|118|15.719|-10.335|
|g1_h7|153|11.832|-7.595|

These are regressions against log(K/marketF) over K/marketF>2, not pointwise derivatives. Market curvature is not uniformly positive. Both linear- and log-strike figures are exported. Call PRICE convexity and implied-VOLATILITY convexity are different statements.

The selected shape candidate changes the four fresh wing slopes relative to baseline by [1.071, 8.61, 3.86, 2.123]. Its fresh shape score changes from 0.444949 to 0.308749. Interpret this alongside price loss, forwards and SPX rather than declaring success from slope alone.

- `g0.25_h7` fresh paired slope changes: [-0.333, 0.773, 0.177, 0.199]; jackknife SE: [1.082, 0.534, 0.302, 0.347].
- `g1_h7` fresh paired slope changes: [1.071, 8.61, 3.86, 2.123]; jackknife SE: [3.377, 2.482, 1.152, 0.474].

Checking uses 32 newly generated joint37D scrambles: eight512-point nets at each of four maturities, exact defensive-mixture likelihood weights, no self-normalization. Selection cannot be revised from these draws. Market quotes are the same in-sample observations. There is no market out-of-sample claim or fresh SPX ensemble.

## Numerical and mathematical checks

- Audit PASS: 16 saved price/parity/SE/objective records; 32 complete joint scrambles regenerated with exact weights; 11 dense independent conditional-factor integrals; 3 direct NumPy old-state Gaussian oracles; 277 archived/current source comparisons.
- Zero-factor full baseline maximum price discrepancy: 1.36e-09.
- Sampled full VIX windows, old-state GH32/time8 to GH48/time16: maximum change 7.68958e-06 VIX points. These are representative checks, not a whole-price convergence certificate.
- Original normalizer GH48 to GH96 sampled relative changes: maximum 5.49e-10.
- Five unit checks pass across the new implementation and prior spike-layer implementation. New checks compare analytic conditional factors with dense Gaussian integration, test exact mean-one centering, compare streamed stock stepping with independent Python, and verify zero-strength nesting of the original kernel.

For the continuous model, conditioning on the independent Z path makes L a bounded continuous time function on a finite horizon. With capped G_base and positive continuous deterministic M_base, conditional stock volatility is bounded, yielding the stock martingale identity. Fixed-time variance moments are finite. This does not guarantee finite variance of stock payoff estimators. Numerical errors and rare-tail sampling remain separate questions. No new Lean certificate or P-measure calibration is claimed.

Negative signed control estimates and missing IVs are retained in the raw tables. A missing IV can also arise from using the market forward for inversion when the model forward differs; it need not mean a negative model price.

- `g0.25_h7` (full): 2 negative native estimates, maximum magnitude 0.330 reported SE. These are signed finite control estimates; no negative true model-price claim follows.
- `g1_h7` (full): 9 negative native estimates, maximum magnitude 0.870 reported SE. These are signed finite control estimates; no negative true model-price claim follows.

Training squared-VIX moments (weighted estimates, not a convergence pass):

|Case|Days|E[VIX²]|Curve target|Eight-net SE|
|---|---:|---:|---:|---:|
|baseline|27|474.820|477.229|1.496|
|baseline|55|539.139|536.027|2.042|
|baseline|118|577.553|580.804|2.057|
|baseline|153|612.534|610.933|1.457|
|g0.25_h7|27|474.489|477.229|2.102|
|g0.25_h7|55|541.385|536.027|3.052|
|g0.25_h7|118|575.798|580.804|2.265|
|g0.25_h7|153|612.947|610.933|1.967|
|g1_h7|27|471.809|477.229|4.324|
|g1_h7|55|547.088|536.027|6.518|
|g1_h7|118|569.796|580.804|5.198|
|g1_h7|153|616.127|610.933|7.654|

Fresh checking squared-VIX moments (weighted estimates, not a convergence pass):

|Case|Days|E[VIX²]|Curve target|Eight-net SE|
|---|---:|---:|---:|---:|
|baseline|27|477.760|477.229|1.007|
|baseline|55|538.668|536.027|1.336|
|baseline|118|580.747|580.804|2.069|
|baseline|153|611.728|610.933|1.229|
|g0.25_h7|27|477.747|477.229|1.039|
|g0.25_h7|55|539.138|536.027|1.706|
|g0.25_h7|118|580.387|580.804|1.710|
|g0.25_h7|153|611.460|610.933|1.202|
|g1_h7|27|477.510|477.229|2.034|
|g1_h7|55|540.318|536.027|3.873|
|g1_h7|118|579.853|580.804|2.547|
|g1_h7|153|609.307|610.933|2.282|

## Evidence and reproduction

Plan: [POSTCAP_VOV_PLAN.md](POSTCAP_VOV_PLAN.md). Plots: [postcap-vov017](visualizations/postcap-vov017/README.md). Run: `runs/m2_postcap_vov_001`. `environment.txt` records the temporary Python3.12 package versions. Reproduce with `PYTHONPATH=matteo/exogenous_reservoir_eq_m2/src` and modules `postcap_vov_study --output NEW_DIRECTORY`, `audit_postcap_vov --run RUN`, `plot_postcap_vov --run RUN --output PLOTS`, `summarize_postcap_vov --run RUN`. All original baseline files are preserved.

Study wall time: 18.39 minutes. No model was automatically adopted; this comparison is a finite grid with the old parameters fixed, not a full enlarged-model calibration.
