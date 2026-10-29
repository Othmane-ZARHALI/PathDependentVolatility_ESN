# Protected nine-level SPX fit: results and durable conclusion

28 October 2026. The user authorized this experiment and explicitly asked to save its conclusion so it is not rediscovered in a later conversation. The complete run is `runs/m2_protected_curve_001`. See [the experiment plan](PROTECTED_CURVE_PLAN.md).

**The nine-level fit improves the specified SPX price objective without sacrificing VIX, but does not solve SPX smile shape. It is not a satisfactory joint calibration.** On fresh 128-step checking, price loss falls by 34.9%, while equal-maturity average native implied-volatility error increases from 4.345 to 4.710 percentage points. Both facts must be retained when describing the result.

## Exact experimental scope

The existing two-speed variance model remains

$$v_t=f+(\xi_t-f)[(1-w)R_t+w\{(1-r)A_t^s+rA_t^f\}].$$

The original M2 dynamics, readout, normalizer, positive spike shift, risk premia, full initial state, both independent stress factors, stress weights and sensitivities, half-lives and floor are fixed at the successful VIX-only values. In particular, w=0.922838 and the volatility floor remains 12%. Only xi indices 0,5,6,9,11,12,13,14,15 are free, with the original floor-to 0.2 variance bounds. All seven VIX-sensitive xi levels stay fixed. No smoothing, new factor, correlation, jump or maturity-specific rule was added.

All 1294 primary SPX quotes at 29,57,119,211,302,666,1030 days enter the unchanged native bid–ask-scaled interval objective, with equal maturity weights. All 1751 eligible quotes are checked and retained in CSVs. All 217 VIX native quotes and four forwards remain protected. Signed estimates and undefined IVs are retained; there is no clipping or smoothing.

## Why VIX stays unchanged

VIX at expiry T depends on conditional variance over[T,T+30 days]. The free curve intervals lie outside all four VIX windows. They do not enter the exact VIX parameter map, and none of the exogenous variance-factor dynamics depends on the stock or on xi. Therefore these changes leave the VIX distribution itself unchanged. This is stronger than merely preserving four forward means.

The audit checks exact array equality of baseline/candidate VIX inputs and every other fixed coordinate. `VIX_protection/invariance.json` records this and hashes inherited original price/forward artifacts. The original checking VIX option loss remains 0.050081, with 153/217 quotes inside bid–ask. The previously measured checking forward errors remain as previously measured: no new population-forward equality certificate is claimed, and inherited VIX CSVs are not presented as new Monte Carlo estimates.

## Fitting and genuinely fresh checking

Fitting uses 512 existing full-M2 search paths and 4096 cheap reference paths, with the terminal-Black-delta auxiliary estimator from the precision study and analytic derivatives. This is an unbiased estimator of the existing discrete model, with sampling and numerical errors assessed separately. It is cheaper than the one-step conditional-value estimator; they must not be conflated. All nine native-price derivative directions pass finite differences. Two bounded starts were tested, each with at most 40 least-squares evaluations. Both terminated by ftol. Training loss changes from 6406.258 to 3665.656. This is a bounded local search, not a global fit or uniqueness result.

Selection was frozen before checking. Fresh checking generates 512 new target-Q paths (256 antithetic pairs in eight batches), including new M2, equity and both stress innovations. Each path carries both candidate models and all seven maturities. The stock 64/128 grids are coupled through the same exact fine OU/stress/equity paths. All curve knots and maturity events are inserted exactly. Independent reference checking uses eight new scrambled OU-bridge nets of 4096 points per expiry, shared between candidates. The old expensive checking cache is not reused. This is independent simulation checking at the same market snapshot, not market out-of-sample testing.

The 12 orthogonal OU modes have zero C/Y/F readout projection and no feedback into retained dynamics, so they are analytically marginalized. The exact retained 24-mode joint OU/equity Gaussian law is preserved and audited, including its Q drift and covariance. This is not a reduced three-projection Markov approximation. State paths are streamed in batches, not stored as an expanded full fine-state cache.

## Price loss versus implied-volatility error

At 64 steps/tradingday, fresh price loss changes 5964.131→3895.393. At 128 steps it changes 5979.323→3895.071. The paired 128-step loss difference is-2084.251; its approximate paired batch-bootstrap 95%interval is[-2495.423,-1690.923]. It excludes zero, supporting improvement of this price metric beyond the observed sampling noise. It excludes fitting variability, normalizer error and stock-grid bias; it is not a general price-error bound.

| Days | Baseline price loss | Protected price loss | Baseline IV MAE (pp) | Protected IV MAE (pp) |
|---:|---:|---:|---:|---:|
| 29 | 5488.66 | 2109.60 | 4.763 | 4.947 |
| 57 | 6830.19 | 3632.75 | 4.301 | 4.211 |
| 119 | 6805.17 | 4818.32 | 3.729 | 4.068 |
| 211 | 10400.56 | 6999.18 | 4.246 | 4.848 |
| 302 | 10674.62 | 8221.86 | 4.533 | 5.256 |
| 666 | 1240.14 | 1135.16 | 4.444 | 4.898 |
| 1030 | 415.93 | 348.62 | 4.397 | 4.742 |


The price objective and IV MAE measure different things. An option's price sensitivity to implied volatility varies with strike and maturity. Reducing the chosen bid–ask-scaled price loss can therefore worsen an average IV error. The fit mostly lowers the SPX curves: it helps the price objective, but the left wing becomes too low and the right wing remains too high relative to the market. The curves retain a much more U-shaped profile than the observed downward skew.

At 128 steps only 8/1294 primary quotes are inside bid–ask, versus 11 for the baseline. All fresh checking native prices are nonnegative in this sample; signed estimates were retained. Median candidate sampling SE is 2.142 price-scale units. Thus this is clearly not a final quote-level calibration despite the statistically supported loss improvement.

[All seven SPX smiles](visualizations/protected-curve028/SPX-protected-smiles.png), [the four inherited VIX smiles](visualizations/protected-curve028/VIX-protected-smiles.png), and [the variance-curve changes](visualizations/protected-curve028/protected-variance-curve.png). Model IVs use the existing percentage-valued inversion routine, and invalid values are flagged in the all-quote CSVs.

## What changed in the curve

| Index | Days | Baseline xi | Protected xi |
|---:|---:|---:|---:|
| 0 | 0–27 | 0.04898962 | 0.02316962 |
| 5 | 85–91.25 | 0.04898962 | 0.01454153 |
| 6 | 91.25–118 | 0.05827299 | 0.01446592 |
| 9 | 148–153 | 0.05827299 | 0.01440000 |
| 11 | 183–211 | 0.05827299 | 0.01440000 |
| 12 | 211–273.75 | 0.05827299 | 0.04831126 |
| 13 | 273.75–302 | 0.05257415 | 0.01448961 |
| 14 | 302–666 | 0.05257415 | 0.04670744 |
| 15 | 666–1030 | 0.05257415 | 0.04766661 |


Several short or intermediate off-window levels approach the frozen floor. Those boundary/near-boundary values are numerical outcomes of this restricted fit, not uniquely identified economic estimates. The original model and this exploratory candidate remain separate.

## Numerical checks and limits

The implementation audit passes: all 9 price derivative directions; independent scalar/NumPy reconstruction of the fresh stock kernel (maximum error 1.11 e-16); all 8 new batch-seed replays; the complete retained Gaussian covariance recurrence over 91040 fine steps (maximum absolute covariance error 2.69 e-12); Q projected means; source fingerprints; input-manifest checks; exact protected-parameter invariance; and positive spike shift. Normalization still uses the frozen GH48/interpolated normalizer, with its earlier sampled interpolation audits. No new global GH-convergence certificate is supplied.

The broad improvement survives the 64→128 grid change. Nevertheless the candidate's median absolute paired price change is 0.220 price-scale units, maximum 2.438; finer-grid convergence is not established. Sampling SEs use antithetic pairs and independent reference-net means. Loss intervals resample the eight expensive and eight cheap batch means separately, with common resamples across candidates. Eight batches give an approximate uncertainty assessment, not a guaranteed tail bound.

## Durable conclusion: do not repeat this experiment as if it were untested

1. Freeing off-window forward-variance levels can improve the original SPX price objective while preserving the successful VIX distribution.
2. In these two bounded starts it does not produce the required SPX downward skew. Price loss improves while average native-IV error worsens. Do not summarize the experiment only by the lower price loss.
3. Further curve-only iterations should have a new concrete hypothesis or new evidence. Repeating the same unrestricted or nine-level optimization is not the recommended next step.
4. A plausible explanation is that the dominant stress variance is independent of the stock shock. It supplies volatility variation without directly linking stress increases to stock falls. This is a mechanism hypothesis, not an established theorem or proof that the entire fixed-factor family is impossible to calibrate.
5. The next proposed structural test is one global stock/stress leverage parameter, with no added stress speed. A constructive driver rotation can preserve the joint variance-factor law and hence VIX; see [STOCK_DRIVER_ROTATION_PROPOSAL.md](STOCK_DRIVER_ROTATION_PROPOSAL.md). That proposal is not implemented here, changes the old stock/M2 correlations by construction, and requires new joint-transition and pricing-control checks.

The original VIX-only artifact and production baseline are unchanged. The candidate is `experimental_model_protected.json`, explicitly not adopted. No new Lean theorem or global impossibility claim is made. Current continuation notes link directly to this result so a future conversation can resume without repeating it.
