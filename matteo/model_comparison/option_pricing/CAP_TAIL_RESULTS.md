# Cap and tail investigation: bounded experiments

Completed 25 October 2026. The selected baseline remains `runs/m2_maturity_curve_fit_001/selected_model.json`. No new joint calibration is adopted. New evidence is isolated in `runs/m2_cap_tail_001/`. Original runs, selected configuration and production kernels are unchanged.

The cap compresses payoff-relevant extreme-state responses and changing it steepens the wing. Its finite support is not the binding obstruction at the quoted strikes. Neither raising the transition nor the tested continuous alternatives solves the smile. New VIX-focused curve fits expose additional flexibility with severe SPX tradeoffs and weakly identified variance pieces.

## Scope and correspondence

All seven SPX dates 29/57/119/211/302/666/1030 and all four VIX dates 27/55/118/153,1511 original OTM quotes and four forwards, cutoff 3, common-call payoff and all 4096 samples at every strike. Market puts use market-forward parity; native model puts use model forwards. OU rates/geometric grids, initial state, correlations and risk-price scenario are fixed. GH32/GH48/time8 and stock 512 remain baseline approximations. Experimental kernels and additional quadrature tables are separate. No jumps or quintic approach was explored.

The new baseline VIX repricing differs from saved selected prices by at most 1.36e-09 VIX option-price points. Independent saved-draw replay verifies the parent selection before using it. Equations, integrability and martingality arguments: [CAP_TAIL_MATHEMATICS.md](CAP_TAIL_MATHEMATICS.md). Designs/literature/cache costs: [STOCHASTIC_VOL_OF_VOL_RESEARCH.md](STOCHASTIC_VOL_OF_VOL_RESEARCH.md).

## Range and weighted tails

Numerically evaluated cap upper bounds are 936.52/949.28/976.28/1003.19, versus highest quoted strikes 100/170/200/200. Those bounds alone do not establish attainability. Saved separate-checking VIX maxima425.17/431.81/543.95/508.56 independently show implemented conditional values exceeding the quoted strikes.

|Days|Highest strike|Proposal hits /4096|Weighted Q probability|Raw call ± scramble SE|Payoff ESS|Top10 payoff share|
|---:|---:|---:|---:|---:|---:|---:|
|27|100|234|0.000423256|0.0092771 ± 0.00036|165.3|10.4%|
|55|170|143|6.96538e-05|0.00193892 ± 0.00012|90.5|19.9%|
|118|200|126|3.88152e-05|0.00125814 ± 8.3e-05|84.6|19.1%|
|153|200|110|3.68207e-05|0.00117303 ± 7.5e-05|64.7|27.2%|

Hits are proposal frequencies, not risk-neutral frequencies. Probability uses mean[w*indicator]; raw call uses mean[w*(VIX-K)+], without self normalization. Payoff ESS and concentration are descriptive. Actual controlled prices retain signed zero-mean corrections. Eight independent scramble means provide sampling SE, excluding deterministic bias, frozen pilot and parameter uncertainty.

At the four window midpoints, only roughly 0.0066–0.0086% of direct-Q probability lies above eta=L, but those states contribute 8.8–10.4% of E[G]. G-weighted mean response h'(eta) is 0.86–0.88 overall. Highest-strike-call-weighted conditional responses at the window midpoint are 0.662/0.597/0.574/0.577. Thus saturation is relevant to wing contributors even though their VIX values need not approach the global upper bound. Marginal diagnostics use eight NEW scrambled 32768-point Gaussian-mixture importance nets per date, with saved seeds and exact p/q weights.

## Fixed-curve cap profiles

Only the tempering function changes here; each candidate rebuilds its own Q normalizer. Curve/floor and original score coefficients stay fixed. These are finite profiles, not optimized families.

|Candidate|Native VIX options|Forwards|Native VIX + forwards|Unavailable call IVs|VIX evaluation seconds|
|---|---:|---:|---:|---:|---:|
|baseline|26.276301|0.890564|27.166865|0|8.23|
|lower_L|37.036268|3.488472|40.524740|0|10.59|
|higher_L|34.462083|0.078745|34.540828|0|9.83|
|soft_k1|34.173899|2.984718|37.158617|0|9.33|
|sharp_k4|27.915860|0.398971|28.314830|0|9.66|
|algebraic|25.937184|0.895556|26.832740|0|9.35|
|logtail025|26.643263|0.509752|27.153015|0|11.13|
|logtail050|28.738746|0.202424|28.941170|0|13.09|
|separate|64.148272|7.150778|71.299051|14|41.43|

The logarithmic continuation is h_original+a*log(1+softplus(k*(eta-L))/k), so G grows only polynomially in the Gaussian score beyond the transition. All its Gaussian multiplier moments exist. Algebraic saturation matches the original transition value/slope and changes its approach to the plateau. Separate-channel caps keep the linear C channel outside individual Y/F caps; the tested 3/2.5 allocation is poor. Their equations are recorded in the mathematics report.

Raising L strengthens tails AND increases M_Q. It therefore reduces normalized variance in many ordinary states. At fixed xi this changes VIX forward and near-ATM fit as well as wings. The L+.75 profile has VIX forwards 19.870/20.632/21.267/21.822 versus baseline 20.227/21.094/21.782/22.346 and market 20.10/20.73/21.19/21.72. The resulting forward improvement does not prevent option-price loss from worsening.

|Days|Market right-wing IV slope|Baseline|L+.75|Log a=.25|Log a=.50|
|---:|---:|---:|---:|---:|---:|
|27|52.38|35.40|42.84|37.95|40.54|
|55|38.31|20.29|26.58|22.61|24.85|
|118|31.23|12.32|17.23|14.18|15.95|
|153|27.32|10.54|14.92|12.19|13.77|

Slope is a descriptive least-squares slope of raw IV against log(K/F), restricted to K/F>2, in percentage points per unit log strike. It is not an asymptotic slope, smoothed price curve or new calibration objective. Curvature, ATM level/MAE and available-IV counts are saved in `visualizations/cap-tail-014/smile-shape.csv`. All variants leave material shape deficits.

At 153 days and K=200, actual training call price is 0.001275 at baseline,0.006352 at L+.75 and0.002362 for log a=.25; market midpoint is 0.141629. Their pointwise sampling SEs are 0.000094/0.000614/0.000161. Tail strengthening is substantial in relative terms but much smaller than the price gap. Small SE is not a guarantee against unsampled tails or deterministic bias.

## VIX-focused fits and actual SPX tradeoffs

Training profiles selected L+.75 and log a=.25 for follow-up before separate checking. Two bounded least-squares fits free 16 curve levels and floor, holding their readout/cap fixed, with no joint safeguards. The fit parametrization enforces xi>=f and xi<=.2; floor volatility remains5–6% for comparability. The 6% endpoint and .2 curve upper bound are implementation choices. Both optimizers report local ftol termination; neither is a global convergence result.

|Candidate|Training native VIX + forward|Training common-call VIX + forward|Separate reused native VIX checking|Training SPX|Training native joint|
|---|---:|---:|---:|---:|---:|
|Selected baseline|27.166865|45.509602|27.523027|43.777476|70.977709|
|higher_L_fit|20.044859|48.202084|20.653505|1656.376316|1676.460018|
|logtail025_fit|18.674677|50.890919|19.124019|318.635148|337.349445|

These are diagnostic VIX-only fits, NOT accepted joint calibrations. Both also worsen the training common-call VIX-plus-forward criterion, despite improving native OTM VIX; the dual-convention safeguard is materially relevant. Both reach the 6% floor endpoint and the xi=.2 bound near day 118. VIX-only data do not identify curve pieces outside their30-day windows; the weak roughness penalty can move those pieces substantially. The higher-L fit, for example, takes the 0–27-day variance from 0.04899 to 0.00367 and damages short SPX. Consequently these SPX tradeoffs mix cap changes, fitted-window changes and poorly identified off-window curve changes; they do not isolate the cap cost or prove joint-fit guards are overly restrictive. The useful finding is that removing joint guards permits further sampled VIX improvement, while joint information is essential to control the curve. A next VIX-only diagnostic should freeze off-window pieces and floor before interpreting SPX costs.

- higher_L_fit: separate-checking native VIX-plus-forward change-6.869522, descriptive paired eight-net jackknife SE 0.342465. Training SPX loss change 1612.599, paired eight-block SE 65.018.
- logtail025_fit: separate-checking native VIX-plus-forward change-8.399008, descriptive paired eight-net jackknife SE 0.194685. Training SPX loss change 274.858, paired eight-block SE 52.174.

Checking streams were used previously and score the SAME market quotes. They are not untouched holdout or market out-of-sample data. SPX checking on separate streams was not repeated. Signed prices, unavailable IVs, shape flags and raw stock-forward diagnostics remain in saved scores/arrays. In particular baseline training day 1030 raw stock-forward ratio is about0.9624 with iid SE 0.00927: no new stock sampling clearance follows.

## Stochastic vol-of-vol follow-up

After cap profiles, two bounded fixed-curve experiments let existing cluster C change the mixture weight between original and stronger-spike responses. Original OU/randomness remain unchanged. Native VIX-plus-forward losses worsen to 87.177224 and117.866292 (spike boost.5 and1). These tests also bend the far wing downward at several dates. No cluster candidate was promoted to a joint fit or a full SPX claim. This rejects those profiles as immediate improvements, not every cluster mechanism.

The linked stochastic-vol-of-vol report explains Ito response strength, the implemented bounded mixture, a dedicated continuous latent OU response state, an alternative latent mixture, identification, boundedness/martingality, randomness and cache reuse. A dedicated state is explicitly proposed rather than implemented. It would require new state simulation, four-dimensional conditional pricing, and full volatility/stock recomputation; one stored float64 state on the current stock cache costs 26.05GiB.

## Numerical checks and unresolved issues

- New readout/stock kernels match independent Python equations and the adapted stock recursion. Tests were run directly because this environment has no pytest module; the initial pytest-launch failure is retained.
- Same-GH direct conditional probes give maximum interpolation relative errors2.23e-10 for cap profiles and1.85e-8 for cluster profiles; independent randomized tensor checks also pass. These are sampled interpolation checks, not global bounds.
- At representative outer states (first, largest VIX and largest weighted highest-strike payoff), all four maturities and baseline/higherL/logtail025 profiles: conditional GH32→48 maximum VIX change 2.26e-06; time order8→16 maximum3.39e-08. Normalizer GH48/64/96 midpoint values are also saved. This does not certify the whole surface or the fitted curves.
- Baseline separate-checking153-day second moment remains606.20794 versus610.93257 with net SE 1.63548. HigherL fixed-curve checking gives605.78268 with SE 1.68207 against the SAME target. Raising the cap does not remove that deficit. Direct marginal moment agreement and the small conditional probes point away from those tested integration errors, while outer-tail/proposal errors remain open.
- Removing the cap is not tested as a solution: selected uncapped first moments are finite, but long-date second multiplier moments are infinite (critical moment power 1.11480). Higher curvatures can also destroy first moments. The continuous stock martingale is not certified merely by a finite normalizer. Bounded and controlled-growth alternatives have separate mathematical arguments in the equation report.
- Independent audit replays every saved price/SE/parity/loss, exact4096 allocation, candidate constraint and captured source snapshot. It does not regenerate original primitives or prove the numerical model equivalent to the continuous SDE. No new Lean or P-compatibility certificate is asserted.

## Decision and continuation

Preserve the selected baseline. The cap is a real finite-state compression mechanism and an adjustable wing-shape lever, but this bounded study does not establish it as the sole cause or exhibit a satisfactory joint remedy. Logarithmic continuation remains a mathematically controlled experimental family; it is not selected after these tests. Avoid an expensive 24-parameter recalibration before a more focused shape design. Next isolate the response-versus-level tradeoff with off-window curve pieces fixed, then compare cluster modulation against a dedicated bounded latent response. Keep jumps and quintic approaches excluded. Do not relax floor/beta constraints or change the existing reservoir.

Plots: [cap location/sharpness](visualizations/cap-tail-014/vix-cap-location.png), [alternative cap shapes](visualizations/cap-tail-014/vix-cap-shapes.png), [weighted tail calls](visualizations/cap-tail-014/weighted-tail-calls.png), [VIX-only fits](visualizations/cap-tail-014/vix-focused-fits.png), [all seven SPX maturities](visualizations/cap-tail-014/spx-all-maturities-iv.png), [separate VIX checking](visualizations/cap-tail-014/vix-separate-checking.png). All have SVG exports; signed-price/error and forward plots are retained.

Reproduction modules: `cap_tail_study`, `cap_tail_diagnostics`, `cap_tail_followup`, `cluster_vov_study`, `cap_tail_refinement`, `audit_cap_tail`, `plot_cap_tail`, `summarize_cap_tail`. Run with `OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=matteo/exogenous_reservoir_eq_m2/src /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.<module>`. Initial driver takes `--output NEW_DIRECTORY`; followup takes `--run`. Auxiliary postprocessors currently name001 explicitly. Do not rerun mutation drivers over existing directories; their mkdir checks protect originals. Full command/config/source/run records are saved.

