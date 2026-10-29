# Stock exposure to future integrated variance — 28 October 2026

Implemented the user-authorized follow-up to TARGETED_SKEW_RESULTS.md. All variance/factor parameters stay fixed; only the stock shock direction changes. No new bank, endogenous factor feedback, readout fitting or constant-theta search. All seven SPX maturities/1751eligible/1294primary quotes and217VIX quotes/four forwards remain. This is an exploratory fixed-model comparison, not an adopted double calibration.

## What the stock now follows

For the same30day horizon L=30/365 at every time/maturity, define

$$F(t,X_t,Z_t)=E_Q[\int_t^{t+L}v_u du\mid\mathcal F_t],\qquad dF_t=\text{drift}dt+\mathbf K_t\cdot d\mathbf B_t.$$

The new stock direction is

$$dS_t/S_t=\sqrt{v_t}\mathbf a_t\cdot d\mathbf B_t,\qquad \mathbf a_t=-\mathbf K_t/\|\mathbf K_t\|.$$

At zero K, use the old mixed theta=.8 loading. The state derivative of F uses all24retained Gaussian modes and their individual decay speeds, plus both stresses, in13independent primitive Brownian coordinates. It is not a three-projection Markov approximation and does not use realized future shocks. Compare the fixed mixed stock loading, previous instantaneous-variance loading, future-variance loading updated4and8times/trading-day. Directions are predictable and unit length between updates. The purpose is to test persistent leverage, not to prove that this choice maximizes finite-wing implied skew.

The original M2 cap/readout, positive spike shift, stress rates/amplitudes/weights, floor,16-level curve, initial state and Q law remain fixed. The entire conditional variance-factor law is unchanged, so VIX options/forwards are identical to the original model. Existing numerical forward gaps (up to.0552points) are inherited; this experiment does not newly establish exact population forward matching. Stock martingality in continuous time is separate from unit Brownian loading.

The last variance level is continued flat beyond day1030 solely to define the common30day signal through day1060. No observed-horizon variance value changes and no new level is fitted. The existing normalizer covers this lookahead.

## Calculation and validation

For each future lag, project the complete centered24mode state with its exponential decays and add the deterministic Q mean. Integrate the independent cluster C dimension of the capped readout derivative using48point Gaussian convolution and an8193knot cubic-Hermite table; integrate feedback/spike Y/F with12point tensor Gaussian quadrature. Stress sensitivities are analytic. Sum the sensitivities on0plus22geometric lag knots from1e-5years to30/365. Piecewise-linear hat weights integrate xi(t+lag)-floor exactly across every deterministic curve break. Time means, normalizers and decayed diffusion matrices are precomputed, not fitted.

Independent probes compare12/24/32Gaussian nodes,22/44lag knots and direct finite differences of independently repriced F in every primitive direction. At60declared Q/tail-amplified probe states, max unit-direction distance GH12vs32 is 0.000514;22vs44lag knots gives 0.0127. These meet the declared exploratory.03direction-distance gate, not a uniform error bound. Full primitive finite-difference sensitivity error is 9.34e-12; exact curve hat-weight integral error is 4.34e-18. The signal's Gaussian old-bank gradient is obtained by differentiating the conditional expectation, not by evaluating the raw readout derivative at the conditional mean.

Missing primitive Brownian increments use the previously audited conditional OU bridge. Frozen baseline/instantaneous full-horizon replay matches earlier stocks to 0 on four declared paths, with original counter stride128. An independent97step scalar recurrence checks both update schedules, both stock grids and antithetic partners to 2.22e-16; unit forward direction error is 2.22e-16. No new Lean or compiled-simulator theorem is claimed.

## Fresh paired results

1024new full-model stock paths,512antithetic pairs/eightindependent batches, coupled128/256stock grids. The4/8update schedules share the identical factors and variance. All choices were fixed before pricing; there is no market-selected winner. The table uses256steps/trading-day; IV MAE equally averages available native IV errors across seven maturities.

| Stock exposure | SPX interval loss | IV MAE,pp | Invalid primary IVs | Median MC SE / quote scale |
|---|---:|---:|---:|---:|
| Fixed exposure | 1909.19 | 3.201 | 1 | 2.07 |
| Instantaneous variance | 1978.79 | 3.163 | 1 | 6.04 |
| Future30days,4updates | 1969.78 | 3.175 | 2 | 6.12 |
| Future30days,8updates | 1941.31 | 3.157 | 3 | 6.10 |

Future30day/8update point-estimate price loss relative change versus fixed exposure: +1.68%(positive means worse). Approximate paired1024resample batch-bootstrap95%loss-change intervals: [{"candidate": 1, "name": "Instantaneous variance", "price_loss_delta_CI95": [-572.5305487200371, 981.3045020440446]}, {"candidate": 2, "name": "Future30days,4updates", "price_loss_delta_CI95": [-556.7310013033479, 810.089215024906]}, {"candidate": 3, "name": "Future30days,8updates", "price_loss_delta_CI95": [-577.7064309209505, 803.3020717036825]}]. Eight batches underpin these sampling intervals; deterministic integration, stock-grid, normalizer and market errors are excluded. Invalid IVs/negative signed price estimates are retained. A lower available-IV error with omitted failed IVs is not a final calibration certificate.

At29days, quoted K/F near.9,1.1 gives IV(right)-IV(left) -6.39pp for fixed exposure, -3.72pp for instantaneous exposure, -4.51pp for future30day/8updates, versus market -14.05pp. All seven maturity panels and per-maturity metrics are saved; assess the slope and both wings, not only the aggregate price loss.

Update-frequency diagnostic4to8: median absolute price change 0.253quote scales, maximum 1.994; median paired SE of the change 0.499scales. Stock128/256grid diagnostics: [{"candidate": 0, "median_grid_change_scales": 0.1493100061200234, "max_grid_change_scales": 1.1570308381595023}, {"candidate": 1, "median_grid_change_scales": 0.3968655236160672, "max_grid_change_scales": 3.8755657872563862}, {"candidate": 2, "median_grid_change_scales": 0.15724673526070293, "max_grid_change_scales": 1.1575141378736822}, {"candidate": 3, "median_grid_change_scales": 0.14413141535652146, "max_grid_change_scales": 1.1931255907434202}]. Nonzero changes remain numerical approximation error, not a separate calibrated model benefit.

## Estimator and limits

Use the same65536path stress-only auxiliary reference and its correct conditional lognormal call prices; that auxiliary stock model has not acquired the new coupling. The expensive correction is target/reference payoff difference minus a deterministic delta times stock difference. Include independent reference-batch uncertainty. No future-stress-dependent delta, sample-forward renormalization, clipping or price smoothing. The reference is reused from earlier research, so the complete estimator is not newly independent even though the1024factor paths are new. This is the same market snapshot, not market out-of-sample validation. Raw downside/upside payoff concentration, stock means/SEs, quote membership, prices, invalid IVs and paired differences remain in the run.

VIX prices/forwards are inherited unchanged, not newly re-estimated here. The unchanged factor/readout law is the invariance argument. No new claim about precise population-forward matching, continuous-stock martingality, optimum leverage direction for option prices or impossibility of fitting the model follows.

## Assessment and memory diagnostic

This particular30day exposure rule does not materially improve the SPX fit. On these paired draws, the8update price loss changes 1909.19to1941.31(+1.68%), while native-IV MAE changes only3.201to3.157pp. Its29day quoted IV gap becomes less negative, -6.39to-4.51pp versus market-14.05pp. All paired price-loss sampling intervals include zero. There is no demonstrated improvement to adopt, and the missing left-wing slope remains visible.

Together with the earlier constant-rotation and instantaneous-gradient tests, this gives no empirical reason to keep adjusting stock-shock direction as the next priority. It does not prove that every possible state-dependent correlation or lookahead kernel fails. The tests have held the variance readout, bank projections and stress distributions fixed.

A separate, exact finite-dimensional Gaussian diagnostic exposes a plausible readout-memory limitation. For the existing feedback projection Y, stationary covariance gives autocorrelation0.170after1day,0.048after7days and0.017after30days. This is a model stationary-memory calculation, not a dated option calibration or historical correlation estimate. Its projection weights are fixed and strongly favour fast modes; the slow OU modes already exist in the bank.

A concrete next hypothesis is to give a slower projection of those existing modes an independent role in the asymmetric readout, while retaining the current fast component and the Zumbach-motivated design. Gaussian state paths remain reusable; new readout normalization, stock paths and VIX conditional expectations would need rebuilding. This readout-memory hypothesis is untested and must not be presented as the diagnosed sole cause of the skew gap. It is more distinct from the failed exposure scans than another constant-theta/global-multiplier variation. No new bank, readout or model is adopted in this experiment.

## What this establishes

The implemented direction targets expected future integrated variance with verified full-state sensitivities and retains the successful VIX law. The response and limits must be read from the paired prices, all-maturity smiles and uncertainty above. This is a diagnostic of one common30day horizon and the existing variance model, not a broad search over horizon functions or readouts. Do not revert to another constant-theta/global-multiplier scan or treat a favorable aggregate metric as a satisfactory joint calibration. Any extension should state which missing maturity/wing mechanism it addresses and keep factor reuse and VIX protection explicit.

Figures: [all seven SPX smiles](visualizations/forward-leverage028/SPX-forward-leverage.png), [unchanged VIX smiles](visualizations/forward-leverage028/VIX-unchanged.png). Detailed run: `runs/m2_forward_leverage_001`. Stock simulation/pricing wall time 163.2seconds, excluding preparation and audits.
