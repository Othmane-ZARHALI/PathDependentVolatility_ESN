# Asymmetric whole-response multiplier: results, 28 October 2026

**Tested, not adopted.** The proposed readout gives modest SPX improvements at stronger settings, but damages VIX smiles and raises stock-pricing uncertainty. The VIX-preserving selection retains alpha=0. This is evidence about this fixed functional shape, this existing factor and this bounded seven-level refit—not a model-family impossibility result.

## Equations and exact scope

Keep all exogenous OU dynamics and the previous mixed stock rotation theta=.8. Use the existing feedback projection

$$X_t=(Y_t-E_Q[Y_t])/s_Y,\qquad s_Y=1.$$

Its instantaneous covariance with the old equity driver is positive; the stock still loads .6 on that driver. Define

$$D_\alpha(X)=\exp[\alpha\log(1+e^{-X})],\qquad H^\alpha_t=\frac{H^0_tD_\alpha(X_t)}{E_Q[H^0_tD_\alpha(X_t)]},\qquad v_t=f+(\xi_t-f)H^\alpha_t.$$

All original M2 parameters, stress rates/amplitudes/weights, floor and nine off-window variance levels are fixed. Alpha=0,.25,.5,1,1.5,2 is shared across every maturity. For each strength, fit only the seven VIX-sensitive variance levels, bounded between the floor and .2, with all four training-estimator forward equalities. Include the original-curve baseline and an alpha=0 curve-refit control. Positive spike shift remains fixed. No endogenous feedback, extra bank, jump or quintic readout is introduced.

## Comparison

The table uses 1024 fresh full-model stock paths on the 128-step/trading-day grid and 4096 fresh full-factor VIX draws per expiry. VIX checking freezes training-fitted levels; it does not force checking forwards back onto target. All 1751 eligible/1294 primary SPX quotes at seven maturities and all 217 VIX quotes/four forwards remain present.

| Alpha | SPX IV MAE, pp | SPX interval loss | VIX interval loss | Invalid primary SPX IVs |
|---:|---:|---:|---:|---:|
| 0 | 3.278 | 1992.1 | 0.047 | 0 |
| 0.25 | 3.227 | 2033.8 | 0.397 | 0 |
| 0.5 | 3.141 | 2038.0 | 2.000 | 0 |
| 1 | 2.956 | 1923.6 | 11.090 | 3 |
| 1.5 | 2.832 | 1896.3 | 25.725 | 6 |
| 2 | 2.931 | 2138.9 | 38.539 | 8 |

Training forward gaps are below 8.45e-10 VIX points for every strength. That is equality of a controlled training estimator, not an exact population-forward guarantee. Fresh forward gaps and both anchor/difference SEs are in `runs/m2_asymmetric_readout_001/fresh_VIX/VIX_forwards_c*.csv`. Signed negative price estimates and unavailable IVs are retained. IV MAE averages available IVs; the missing counts matter. The original unchanged-curve alpha0 baseline on these new stock draws has loss 1988.93 and IV MAE 3.277 pp, so comparisons should use this experiment's paired alpha0 controls, not a different experiment's Monte Carlo estimate.

The predeclared selection allowed VIX loss no more than alpha0-refit loss+.1 before minimizing replay SPX IV MAE. No positive tested strength met that VIX screen; candidate1, alpha0, was selected before fresh stocks. The fresh full-grid curves are exploratory response evidence, not independent validation of a newly selected winner. All bounded VIX fits completed with forward feasibility; several curve levels reach their bounds. Wider readout/stress/curve optimization was not attempted.

## Why the normalization matters

Since the stress factors are independent of the old bank,

$$E[H^0D]=(1-w)E[RD]+wE[D].$$

The original R already favours negative feedback states. D favours those states too. Thus E[RD] grows faster than E[D], and whole-product normalization changes the average M2/stress balance despite w staying fixed. At day153 the stress fraction of mean excess variance is 92.3% at alpha0, 85.6% at alpha.5, 72.9% at alpha1, 53.2% at alpha1.5 and 32.1% at alpha2. This mechanism is directly computed under the Gaussian law; attributing the entire VIX loss to it would require an additional controlled experiment.

The volatility reaction is stronger on bad-news states, but the normalization also shifts mean variance away from the stress mechanism that produced the successful VIX fit. Merely matching the four forwards does not recover the VIX distribution or smiles.

## Numerical evidence and limits

The stock sampler rebuilds candidate variances and stock paths from the same factor stream once per batch for every strength. Alpha0 full replay matches previous rotated stocks to 5.6e-14. Independent NumPy direct tensor GH versus the C-convolution/YF implementation agrees to 3.68e-11 scaled error; a scalar stock recursion agrees to 1.11e-16; all seven forward derivatives agree within 1.4e-09. Normalizer GH48 versus GH64 discrepancy is at most 3.08e-08 at declared probes. Midpoint time interpolation error is about1.1e-7. These are sampled numerical checks, not rigorous global bounds.

Stock uncertainty uses antithetic pairs and eight independent batches; auxiliary reference uses the previous65536 paths and is shared across screening/checking. Median stock SE rises from about1.66 quote scales at alpha0 to7.88 at alpha1.5; strong-readout numerical precision is not sufficient for a final calibration claim. Coupled64/128 grid changes, paired bootstrap loss intervals, stock means and invalid prices are saved. Eight batches underpin approximate uncertainty, not a certified tail bound.

VIX training uses1024 randomized rows from the earlier65536 state/stress stream and a controlled old-model anchor. Follow-up VIX checking uses genuinely fresh joint38-dimensional Gaussian-mixture draws (eight independent scrambled Sobol batches), with exact likelihood weights and no sample-mass normalization; all strengths share those draws. The old anchor price estimate is still shared, and its uncertainty is retained separately and combined with the independent new correction. Weighted-tail ESS/concentration are saved for every strength and maturity. GH32/time8 checking and matched-state GH48/time12 probes diagnose deterministic VIX error. The same market snapshot is used throughout.

Existing Lean algebraic identities were not changed. The new quadrature/sampler correspondences are numerical audits; neither continuous stock martingality nor a finite-wing-skew response theorem is claimed. No production baseline was overwritten.

## Separating readout and curve effects

The 29-day improvement in the combined fit must not be attributed to the multiplier alone. A four-way alpha1.5 ablation replays the exact same fresh stock Gaussian streams and uses the same fresh VIX states:

| Variant | SPX IV MAE, pp | 29-day IV(right)-IV(left), pp | VIX loss | Invalid primary SPX IVs |
|---|---:|---:|---:|---:|
| original | 3.277 | -5.86 | 0.050 | 0 |
| multiplier only | 2.940 | -6.05 | 21.278 | 12 |
| curve only | 3.528 | -8.51 | 19.400 | 0 |
| multiplier and curve | 2.832 | -13.44 | 25.725 | 6 |

The 29-day gap compares actual quotes nearest K/F=.9 and1.1; its market value is -14.05pp. At fixed curve, the multiplier changes this gap only from -5.86 to -6.05pp. Curve adjustment alone gives -8.51pp, while both give -13.44pp. The combined effect is nonlinear: neither isolated modification reproduces it. The first VIX-sensitive interval spans days27–29, so its refit also directly affects the 29-day SPX expiry; the interpretation must acknowledge this calendar interaction. No maturity-specific model parameter was introduced. The same parameters apply to all seven expiries, where shape problems remain.

The multiplier-only model has a maximum fresh VIX-forward gap of2.74points; combined readout/curve brings this down to.337points without restoring the VIX smiles. The curve-only ablation is a diagnostic control, not a forward-constrained candidate in its own right. Full-candidate stock replay in the ablation agrees exactly with the main fresh comparison. This is a post-selection mechanism check, not an independent model-selection sample.

All approximate paired price-loss confidence intervals for positive strengths versus the alpha0-refit control contain zero. Strong-readout stock sampling and timestep errors remain material: for alpha1.5 the median absolute64/128 price change is1.50quote scales, with maximum11.23. Thus the observed available-IV improvements and quoted skew gaps are finite-grid numerical diagnostics, not accuracy-certified market fits. The matched-state VIX GH32/time8 versus GH48/time12 probe changes VIX by at most8.76e-6points across the declared128-state probes per expiry; this does not resolve rare-tail sampling.

The tested readout/factor combination therefore provides no satisfactory joint calibration. Its useful short-maturity interaction is saved for future work; it does not justify adopting the model or concluding that asymmetric readouts in general fail. See [four-way ablation](visualizations/asymmetric-readout028/SPX-readout-curve-ablation.png).

## Conclusion to retain

The first whole-response asymmetric multiplier does not solve the joint-fit problem in this bounded trial. Strong settings modestly improve SPX while spoiling VIX and making stock estimation harder. Do not repeat the same strength scan or claim that a lower available-IV MAE alone gives an acceptable double calibration.

A possible next refinement is to keep each component's mean variance allocation fixed when adding asymmetry, using the correct joint normalization inside the M2 component. That is a separate proposed model, not tested or adopted here. A different Gaussian factor/timescale is another hypothesis; this experiment does not establish that all exogenous readouts fail.

Figures: [SPX smiles](visualizations/asymmetric-readout028/SPX-asymmetric-smiles.png), [VIX smiles](visualizations/asymmetric-readout028/VIX-asymmetric-smiles.png), [response and variance allocation](visualizations/asymmetric-readout028/readout-tradeoff.png). Detailed results: `runs/m2_asymmetric_readout_001`.
