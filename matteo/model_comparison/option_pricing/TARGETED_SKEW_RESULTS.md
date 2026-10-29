# Targeted SPX skew experiments — 28 October 2026

**Both constructions implemented and checked; neither adopted.** Adaptive stock exposure preserves the complete VIX law but gives only modest SPX changes in these tests. Conditional stress-regime normalization reduces the VIX damage of the earlier global multiplier, but strong settings still spoil VIX smiles. Neither yields a satisfactory fit across the seven SPX maturities. These are bounded numerical experiments, not a model-family impossibility result.

## Fixed model and scope

All original M2 readout parameters, positive spike shift, exogenous OU dynamics, stress speeds/amplitudes/weights and floor remain fixed. No new bank, jumps, quintic, endogenous OU feedback or new constant-theta sweep. The previous protected16-level curve is the original baseline. The stock driver for readout-only comparisons remains the mixed theta=.8 rotation. All1751eligible/1294primarySPX quotes, sevenSPX maturities,217VIX quotes and fourVIX forwards remain. Prices are signed; invalid IVs and negative estimates are retained.

## Experiment1: adaptive stock exposure

In independent primitive Brownian coordinates write

$$dv_t=\mu_tdt+\mathbf k_t\cdot d\mathbf B_t,\qquad dS_t/S_t=\sqrt{v_t}\mathbf a_t\cdot d\mathbf B_t.$$

Use unit stock direction opposite the current variance sensitivity. For xi>floor,

$$\mathbf a_t=-\mathbf k_t/\|\mathbf k_t\|.$$

This minimizes instantaneous stock/variance covariance among unit directions. It does not maximize a finite-maturity option-skew slope. All11relevant old primitive Brownian directions and two stress directions are included. The original M2 cap derivative and instantaneous projection diffusion coefficients are used, not finite-step normalized OU innovations.

Implementation boundary convention: the common positive xi-floor factor cancels in normalization, so the code uses the gradient of H0. At xi=floor, the continuous limiting H0 direction is retained; all directions have zero instantaneous variance covariance there. If H0's gradient vanishes, use the old mixed stock loading. This refines the initial plan's generic zero-gradient fallback and is recorded explicitly in extended_audit.json.

Only the stock equation changes. The full Gaussian factor law, variance paths/readout and conditional future variance law stay fixed; hence all VIX prices and forwards are identical to the original model. The inherited checking forward gap up to.0552VIXpoints and old estimator error are unchanged, not newly eliminated. No exact market-forward equality or continuous-stock-martingale certificate follows from unit loading norm.

Existing retained OU streams are replayed unchanged. Missing old-bank Brownian increments are sampled conditional on the original Gaussian innovations with the full covariance, plus independent residuals. Stock loading and variance are held at each stock-grid start, including coarse-grid controls. No future-state loading enters a coarse step.

## Experiment2: stress-regime-normalized readout

Let a=1-w, A=(1-r)As+rAf, H0=aR+wA and D_alpha=exp(alpha softplus(-X)), with existing X=(Y-E_QY)/stationary_sd(Y). Set mRD(t)=E_Q[R_tD_t],mD(t)=E_Q[D_t] and

$$H_t^{new}=H_t^0D_t\frac{a+wA_t}{a\,m_{RD}(t)+wA_t m_D(t)}.$$

Then, using independence of the old bank and stress factors,

$$E_Q[H_t^{new}\mid Z_t^s,Z_t^f]=a+wA_t,\qquad E_Q[H_t^{new}]=1.$$

This fixes the average total variance response within each stress state. It does not preserve the internal mean allocation to the two algebraic terms, the complete conditional future variance law, VIX distribution or VIX forwards. VIX is repriced. Alpha0 recovers H0.

For future-time pricing define positive stress functions

$$F_1(A)=\frac{a(a+wA)}{a m_{RD}+wA m_D},\quad F_2(A)=\frac{wA(a+wA)}{a m_{RD}+wA m_D}.$$

Conditional E[Hnew] equals E[RD|old state]E[F1(A)|stress state]+E[D|old state]E[F2(A)|stress state]. This uses correct joint old-bank product moments, plus two-dimensional stress Gaussian quadrature. Products of separately normalized old factors are not substituted.

Alpha0,.5,1,1.5,2 is shared across all maturities. Each positive strength is evaluated on the original curve; a separate fit adjusts only seven VIX-sensitive variance levels with four training-estimator forward equalities, leaving nine off-window levels fixed. The alpha0-refit control is included. Every bounded fit is forward-feasible; maximum training forward gap is 3.36e-10VIXpoints. This is sampled-estimator equality, not population equality. Checking does not re-enforce those equalities. The VIX-loss-preserving preselection retains alpha0; no positive tested strength met the declared baseline+.1loss screen.

## Finer-grid results

The following SPX numbers use2048fresh paths and256steps/trading-day, with paired128/256 grids. The main experiment uses512replay and1024fresh stock paths, with paired64/128 grids. Alpha1.5 was frozen as an illustrative diagnostic before the finer-grid run; no winner was reselected from it.

| Variant | SPX interval loss | SPX IV MAE,pp | Invalid primary IVs | VIX interval loss |
|---|---:|---:|---:|---:|
| Original fixed loading | 1869.58 | 3.203 | 0 | 0.0501 |
| Adaptive stock loading | 1821.00 | 3.040 | 4 | 0.0501 |
| Regime α1.5, fixed curve | 1705.08 | 2.960 | 4 | 7.6842 |
| Regime α1.5, seven-level refit | 1560.53 | 2.897 | 2 | 6.1592 |

SPX IV MAE averages available IVs equally across maturities; unavailable IVs are counted. VIX numbers use the same separately generated4096full-factor states per expiry previously used by the asymmetric-readout study, plus the old65536path controlled anchor. These are previously examined research draws, not a newly independent market sample. The original/adaptive VIX entries are the same model law, not separate numerical fits. Conditional-normalization alpha1.5-refit has checking forward gap 0.0784points; its alpha0-refit control has VIX loss0.0469 and checking gap0.0024points.

At29days, IV(right)-IV(left) at quoted strikes nearest K/F=.9,1.1 is -5.22pp for the baseline, -4.90pp with adaptive stock exposure and -6.29pp for the alpha1.5-refit readout; market is -14.05pp. The missing short-end slope remains. The dramatic29day improvement of the earlier global-multiplier/curve combination is not reproduced here.

Compared with the earlier whole-response alpha1.5 multiplier's VIX loss25.725, the new alpha1.5-regime refit yields6.159 on the same VIX checking states/anchor and quotes. This supports the value of conditional normalization for limiting VIX damage, not adequacy of the resulting smile or a controlled comparison with every calibration parameter fixed.

## Uncertainty and implementation evidence

Eight antithetic-pair batches underpin stock uncertainty. The unchanged stress-only auxiliary reference uses the previous65536paths, its correct conditional lognormal pricing and a deterministic delta stock-difference control. That reference has no new multiplier/adaptive loading, so its conditional formula stays valid. It is shared across screening/checking; no future-stress-dependent delta is used. Price differences are paired on common factors. For the2048path finer-grid check, all reported approximate paired95%price-loss change intervals include zero. Available-IV point estimates are not accuracy certificates or proof of improved market fit.

At128/256 refinement, median absolute grid changes in quote scales are 0.090(baseline),0.230(adaptive),0.866(regime fixed) and0.831(regime refit). Sampling and strong-readout grid errors remain unresolved. Stock means/SEs, signed prices, native-IV failures, paired intervals and raw likelihood-weighted payoff concentration are saved; no mean renormalization, clipping or smoothing is applied.

Checks pass: original baseline stock replay 5.6e-14; independent scalar primitive/stock recursion 0; unit direction norm 4.44e-16; readout-gradient finite differences 6.58e-11; old-bank Brownian/OU standardized covariance error 1.08e-11; independent stress quadrature 1.84e-15. The conditional-mean identity holds algebraically, and interpolated means versus independentGH64 at declared probes agree to 3.21e-08; seven new forward derivatives agree within 1.3e-09. MomentGH48vs64 error at declared probes is 4.46e-08. Matched128state/expiry GH32/stress32/time8 versus GH48/stress48/time12 probes change VIX by at most 1.48e-05points. These are finite numerical checks, not uniform error bounds.

No new Lean theorem was introduced. Existing coefficient/normalization identities retain their original premises. Numerical audits do not certify the compiled simulator, VIX estimator convergence or continuous-time stock martingality. All original source fingerprints are checked at finalization; no prior production or experimental baseline is overwritten.

## Conclusion and next research question

The first construction is a valid way to change SPX dependence while preserving VIX, but maximum negative instantaneous covariance does not deliver the required SPX skew in this diagnostic. The second construction preserves more of the successful stress response and reduces VIX damage, but does not produce an acceptable double calibration in this bounded search.

Do not repeat the same constant-theta or global-multiplier strength scan. A more informative next question is which existing shocks predict persistently higher future integrated variance, rather than only a large instantaneous variance change. Aligning stock exposure with that future-variance sensitivity is a distinct, untested refinement that can still leave the factor/VIX law fixed. Whether short-lived feedback noise explains the present outcome needs that diagnostic; it is not established here. Broader stress/readout refitting and additional banks were not tested, so no model-family impossibility follows.

Figures: [finer SPX comparison](visualizations/targeted-skew028/SPX-targeted-256.png), [main SPX comparison](visualizations/targeted-skew028/SPX-targeted-128.png), [VIX comparison](visualizations/targeted-skew028/VIX-targeted.png), [readout response](visualizations/targeted-skew028/regime-response.png). Detailed run: `runs/m2_targeted_skew_001`.
