## 29 October 2026 — retain the endpoint/path lesson and its empirical limitation

> VIX options constrain where VIX can finish; SPX options also care about the volatility experienced along the way.

This may be useful for future modelling, but it is not yet established as fundamental to M2. The results of the approaches tested so far are not really encouraging: the missing SPX downside-volatility level remains, and no satisfactory joint calibration has been demonstrated. Preserve this distinction between a useful modelling insight, mathematical feasibility and empirical success.

See the equations, qualifications and experiment links in [the research-ledger lesson](RESEARCH_LEDGER.md#29-october-2026--lesson-to-retain-vix-endpoints-and-the-spx-volatility-path). The markets still share one variance process; finite VIX quotes constrain expiry distributions, while SPX also depends on variance along the stock path and its interaction with shocks.

## 29 October 2026 — existing stress variance as correction budget tested

Outcome clarification: sampled new gain1.20% has paired95% loss-change interval[-120.12,76.59], includes0; new-minus-previous interval[-46.34,128.55] also includes0. No improvement is established.29d left IV18.97->19.26% vsmarket26.90%; missing downside level persists. Day7 state probes exceed old absolutef/2limit, so lack of numerical success is not evidence that the budget was never enlarged.

Read [AVAILABLE_VARIANCE_RESULTS.md](AVAILABLE_VARIANCE_RESULTS.md), run `runs/m2_available_variance_001`, plots `visualizations/available-variance029`. Potential P_j=c_j(t)beta_j(t)erf(U/sqrt2)A_j exp[-lambda(Z_j-b_j Q_j(t))^2], j=slow/fast stress; h_j=(partial_t+L)P_j, v=v0+eta sumh. Include full time/stress generator terms; simple multiplication of old h by variance is not protected. Clock zeros at all16curve breakpoints (includes quoted VIX endpoints) remove curve-coefficient jump terms. Analytic |h_j|<=c_j A_j implies v>=ORIGINALf=.0144 for|eta|<=1; can exceed old fixed f/2 budget. OU law, oldM2/cap/Qnormalizer/positive spike shift, stresses, curve and stock driver fixed. Odd independent U preserves ideal mean variance. Four quoted VIX values protected state by state, baseline errors inherited; no exact market-forward claim, no unquoted/full-process protection. This is a whole stress-scaled construction, not a pure amplitude ablation versus prior rich correction; calendar/damping/floor differ too.

512pilotpaths/38configs, selected half-life3.329d/lambda0.125/eta1.0 frozen before1024fresh targets. All1751eligible/1294primary/sevenSPXmaturities,217VIX/fourdates retained. Fresh256step losses baseline1948.66/previousrich1888.60/new1925.31, gain1.20%vsbaseline/-1.94%vsprevious.29d gaps-5.72/-6.69/-6.14pp vsmarket-14.05; new leftIV19.26%vsmarket26.90%. Paired intervals, stock128/256, conditional discrete sums, tail/budget probes saved. Compiled full formula, independent full-potential generator differences, conditional analytic/GH oracle, old replay and scalar stock recurrence audited. No new Lean proof, uniform numerical bound or martingale theorem. Auxiliary65536paths and market snapshot reused. No candidate adopted or production/old artifacts changed.

## 29 October 2026 — richer protected variance shapes tested

Read [RICH_PROTECTED_VARIANCE_RESULTS.md](RICH_PROTECTED_VARIANCE_RESULTS.md), run `runs/m2_rich_protected_variance_001`, figures `visualizations/rich-protected-variance029`. Same frozen M2/OU bank/stresses/stock theta=.8/curve. Replace first product-envelope/global scale by interval sine-squared potentials b(t)phi(U), normalized locally by K=kappa*C/2+sqrt(A^2+(kappa*C)^2/4). At all quoted window endpoints b=b'=0, so jumps of the interval scale add no jump term. Tanh/erf shapes, widths1/2, existing modes7/6/5, both signs. Global |h|<=1 and L1 coefficient budget retain the same v>=f/2=.0072; no clipping or new bank. Conditional Dynkin still protects all217VIX quotes/four quoted dates state by state, inherits baseline errors, not exact market forwards. No unquoted/full-process invariance.

512pilotpaths/26choices; frozen best new choice {'mode': 1, 'width': 1.0, 'kind': 1, 'name': 'local 9.5d erf width1'};1024fresh checking paths and reused65536path auxiliary. Fresh256step SPX losses baseline1945.52,previous1909.39,richer1884.31; richer gain3.15%vsbaseline/1.31%vsprevious. IVMAE baseline3.209/previous3.162/richer3.137pp.29day gaps-5.36/-5.75/-5.85pp vsmarket-14.05. The richer shape improves aggregate SPX fit in this checking sample, but the observed downward skew remains far from market. All1751eligible/1294primary/sevenSPXmaturities, signed estimates/invalid IVs retained. Paired intervals, grid sensitivity and tails saved. Compiled formula/conditional cancellation, old-model replay, independent scalar recurrence and finite-difference generator checked. No new Lean proof or full positivity certificate from numerical probes alone; positivity rests on the analytic bound in the report.

No candidate adopted; production model and old artifacts unchanged. Preserve this failed-or-partial shape search rather than repeating it. The following entries are historical.

## 29 October 2026 — VIX-protected additive variance correction tested

Read [PROTECTED_VARIANCE_RESULTS.md](PROTECTED_VARIANCE_RESULTS.md), run `runs/m2_protected_variance_001`, figures `visualizations/protected-variance029`. Retain full M2, stress mixture, curve, cap, OU dynamics and theta=.8. Add v=v0+(f/2)sum(c_j h_j), h_j=(partial_t+L_j)[a(t)tanh(U_j)/K_j], common envelope zero at all eight quoted VIX-window endpoints. Three existing modes7/6/5 have half-lives3.329/9.481/27.004days. L1(c)<=1 proves v>=f/2=.0072 without clipping; original exact floor f=.0144 is not retained. Centered odd correction preserves ideal mean variance. Conditional Dynkin identity protects all quoted-date VIX values/strikes/forwards state by state; inherits baseline VIXloss.05008/maxforwardgap.05517, does not newly match market forwards exactly, and does not protect unquoted VIX dates/full process law.

512pilotpaths, baseline+12axis configurations, selection frozen at27day mode/c=+1 before1024fresh paths. All1751eligible/1294primary SPX/sevenmaturities,217VIX/fourdates, signed estimates retained. Fresh256step price loss1841.07->1801.18 (2.17% improvement), IVMAE3.169->3.121pp;29day skew-6.00->-6.59pp versusmarket-14.05. Paired bootstrap loss-changeCI[-51.17,-28.63], excludes deterministic/model error. Compiled formula, all-mode conditional cancellation, zero-correction replay and independent nonzero scalar recurrence audited; stock128/256refinement and discrete-window audit saved. No new Lean theorem. Fresh targets use reused65536path auxiliary, same market snapshot.

Conclusion: modest SPX improvement with explicit quoted VIX protection, missing skew remains. No candidate adopted or production/old run changes. Do not repeat this exact test as untested; investigate positivity-budget efficiency and richer protected bases before new banks. New amplitudes require a fresh global positivity argument. The following proposal entry is historical and has now been tested.

## 29 October 2026 — proposed M2 variance correction with quoted VIX protection

User preference: retain as much of current M2 as possible while separating volatility-path controls from VIX expiry distributions. This entry is a mathematical design and finite Gaussian/quadrature audit, not a fitted/adopted model or an implementation of the VDV paper.

Keep the full exogenous Gaussian OU bank and pricing filtration, old capped/Q-normalized M2 function R with positive spike shift, independent slow/fast stresses, baseline mixture and deterministic forward-variance curve. Write the resulting base variance as v0(t,Z_t). Add v_new=v0+epsilon*h(t,Z_t). The new term acts on variance after the existing readout; it does not feed returns or volatility into the OU equations. Existing Gaussian paths remain reusable.

For each quoted VIX expiry Ti, seek E_Ti[integral_Ti^(Ti+Delta) h(s,Z_s)ds]=0, with Delta30/365. This preserves model-derived VIX at that date state by state, hence all VIX option strikes and the forward there, assuming the same pricing filtration and integrability. It protects the model values, including their existing fitting/estimation error; it does not newly certify exact market/population forwards. It is stronger than matching four means or preserving unconditional E[v]. It does not preserve VIX at every unquoted date or allow arbitrary changes of the VIX copula at the quoted dates.

An explicit construction uses the unchanged factor generator L. Let q(t,z) be smooth with q(Ti,z)=q(Ti+Delta,z)=0 for all states and all four quoted expiries. Define

$$h=(\partial_t+\mathcal L)q.$$

Dynkin's conditional identity gives

$$E_{Ti}\int_{Ti}^{Ti+\Delta}h(s,Z_s)ds=E_{Ti}[q(Ti+\Delta,Z_{Ti+\Delta})]-q(Ti,Z_{Ti})=0.$$

This can change the distribution and timing of realized integrated variance and its interaction with stock shocks even though the conditional expected VIX-window total is unchanged. No SPX skew improvement follows automatically. The correction is not the old positive multiplier or averaging filter; its compensation is conditional and window-specific.

For a concrete scalar design use an existing feedback mode, centered by its dated deterministic Q mean, dU=-kappa Udt+sqrt(2kappa)dB, U0=0. Set q=a(t)tanh(U)/K. The common smooth envelope a(t)=product_b tanh((t-b)/tau)^2 vanishes at all quoted window endpoints b. Tau10days is a generic common width, not a special factor for one maturity. All quotes are treated by the same condition. Here

$$h=\{a'(t)\tanh U+a(t)[-\kappa U\,\mathrm{sech}^2U-2\kappa\tanh U\,\mathrm{sech}^2U]\}/K.$$

One analytic global bound is K=1/tau+kappa[1/2+4/(3sqrt3)], giving |h|<=1. Product-envelope derivative bound |a'|<=1/tau follows from |f'|<=(1-f^2)/tau for f=tanh^2 and the product identity (1-x^2*y^2)-[(1-x^2)y+x(1-y^2)]=(1-x)(1-y)(1-xy)>=0. Centered Gaussian symmetry makes E[q]=E[h]=0, so the ideal forward variance mean E[v_new]=xi is preserved. Multiple corrections can use distinct existing modes, with a sum-of-amplitudes bound.

Positivity is a real constraint: v0>=f and |epsilon|<=f/2 imply v_new>=f/2>0. Thus this simple construction retains a positive floor but does not automatically retain the original exact hard lower bound f. Preserving that exact bound requires additional support/shape restrictions; in particular where xi=f, nonnegative variance above f with the same mean forces v=f almost surely. Do not silently claim that every original floor interpretation is unchanged.

Design audit `m2_protected_variance_design.py`, output `runs/m2_protected_variance_design_001/design_audit.json`: existing mode6 has half-life9.481468days. Fifty conditional-window probes (four expiries, five conditioning states, GH64/96) have maximum integral residual3.20e-11; zero unconditional-mean probes<1.6e-17; nonzero observed h magnitude.48966 under the analytic bound1. A deliberately unquoted day40window has nonzero conditional integral, confirming the finite protection scope. These are finite numerical correspondence checks, not a new Lean proof, full compiled-model audit, option calibration or population-forward certificate.

Calendar dependence is a calibration restriction and a trade-off. The envelope protects all current quoted VIX dates, not a universally prescribed full VIX surface. It is an initial diagnostic more restrictive than preserving only marginals. A less calendar-dependent extension could fit generic time/state bases under explicit VIX price constraints, but that numerical relaxation would require new protection audits. Do not present it as already implemented.

Recommended next test, not yet run: keep the VIX-fitted M2/stress baseline frozen, use one or two such bounded corrections from existing Gaussian modes, independently check the new conditional-exposure cancellation in pricing code, then fit correction amplitudes to all seven SPX maturities under positivity bounds. Reprice every SPX path; do not resimulate parameter-independent OU paths. Keep all quote memberships and signed estimates. The protected VIX prices may be reused only after the algebra-to-code bridge is checked. This is a proposed M2-preserving extension; production artifacts remain unchanged.

## 29 October 2026 — VIX first inverse construction studied

Source: Nicola F. Zaugg and Lech A. Grzelak, [The VIX-Derived Volatility Model: A VIX-first Joint SPX-VIX Framework](https://arxiv.org/abs/2608.01479), v1. Local original and prepared page images: `literature/vdv_2608_01479`. Relevant PDF pages: 4–8 (model and backward relation), 10–13 (polynomial/general coupling), 15–17 (marginal-preserving dynamics changes), 18–20 (market example), 22–25 (local volatility and proofs). The source identifies its license as CC BY4.0. The full27-page paper was prepared locally; equations and numerical coefficients were cross-checked with the original PDF and page19/20 images.

The new hypothesis is an inverse construction: start from an explicitly calibrated VIX process X, then solve for instantaneous SPX variance. With X=VIX/100, Delta=30/365,

$$dX_t=b(t,X_t)dt+a(t,X_t)dB_t,\qquad v_t=\psi(t,X_t),\qquad dS_t/S_t=\sqrt{v_t}\,dW_t.$$

The rolling consistency requirement is

$$\Delta X_t^2=E_t\int_t^{t+\Delta}\psi(u,X_u)du.$$

Under the stated Markov/regularity assumptions, the source's backward relation is

$$\psi(t,x)=E[\psi(t+\Delta,X_{t+\Delta})\mid X_t=x]-\Delta[2xb(t,x)+a(t,x)^2].$$

A compatible terminal strip and nonnegative variance are additional requirements. The general local-volatility implementation uses numerical projection/backward PDEs; the consistency guarantee is therefore conditional on admissibility and numerical error, rather than automatic from naming X a VIX process.

The paper's key separation is that changing mean reversion kappa requires recalibrating VIX local volatility and its deterministic mean curve to retain the prescribed VIX expiry marginals. European VIX prices depend on those marginals; SPX depends additionally on temporal dynamics. This differs from keeping the entire factor law fixed and changing only stock coupling. It is not another fixed-sigma mean-reversion sweep. The reported experiment uses six maturities and rho=-.95; the SPX call wing retains discrepancies (PDFpage20). Its market snapshot and quote selection differ from ours, so its figures are evidence for that implementation, not a replication on our data.

Relation to the supplied Guyon slides: their R2=h+K*v averages past variance, with v=(1-chi)nu(t)g(Z)+chiR2 and Gaussian OU Z. VIX squared averages conditional expected future variance. Both involve integral operators, but the operators, information sets and calibration directions differ. Our previous normalized positive filtering experiment tested the forward construction, not this VIX-first inverse construction. Do not describe VDV or all hysteretic models as already tested or rejected by that filter result.

Independent positivity checks: the rounded terminal polynomial displayed on PDFpage19, 1.16x^2-4.19x-101.37, has positive root11.32704764 in the example's VIX-point units and is negative below that value on the positive axis. This is a sign check of the displayed terminal polynomial; the actual reachable support, extrapolation and whether author code floors it remain unverified. The text permits post-hoc flooring (Section3). If an exactly consistent psi is replaced by psi+=max(psi,0), implied squared VIX increases by E_t[integral(-psi)+]/Delta. A floor therefore needs a new rolling-consistency audit; it cannot silently be regarded as preserving the original identity.

A separate illustrative Brennan-Schwartz calculation, not a replication of the authors' fitted local-vol model, uses kappa12/theta.2/sigma1 and Delta30/365. Its exact polynomial rolling-operator inversion gives coefficients[.00585282,-.28526228,2.226658], algebraic residual below3e-17, but minimum-.00328358 on an interval between X=.02565 and.10246. A positive GBM control satisfies its rolling identity to numerical precision. This demonstrates that solving the algebraic inverse equation does not alone establish nonnegative variance; it does not prove impossibility of every time-dependent coupling. Code `vdv_consistency_diagnostic.py`; exact parameters, independent matrix-quadrature checks and results in `diagnostics/algebraic_audit.json`. No flooring, market calibration, new Lean result or production change.

Architecture implication for our fixed-bank preference: the source's local-volatility VIX paths are parameter-dependent; Brownian-draw reuse is not complete OU-state reuse. An alternative to investigate is X_t=F_t^{-1}(Phi(U_t)), with fixed-bank Gaussian U_t standardized using its dated Q mean/variance and target VIX CDF F_t. For t>0 this gives the target marginal exactly in the ideal model. Initialization/continuity, the entire interpolated CDF family, tails/second moments, conditional rolling consistency and nonnegative psi must still be checked. A one-factor Gaussian prototype keeps the conditional state Markov. Mixtures of different OU speeds require the full multivariate state, not an assumption that the scalar mixture is Markov. Such a construction is our proposal, not the authors' model, and it need not retain the current M2 readout unchanged.

A dataset-specific limitation is important: current SPX maturities extend to1,030days while quoted VIX maturities end at153days. A VIX-first model needs a documented long-horizon extension, rather than claiming those unobserved VIX marginals are market-identified. The paper's roughly six-month comparison does not by itself validate our longer horizon.

Next research priority: freeze/check VIX marginals including forwards and second moments; test a small VIX-first inverse prototype for positive coupling and pointwise rolling consistency without silent flooring; only then investigate SPX fitting while rechecking VIX prices. Retain every existing market quote/maturity and simulate variance through max(last SPX expiry,last VIX expiry+30days). Do not launch another memory-filter or constant-correlation grid by default. No VDV calibration has yet been implemented or adopted.

## 29 October 2026 — normalized linear-memory filters tested

Read [MEMORY_FILTER_RESULTS.md](MEMORY_FILTER_RESULTS.md); corrected run `runs/m2_memory_filter_001`, figures `visualizations/memory-filter029`. User authorized the Guyon-slide-inspired linear memory experiment. No nonlinear GLE or actual-return feedback was implemented. Gaussian OU dynamics, fast capped/Q-normalized M2 with positive spike shift, two stress factors, stock theta=.8 and floor stay fixed. For J=R,As,Af, dMJ=a(J-MJ)dt with MJ0=1. Test M2-only filtering and whole-multiplier filtering with half-lives3/10/30days, shares.25/.5, uniform across maturities. Variance=f+(xi-f)H; unit means preserve E[v]=xi algebraically, not exact population VIX forwards. This filters nonlinear signal history rather than slower Gaussian states. All downstream prices are rebuilt with shared factors.

Pilot26configurations: baseline and12nonzero fixed-curve controls plus13refits of the same seven active curve levels under four training-forward equalities; nine off-window levels frozen.512pilotstock paths,512weighted terminal-plus-history VIX samples/expiry. Selection M2-only3days/share.25, frozen before1024new stock paths and1024new joint VIX terminal states/conditional histories per expiry. Conditional VIX includes realized memory states, not just current OU projections. Keep1751eligible/1294primarySPX/sevenmaturities and217VIX/fourforwards, signed estimates/no clipping.

Fresh256step results: baseline SPX loss1903.59/IVMAE3.238pp/VIXloss0.05008; M2 fixed1899.76/3.236/0.05181; M2 refit1901.94/3.237/0.04882; curve alone1905.77/3.239/0.04579; whole3days/share.25refit1912.77/3.241/0.05125. M2refit price-loss reduction0.087%, approximate paired eight-batch bootstrap loss-change interval[-18.28,15.92], includes0.29day right-minus-left IV gap baseline-5.57,M2refit-5.56,whole-5.29pp versusmarket-14.05. No meaningful SPX improvement; missing skew persists. Fresh M2refit max forward gap0.00762points; all four training equalities remain sample constraints.

History conditioning numerical issue was caught and corrected before final conclusions: direct inversion was unreliable; a whitened cutoff1e-14 pilot was rejected for projected history covariance error.02146 and preserved as `m2_memory_filter_001.conditioning_diagnostic001`. Corrected cutoff1e-11 retains21/21/22/22terminal directions, covariance-gated below2e-6; measured maximum5.53e-7. Actual terminal projection error<9.5e-9. Stable cutoff1e-10vs1e-11 gives identical matched-history probes. Scalar stock/filter recurrence and zero-share full-horizon replay errors0; independent scalar bridge-history recurrence2.22e-15. No Lean certification/uniform pricing bound.

Historical64/128 and futureGH32/time8/32lag versus48/time12/64lag refinements have large unweighted tail maxima. Do not interpret these without likelihood weights: maximum weighted absolute VIX change.000956points and native option-price change.001277quote scales in the declared probes. All extreme states retained. Stock grids128/256, paired SE/bootstrap intervals, raw/weighted payoff concentration, memory means/SEs and invalid IVs are saved. Reused65536-path conditional-lognormal auxiliary reference and same market snapshot limit independence claims.

Durable conclusion: positive linear averaging preserves mean variance and adds persistence, but it weakens the filtered signal's direct Brownian response. These tested filters do not supply the missing SPX skew. Longer whole-multiplier filtering also damages the VIX training fit. Do not repeat this same filter grid as untested; this is not a no-go theorem for nonlinear hysteresis or other exogenous Gaussian readouts. No candidate adopted; production model and all previous experiment artifacts remain unchanged.

## 28 October 2026 — slower existing-bank downside signals tested

Read [SLOW_READOUT_RESULTS.md](SLOW_READOUT_RESULTS.md); runs `runs/m2_slow_readout_001` and `runs/m2_slow_readout_002`, figures `visualizations/slow-readout028`. The previously proposed slow-readout hypothesis is now tested. Existing feedback modes nearest log-half-life targets2/8/32days have actual half-lives3.329/9.481/27.004days. Original fast M2 function and positive spike shift retained; no new Gaussian bank or endogenous return/volatility feedback. Stock theta=.8 fixed. One slow centered unit-stationary-SD mode at a time enters D=exp(alpha softplus(-X)). Conditional-stress normalization H=(aR+wA)D(a+wA)/(a E[RD]+wA E[D]) preserves E[H|stress]=a+wA, not the VIX law. R and D are correlated; product moments use their full joint Gaussian law.

Exploratory alpha.5/1 then reactive.125/.25 grids, fixed curve and same seven active curve levels refitted under four training-forward equalities; nine off-window levels frozen.512paired pilot stock paths and512training VIX states/expiry; previously examined4096checking states. Frozen weak selection9.481days/alpha.125 before1024new stock paths and2048new VIX states/expiry. Five fresh controls include weak fixed curve, weak refit, curve-only and strong3.329day/alpha.5refit. All1751eligible/1294primarySPX/sevenexpiries and217VIX/fourforwards retained, signed prices/no smoothing/clipping.

Fine256step results: baseline SPX loss1866.88/IVMAE3.176pp/VIXloss0.050; weak fixed1816.32/3.130/0.089; weak refit1811.15/3.128/0.085; curve alone1862.10/3.174/0.041; strong refit1749.73/3.030/0.553. Weak refit improves SPX price loss2.99%; approximate paired eight-batch bootstrap loss-change interval[-95.71,-7.82], excluding deterministic error. Strong contrast interval includes0.29day right-minus-left IV gap baseline-5.98,weak-6.17,strong-7.05pp versusmarket-14.05. Missing skew persists. Weak refit checking forward gap max0.0174points: training equality is not exact population matching.

Alpha-zero old-kernel replay2.22e-16; independent conditional tensor product-moment checks<6e-7relative; GH48/64normalizers<6.66e-7; independent97step scalar five-candidate stock recurrence0error. Coupled128/256stock median changes.117/.135/.128/.118/.183quote scales; matched64state/expiry VIX GH32/time8 versus48/time12 max7.32e-6points. Numerical probes only, no Lean certification/uniform tail bound. Stock/VIX likelihood-weighted payoff concentration, SEs, stock means and invalid IVs retained.

Durable conclusion: this simple slow one-sided multiplier provides a small SPX response but no satisfactory double fit. Stronger changes damage VIX despite conditional-mean normalization. Do not repeat this same strength scan or describe the slow-readout hypothesis as untested. This is not a no-go theorem for other slower projections/readout structures. Production model unchanged; no candidate adopted. Keep fast Zumbach-motivated component in any subsequent design and rebuild downstream prices whenever a readout changes.

## 28 October 2026 — future-integrated-variance stock exposure tested

Read [FORWARD_LEVERAGE_RESULTS.md](FORWARD_LEVERAGE_RESULTS.md); run `runs/m2_forward_leverage_001`, figures `visualizations/forward-leverage028`. User authorized the proposed persistent-leverage test. Stock follows negative Brownian sensitivity of conditional next30day integrated variance; same horizon for every maturity. Existing full24mode Gaussian bank, variance readout/positive spike shift, stress parameters, floor and16-level curve stay fixed. VIX distribution/options/forwards unchanged, including inherited numerical error. No endogenous factor feedback, new bank, readout fit or constanttheta sweep.

Full-state conditional Gaussian gradient: independent C convolution GH48, tensor Y/F GH12,0plus22geometric lag knots; hat weights integrate xi-floor exactly across curve breaks. Fixed last xi continued only for signal through1060days. Direction updated4or8times/tradingday on shared factors; paired128/256stock grids. GH12vs32 max probe direction difference.000514;22vs44lagknots.01269; all13primitive finite differences9.3e-12; scalar97step stock recursion/unitnorm2.2e-16; frozen baseline/instant full replay0. No new Lean or continuous-martingale theorem.

1024freshstock paths/eightantithetic batches, prior65536auxiliary reference reused. At256steps: fixed price loss1909.19/IVMAE3.201pp; instantaneous1978.79/3.163; future4updates1969.78/3.175; future8updates1941.31/3.157. Invalid primary IVs1/1/2/3 retained. Future8 loss is1.68%higher; all approximate paired loss-difference intervals include0.29day quoted right-minus-left IV gap nearK/F=.9,1.1 fixed-6.39pp,future8-4.51,market-14.05. No demonstrated fit improvement; missing skew persists. VIX217quotes/four forwards inherited unchanged; no new exact population-forward claim. All1751eligible/1294primarySPX/sevenmaturities retained, signed prices/no clipping. Update4to8 median.253quote scales versus paired SE.499; future8 stock128/256 median.144scales. Final sampling/tail accuracy unresolved.

Durable conclusion: do not repeat this same30day stock-direction rule, constanttheta grids or the global multiplier scan. This is not a no-go theorem for all correlations/horizons. Next untested hypothesis: change readout memory using slower projections of the existing OU modes, preserving the fast Zumbach-motivated component. Exact stationary Gaussian memory diagnostic of current feedbackY gives correlation.170after1day,.0479after7days,.0173after30days; this suggests a target but does not prove the cause of SPX misfit. No new projection/readout/bank adopted; baseline artifacts preserved.

## 28 October 2026 — adaptive stock exposure and stress-regime normalization tested

Read [TARGETED_SKEW_RESULTS.md](TARGETED_SKEW_RESULTS.md); run `runs/m2_targeted_skew_001`, figures `visualizations/targeted-skew028`. Both user-authorized constructions implemented separately. Exogenous OU dynamics/banks, original M2 parameters/positive spike shift, stress parameters and floor fixed; all sevenSPX maturities/1751eligible/1294primary and217VIX quotes/four forwards retained. No new constant-theta sweep, endogenous OU feedback or production adoption.

Adaptive stock loading is opposite current H0 diffusion sensitivity in all13relevant independent Brownian directions. For xi>floor this maximizes negative instantaneous logstock/variance covariance; at xi=floor use the continuous limiting H0 direction, and at zero H0 gradient fall back to the original rotation. All factor/variance/VIX laws unchanged. Correct missing old-bank Brownian increments are sampled conditional on retained OU innovations with covariance error1.1e-11; original OU streams replay unchanged.

Regime readout Hnew=H0 D[(1-w)+wA]/[(1-w)E(RD)+wA E(D)], with existing feedback X and alpha0,.5,1,1.5,2. It preserves E[Hnew|stress states]=(1-w)+wA, not the full VIX law or separate internal component budgets. Compare fixed curve and seven-level forward-constrained refits; every training forward equality feasible(maxgap3.4e-10), but no positive strength meets the declared VIX baseline+.1loss screen. Alpha0 retained before fresh stock checking.

Main512replay/1024fresh stock paths; extra2048fresh paths on paired128/256grids. At256steps, SPX price loss/IVMAE: baseline1869.58/3.203pp; adaptive1821.00/3.040; regimealpha1.5fixed1705.08/2.960; regimealpha1.5refit1560.53/2.897. Invalid primary IVs0/4/4/2; all approximate finer-check paired loss-change intervals include0. Median128/256 grid changes.090/.230/.866/.831quote scales. 29day quoted right-minus-left IV gap nearK/F=.9,1.1: baseline-5.22pp,adaptive-4.90,regime1.5refit-6.29,market-14.05; required slope still missing.

VIX baseline/adaptive loss.05008identical law. Regimealpha1.5fixed loss7.684; refit6.159 with checkingmaxforwardgap.0784points, versus earlier globalmultiplierrefit25.725. Same4096VIX checking states and old65536anchor; research stream reused, not newly independent market validation. All signs/IVfailures/likelihood tails retained. No satisfactory double calibration or family-impossibility conclusion. Independent scalar,gradient,unitnorm,covariance,stressconditionalmean,7forwardderivative and quadrature checks pass; no new Lean/continuousmartingale theorem. Stock sampling/grid accuracy unresolved.

Do not repeat constanttheta/globalmultiplier scans. Next untested question: target shocks that raise future integrated variance persistently, rather than maximizing only instantaneous covariance. Keep OU factor-law reuse and VIX preservation when studying that refinement. Additional banks and broad stress/readout refitting remain untested.

## 28 October 2026 — asymmetric whole-response readout trial completed

User rejected further theta scans and authorized the asymmetric readout. Read [ASYMMETRIC_READOUT_RESULTS.md](ASYMMETRIC_READOUT_RESULTS.md); run `runs/m2_asymmetric_readout_001`, figures `visualizations/asymmetric-readout028`. New D_alpha=exp(alpha softplus(-X)), X=(Y-E_QY)/stationary_sd(Y), uses the existing equity-linked feedback bank. Normalize the whole H0*D product correctly. Theta=.8 and all Gaussian dynamics/readout/stress parameters/floor stay fixed; no endogenous feedback or new bank. Test alpha0,.25,.5,1,1.5,2; seven VIX-sensitive variance levels refit per strength, nine off-window levels frozen. All four training controlled-forward equalities hold within8.5e-10, not exact population-forward matching.

Fresh1024stock paths and fresh4096full-factor VIX draws/expiry: alpha0-refit SPX IVMAE3.278pp, VIX loss.0468; alpha1.5 combined IVMAE2.832pp, VIX loss25.725; six primary signed prices negative/invalid IVs and median stock SE7.88quote scales. Whole-product normalization changes mean stress excess-variance share from92.3% to53.2% at day153, despite unchanged nominal w. The VIX-preserving preselection retains alpha0. No model adopted or original baseline overwritten. All seven SPX maturities/all1751eligible/1294primary and all217VIX quotes/four forwards retained; no clipping/smoothing.

Important four-way paired ablation: at alpha1.5,29-day IV(right)-IV(left) near K/F=.9,1.1: original -5.86pp; multiplier only -6.05; curve only -8.51; both -13.44 versus market -14.05. Strong short-end improvement needs both changes and their interaction; the refitted27–29day interval also affects the29daySPX expiry. Do not attribute the combined result to the multiplier alone. This is no new maturity-specific model parameter. Multiplier-only VIX loss21.278 and maxforwardgap2.74points; combined freshgap.337 despite training equality. Other SPX maturity shapes remain mismatched.

Checks: independent direct-tensor/C-convolution3.7e-11scaled; scalarstock1.1e-16; zero-alpha fullreplay5.6e-14;7forwardderivatives1.4e-9; normalizerGH48vs64 max3.1e-8; VIX GH/time matched probe max8.8e-6points. Strong-stock grid error unresolved (alpha1.5 median1.50,max11.23quote scales for64/128); approximate paired price-loss intervals include zero. Stocks/auxiliary reference and VIX anchor uncertainty remain explicit. Fresh VIX factor draws use full38-dimensional defensive-mixture weights. No new Lean or continuous-martingale theorem.

Durable conclusion: this bounded whole-response multiplier on the existing feedback factor does not solve double calibration; do not repeat the same alpha/theta scan. A possible next hypothesis is component-specific joint normalization that preserves mean M2/stress allocations, or a different equity-linked memory factor. Neither has been tested/adopted here. Full stress/readout refitting was outside this isolated trial, so no model-family impossibility claim follows.

## 28 October 2026 — fixed-factor leverage scan completed

Read [LEVERAGE_SCAN_RESULTS.md](LEVERAGE_SCAN_RESULTS.md); run `runs/m2_leverage_scan_001`, figures `visualizations/leverage028`. Tested theta=0,.2,.4,.6,.8,.95 with mixed, slow-only and fast-only stock/stress driver rotations. All variance parameters, protected nine-level curve and exogenous OU factor dynamics/law are fixed; VIX law unchanged. No endogenous volatility or return feedback. Four sufficient statistics per path/expiry/grid make all theta/direction evaluations reusable.

Mixed theta=.8 was selected on replay price loss before fresh expensive checking. Fresh 1024 paths: fine-grid equal-maturity price loss 3788.477 -> 1861.376 (-50.9%); native IV MAE 4.493 -> 3.199 percentage points; inside bid-ask 7 -> 18 /1294; no primary invalid IVs or negative prices in this fresh comparison. Approximate paired bootstrap loss-change interval [-2182.8,-1622.7]. Mixed beats fast-only/slow-only. Theta=.95 looks better on the fresh full grid but was not the preselected winner. Short-maturity left wing is still too low: 29-day IV(right)-IV(left), near K/F=.9,1.1, moves -1.52 -> -5.58 pp versus market -14.05 pp. This is a material shape improvement, not a satisfactory double calibration or adopted model.

Eight independent expensive batches, antithetic pairs; 65536 auxiliary paths and eight independent reference batches. Cheap reference reused between screening/checking, with its uncertainty included; not a wholly independent checking estimator or market OOS. Median fresh selected MC SE is 1.36 quote scales; nested 64/128 changes remain material. Correct OU/Brownian joint sampling and rederived deterministic-delta conditional-reference control are implemented. Full theta=0 replay, independent scalar target/reference statistics and covariance audits pass. Existing Lean algebra remains scoped; no new continuous-time martingale/VIX/simulator theorem. Next: independent reference/tail precision in the promising mixed-theta range, then revisit remaining left-wing shape and protected curve/readout fitting. Keep existing baselines intact.

## 28 October 2026 — exogenous-only skew research and fast-scan Lean bridge

User explicitly reiterated: no volatility/return-dependent OU dynamics or endogenous feedback; retain fixed-bank reuse and the existing Zumbach-motivated design. Read [SPX_SKEW_RESEARCH_AND_FAST_ABLATIONS.md](SPX_SKEW_RESEARCH_AND_FAST_ABLATIONS.md). Othmane's strong leverage and asymmetric return-feedback readout explain his SPX skew, but that endogenous bank is outside scope; his 19 October report retains VIX wing and estimator-validation limitations. Literature supports exogenous Gaussian-factor alternatives and targeted leverage studies. No rotation pricing engine, new calibration or response sweep implemented in this review.

Proposed next screen: keep variance/factor dynamics fixed; scan a single negative stock/stress driver rotation, plus fast-only versus slow-only leverage directions. This preserves the entire variance-factor/VIX law. Cache I=sum(v dt), U=sum(sqrt(v) dW 0), V=sum(sqrt(v) dBstress); terminal log-stock equals -I/2+sqrt(1-theta^2)U-theta V. Separate Bs/Bf integrals allow direction scans with four numbers per path/expiry. No OU resimulation is needed for these fixed-variance scans. Readout/xi changes require rebuilding downstream sums; OU dynamics/rate changes require new state construction. Joint OU/Brownian increments and control identities must be checked/rederived before pricing.

New formalization/LeverageRotation.lean: coefficient normalization, sign of the proposed stress covariance contribution, and exact three-statistic reconstruction of the declared discrete log-stock recursion are Lean-checked in leverage_audits/attempt_001 with warnings-as-errors, no sorryAx/custom axioms. This is a restricted algebraic bridge, not a Brownian-process theorem, full VIX-invariance formalization, simulator certificate or market fit. Keep numerical response/ablation evidence separate. Existing protected-fit conclusion remains binding: price loss improves while IV shape does not; do not repeat curve-only fitting without a new hypothesis.

## 28 October 2026 — protected nine-level fit completed; durable conclusion saved

User explicitly asked to save this conclusion. Read [PROTECTED_CURVE_CONCLUSION.md](PROTECTED_CURVE_CONCLUSION.md) and [PROTECTED_CURVE_RESULTS.md](PROTECTED_CURVE_RESULTS.md) before repeating a curve-only experiment. Run `runs/m2_protected_curve_001`, figures `visualizations/protected-curve028`. Only nine off-window xi levels changed; full VIX law and original successful VIX artifact preserved. All 1294 primary/seven SPX maturities fitted, all 1751 eligible checked; all 217 VIX quotes/four forwards retained. Revised numerical policy explicitly authorized by user: all-quote precision screen is final acceptance, not a prerequisite for a bounded exploratory fit.

Training price loss 6406.258→3665.656. Genuinely fresh 512-path Q check,256 antithetic pairs/eightbatches, coupled 64/128 stock grids and new reference scrambles:128-step price loss 5979.323→3895.071(34.9% decrease), paired approximate 95%loss-difference interval[-2495.4,-1690.9]. BUT equal-maturity native-IV MAE 4.345→4.710 pp, inside bid–ask 11→8/1294; downward skew still missing. This is not a satisfactory joint calibration. Full-M2/equity/stress checking draws are new, not a reused old-state cache; same market snapshot, not market OOS. Exact 24-mode retained Gaussian marginal used;12 irrelevant orthogonal modes marginalized, not a three-state Markov approximation. Independent fresh-kernel, Gaussian-law,9 gradient,8 seed, input/source and VIX-invariance checks PASS. Final sampling/grid/normalizer gates unresolved. Candidate saved separately, not adopted; no Lean/infeasibility theorem.

Next proposed structural experiment: one global stock/stress leverage parameter, no new stress speed. [STOCK_DRIVER_ROTATION_PROPOSAL.md](STOCK_DRIVER_ROTATION_PROPOSAL.md) derives a valid rotation preserving the joint variance-factor law and hence VIX; it also weakens old stock/M2 correlations and invalidates unchanged reuse of the current independent-stress pricing controls. Proposal only, not implemented or calibrated. Do not silently treat a lower price loss as improved IV shape or repeat the same curve-only search without a new hypothesis. Earlier precision note's optimizer prohibition was superseded by the later explicit user authorization recorded in PROTECTED_CURVE_PLAN.md.

## 28 October 2026 — protected two-speed SPX precision phase

User authorized improving SPX precision with the successful two-speed VIX-only model fixed, then a protected fit of only nine off-window xi levels, with stress/equity correlation a later possible experiment. All model parameters remain unchanged in this phase: no joint candidate adoption, correlation, factor, readout or floor change. All1751eligible SPX quotes repriced,1294primary at all seven retained maturities; all217VIX quotes/four forwards remain protected. See [TWO_SPEED_PRECISION_RESULTS.md](TWO_SPEED_PRECISION_RESULTS.md), [TWO_SPEED_PRECISION_PLAN.md](TWO_SPEED_PRECISION_PLAN.md), `runs/m2_two_speed_precision_002`, `runs/m2_two_speed_precision_conditional_002`, and `visualizations/two-speed-precision027`.

Implemented an independent-stress auxiliary stock, likelihood-weighted differences, terminal/dynamic hedge controls, then an exact discrete one-step conditional Black-value martingale correction and independently scrambled OU-bridge reference integration. Checking median SE/scale improves about8–29x, from4.68/15.03/12.17/12.76/13.08/5.74/5.01 to0.579/0.550/0.903/1.153/1.239/0.435/0.174. Only1/1294primary quotes passes the existing sampling screen. The nine-level price optimizer was therefore NOT launched. Broad frozen-model SPX smile-shape mismatch remains visible; it is not a model-family infeasibility theorem.

Remaining bottleneck is the512-path full-M2 correction. More cheap stress paths alone cannot clear the gate. Pilot equal-component budget projections give median expensive counts about1.8k–93k across maturities and worst quote about1.2million, excluding timestep/normalizer errors; these are not guaranteed budgets. Next implement streaming fresh full-M2 corrections with adaptive/multilevel allocation and coupled stock-grid checks, rather than expanding the full cache blindly. Then free ONLY indices0,5,6,9,11,12,13,14,15; keep floor, both stresses and VIX-sensitive indices1,2,3,4,7,8,10 fixed. The exact VIX parameter map is invariant to those nine changes. Existing forward accuracy is preserved, not upgraded to population equality. Independent target/reference/conditional/RNG/covariance/full-horizon derivative and finite-law control checks PASS. Checking streams reused; not market OOS, full convergence clearance or a Lean proof. Current chat can continue; these notes are the restart handoff if a new context is chosen.

## 28 October 2026 — two-speed joint SPX/VIX refit completed

User resumed SPX and chose the selected two-speed stress model. Restored seven SPX maturities29/57/119/211/302/666/1030, cutoff3(1294primary;1751eligible including diagnostic wings), all217VIXquotes/fourforwards. Freed all16xi levels, floor and six existing stress parameters(23total); original M2 readout/premia/dynamics fixed. See [TWO_SPEED_JOINT_RESULTS.md](TWO_SPEED_JOINT_RESULTS.md), [TWO_SPEED_JOINT_PLAN.md](TWO_SPEED_JOINT_PLAN.md), `runs/m2_two_speed_joint_001`, `visualizations/two-speed-joint026`. Checking SPX loss5467.10→618.74, VIX0.05008→37.9711. Joint stress weight collapses0.92284→0.01943; slow branch effectively inactive. Forward calibration gap1.24e-8, checking max0.01782. Unsatisfactory overall: VIX fit badly worsens, all SPX sampling gates fail. Coupled64/128/512checks retain broad trade-off and unresolved timestep/sampling error. Source/input/56seed/23derivative/NumPy checks pass. This is one bounded local staged fit, both stages iteration-limited, with old readout frozen; no structural impossibility claim. Two-speed VIX-only and joint candidates preserved separately; no production baseline replacement. SPX scope is active for research; do not silently resume the earlier pause. Improve SPX numerical precision before further calibration conclusions.

## 28 October 2026 — fixed geometric four-speed stress bank finalized

Authorized bank experiment completed with half-lives7,28,112,448days fixed before fitting (ratio4, calendar-based endpoints). Four independent normalized stress responses plus frozen M2; globally shared positive weights and amplitudes. All217nativequotes/fourforwards, seven free active xi levels and fitted floor retained. Checking option loss 0.050081→0.062160; forward gaps/SE and all-maturity wings are in [STRESS_BANK_RESULTS.md](STRESS_BANK_RESULTS.md). See [STRESS_BANK_PLAN.md](STRESS_BANK_PLAN.md), `runs/m2_stress_bank_002`, `visualizations/stress-bank025`.16 fitted parameters versus14 previously; this is not a nested comparison. Cached terminal coordinates reused with new additional nodes; checking partly reused and not market OOS. Source/likelihood/scalar audits PASS; no population-forward certificate, SPX work, baseline replacement or Lean theorem.

## 27 October 2026 — shared fast/slow stress, all maturities jointly

The user explicitly requested one common model rather than special treatment of 27 days. The proposed isolated first-expiry diagnostic is superseded as the next fitting step. Completed the joint two-speed pilot: see [TWO_SPEED_RESULTS.md](TWO_SPEED_RESULTS.md), [TWO_SPEED_PLAN.md](TWO_SPEED_PLAN.md), `runs/m2_two_speed_001`, and `visualizations/two-speed024`. Checking option loss 1.725555 → 0.050081; all 217 quotes/four forwards retained. Four training forward equality constraints, seven free active xi levels and fitted floor. Checking gaps and SEs are reported, not declared exact. Reused old/slow checking states plus new fast normals; not wholly fresh. Previous continuous Beta measure not simultaneously refitted; two-speed extension tests the empirically motivated simpler mixture. Numerical audit PASS; SPX paused and both baselines unchanged.

## 27 October 2026 — teaching notes and first VIX maturity diagnosis

Saved CONTINUUM_EMPIRICAL_REDUCTION.md: equations, Bergomi attribution, the Beta-to-two-atom limit, bounded empirical reduction and why the general integral remains useful for future market snapshots. The reduction concerns mixture/readout complexity, not stochastic state dimension; only one stored market snapshot was tested.

Saved VIX_27_DAY_DIAGNOSIS.md with reproducible read-only diagnostics in runs/m2_continuum_reduction_note_001. The27day right wing is too low while central calls are too high. In the matching boundary reference, the strike100 payoff-weighted stress variance share is98.71%; stress dominance and shared response parameters are a leading shape mechanism. Exact-forward necessary second moment468.99933 exceeds selected467.44959. A27day-only constrained calibration is a proposed discriminating experiment, not an authorized scope change or a completed fit. All217quotes/fourforwards, the working candidate and SPX pause are unchanged.

## 27 October 2026 — continuous measure, forward constraints and multilevel repair

Completed VIX-only continuous M2/stress Beta-measure experiment003, after rejecting002's noisy55-day forward calibration. The model is unchanged: normalized G^(1-u) times the normalized OU stress exponential of strength b*u, integrated against one Beta measure. Four forward equalities use32768-path boundary plus8192-path paired difference; checking uses65536 NEW boundary paths plus16384 reused difference paths. Do not claim wholly fresh validation or exact population equality.

Best bounded fit remains near kappa=.01; a broader measure is retained as a training-selected diagnostic, with b=5 at its bound. No clear improvement established over the two-atom reference. All217 quotes, four maturities, seven free active xi levels and fitted floor preserved. SPX paused; original baseline and previous working VIX candidate unchanged. See CONTINUUM_RESPONSE_RESULTS.md, runs/m2_continuum_response_003 and visualizations/continuum-response023.

# Copy-and-paste prompt for a new conversation

**27 October: tighter VIX forwards plus fitted-floor experiment completed.**
See [FORWARD_FLOOR_RESULTS.md](FORWARD_FLOOR_RESULTS.md),
[FORWARD_FLOOR_PLAN.md](FORWARD_FLOOR_PLAN.md), run `runs/m2_forward_floor_001`,
plots `visualizations/forward-floor021`. All16xi intervals retained, seven active
levels independently fitted in both comparisons; SPX remains paused. Four training
forward equality constraints hold within1.9e-10VIXpoint. Fresh32768-point/expiry
comparison: previous mixture loss10.2653; tight forwards/fixed floor12.7780;
tight forwards/free floor1.7830. Largest fresh forward gap falls0.9001->0.02694,
with sampling SEs reported. Releasing the floor chooses12%volatility, at the
search ceiling, b3.1710,w0.66150,half-life104.215days. Best tested, not a global
optimum or settled floor estimate.55-day lowest put is close to market,27-day
left end improves; lowest long-dated puts become too cheap. Wing IV MAE8.605->8.445
on average;27-day far-right level worsens.84/217quotes within bid–ask.
Five small negative fresh native estimates and missing IVs remain, unaltered;
raw nonnegative put-payoff diagnostics are saved separately. No clipping/smoothing.
Six independent price replays,11Jacobian directions,32fresh seed regenerations,
exact budgets,floor support,291source/fiveinput hashes pass. This is numerical
verification, not population-forward equality or complete calibration. Current
working VIX candidate is the free-floor experimental model in this run; original
selected model unchanged. No SPX computation, new component or Lean theorem.
Historical checkpoints below keep their original context.

**26 October: point 2, two-response VIX-only comparison completed.**
See [MIXED_RESPONSE_RESULTS.md](MIXED_RESPONSE_RESULTS.md),
[MIXED_RESPONSE_PLAN.md](MIXED_RESPONSE_PLAN.md), run `runs/m2_mixed_response_001`,
figures `visualizations/mixed-response020`. Both product and additive mean-one
mixtures fitted seven independent VIX-sensitive levels of the original16-piece
curve. The nine off-window entries remain unchanged; all217quotes/four forwards
remain, SPX paused, no clipping/smoothing. Native-price winners and descriptive
shape selections were frozen before new16384-point/expiry numerical checking.
The additive M2/stress price fit improves fresh native loss20.0052->10.2960 versus
the prior fitted single factor (paired change-9.7093, jackknife SE0.2195), wing
IV MAE30.81->8.51points. Long-dated slopes now have enough steepness; short-dated
slopes overshoot.49/217quotes inside bid–ask, residual forward/ATM errors remain.
The additive fit gives99.57% weight to the new component, b at5bound, half-life
6.588days; it largely replaces the old response in this exploratory readout.
Two short variance pieces reach0.20. No acceptable full-smile/joint fit or model
adoption is claimed.72implementation checks and independent14price replays,
32fresh scramble regenerations, exact budgets,285source/fiveinput hashes pass.
Sampled GH/time VIX refinement max0.000118point for additive winner,0.0332 across
all candidates' extreme probes. Same market targets; no market OOS or Lean claim.
Preserve baseline, scope and historical records. Further timescale/readout changes
remain proposals; no second response timescale has been implemented in this run.

**26 October: VIX-only free-curve comparison completed.**
See [VIX_ONLY_CURVE_RESULTS.md](VIX_ONLY_CURVE_RESULTS.md),
[VIX_ONLY_CURVE_PLAN.md](VIX_ONLY_CURVE_PLAN.md),
`runs/m2_vix_only_curve_fit_002`, figures `visualizations/vix-only-curve019`.
The original curve has **16 levels**, already fitted in the earlier maturity-curve
study; recent factor trials had frozen them. Seven pieces affect the four VIX
averaging windows and now vary independently. The nine off-window pieces retain
their old values because these VIX targets have zero sensitivity to them.
No reversion to three levels or four broad scalers. SPX remains paused.
Fresh native VIX+forward losses: baseline27.3491, curve-only20.5245,
fixed strong factor+curve20.9500, fitted factor+curve19.8086. Additional fitted-factor
gain versus curve-only is -0.7159 with paired jackknife SE0.5476. Wings remain too
flat; one short curve piece hits its0.20 bound and fitted half-life hits365days.
All217native quotes and four forwards retained; no smoothing/clipping. Only21quotes
inside bid–ask for the joint fit; no acceptable full-smile calibration or model
adoption claimed. Baseline preserved. Eight candidate/stream price replays,
nine Jacobian directions, exact budgets,282source/fiveinput hashes pass; sampled
GH/time refinement differs by less than4.6e-6VIXpoints. Numerical checks use the
same market quotes. No SPX pricing, new bank/readout or new Lean theorem.
Historical checkpoints below retain their original context.

**26 October: current phase is VIX-only research and calibration planning.**
SPX fitting/repricing is temporarily deferred by explicit user direction.
Read the new section at the top of `STOCHASTIC_VOL_OF_VOL_RESEARCH.md` and
`runs/m2_vix_design_research_001`. All three Othmane prepared PDFs and the14October
email figures are readable; no original-PDF request is needed. The later cs_doc
still reports too-flat2–5month VIX wings, fits a narrower strike range and uses
an approximate VIX estimator with a nested check pending. Its Gamma modulates
bank noise, unlike our final readout multiplier; no matching-code repricing is claimed.

New model-independent finding: exact-rational static-payoff certificates show
that the frozen xi curve's E[VIX²] budget is too low to fit all native VIX
bid–ask quotes and forwards at55/118/153days. This persists with forward tolerance
plus/minus0.5points. `moment_certificate_audit.json` passes12 global payoff/price
checks; it is exact computation for archived decimal data, not Lean or a market
arbitrage claim. Necessary exact-forward lower bounds are556.673/648.626/695.295,
versus old budgets536.027/580.804/610.933. The27day bound does not exclude baseline.

Next comparison should allow VIX-window xi levels and shape parameters to move
together. First re-evaluate the current after-cap layer fairly; then prefer a
mixed-response readout/independent mixed-standard-Bergomi benchmark. A continuous
Gamma-driven new bank is the more direct Othmane analogue but needs a new
conditional pricer. No new reservoir calibration, SPX run, jump adoption or
selected-model replacement has been performed in this research turn.
All217VIX quotes and four forwards remain; old SPX membership/cutoff3 are retained
for the later return. Earlier checkpoints below remain historical.

**26 October checkpoint: independent amplitude AFTER the original cap completed.**
See [POSTCAP_VOV_RESULTS.md](POSTCAP_VOV_RESULTS.md),
[POSTCAP_VOV_PLAN.md](POSTCAP_VOV_PLAN.md), run `runs/m2_postcap_vov_001`,
and figures `visualizations/postcap-vov017`. This is the authorized first
alternative; no mixture was implemented. One independent Q OU exponential
L=exp(gamma*Z-gamma^2*Var(Z)/2) multiplies G_base/M_base after the old cap.
The old normalizer, reservoir, coefficients, all16 xi levels and floor remain fixed.
Nine nonzero strength/persistence cases plus baseline have full VIX diagnostics.
Price selection gamma=.25/7d and shape selection gamma=1/7d were frozen before
full SPX repricing (all seven dates) and32 fresh joint37D VIX scrambles.
The strong case improves fresh measured wing slope at all four dates and lowers
the descriptive slope/curvature error0.44495->0.30875, but remains too flat.
Fresh native VIX+forward loss26.8525->29.7000; SPX loss43.7775->121.1470.
The mild case's training joint-loss reduction70.9777->69.8756 is smaller than
its paired jackknife sensitivity (change -1.1021, SE3.9335); fresh VIX price loss
does not improve. No accepted joint fit or new model adoption follows.
Five unit checks and the saved-price/source/seed/independent-integral audit pass.
Signed estimates and missing IVs are retained. No population convergence,
market OOS, new P calibration or Lean certificate is claimed.
Keep `m2_maturity_curve_fit_001/selected_model.json`; no experiment remains running.
The earlier checkpoints below are historical, not an instruction to rerun them.

**COMPLETED — dedicated stochastic spike-sensitivity trial, 25 October:**
Read `STOCHASTIC_VOV_RESULTS.md` and `STOCHASTIC_VOV_PLAN.md`.
Run `runs/m2_stochastic_vov_001`; figures `visualizations/stochastic-vov016`.
Selected baseline remains `m2_maturity_curve_fit_001/selected_model.json`.

One added independent Q OU state modulates qJ by exp(s*tanh(Z)); Z0=0,
stationary variance1. Strength ranges 1.5x/2x and half-lives7/30/90 days,
crossed with original cap and logarithmic continuation a=.5. All original
OU rates/grids/state/correlations/premiums/readout coefficients, all16 xi
levels and floor .00264 are exactly fixed. No curve refit, jumps or quintic.
All four VIX dates/four forwards and4096 samples per strike remain. All12
layer cases and both zero-strength references have full VIX diagnostics.

Wing-level training selection picks2x/90d for BOTH caps, before new SPX
results or checking. Full actual SPX repricing covers all seven dates for
both selections and both cap references (four full comparisons,1511 OTM
quotes each, cutoff3/common-call/native model-forward parity preserved).
No full-SPX search over the entire12-point grid is claimed.

Main finding: closer wing IV LEVELS do not mean greater STEEPNESS. On NEW
joint37D scrambles, adding the selected state FLATTENS measured right-wing
slope at ALL FOUR VIX maturities relative to each same-cap no-layer reference.
New-stream wing MAE33.841 ->15.072/11.733, but near-ATM12.670 ->18.994/18.555.
Native VIX+forward loss27.597 ->272.754/482.047. Training SPX loss43.777
->561.511/1376.764; no-layer slow-growth reference112.992. Paired SPX changes
are positive in every leave-one-block comparison for both selected cases.
The frozen-curve grid does not give the desired joint improvement. Keep baseline;
this does not reject a different interaction or a recalibrated stochastic model.

Normalizer M at153d rises2.408x/3.074x for selected layers vs same-cap refs;
stronger raw response plus stronger normalization redistributes variance.
Targets E[VIX^2] stay identical because xi is fixed. New153d baseline moment
606.410 vs610.933, SE2.466; selected original-cap state608.135, SE12.257.
These are numerical diagnostics, not a convergence pass or resolution of the
historical second-moment discrepancy. No fresh SPX ensemble is claimed.

Three tests pass; audit replays22 price/SE/parity/loss records,270 archived
sources and all32 NEW complete joint mixture scrambles. Eight independent
4D oracles have max relative error5.63e-10. Representative normalizer GH48->96
changes <=5.25e-8, extra-state GH and time refinements are recorded separately.
Original GH32/GH48 and time8/stock512 mesh retained; extra Z usesGH128/GH192.
Exact OU stock updates stream reproducible shared innovations without an
extra all-time state cache. OU primitives reused; volatility/stock recomputed.
Nine signed SPX302 wing estimates are negative, all <.82 reported SE below0.
Separate VIX IV failures reflect market-forward inversion; no clipping/smoothing.
Eleven inspected figures have22 PNG/SVG exports. Baseline full price replay
max difference1.37e-9. No new accepted fit, Lean/P/global-optimum certificate.
All numerical jobs finished; finalization includes audit/environment/hashes.

**COMPLETED — cleaner cap comparison, 25 October:**
Read `CAP_CONTROLLED_RESULTS.md` and `CAP_CONTROLLED_PLAN.md`.
Run `runs/m2_cap_controlled_001`; figures `visualizations/cap-controlled015`.
Selected baseline stays `m2_maturity_curve_fit_001/selected_model.json`.
No new accepted joint fit or dedicated latent layer. Jumps/quintic excluded.

Floor .00264 and every off-window curve piece are EXACTLY fixed. Four broad
VIX-window adjustments average overlaps, with maximum15% excess-variance change;
selected changes are only about1--3%. Same OU rates/geometric grids/state/
correlations/Q scenario/score coefficients, beta_J>0,1511 original OTM quotes,
all seven SPX/four VIX maturities/four forwards, cutoff3/common calls and4096
paths per strike. No clipping/smoothing. Nine cap shapes: L increases0/.375/.75
crossed with logarithmic continuation a=0/.25/.5, k=2. Local fits hit the
common-call safeguard45.509602; several also hit +.5-point forward bounds.
Those screens are implementation choices; no family/global-optimum conclusion.

TRAINING selection frozen before NEW draws: shape03 (L+.375,a=0) and shape02
(L unchanged,a=.5). Wing IV MAE33.646->27.717/26.664; near-ATM12.389->14.921/15.180.
Actual full seven-maturity SPX loss43.777476->87.400641/98.795211. Paired SPX
loss-change SE5.959/6.734; all eight leave-one-block comparisons worsen.
Native joint70.977709->114.299387/125.738973. Keep the unchanged baseline.

Genuinely NEW8x512 VIX scrambles per date, seed base2026102511, same proposal/
frozen pilots/SAME in-sample market quotes. Native VIX+forward baseline26.986598;
shape03/02 give26.845187/26.914668: changes-.141411/-.071930, paired eight-net
SE.203733/.228777. Aggregate gains are not precise. New wing MAE33.866->
27.963/26.919; near-ATM12.590->15.183/15.435. No market OOS or fresh SPX check.
New baseline VIX153 moment608.632 vs610.933, SE2.294 (~1SE); historical deficit
is preserved, not cleared. Baseline139/1511 inside bid/ask vs94/87.

Three tests pass. Audit replays15 scores,139 fit trials,270 frozen-source
snapshots and32 complete new scramble/state/likelihood regenerations, with
seed disjointness against recorded historical pilots/production/diagnostics.
Nine PNG/SVG plots cover all agreed maturities, signed errors, forwards and tails.
Original fresh-cache setup failed BEFORE fresh pricing, after fits/SPX/draws
completed. Current driver is corrected. `resume_cap_controlled.py` retained all
fits/SPX/fresh draws and resumed only VIX pricing. Original failure/sources and
separate correction snapshots are preserved. All numerical jobs are finished.
These tested gentler brakes improve wings but do not give the desired joint fit.
Further continuous response-layer work remains proposed, not automatically run.

**COMPLETED — cap/tail profiles and continuous response diagnostics, 25 October:**
Read `CAP_TAIL_RESULTS.md`, `CAP_TAIL_MATHEMATICS.md` and
`STOCHASTIC_VOL_OF_VOL_RESEARCH.md`. New run `runs/m2_cap_tail_001`;
plots `visualizations/cap-tail-014`. Selected baseline remains
`m2_maturity_curve_fit_001/selected_model.json`; no new joint fit is adopted.
Current user priorities override archived queues: cap first, then continuous
stochastic vol-of-vol; exclude jumps and quintic proposals. Existing OU rates,
geometric grids, initial state, correlations and Q scenario remain fixed.
All seven SPX/four VIX dates, original 1511 OTM inputs/four forwards, cutoff 3,
common calls, market/model-forward parity and all 4096 samples remain active.

Nine fixed-curve VIX profiles (baseline, L +/-0.75, k=1/4, algebraic cap,
logarithmic continuations a=.25/.5 and separate channel caps) rebuild Q M.
Baseline prices replay within 1.37e-9. Raising L steepens the right wing but
worsens native VIX+forward loss 27.166865 -> 34.540828; log a=.25 gives
27.153015, algebraic 26.832740. None solves finite-strike shape. Saved checking
maxima exceed quoted strikes; likelihood-weighted probability/payoff tables
show why raw hits do not establish adequate tail prices. Upper cap bounds
about 937--1003 are not themselves attainability proofs. Payoff-relevant
conditional midpoint responses are compressed to about .57--.66.

Two explicitly VIX-focused curve/floor fits improve native VIX+forward training
to 20.044859 / 18.674677, reused separate checking to 20.653505 / 19.124019,
but actual seven-maturity training SPX losses rise 43.777476 -> 1656.376316 /
318.635148; common-call VIX criteria also worsen. NOT accepted joint fits.
Both hit implementation floor/curve bounds; VIX-only data leave off-window
variance pieces weakly identified, so SPX costs cannot be attributed to cap alone.
Next focused VIX diagnostics should freeze off-window pieces and floor rather
than repeat an unconstrained full-curve VIX-only search.

After cap profiles, two existing-cluster bounded spike-response mixtures worsen
VIX loss to 87.177224 / 117.866292. No SPX or joint success claimed for these.
Dedicated continuous latent response state is mathematically specified, with
parameters/randomness/cache/cost/integrability/martingality discussion and primary
literature; it is not implemented. Current nonlinear M2 already has state-dependent
Ito vol-of-vol; the proposal adds control of that response and persistence.

Selected uncapped first Gaussian multiplier moments actually remain finite
(stationary quadratic eigenvalue .897020488), while second moments diverge at
long dates (critical power 1.114801739). Do not claim uncapping necessarily makes
the selected M infinite, or that first-moment existence proves stock martingality.
GH32->48 and time8->16 representative extreme-state changes are at most
2.27e-6 / 3.40e-8 VIX points; no surface convergence certificate. Baseline VIX153
second-moment deficit 606.208 vs 610.933, SE1.635, persists; higher cap does not
clear it. Saved-draw audit passes 17 score records and 259 frozen-source checks.
Independent Python equation/quadrature/stock tests pass (direct runner; pytest
unavailable). Eleven PNGs plus SVGs cover all agreed maturities and diagnostics.
All numerical jobs are finished. Preserve baseline/historical sources, failures,
signed prices/unavailable IVs; no new Lean or P-compatibility certificate.

Copy the text below the separator into a new chat with this repository open.
The repository records, rather than the old chat history, are the evidence.

**COMPLETED — maturity-aligned curve calibration, 25 October:**
`runs/m2_maturity_curve_fit_001` completed after the user pause/resumption.
Report `MATURITY_CURVE_FIT_RESULTS.md`; plots `visualizations/maturity-curve-013`.
Sixteen curve pieces retain old 91.25/273.75-day boundaries and add all retained
SPX expiries and VIX expiry/window ends; previous curve exactly nested. 24 free
parameters across the exploration. Fixed OU dynamics/rates/geometric grids,
Q scenario/correlations/initial state; original OTM scope/cutoff3/seven SPX/four
VIX maturities, all4096 paths/common calls/GH32/GH48/time8/stock512 unchanged.
Hard xi_i>=floor, floor volatility 5--6%, beta_J>0; no xi monotonicity restriction.

Three full proposals, one accepted. Accepted step changes curve/floor only;
subsequent all24-free proposals improve VIX further but fail full joint safeguards.
Selected readout/cap stay unchanged despite being searched. Actual training:
canonical joint95.117702->89.320446; native joint84.982024->70.977709;
SPX49.518173->43.777476; native VIX+forwards35.430506->27.166865 (-23.32%).
Forward component alone worsens0.838006->0.890564; VIX smile still poor at high
strikes, with nonuniform near-ATM effects. Do not call this a satisfactory fit.
Separate reused VIX checking (selection frozen first): native VIX+forwards
35.862036->27.523027 (-23.25%), paired eight-net change-8.339009,
jackknife SE0.191220. Common-call change-0.601058 has SE0.781157, so its small
gain is not precise. Checking streams were used before: NOT pristine holdout or
market OOS. No SPX separate checking repricing in this run.

Selected floor0.00264 =>5.13809% minimum instantaneous volatility;
beta_J0.369729973>0. No negative native prices or unavailable common-call IVs
in training / VIX checking. SPX training still has1 monotonicity/66 convexity
flags; VIX has0/0. Raw signed values remain. Selected VIX153 second moment
606.208 vs target610.933 (SE1.635): rare-tail/numerical issues are not cleared.
Ten focused tests pass; audit replays70 local trials, all full acceptance decisions,
prices/SE/parity/loss/policy and74 sources/inputs. Eight PNG plots visually reviewed,
16 PNG/SVG hashes verified (`validation/plot_review.json`). Source snapshots and
original pause record preserved; continuation restores completed baseline and
reproduces saved proposal0 parameters, without new OU cache generation.
`selected_model.json` is the exploratory selected configuration. Not global
optimality, identified curve, price convergence or a new Lean/P certificate.
All processes are finished. Preserve prior runs and full agreed scope. Do not
restart archived searches automatically; remaining questions concern VIX smile
shape, forward tradeoffs, local-model safeguards and SPX estimator accuracy.

**25 October — explicit forward-variance hard constraint:**
Every new curve level must satisfy xi_i>=variance_floor, with the existing
variance_floor>=0.05^2 policy. `forward_variance_constraints.py` now supplies
this guard for any number of pieces; `calibration_policy.validate` calls it
for the current three pieces. Nonmonotonic levels and equality are allowed.
The fit launcher snapshots the new dependency. Fourteen focused tests pass,
including rejection of a below-floor value in each of14 positions and the
current calibration's relational check independently of its tighter box bounds.
Evidence: `runs/m2_forward_variance_constraint_001/validation.json` and tests.log.
No recalibration or maturity-grid expansion was launched. Historical results
and frozen sources are unchanged. Current policy/launcher source revisions are
intentional; replay historical source hashes against the preserved snapshots.

**Latest completed work — 5% volatility minimum and point1 acceleration:**
User requires minimum annualized volatility around5%, then acceleration point1.
For ALL NEW calibrations use `calibration_policy.py`: variance_floor in
[0.0025,0.0036] (5--6% floor volatility), beta_J>0; all eleven remain free.
The6% upper endpoint is the implementation choice keeping the floor free in
its stated reasonable range. Historical bounds/entrypoints remain reproducible.
Since xi>=floor and G/M>0, instantaneous sqrt(v)>=sqrt(floor)>=5%.
This is not an option-IV bound. Selected floor is0.00254, volatility5.039841%.

Completed pilot `runs/m2_joint_local_fit_001`, report
[JOINT_LOCAL_ACCELERATION_RESULTS.md](JOINT_LOCAL_ACCELERATION_RESULTS.md), plan
`POINT1_ACCELERATION_PLAN.md`, figures `visualizations/joint-local-012`.
`joint_price_local.py` approximates BOTH SPX and VIX prices/forwards at a true
anchor; native puts include MODEL-forward derivatives. `joint_local_search.py`
uses normalized SLSQP, a small radius-squared reserve and the unchanged strict
full-price guards. Accepted full checks refresh all eleven Jacobian columns.
Ten new/existing tests pass;119 local trial scores, saved full means/SE/parity/
losses/acceptance decisions,68 sources/inputs and selected-policy/config replay.
All six PNG figures visually reviewed; PNG/SVG output hashes verify.

Measured119 optimizer evaluations in0.3236s total (about2.54ms/evaluation), versus
historical33.35s nonlinear-VIX trial median. Approximate versus nonlinear work:
NOT an end-to-end speed factor or identical trial accuracy. Full baseline363s,
two full checked proposals377s/379s, whole pilot1143s. Three full pricing calls
still dominate. Both proposals accepted: common joint95.14831 ->95.11770;
native joint85.03863 ->84.98202; native VIX+forward35.51359 ->35.43051. Gains are
small. SPX and forward components individually worsen slightly; all four guarded
joint/VIX criteria improve. Selected beta_J0.36972997 remains positive.

Same1511 OTM inputs/four forwards/cutoff3/maturities, all4096 paths, fixed OU/Q/
initial state/geometric grids and GH32/GH48/time8/stock512. No coarsening or block
optimization. No new independent checking priced: overlays use actual TRAINING
prices. Selected training has0 negative native estimates and2 VIX55 missing IVs;
this does NOT clear the prior fresh-stream SPX negative-price failures or VIX
wing/forward misfit. No convergence/global optimality/P-compatibility claim.

Pilot driver `fit_joint_local.py` preserves its prior003 start for reproducibility.
Future continuation can restore THIS run's selected full-price/Jacobian checkpoint
with `parity_checkpoint.load_result` and build `JointPriceLocal` directly, avoiding
another identical baseline simulation. Preserve old runs and immutable sources.

LATEST COMPLETED PARAMETER ADJUSTMENT (24 October):
`runs/m2_parity_parameter_fit_003` is complete. Read
[PARITY_PARAMETER_FIT_RESULTS.md](PARITY_PARAMETER_FIT_RESULTS.md) and
[all-maturity plots](visualizations/parity-parameter-fit-011/README.md).
All eleven agreed parameters were free, beta_J>0; same OU/Q/initial state/grids,
1511 original OTM quotes/four forwards, cutoff3/retained maturities, GH32/GH48/
time8/stock512, common-call simulation and uniform4096 paths. One full true step
accepted; two subsequent20-evaluation local SLSQP cycles found no reserve-feasible
proposal. Budget termination is NOT convergence or evidence of model infeasibility.
The final candidate has beta_J0.34972997, cap_sharpness2 and floor0.0025 (the latter
two at the existing bounds). Starting point remains the previous positive-shift fit.

Training common-call joint98.1192 ->95.1483; native joint85.5281 ->85.0386.
Parameters frozen BEFORE checking. NEW independent4096 samples, same frozen
training-pilot control functions: common-call joint96.4690 ->92.2402; native joint
82.8286 ->81.2319. Original-OTM VIX-plus-forward36.1043 ->35.9254 (only0.50%);
paired eight-block loss-change SE0.0454 versus change-0.1789, descriptive only.
High-strike VIX IV errors remain large/slightly worse. Forward errors after:
[+0.3644,-0.5400,+0.6495,+0.0754], much larger than sampling SE at most dates.
This is a small numerical improvement, NOT an acceptable final smile fit.

Fresh checking exposes8 negative selected native SPX prices (57/302 days), hence
8 missing SPX IVs; two VIX55 IVs also unavailable due model/market-forward mismatch.
Do NOT carry the old ensemble's no-negative/fully-defined-SPX-IV result over to
new draws. No clipping/smoothing. Means/SE/parity/loss/acceptance and source/seed/
frozen-control replay passed; five focused fitting tests and production derivative
screens passed. Twelve representative selected-state GH32/GH48 VIX probes differ
by at most7.9e-7 points; NOT whole-surface/time/stock convergence. All3250 frozen input files verify byte-for-byte; all six PNG figures were visually
reviewed and PNG/SVG output hashes verify. Records are saved alongside the run.

Run001's first true canonical-improving proposal worsened native VIX and was
rejected. Run002 supplies the validated full baseline checkpoint and a local warm
point; both were explicitly stopped before checking, not adopted calibrations.
Run003 fully reprices/accepts that warm point; subsequent native-OTM-VIX-priority
local guidance guards BOTH joint/BOTH VIX criteria. Reserve constraints are extra
numerical margins, not a user requirement; failure to meet them does not exclude
smaller improving directions. Review solver scaling/feasible starts and SPX
estimator stability before treating stalled local progress as a family limitation.
Preserve all run/prototype/checking records and fixed-source manifests.

LATEST COMPLETED ESTIMATOR EXPERIMENT (24 October): read PARITY_UNIFORM_RESULTS.md
and visualizations/parity-uniform-010/README.md FIRST. The user's final decision
is one CALL payoff using existing OTM inputs converted by parity, and ALL4096
paths at EVERY strike. Adaptive-prefix allocation is superseded. Run
m2_parity_uniform_001 is complete: fixed old/selected parameter comparison, not
another calibration. Same 1511 options/four forwards/all 11 maturities/cutoff3;
OU dynamics/grids/Q prices/initial state and GH32/GH48/time8/stock512 remain.
VIX has 8 independent 512-point scrambles. Four old SPX57 gaps disappear, all
selected native OTM prices are positive; two VIX55 common-call IVs remain below
market-forward intrinsic because model E[VIX] is 0.5167 below market. Preserve
both canonical and native OTM objectives and actual model forward outputs.

Same-estimator common-call losses 1941.75/98.12 and native OTM 1724.67/85.53
(old/selected) are NOT comparable to old 512 controlled losses. A rare old SPX302
stock ratio 56.63 (selected 5.58) dominates the huge baseline wing score; paired
loss-change jackknife SE 1687 vs reduction 1844 prevents precise improvement claims.
Only 100/1511 native selected prices fit bid/ask; 1 SPX order/72 convexity failures
remain and many SEs are several quote scales. VIX shape checks pass on this sample.
Fourteen tests, source/price/SE/parity/loss replay, 73 source checks and 1674 immutable
array checks pass. All jobs complete. Cache/pilot paths retain adaptive names but
active allocation is uniform4096. Queued 001/002 were interrupted before pricing.

reusable_ou_initial.py tests exact exp(-Lambda*t)*(newX0-oldX0) translation and
unchanged equity increments/importance weights/scores/components. Recompute G,
normalizer, stock and VIX after an initial-state change. Full initial OU states
are needed; the project's minimum-Mahalanobis mapping from C/Y/F is available.
This capability does not add initial state to the agreed 11 parameters. Future
calibration should retain those 11/beta_J>0, new common-call convention/all 4096,
report original OTM and forward errors, and preserve rare-tail uncertainty.
Do not resume old raw-fit/Q-scenario queues or claim another fit already ran.

LATEST COMPLETED EXPLORATION: read CONTROLLED_READOUT_EXPLORATION_RESULTS.md
and visualizations/controlled-readout-fit-009/README.md FIRST. The user asked
for actual broader exploration, not stopping at reproduction of the old fit.
`runs/m2_controlled_readout_fit_002` completes all eleven agreed free parameters
at fixed OU grids/correlations/projection powers/Q prices/C0=1,Y0=F0=0, same
full OTM cutoff-3 market scope. Three starts/49 evaluations, restored 512-path
SPX searching/checking and 8x64 VIX, conditional GH32/normalizer GH48/time8.
Continuous zero-mean controls avoid donor-switch jumps; the same GH32 integral
is accelerated and the selected checking prices also replay by direct GH32.
The direct startup run 001 was interrupted and retained; cache 001 is complete.
Search loss 256.33→83.53, checking 305.21→103.34; all aggregate families improve
and all eight leave-one-block changes favor selection. Beta_J=0.29972997 is
interior; floor 0.00248112 and cap sharpness 2.02455 approach upper/lower bounds.
Only 114/1511 prices are inside bid/ask; four checking estimates are negative,
two raw-payoff contracts have no hits, long SPX and VIX153 worsen, and all local
runs reached evaluation limits. It remains exploratory; P compatibility and
price convergence are unproved. Fourteen focused tests and independent source,
immutable-cache and saved-score/SE replay pass. All numerical jobs are complete.
`selected_model.json` in the run is the self-contained positive-shift candidate;
`frozen_stage3()` still loads the historical zero-shift anchor. Preserve both.
Next bounded questions: profiles of newly hit floor/cap-sharpness bounds and
candidate-appropriate continuous controls, with fair same-estimator comparisons.
Neither extension is launched. The older raw-fit/Q-scenario queue below is
historical and must not override this completed checkpoint or the user's scope.

Latest steering after the broad fit: start from the previous all-maturity
parameters (only beta_J changed to its required positive lower bound 0.01),
retain that feasible start as an incumbent, and repair strike jitter while
keeping useful variance reduction. Read STRIKE_CONTINUITY_RESULTS.md first.
The saved-draw continuous-control experiment removes all observed strike-order
violations at both 512-path ensembles and preserves much of the precision,
but still has convexity failures. Its affine comparison has negative/order
failures. Neither is a fully qualified pricing repair. Final retained diagnostic
is runs/m2_continuous_strike_controls_003; no new calibration was launched.
Do not restart the low-budget raw-estimator fit as though this follow-up had
established comparable pricing accuracy or solved every shape issue.

LATEST COMPLETED ATTEMPT, 24 October: read
`FULL_READOUT_FIT_RESULTS.md` and `visualizations/full-readout-fit-008/README.md`
alongside the ledger before following older queue instructions below. The user
confirmed eleven free parameters with unchanged OU architecture, Q prices and
initial state: five score coefficients, two cap parameters, three curve levels
and variance floor. `runs/m2_full_readout_fit_001` is complete. Positive beta_J
selected 0.01 at the lower bound. A new raw positive importance estimator omits
strike-dependent controls; old parameters were repriced with it. Checking loss
1044.10→906.61 is not a reliable gain: descriptive jackknife SE 876.70, one
leave-one-block reversal, worse VIX options/forwards and only 61/1511 interval
hits. Zero-payoff tails and quadrature/normalizer accuracy remain unresolved.
All eleven maturity panels and parameter/forward/price-error plots are saved.
Six focused tests and independent saved-score/source/seed/monotonicity replay
pass; no new formal certificate or P compatibility claim exists. Keep the new
candidate exploratory, preserve all historical runs, and do not repeat this
completed search or automatically launch the older Q-scenario queue.

---

Please continue the important M2 SPX/VIX option-fitting investigation in:

`/Users/matteogambara/Desktop/code/PathDependentVolatility_ESN`

We have spent substantial time diagnosing numerical instability before fitting.
I want to preserve that knowledge and move toward a defensible market fit,
with bounded experiments and explicit decisions rather than an endless sequence
of variance-reduction ideas. Please begin actual work after reading the records.

First check the branch, git status and repository instructions. Expected branch:
`codex/m2-option-pricing-benchmark`. Committed handoff base is `df27b4a`;
the current joint-fit/Lean additions may be uncommitted. Numerical sources are
captured in the completed curve/state/refinement run JSONs. Current Lean
proof/compiler/axiom evidence is under `formalization/audits/attempt_006/`;
attempt 003 remains the historical ten-theorem checkpoint.
Do not revert unrelated user changes:
there was a pre-existing edit in
`matteo/exogenous_reservoir_eq_m2/experiments/physical_stability_001/RESULTS.md`.

Read these files, in order:

1. `matteo/model_comparison/option_pricing/RESEARCH_LEDGER.md` — consolidated
   attempt index, findings, rejected ideas, preserved decisions and current queue.
2. `matteo/model_comparison/option_pricing/active_scope.json`.
3. `matteo/model_comparison/option_pricing/ALL_MATURITY_FIT_RESULTS.md` and
   `ALL_MATURITY_FIT_PLAN.md` — latest user-authorized all-observation scope,
   actual recalibration and all-seven-maturity model pictures. Then
   `JOINT_STATE_REFINEMENT_RESULTS.md` for the preceding completed checkpoint;
   then `JOINT_CURVE_FIT_RESULTS.md` for the preceding curve-only checkpoint.
4. `matteo/model_comparison/option_pricing/formalization/CORRESPONDENCE.md` —
   thirteen checked Lean theorems, numerical mappings and unproved bridges.
5. `matteo/model_comparison/option_pricing/CURVE_PROFILE_RESULTS.md` and
   `NUMERICAL_STABILITY_RESULTS.md` — why the apparent fit gain reversed.
6. `matteo/model_comparison/option_pricing/CALIBRATION_PROGRESS.md` and
   `matteo/exogenous_reservoir_eq_m2/documentation/CALIBRATION_PLAN.md` —
   common objective, input/state treatment, fitting and evaluation protocol.

Consult the ledger's linked reports and AUDIT.md as needed. Older “next step”
paragraphs are historical when later experiments already tested those ideas.

Preserve these decisions:

- LATEST parameter constraint, 24 October: enforce **beta_J > 0**, where
  `beta_J = spike_shift` enters `q_J*(F-beta_J)^2`. Keep the spike channel
  active (`q_J > 0`). The current loader defaults to beta_J=0 and the completed
  runs used zero; those remain historical baselines and do not satisfy this
  new requirement. The later full-readout attempt selected 0.01 but did not
  establish a reliable joint fit. Before future fitting,
  explicitly configure and record a positive shift or a search domain with a
  strictly positive lower bound, and reject nonpositive candidate shifts.
  The old illustrative [0,0.5] bounds no longer authorize zero. Rebuild the
  readout-dependent normalizer, stock paths and VIX exposures. Older directions
  to freeze the readout during Q-scenario comparisons apply only after a
  compliant positive-shift baseline is established. Preserve all archived
  zero-shift results; this documentation change launches no numerical job.

- LATEST user update, 23 October: include **29/57/119/211/302/666/1030 days in
  SPX recalibration and model plots**, cutoff 3, excluding only 1/1401/1765/2129.
  The user explicitly says no out-of-sample results and confirms recalibration
  is needed. Use all eligible original quote splits in this feasibility run:
  **1294 SPX options, 217 VIX options, all four VIX forwards**. Old reserved/
  diagnostic maturity and strike splits are historical metadata, not current
  exclusions. No model-only-for-four or market-only extra panels in new figures.
  Preserve raw labels/data and all historical results. `active_scope.json`
  records the new scope and points to the saved previous scope bytes/hash.

- Keep all four VIX maturities 27/55/118/153 days, now all 217 in-sample options and
  four forwards. Earlier two-maturity exclusions were superseded explicitly.
- SPX central priority: `abs(log(K/F(T))) <= 3*sigma_ATM_proxy(T)*sqrt(T)`.
  Cutoff 3 is provisional and user-selected: now 1294/1751 eligible SPX quotes
  at all seven dates, with all-observation ATM proxies. Historical 707/980
  training-only counts remain the prior scope. Keep wings as secondary diagnostics; representative
  contracts do not replace the full priority population.
- Freeze the current Stage-3 candidate 18/scenario 0, initial state and selected
  nuisance curve during numerical comparisons. Original vol7/ridge5 remains
  the SPX reference. Explicitly record any later fitting parameter changes.
- Recent SPX diagnostics use 512 steps per trading day. Sampling, stock timestep,
  conditional/time quadrature and normalizer errors are separate gates. Exact OU
  updates and passing unit tests do not certify option-price convergence.
- Retain negative controlled estimates and zero-hit/tail failures without clipping.
  Record inconsistencies between code and reports in AUDIT.md.
- Current feasibility run uses all retained market observations in sample;
  preserve original market split labels as metadata. Fresh simulation seeds
  check numerical stability on those SAME market observations. No market
  out-of-sample interpretation or reliable final model ranking exists yet.

LATEST COMPLETED RUN: `runs/m2_all_maturity_fit_001` prices all seven dates and
actually recalibrates all three curve levels on the full 1511-option in-sample
objective at fixed C0=1/Y0=F0=0, readout/structure/Q scenario 0. Three starts
converge in five reported evaluations each. Selected exact levels are
**[0.057620550549099724,0.07134085119966904,0.06321057443478936]**. The long
square-root variance input rises 18.89%→25.14%, substantially improving
666/1030-day SPX. Searching loss 300.94→253.44; separate numerical ensemble
on the SAME market quotes gives **351.36→308.64 (12.16%)**, with SPX/VIX/forward
components all improving. Only **103/1511** prices lie inside bid/ask and all
SPX precision screens still fail; no qualified fit or global capacity conclusion.
All four score records, full membership, source hashes, sample means/SEs and
new-seed metadata pass independent replay. An audit reader mistook raw VIX
samples for controlled price arrays; it was corrected without changing any
pricing/optimization output, and its failure log is preserved. All numerical
jobs are complete; do not repeat them. All 182 Python tests and thirteen Lean
identities (audit 006) pass. Final figures and their reading guide are in
`visualizations/all-maturity-fit-007/`, with labelled searching previews retained.
Next bounded fit: declared Q-risk-price scenarios one at a time, using this new
curve and the full 1294-SPX/217-VIX in-sample population. The older four-maturity
baseline/175-VIX/market-holdout descriptions below are historical, not current.
Older strict current-source audits see the deliberate active_scope version
change; old archived source hashes and historical outcomes are not rewritten.

Latest result: the user requested actual exploratory joint fitting now, rather
than requiring more standalone numerical experiments, and explicitly requested
Lean checks. This superseded the unlaunched nine-point/two-level/24-SPX proposal.
`m2_joint_curve_fit_002` optimizes all three declared curve levels on all 707
cutoff-3 SPX options, all 175 VIX options and four forwards. Readout/state/Q are
fixed. Other 273 SPX training quotes remain secondary diagnostics; held-out
market observations were not opened. Eight starts converge near
`[0.0600941554,0.0748505719,0.0383634764]`, using 104 curve evaluations.
Searching loss falls 637.32→333.68; all three fresh checks favor the selection.
Pooled fresh loss falls 597.90→347.10 (41.95%), with SPX 529.20→249.04, VIX
options 67.88→84.91 and forwards 0.82→13.14. Only 50/882 objective options
are inside bid/ask. VIX forward errors are +1.70/+0.73/+2.48/+1.89 points.
This is a persistent sampled improvement and a poor cross-market compromise,
not a qualified market fit. Only 27/882 pooled nominal SE screens pass, all VIX;
SPX sampling, stock-grid and VIX quadrature/normalizer gates remain unresolved.

The fit uses unchanged vol7 sampling with fixed anchor ridge5 controls, sharing
nearest same-side representative hedge gains across all strikes. It reprices
the stock for every curve and reuses only compatible G/M0 and primitive shocks.
VIX uses the unchanged targeted mixture/pooled28 Gaussian controls with fresh
exposures. The exploratory loss retains signed numerical estimates; production
validation remains strict. Two poor starts have 4/7 negative primary estimates,
and secondary-wing negatives remain in the full-training diagnostic files.
Run001 startup and failed Lean development attempts remain retained.

All 7056 searching/fresh objective quote scores replay; 58 captured source/input
hashes are verified, and 200 pilot/search/check seeds are distinct and disjoint
from the reused VIX control pilots. This does not mean all random streams were
independently regenerated by the score auditor. Ten Lean theorems compile with
warnings treated as errors and no sorryAx/custom axioms; 170 Python tests pass.
Lean checks Gaussian one-step means, fixed-measure control reuse and relevant
finite/affine/loss identities; it does not certify Python or price convergence.

Live continuation after this checkpoint: state run 002 completed seven profiles,
selecting C0=1 with search loss 288.8213 versus 333.6764 at zero, then failed
before a fresh score on a final stock-chunk shape mismatch. Five-piece run 002
completed its search (288.8213 to 284.8048), then stopped waiting for the parent.
Failures and selections are retained. Corrected `runs/m2_joint_state_fit_003`
and `runs/m2_curve_refinement_003` completed these frozen selections' fresh
checks; neither search was repeated. The first startup state run was interrupted before optimization,
and five-piece startup 001 failed on relative provenance paths before pricing.
The five-piece refinement adds date-based knots near days 57/151, at the
fixed selected state. It reuses the state run's NEW checking draws, so the two
comparisons are paired/correlated rather than independent validation studies.
See `JOINT_STATE_REFINEMENT_RESULTS.md` and `CURVE_REFINEMENT_PLAN.md`. Check
status/failure/run JSON before starting anything; all current numerical jobs and
audits are complete. Do not repeat searches or rewrite historical inputs. Thirteen Lean
identities pass current audit 005; 179 Python tests pass. The model-forward IV
diagnostic `m2_vix_forward_convention_001` confirms the display gap is largely
forward-driven; it changes no prices or fitting residuals. An independent
market-forward inference cannot be obtained from premiums reconstructed from
the supplied forward and raw IVs. Diagnostic 001 used broader raw VIX pairs for
market parity checks, exposing 18 validation and 17 final-test strike keys;
no held-out model performance score or parameter selection used them. Their
key-only list is in `m2_vix_forward_convention_002/historical_parity_exposure.csv`.
Future evaluation on those keys is quote-exposed, not wholly blind. Diagnostic
002 and current run diagnostics use exact training option keys (no matched pairs).
No qualified new market fit exists yet.

The completed restricted-state/curve refit selects C0=1/Y0=F0=0 and improves
all three fresh trials. Pooled loss 355.46→316.93 improves SPX/VIX-option/forward
components, but only 41/882 estimates are inside bid/ask, all SPX precision
screens fail, and forwards are +1.7115/+0.6619/+2.2645/+1.7582 points too high.
The five-piece refinement adds only 0.64% (316.93→314.90), with one reversal,
descriptive uncertainty larger than the gain, worse forwards and 38 interval
hits. Keep the selected-state THREE-piece curve as the working reference;
retain five-piece evidence without adoption. Do not add further time pieces now.
The new plots/reading guide are `visualizations/joint-state-curve-005/README.md`;
the curve/window explanation is `visualizations/variance-pieces-006/README.md`.

Next bounded fit, now using the all-maturity curve and full in-sample population:
the already declared unspanned Q-risk-price scenarios, one
at a time, starting with feedback pair (+0.25,0)/(-0.25,0) against scenario 0;
only then consider spike pair (0,+0.25)/(0,-0.25). Keep the selected state,
equity premium, readout/architecture fixed and refit the three curve levels.
Declare the finite budget before launch, freeze before NEW checks, and rebuild
candidate normalization/state mapping/stock/VIX exposures. Revisit primitive-law,
control and importance-likelihood correspondence, plus P diagnostics after risk
changes. The current G/M0/exposures cannot be reused unchanged. This next fit
has not been launched. Do not free state/readout/Q risk prices simultaneously
or automatically reopen variance-reduction searches.
The curve-only result is not proof of M2 infeasibility. Keep numerical accuracy,
physical compatibility and eventual model comparison as separate assessments.
Retain all failed points and do not normalize market errors by model MC SEs.
The pending VIX payoff-shaped Gaussian tail-control idea is recorded; it was
not launched and should not displace the current SPX/fitting thread by default.

Use fresh non-overwriting output directories, save source/configuration hashes
and seeds, and update RESEARCH_LEDGER.md with each meaningful outcome. Check
whether the lean-numerics-bridge skill is available and applicable; distinguish
mathematical arguments, implementation checks and numerical evidence. The user
explicitly requested Lean support; maintain the scoped proofs/correspondence
and axiom audits without describing the full simulator or fit as certified.

Previous Python environment: `/private/tmp/pdv-options-venv/bin/python`, with
`OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1` and package source at
`matteo/exogenous_reservoir_eq_m2/src`. Verify it still exists. Keep me informed
with concise updates and explain important conclusions at undergraduate level.
Ask only when a missing decision materially blocks progress.
