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

# Recorded next steps — 23 October 2026

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

**25 October checkpoint: dedicated stochastic sensitivity trial completed.**
See `STOCHASTIC_VOV_RESULTS.md` and `visualizations/stochastic-vov016`.
The prescribed12-point two-cap grid improves wing IV levels but the selected
added state flattens fresh measured wing slopes under both caps, while body
and SPX fits worsen. All original coefficients/xi/floor were fixed; this was
an exploratory mechanism test, not full recalibration. Keep the selected baseline.
No further experiment is running. A next design should distinguish actual
slope/curvature improvement from lifting the whole smile, and justify its
interaction and cross-market calibration before expensive production searches.
The following older queues remain historical context, not automatic authorization.

For a consolidated attempt index and current ordered return-to-fitting queue,
read [RESEARCH_LEDGER.md](RESEARCH_LEDGER.md). The conversation handoff is
[CONTINUE_FITTING_PROMPT.md](CONTINUE_FITTING_PROMPT.md). Older hypotheses below
remain historical when a later checkpoint completed them.

## Latest user request: recalibrate all retained maturities in sample

Read [ALL_MATURITY_FIT_RESULTS.md](ALL_MATURITY_FIT_RESULTS.md) and its captured
plan. `m2_all_maturity_fit_001` is complete. It fits 1294 cutoff-3 SPX options at
29/57/119/211/302/666/1030 days and all 217 VIX options/four forwards. Exclude
only SPX 1/1401/1765/2129. All old market split labels are metadata in this
feasibility exercise; there is no market out-of-sample result. The user confirms
that actual recalibration is required. Refit the three curve levels at the
selected state, fixed structure/readout/Q, and show previous and recalibrated
parameters on every panel. Numerical checking uses the SAME market targets.
The selected levels are [0.05762055,0.07134085,0.06321057]. Separate numerical
loss improves 351.36→308.64, with substantial long-SPX improvement but only
103/1511 interval hits and persistent short-SPX/VIX shape errors. All four score
records and membership/means/SEs replay. The audit-reader raw-VIX shape error
was fixed without repricing; its failure log is retained. No jobs remain running.
This expanded curve is the new working reference for the Q-scenario next step
below. Keep the full in-sample population; do not revert to four SPX dates or
training-only quotes. Old source labels/data/results are retained.

## Preceding state/five-piece checkpoint; use the expanded baseline above next

[JOINT_STATE_REFINEMENT_RESULTS.md](JOINT_STATE_REFINEMENT_RESULTS.md) records
the complete seven-state and nested five-piece investigations. Runs 003 resume
frozen searches after retained implementation failures, without reoptimization.
State/three-piece refitting improves pooled loss 355.46→316.93 on new shared
draws, with three favorable trials but only 41/882 prices inside bid/ask.
Five pieces improve only to 314.90 (0.64%), with one trial reversal, descriptive
uncertainty exceeding the gain, worse forwards and 38/882 interval hits.
Retain the C0=1/Y0=F0=0 three-piece working reference; do not add more pieces.

Next bounded fit: declared unspanned-risk-price scenarios, feedback pair first,
one at a time, with selected state/readout/architecture fixed and three curve
levels refitted. Rebuild dated normalization/state mapping/stock/VIX exposures,
review controls and sampling-law correspondence, revisit P diagnostics and use
NEW numerical checking streams. This next fit has not been launched. Thirteen
Lean identities and 179 Python tests pass; all SPX precision screens and
deterministic bias qualification remain unresolved.

Model-own same-expiry forward IV is a display diagnostic; it changes no prices
or forward residuals. Historical raw-parity diagnostic exposed 18 validation
and 17 final-test VIX strike keys to a market consistency calculation; no
model scoring/selection used them. Retain the exposure list and label future
evaluation at those keys as quote-exposed. Current diagnostics restrict to
training option keys. See AUDIT.md. No held-out model performance was scored.

## Earlier curve-only checkpoint; its restricted-state step is now complete

The user subsequently asked for actual joint fitting now and Lean checks,
superseding the unlaunched nine-point/24-SPX proposal below.
[JOINT_CURVE_FIT_RESULTS.md](JOINT_CURVE_FIT_RESULTS.md) records optimization
of three curve levels on all 707 cutoff-3 SPX options, 175 VIX options and four
forwards, followed by three fresh numerical checks. Pooled loss improves
597.90→347.10, but only 50/882 options are inside bid/ask and VIX forwards
worsen. Retain a provisional candidate, not a qualified market fit.

Ten Lean identities pass compiler/axiom audits; 170 Python tests pass.
[CORRESPONDENCE.md](formalization/CORRESPONDENCE.md) states exactly what is
proved and what remains outside Lean. Accuracy gates remain separate.

Next test the declared restricted C/Y/F initial-state coordinates, with the
existing prior penalty and Q scenario 0, updating candidate normalization and
state paths. Preserve the priority objective/scales and fresh checking. Do
not reuse the current G/M0/exposure cache after changing state or restart a
standalone variance-reduction feature search. The older curve-profile proposal
below is historical, not a pending task.

The user's subsequent visual review requests undergraduate-level explanations
and a later explicit review of maturity/strike coverage. See the
[complete maturity and variance context](visualizations/joint-curve-context-003/README.md).
The existing four-of-eleven SPX maturity choice is a pilot split, not an
optimized liquidity criterion; primary OTM selection likewise does not verify
executed-volume liquidity. Preserve the current cutoff-3 objective until a new
selection rule is explicitly adopted. The 119-day market curve is now visually
exposed, remains excluded from fitting and has no model score. Additional SPX
expiries still need model repricing before they can support model comparisons.

## Current decision: conditional controls completed, reference retained

[SPX_CONDITIONAL_CONTROL_RESULTS.md](SPX_CONDITIONAL_CONTROL_RESULTS.md) completes
the previously pending extension. Conditional existing controls and a valid
retained-shock proposal pass their Gaussian/Black/control oracles, but lose
whole-panel variance-cost efficiency in both fresh trials (0.7483/0.8541).
Tiny conditional-tail nominal SE passes are not tail clearance. Preserve all
negative reference estimates, frozen pilots, source hashes and failed results.
Keep original vol7/ridge5; do not begin another feature search by default.

[FITTING_REENTRY_DECISION.md](FITTING_REENTRY_DECISION.md) records the separate
accuracy gates and the concrete next fitting experiment: a nine-point local
profile of two curve levels, fixed readout/state/Q, on 24 explicitly declared
SPX representatives plus all 175 VIX options/four forwards. It has not run.
The provisional subset never replaces the 707 cutoff-3 priority population.
No qualified joint fit, feasibility decision or model ranking is ready.

## Historical conditional prototype handoff (extension now completed)

The [partial conditional investigation](SPX_CONDITIONAL_PRICING_RESULTS.md)
now establishes the Gaussian construction and a positive raw-payoff prototype.
Only 0.636685% of frozen-model equity variance is independent of the OU drivers;
the raw target-Q market screen shows a modest median variance gain 1.025 and
no convergence clearance. It does not replace original vol7/ridge5.

Next derive and test conditional expectations of the **existing** stock,
reference-payoff and adapted delta controls. A small raw-payoff gain does not
rule out a larger controlled-residual gain. The conditional gain generally is
not zero when retained shocks include correlated equity; use the Gaussian
expectation identities recorded in the report. Derive a new importance proposal
on retained shocks before pricing with weights; old full-path weights cannot
be reused after integration. Freeze pilots/checking before any market comparison.
This extension is pending, not implemented or validated. Keep cutoff 3/all 707
priority quotes, secondary wings and all four VIX maturities, with independent
stock-timestep and normalizer refinement gates still pending.
Older entries below retain the prior hypothesis and experiment history.

## Latest sampling result and next structural candidate

The [fixed +3 equity-tilt study](SPX_CALL_TILT_RESULTS.md) completes pilot-only
component diagnosis and a fresh, paired comparison with doubled pilot budgets.
Outer-call coverage improves at every maturity, but no robust full-panel
variance-cost gain follows; the original sampler/ridge5 remains reference.
Both proposal and pilot-size comparisons retain failures and mean warnings.
Cutoff 3/all 707 priority quotes and all four VIX maturities remain active.

Next investigate conditional SPX Monte Carlo before adding further mixture
magnitudes or control features: derive the joint-Gaussian equity/OU transition
conditioning in the **full** state, identify which equity noise can be integrated
analytically under the existing discrete stock scheme, and audit numerical rank
and covariance correspondence. Validate Gaussian/constant-volatility and
same-model checks before market experiments. Derive marginal likelihoods if
importance sampling is retained; full-path weights cannot be reused blindly
after integrating noise out. This is a pending candidate, not completed pricing
or a positivity guarantee. Keep stock-timestep/normalizer refinement gates
separate. Older sampling queues below are historical where superseded.

## Latest completed experiment and immediate sampling work

The user-approved [24-contract cutoff-3 panel](SPX_CENTRAL_PANEL_RESULTS.md) now
has independent pilot/checking results at all four SPX training maturities.
No negative estimates were observed, but 0/48 sampling screens pass and the
linear controls lose full-panel variance-cost efficiency in all eight checks.
Long outer-call pilots contain only one positive payoff. Keep the baseline as
reference and cutoff 3/all 707 training quotes as the priority population.

Next inspect saved pilot residual/payoff contributions by side, distance band
and proposal component. Freeze a better-covered pilot and/or horizon-specific
sampling comparison before new checking; no automatic multi-million-path run
is authorized by unstable plug-in budget forecasts. If old checking results
inform a new design, mark them historical/design evidence and use fresh checks.
Keep coupled 256/512 and independent normalizer refinement pending; neither was
performed by this direct 512-step sampling run. All wings/VIX maturities stay recorded.
The earlier queue below is historical where superseded by this checkpoint.

## Current membership rule: maturity-scaled central band

User requested maturity adaptation and meaningful coverage. Prioritize
`abs(log(K/F(T))) <= 3 * sigma_ATM_proxy(T) * sqrt(T)`, with market training-only
nearest-ATM-IV anchors. This supersedes the fixed 5% band for priority membership.
User selected cutoff **3 provisionally**, superseding cutoff 2.
It includes **707/980 training quotes (72.1%)**, over half at every training
expiry (29/57/211/302 days). The four earlier near-ATM contracts remain diagnostic
representatives, not the entire priority set. All quote and VIX scopes remain.
See [coverage and denominator correction](SPX_SCALED_MONEYNESS_COVERAGE.md).

Next broaden representative strike checks within the central band, keeping
sampling, stock timestep and normalizer gates separate. Longer SPX maturities
need their own convergence/cost checks; no new pricing accuracy is established
by the coverage count. Older priority text below is historical where superseded.

## Current priority: near-ATM SPX (user confirmed)

User: **"Yes, I would also focus more on near-ATM options for now."**
This supersedes wing-control stabilization as the immediate research priority.
It does not remove quotes, maturities or historical warnings from the scope.
`active_scope.json` records the four previously declared training contracts:
29-day 6755 put / 7030 call, and 57-day 6775 put / 7050 call (about 2% OTM).
All four VIX maturities remain active; wings remain secondary diagnostics.

Next bounded work:

1. Assess near-ATM residual variance, control-mean uncertainty and pilot stability
   using saved data, then freeze the baseline/linear-control comparison and a
   sample budget before fresh checking. Do not select gains on checking results.
2. Check fresh near-ATM sampling accuracy and actual runtime at 512 steps/day.
   Positive estimates and roughly 1% SE/price are encouraging but do not meet
   the existing quote-relative target (SE <= scale/8).
3. Once sampling permits interpretation, assess coupled 256/512 stock-price
   changes on these contracts and refine further if needed. Independently
   check normalizer quadrature; control variance reduction does not clear bias.

Retain signed negative controlled wing estimates as unusable-price diagnostics,
not clipped prices. No model/curve/state refit or new broad simulation was run
for this priority update. See [the moneyness audit](SPX_MONEYNESS_POSITIVITY_AUDIT.md)
for evidence and limits. The older queue below is retained as research history;
this section determines the current priority.

User direction: record the VIX tail-control next step, move forward for now,
and keep all maturities. The working interpretation is to return to SPX
numerical diagnostics while leaving the VIX extension pending.

The user subsequently confirmed: **"Yes, let's focus now on that for a while."**
SPX sampling, timestep and normalization diagnostics are therefore the current
explicit focus. The first bounded path-sampling experiment is recorded in
[SPX_PATH_IMPORTANCE_RESULTS.md](SPX_PATH_IMPORTANCE_RESULTS.md); the VIX queue
and full maturity scope below remain in force.

SPX follow-up: [SPX_VOLATILITY_SAMPLING_RESULTS.md](SPX_VOLATILITY_SAMPLING_RESULTS.md)
records independently selected volatility-path sampling that improves fresh
variance in 20/24 cases, while no quote-relative sampling screen passes.
Additional Gaussian-score controls help only modestly. Next prioritize paired
time corrections: standalone level costs, independent correction-pilot fits,
and a bounded matched-cost MLMC comparison. The 128->256 warnings, pilot
instability and tail precision remain open; no model refit is cleared.

The [first MLMC prototype](SPX_MLMC_RESULTS.md) now checks independent
64/128/256/512 levels with frozen allocation and correction-specific controls.
It loses on the central variance-cost objective in both fresh trials and
**does not replace direct vol7 simulation**. Sparse top-level tails and
256->512 warnings remain. Before another allocation study, investigate better
coarse/fine Gaussian coupling on tractable oracles; conditional bridge
reflection is a candidate only, with target-law/likelihood/adaptedness checks
still to be established. Do not interpret the completed prototype as a
completed improved coupling or continuous-model accuracy clearance.

The [endpoint reflection study](SPX_BRIDGE_COUPLING_RESULTS.md) has now completed
the algebraic/Gaussian/Black and independent checking step. Endpoint and target
law checks pass, but all four central variance-cost comparisons lose (gains
0.67–0.77). Direct vol7 and ordinary corrections remain the references.
The same-coefficient averaged residual retains 86–96% of controlled variance.

**Next SPX hypothesis, pending:** cheap projected-OU innovation martingale
controls, including centered quadratic innovations with known conditional
covariance, as complements to the existing stock hedges. Preserve full-state
conditional means; the three readout projections alone are not Markov. Prove
the zero-mean/likelihood correspondence, validate Gaussian/Black and forward
oracles, train on separate pilots, then check direct-price and correction
variance/cost separately. This is not an implemented control or a bias fix;
512-step, tail/pilot and normalization gates stay open. VIX scope is unchanged.

The [full OU-innovation study](SPX_INNOVATION_CONTROL_RESULTS.md) now implements
and independently checks that hypothesis. Central direct-price variance gains
around 1.8–1.9 are nearly consumed by a 1.8 runtime multiplier; corrections
lose on cost. Repeated far-call paired mean shifts, contract-control mean
warnings and pilot sensitivity remain unresolved, despite Gaussian/covariance
and bare-innovation checks passing. The full-nine dictionary remains experimental.

**Next candidate:** an actually cheaper three-linear-control kernel. Old-pilot
two-fold ablations retain most direct-price benefit, but are post-hoc design
evidence only. Avoid unused square/mixed calculations; profile standalone
pricing; freeze before fresh checking and retain every representative strike.
Audit contract-specific control tails/paired mean shifts as well as bare
innovation sums. Keep correction performance and timestep bias separate;
baseline controls remain the reference and no market refit is cleared.

## Scope retained

[active_scope.json](active_scope.json) remains authoritative and unchanged:
**27/55/118/153-day VIX maturities, all 175 training options and four forwards**.
No VIX maturity is deferred or removed. The existing SPX maturity/split scope
also remains unchanged; representative pilot contracts are diagnostics, not
an exclusion of the other targets. Held-out market observations stay reserved.

## Pending VIX extension — recorded, not launched

Screen payoff-shaped Gaussian tail controls with analytically known means,
using the existing nine-control/ridge design and pooled 28-pilot fit as the
reference. Learn and compare candidates on design-only pilot replications,
freeze the dictionary and coefficients, then validate on new shared checking
paths at all four maturities. Retain raw/legacy estimates, payoff-hit counts,
negative adjusted estimates and paired sample-count checks. Keep conditional
quadrature, time integration and normalization checks separate.

This extension remains pending for a later bounded experiment. It is not a
completed result, an automatic background job or a convergence clearance.
The current checkpoint is
[VIX_GAUSSIAN_STABILITY_RESULTS.md](VIX_GAUSSIAN_STABILITY_RESULTS.md): all four
pooled forward SE screens pass, but the 153-day pass is borderline and not
supported by every pilot fit; only 45/175 pooled option SE screens pass, and
118-day second-moment warnings remain recorded. No market refit is authorized
by those numerical screens alone.

## SPX work to advance next

Resume the unresolved sampling, stock-integration and normalization work in
[NUMERICAL_STABILITY_RESULTS.md](NUMERICAL_STABILITY_RESULTS.md):

1. Design a small path-level Gaussian importance-sampling pilot against the
   current adapted hedging-control baseline. SPX depends on the whole stock/
   volatility path; the VIX terminal-state proposal cannot simply be reused.
   Preserve the Q model, use exact path likelihood weights and separate pilot
   learning from checking. Check martingale/reference expectations and rare
   option payoffs with tractable Gaussian/Black oracles before market runs.
2. Keep coupled 64/128/256 steps-per-trading-day comparisons. Exact OU state
   transitions do not remove nonlinear stock-integration bias; the old 2/8-step
   diagnostics remain coarse historical evidence.
3. Measure sampling error, paired-grid correction variance/cost and normalizer
   refinement separately. Preserve zero-hit and unsuccessful cases; do not
   infer accurate tail prices from a zero empirical SE or fit a new market
   model before the numerical gates are resolved.

The latest source review confirms that the SPX runner already supplies a
12-contract diagnostic, independent pilot coefficients, adapted hedging controls
and coupled 64/128/256 grids. Its recorded results do not qualify the full
surface. No new SPX simulation has been launched by this queue update.
## Linear-only follow-up completed; next diagnostic

The [direct512 linear-only experiment](SPX_LINEAR_INNOVATION_RESULTS.md) retains
a modest central variance-cost gain (1.14x in both fresh trials), but far-put
outlier leverage and 30.52x half-pilot SE inflation prevent general adoption.
All primary sampling screens still fail. Keep the baseline as reference.
Next examine saved-data coefficient/scaling and tail-contribution diagnostics;
any shrinkage/stability rule is exploratory until frozen using pilot data and
validated on new checking samples. Do not select a rule on existing checking
performance. Stock-timestep and normalization accuracy remain separate open gates.

## Moneyness and positivity follow-up (22 October 2026)

[SPX_MONEYNESS_POSITIVITY_AUDIT.md](SPX_MONEYNESS_POSITIVITY_AUDIT.md) checks saved
linear-control results. All 12 representatives are technically OTM; the four
about 2% from the forward are positive, with only 0.993–1.026x half-pilot SE
ratios. All eight negative estimates across pooled/half-pilot methods are puts
at K/F about 0.493–0.499. Near-ATM SE/price is 0.95–1.43%, but quote-scale screens
still fail. These facts support a central research focus, not convergence or
ITM/full-surface claims. Negative controlled estimates are unusable as quotes
and retained as signed diagnostics; no clipping or checking-selected fallback.
Raw positive-weight payoffs are nonnegative by construction but can be noisy.
All maturity/quote scopes are retained. Reproduction inputs and hashes are saved
in `runs/m2_spx_moneyness_audit_001`.
