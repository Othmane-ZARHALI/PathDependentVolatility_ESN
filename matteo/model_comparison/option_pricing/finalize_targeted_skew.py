"""Durable report, uncertainty/tail records and integrity manifest."""
import hashlib,json,shutil
from pathlib import Path
import numpy as np
import pandas as pd
from .targeted_skew import OUT,HERE,ALPHA,VIX,arr,dump
from .asymmetric_readout import OUT as OLD
from .vix_only_curve_fit import DAYS,ACTIVE

def main():
    plot=HERE/'visualizations/targeted-skew028';summary=json.loads((plot/'summary.json').read_text());config=json.loads((OUT/'candidates.json').read_text());audit=json.loads((OUT/'implementation_audit.json').read_text());extended=json.loads((OUT/'extended_audit.json').read_text());bridge=json.loads((OUT/'primitive_bridge_audit.json').read_text());ledger=json.loads((OUT/'VIX_fit_ledger.json').read_text());probe=json.loads((OUT/'quadrature_probe.json').read_text());momerr=json.loads((OUT/'moment_order_audit.json').read_text())['maximum_GH48_vs64_relative_error'];main=json.loads((OUT/'run.json').read_text());refinement=json.loads((OUT/'grid_refinement/run.json').read_text());fine=summary['refined'];metrics=pd.read_csv(plot/'SPX_maturity_metrics.csv');m29=metrics[(metrics.stage=='refined256')&(metrics.day==29)].set_index('candidate');grids=summary['grid'];ci=summary['paired'];tail=[]
    # Raw positive payoffs diagnose concentration; never repair signed prices.
    q=pd.read_csv(OUT/'SPX_membership.csv');stock=np.load(OUT/'grid_refinement/stocks.npz')['stocks'];selected=[0,1,4,9];days=sorted(q.expiry_days.unique())
    for j,c in enumerate(selected):
        for e,day in enumerate(days):
            g=q[(q.expiry_days==day)&q.priority_selected];S=stock[:,e,j,1]
            for side in ['put','call']:
                gg=g[g.option_type==side]
                if len(gg)==0:continue
                K=float(gg.K_over_F.min() if side=='put' else gg.K_over_F.max());pay=np.maximum(K-S,0) if side=='put' else np.maximum(S-K,0);total=pay.sum();tail.append(dict(instrument='SPX',candidate=c,day=int(day),side=side,K_over_F=K,raw_payoff_mean=float(pay.mean()),payoff_concentration_ESS=float(total**2/(pay@pay)) if total else 0.,top10_share=float(np.sort(pay)[-10:].sum()/total) if total else 0.,scope='unweighted target-Q raw payoff concentration; antithetic dependence means this is not an independent sample count'))
    H={day:np.load(OUT/f'checking_exposures/d{day}.npz')['H'] for day in DAYS};draws={day:dict(np.load(OLD/f'fresh_VIX_draws/d{day}.npz')) for day in DAYS};x=arr(config['models'])[0];qv=pd.read_csv(OUT/'VIX_membership.csv');f=pd.read_csv(OUT/'forward_targets.csv');v=VIX(H,draws,x,qv,f)
    for c in [0,4,9]:
        idx=config['mean_indices'][c];r=v.predict(idx,arr(config['models'])[c,ACTIVE])
        for day,a in r['raw'].items():
            K=float(qv[qv.expiry_days==day].strike.max());pay=a['weights']*np.maximum(a['V']-K,0);total=pay.sum();tail.append(dict(instrument='VIX',candidate=c,day=int(day),strike=K,weighted_payoff_mean=float(pay.mean()),mass=float(a['weights'].mean()),payoff_concentration_ESS=float(total**2/(pay@pay)) if total else 0.,top10_share=float(np.sort(pay)[-10:].sum()/total) if total else 0.,scope='likelihood-weighted actual VIX payoff; no sample mass normalization'))
    dump(OUT/'tail_diagnostics.json',tail)
    labels=['Original fixed loading','Adaptive stock loading','Regime α1.5, fixed curve','Regime α1.5, seven-level refit'];vix={r['candidate']:r for r in main['VIX_scores']};vix[1]=vix[0];table=[]
    for r,label in zip(fine,labels):
        c=r['original_candidate'];table.append(f"| {label} | {r['loss']:.2f} | {r['IV_MAE_pp']:.3f} | {r['invalid_IV']} | {vix[c]['option_loss']:.4f} |")
    text=f'''# Targeted SPX skew experiments — 28 October 2026

**Both constructions implemented and checked; neither adopted.** Adaptive stock exposure preserves the complete VIX law but gives only modest SPX changes in these tests. Conditional stress-regime normalization reduces the VIX damage of the earlier global multiplier, but strong settings still spoil VIX smiles. Neither yields a satisfactory fit across the seven SPX maturities. These are bounded numerical experiments, not a model-family impossibility result.

## Fixed model and scope

All original M2 readout parameters, positive spike shift, exogenous OU dynamics, stress speeds/amplitudes/weights and floor remain fixed. No new bank, jumps, quintic, endogenous OU feedback or new constant-theta sweep. The previous protected16-level curve is the original baseline. The stock driver for readout-only comparisons remains the mixed theta=.8 rotation. All1751eligible/1294primarySPX quotes, sevenSPX maturities,217VIX quotes and fourVIX forwards remain. Prices are signed; invalid IVs and negative estimates are retained.

## Experiment1: adaptive stock exposure

In independent primitive Brownian coordinates write

$$dv_t=\\mu_tdt+\\mathbf k_t\\cdot d\\mathbf B_t,\\qquad dS_t/S_t=\\sqrt{{v_t}}\\mathbf a_t\\cdot d\\mathbf B_t.$$

Use unit stock direction opposite the current variance sensitivity. For xi>floor,

$$\\mathbf a_t=-\\mathbf k_t/\\|\\mathbf k_t\\|.$$

This minimizes instantaneous stock/variance covariance among unit directions. It does not maximize a finite-maturity option-skew slope. All11relevant old primitive Brownian directions and two stress directions are included. The original M2 cap derivative and instantaneous projection diffusion coefficients are used, not finite-step normalized OU innovations.

Implementation boundary convention: the common positive xi-floor factor cancels in normalization, so the code uses the gradient of H0. At xi=floor, the continuous limiting H0 direction is retained; all directions have zero instantaneous variance covariance there. If H0's gradient vanishes, use the old mixed stock loading. This refines the initial plan's generic zero-gradient fallback and is recorded explicitly in extended_audit.json.

Only the stock equation changes. The full Gaussian factor law, variance paths/readout and conditional future variance law stay fixed; hence all VIX prices and forwards are identical to the original model. The inherited checking forward gap up to.0552VIXpoints and old estimator error are unchanged, not newly eliminated. No exact market-forward equality or continuous-stock-martingale certificate follows from unit loading norm.

Existing retained OU streams are replayed unchanged. Missing old-bank Brownian increments are sampled conditional on the original Gaussian innovations with the full covariance, plus independent residuals. Stock loading and variance are held at each stock-grid start, including coarse-grid controls. No future-state loading enters a coarse step.

## Experiment2: stress-regime-normalized readout

Let a=1-w, A=(1-r)As+rAf, H0=aR+wA and D_alpha=exp(alpha softplus(-X)), with existing X=(Y-E_QY)/stationary_sd(Y). Set mRD(t)=E_Q[R_tD_t],mD(t)=E_Q[D_t] and

$$H_t^{{new}}=H_t^0D_t\\frac{{a+wA_t}}{{a\\,m_{{RD}}(t)+wA_t m_D(t)}}.$$

Then, using independence of the old bank and stress factors,

$$E_Q[H_t^{{new}}\\mid Z_t^s,Z_t^f]=a+wA_t,\\qquad E_Q[H_t^{{new}}]=1.$$

This fixes the average total variance response within each stress state. It does not preserve the internal mean allocation to the two algebraic terms, the complete conditional future variance law, VIX distribution or VIX forwards. VIX is repriced. Alpha0 recovers H0.

For future-time pricing define positive stress functions

$$F_1(A)=\\frac{{a(a+wA)}}{{a m_{{RD}}+wA m_D}},\\quad F_2(A)=\\frac{{wA(a+wA)}}{{a m_{{RD}}+wA m_D}}.$$

Conditional E[Hnew] equals E[RD|old state]E[F1(A)|stress state]+E[D|old state]E[F2(A)|stress state]. This uses correct joint old-bank product moments, plus two-dimensional stress Gaussian quadrature. Products of separately normalized old factors are not substituted.

Alpha0,.5,1,1.5,2 is shared across all maturities. Each positive strength is evaluated on the original curve; a separate fit adjusts only seven VIX-sensitive variance levels with four training-estimator forward equalities, leaving nine off-window levels fixed. The alpha0-refit control is included. Every bounded fit is forward-feasible; maximum training forward gap is {max(r['max_forward_gap'] for r in ledger):.3g}VIXpoints. This is sampled-estimator equality, not population equality. Checking does not re-enforce those equalities. The VIX-loss-preserving preselection retains alpha0; no positive tested strength met the declared baseline+.1loss screen.

## Finer-grid results

The following SPX numbers use2048fresh paths and256steps/trading-day, with paired128/256 grids. The main experiment uses512replay and1024fresh stock paths, with paired64/128 grids. Alpha1.5 was frozen as an illustrative diagnostic before the finer-grid run; no winner was reselected from it.

| Variant | SPX interval loss | SPX IV MAE,pp | Invalid primary IVs | VIX interval loss |
|---|---:|---:|---:|---:|
{chr(10).join(table)}

SPX IV MAE averages available IVs equally across maturities; unavailable IVs are counted. VIX numbers use the same separately generated4096full-factor states per expiry previously used by the asymmetric-readout study, plus the old65536path controlled anchor. These are previously examined research draws, not a newly independent market sample. The original/adaptive VIX entries are the same model law, not separate numerical fits. Conditional-normalization alpha1.5-refit has checking forward gap {vix[9]['max_forward_gap']:.4f}points; its alpha0-refit control has VIX loss{vix[6]['option_loss']:.4f} and checking gap{vix[6]['max_forward_gap']:.4f}points.

At29days, IV(right)-IV(left) at quoted strikes nearest K/F=.9,1.1 is {m29.loc[0,'model_right_minus_left_IV_pp']:.2f}pp for the baseline, {m29.loc[1,'model_right_minus_left_IV_pp']:.2f}pp with adaptive stock exposure and {m29.loc[9,'model_right_minus_left_IV_pp']:.2f}pp for the alpha1.5-refit readout; market is {m29.loc[0,'market_right_minus_left_IV_pp']:.2f}pp. The missing short-end slope remains. The dramatic29day improvement of the earlier global-multiplier/curve combination is not reproduced here.

Compared with the earlier whole-response alpha1.5 multiplier's VIX loss25.725, the new alpha1.5-regime refit yields6.159 on the same VIX checking states/anchor and quotes. This supports the value of conditional normalization for limiting VIX damage, not adequacy of the resulting smile or a controlled comparison with every calibration parameter fixed.

## Uncertainty and implementation evidence

Eight antithetic-pair batches underpin stock uncertainty. The unchanged stress-only auxiliary reference uses the previous65536paths, its correct conditional lognormal pricing and a deterministic delta stock-difference control. That reference has no new multiplier/adaptive loading, so its conditional formula stays valid. It is shared across screening/checking; no future-stress-dependent delta is used. Price differences are paired on common factors. For the2048path finer-grid check, all reported approximate paired95%price-loss change intervals include zero. Available-IV point estimates are not accuracy certificates or proof of improved market fit.

At128/256 refinement, median absolute grid changes in quote scales are {grids[4]['median_grid_change_scales']:.3f}(baseline),{grids[5]['median_grid_change_scales']:.3f}(adaptive),{grids[6]['median_grid_change_scales']:.3f}(regime fixed) and{grids[7]['median_grid_change_scales']:.3f}(regime refit). Sampling and strong-readout grid errors remain unresolved. Stock means/SEs, signed prices, native-IV failures, paired intervals and raw likelihood-weighted payoff concentration are saved; no mean renormalization, clipping or smoothing is applied.

Checks pass: original baseline stock replay {json.loads((OUT/'baseline_replay.json').read_text())['maximum_stock_error']:.3g}; independent scalar primitive/stock recursion {audit['scalar_stock_recursion_max_error']:.3g}; unit direction norm {audit['unit_stock_loading_max_error']:.3g}; readout-gradient finite differences {audit['readout_gradient_max_scaled_error']:.3g}; old-bank Brownian/OU standardized covariance error {bridge['maximum_standardized_Brownian_OU_cross_error']:.3g}; independent stress quadrature {audit['stress_quadrature_max_scaled_error']:.3g}. The conditional-mean identity holds algebraically, and interpolated means versus independentGH64 at declared probes agree to {max(r['maximum_relative_conditional_mean_error'] for r in extended['conditional_mean_interpolation_checks']):.3g}; seven new forward derivatives agree within {extended['new_regime_forward_derivative_max_absolute_error']:.3g}. MomentGH48vs64 error at declared probes is {momerr:.3g}. Matched128state/expiry GH32/stress32/time8 versus GH48/stress48/time12 probes change VIX by at most {max(max(r['maximum_VIX_change_points']) for r in probe):.3g}points. These are finite numerical checks, not uniform error bounds.

No new Lean theorem was introduced. Existing coefficient/normalization identities retain their original premises. Numerical audits do not certify the compiled simulator, VIX estimator convergence or continuous-time stock martingality. All original source fingerprints are checked at finalization; no prior production or experimental baseline is overwritten.

## Conclusion and next research question

The first construction is a valid way to change SPX dependence while preserving VIX, but maximum negative instantaneous covariance does not deliver the required SPX skew in this diagnostic. The second construction preserves more of the successful stress response and reduces VIX damage, but does not produce an acceptable double calibration in this bounded search.

Do not repeat the same constant-theta or global-multiplier strength scan. A more informative next question is which existing shocks predict persistently higher future integrated variance, rather than only a large instantaneous variance change. Aligning stock exposure with that future-variance sensitivity is a distinct, untested refinement that can still leave the factor/VIX law fixed. Whether short-lived feedback noise explains the present outcome needs that diagnostic; it is not established here. Broader stress/readout refitting and additional banks were not tested, so no model-family impossibility follows.

Figures: [finer SPX comparison](visualizations/targeted-skew028/SPX-targeted-256.png), [main SPX comparison](visualizations/targeted-skew028/SPX-targeted-128.png), [VIX comparison](visualizations/targeted-skew028/VIX-targeted.png), [readout response](visualizations/targeted-skew028/regime-response.png). Detailed run: `runs/m2_targeted_skew_001`.
'''
    (HERE/'TARGETED_SKEW_RESULTS.md').write_text(text)
    note='''## 28 October 2026 — adaptive stock exposure and stress-regime normalization tested

Read [TARGETED_SKEW_RESULTS.md](TARGETED_SKEW_RESULTS.md); run `runs/m2_targeted_skew_001`, figures `visualizations/targeted-skew028`. Both user-authorized constructions implemented separately. Exogenous OU dynamics/banks, original M2 parameters/positive spike shift, stress parameters and floor fixed; all sevenSPX maturities/1751eligible/1294primary and217VIX quotes/four forwards retained. No new constant-theta sweep, endogenous OU feedback or production adoption.

Adaptive stock loading is opposite current H0 diffusion sensitivity in all13relevant independent Brownian directions. For xi>floor this maximizes negative instantaneous logstock/variance covariance; at xi=floor use the continuous limiting H0 direction, and at zero H0 gradient fall back to the original rotation. All factor/variance/VIX laws unchanged. Correct missing old-bank Brownian increments are sampled conditional on retained OU innovations with covariance error1.1e-11; original OU streams replay unchanged.

Regime readout Hnew=H0 D[(1-w)+wA]/[(1-w)E(RD)+wA E(D)], with existing feedback X and alpha0,.5,1,1.5,2. It preserves E[Hnew|stress states]=(1-w)+wA, not the full VIX law or separate internal component budgets. Compare fixed curve and seven-level forward-constrained refits; every training forward equality feasible(maxgap3.4e-10), but no positive strength meets the declared VIX baseline+.1loss screen. Alpha0 retained before fresh stock checking.

Main512replay/1024fresh stock paths; extra2048fresh paths on paired128/256grids. At256steps, SPX price loss/IVMAE: baseline1869.58/3.203pp; adaptive1821.00/3.040; regimealpha1.5fixed1705.08/2.960; regimealpha1.5refit1560.53/2.897. Invalid primary IVs0/4/4/2; all approximate finer-check paired loss-change intervals include0. Median128/256 grid changes.090/.230/.866/.831quote scales. 29day quoted right-minus-left IV gap nearK/F=.9,1.1: baseline-5.22pp,adaptive-4.90,regime1.5refit-6.29,market-14.05; required slope still missing.

VIX baseline/adaptive loss.05008identical law. Regimealpha1.5fixed loss7.684; refit6.159 with checkingmaxforwardgap.0784points, versus earlier globalmultiplierrefit25.725. Same4096VIX checking states and old65536anchor; research stream reused, not newly independent market validation. All signs/IVfailures/likelihood tails retained. No satisfactory double calibration or family-impossibility conclusion. Independent scalar,gradient,unitnorm,covariance,stressconditionalmean,7forwardderivative and quadrature checks pass; no new Lean/continuousmartingale theorem. Stock sampling/grid accuracy unresolved.

Do not repeat constanttheta/globalmultiplier scans. Next untested question: target shocks that raise future integrated variance persistently, rather than maximizing only instantaneous covariance. Keep OU factor-law reuse and VIX preservation when studying that refinement. Additional banks and broad stress/readout refitting remain untested.

'''
    for name in ['NEXT_STEPS.md','RESEARCH_LEDGER.md','CONTINUE_FITTING_PROMPT.md']:
        path=HERE/name;path.write_text(note+path.read_text())
    path=HERE/'formalization/CORRESPONDENCE.md';s=path.read_text();addition='''## 28 October 2026: two targeted skew numerical implementations

`targeted_skew.cpp` adds13-coordinate state-dependent unit stock exposure and a stress-regime-normalized variance readout. Scalar recursion, readout gradient, primitive/OU covariance, unit norm, stress quadrature, interpolated conditional means and new forward derivatives have numerical audits in `../runs/m2_targeted_skew_001`. No new Lean theorem is claimed. The variance/VIX-law invariance argument applies only to the stock-loading experiment; the new readout is repriced. The conditional-mean cancellation is not preservation of the full future conditional variance law. Maximizing instantaneous covariance is not a theorem of optimal finite-maturity implied skew. Normalized H0 sensitivity uses its continuous extension at xi=floor; the original rotation is the fallback only at zero H0 gradient. Existing formal sources/audits remain unchanged.

''';path.write_text(s.replace('# Lean checks for the joint curve fit\n\n','# Lean checks for the joint curve fit\n\n'+addition,1))
    sources=['audit_targeted_skew.py','extended_audit_targeted_skew.py','check_targeted_quadrature.py','refine_targeted_skew.py','plot_targeted_skew.py','finalize_targeted_skew.py','TARGETED_SKEW_RESULTS.md']
    for name in sources:shutil.copy2(HERE/name,OUT/'source'/name)
    design=json.loads((OUT/'design.json').read_text())
    for file,h in design['source_hashes'].items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==h,file
    prior_counts={}
    for name in ['m2_leverage_scan_001','m2_asymmetric_readout_001']:
        folder=OUT.parent/name;manifest=json.loads((folder/'artifact_manifest.json').read_text())
        for file,h in manifest.items():assert hashlib.sha256((folder/file).read_bytes()).hexdigest()==h,(name,file)
        prior_counts[name]=len(manifest)
    dump(OUT/'finalization.json',dict(status='completed both mechanisms, two checking grids, numerical audits, uncertainty and tails, plots and durable report',adopted=False,authoritative_report='TARGETED_SKEW_RESULTS.md',prior_artifact_hash_counts_verified=prior_counts,initial_source_fingerprints_unchanged=True,wall_seconds_main=main['wall_seconds'],wall_seconds_refinement=refinement['wall_seconds'],visualization_hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in plot.iterdir() if p.is_file()}))
    files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='artifact_manifest.json');manifest={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files};dump(OUT/'artifact_manifest.json',manifest)
    for file,h in manifest.items():assert hashlib.sha256((OUT/file).read_bytes()).hexdigest()==h
    print('TARGETED FINALIZATION PASS',len(manifest),'artifact hashes; prior',prior_counts,flush=True)

if __name__=='__main__':main()
