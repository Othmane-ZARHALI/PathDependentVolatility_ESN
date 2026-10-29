"""Report the successful and unsuccessful readout strengths without clipping."""
import hashlib,json,shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .asymmetric_readout import OUT,HERE,ALPHA,dump
from .protected_curve_check import loss
from .plot_joint_curve_fit import implied_vol

def main():
    out=HERE/'visualizations/asymmetric-readout028';out.mkdir(exist_ok=False,parents=True);s=pd.read_csv(OUT/'fresh/SPX_scores.csv');s=s[s.steps==128];vs=json.loads((OUT/'fresh_VIX/scores.json').read_text());ledger=json.loads((OUT/'VIX_fit_ledger.json').read_text());audit=json.loads((OUT/'implementation_audit.json').read_text());run=json.loads((OUT/'run.json').read_text());frames={i:pd.read_csv(OUT/f'fresh/SPX_prices_c{i}_128.csv') for i in range(7)};days=sorted(frames[0].expiry_days.unique());colors=plt.cm.viridis(np.linspace(.05,.95,6));metrics=[]
    def save(fig,name):
        fig.tight_layout()
        for ext in ['png','svg']:fig.savefig(out/f'{name}.{ext}',dpi=170)
        plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(16,8))
    for ax,day in zip(axes.flat,days):
        g=frames[0][(frames[0].expiry_days==day)&frames[0].priority_selected];ax.scatter(g.K_over_F,g.market_IV,color='black',s=8,alpha=.5,label='Market')
        for c,color in zip(range(1,7),colors):
            a=frames[c][(frames[c].expiry_days==day)&frames[c].priority_selected];alpha=float(s[s.candidate==c].alpha.iloc[0]);ax.scatter(a.K_over_F,a.model_IV,color=color,s=7,alpha=.7,label=f'α={alpha:g}');lo=a.loc[(a.K_over_F-.9).abs().idxmin()];hi=a.loc[(a.K_over_F-1.1).abs().idxmin()];metrics.append(dict(candidate=c,alpha=alpha,day=int(day),IV_MAE_pp=float(np.nanmean(abs(a.model_IV-a.market_IV))),invalid_IV=int(a.model_IV.isna().sum()),stock_mean=float(a.stock_mean.iloc[0]),stock_mean_SE=float(a.stock_mean_SE.iloc[0]),model_right_minus_left_IV_pp=float(hi.model_IV-lo.model_IV),market_right_minus_left_IV_pp=float(hi.market_IV-lo.market_IV)))
        ax.set_title(f'{day} days');ax.set_xlabel('Strike / forward');ax.set_ylabel('Implied volatility (%)');ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=7);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,'Whole-response downside multiplier\nExisting feedback factor\nθ=.8, OU dynamics fixed\n\nEach α: seven VIX-sensitive levels fit\nFour training forward equalities\n1024 fresh stock paths\nSame auxiliary reference sample\n\nSigned prices and invalid IVs\nretained in all CSV files',va='top');fig.suptitle('SPX: modest gains from stronger asymmetry, with growing sampling error');save(fig,'SPX-asymmetric-smiles');pd.DataFrame(metrics).to_csv(out/'SPX_maturity_metrics.csv',index=False)
    fig,axes=plt.subplots(1,4,figsize=(16,4.4));vixIV=[]
    for ax,day in zip(axes,DAYS):
        for c,color in zip(range(1,7),colors):
            a=pd.read_csv(OUT/f'fresh_VIX/VIX_prices_c{c}.csv');a=a[a.expiry_days==day].copy();a['model_IV']=[implied_vol(r.model_price_U,r.forward,r.strike,r.T_years,r.option_type) for r in a.itertuples()];a['market_IV']=[implied_vol(r.price_mid_U,r.forward,r.strike,r.T_years,r.option_type) for r in a.itertuples()];alpha=float(s[s.candidate==c].alpha.iloc[0]);ax.plot(a.strike,a.model_IV,'.-',lw=1,color=color,label=f'α={alpha:g}');a['candidate']=c;vixIV.append(a)
            if c==1:ax.scatter(a.strike,a.market_IV,s=10,color='black',label='Market')
        ax.set_title(f'{day} days');ax.set_xlabel('VIX strike');ax.set_ylabel('Implied volatility (%)');ax.grid(alpha=.15)
    axes[0].legend(fontsize=7);fig.suptitle('VIX: forwards constrained during training, smiles deteriorate at larger α\nFresh factor draws; IVs use the supplied market forward, forward gaps reported separately');save(fig,'VIX-asymmetric-smiles');pd.concat(vixIV,ignore_index=True).to_csv(out/'VIX_native_IV.csv',index=False)
    fig,axes=plt.subplots(1,3,figsize=(14,4));a=s[s.candidate>0];axes[0].plot(a.alpha,a.IV_MAE_pp,'.-',color='#177c92');axes[0].set_ylabel('SPX IV MAE, percentage points');axes[1].plot(ALPHA,[z['option_loss'] for z in vs[1:]],'.-',color='#b76030');axes[1].set_ylabel('Fresh VIX price interval loss');share=audit['stress_mean_budget_share'][1];axes[2].plot(ALPHA,np.array(share['mean_stress_share'])*100,'.-',color='#75549b');axes[2].set_ylabel('Stress share of mean excess variance (%)')
    for ax in axes:ax.set_xlabel('Multiplier strength α');ax.grid(alpha=.15)
    fig.suptitle('The multiplier changes the balance between M2 and stress after normalization');save(fig,'readout-tradeoff')
    # Paired fixed-candidate loss uncertainty; cheap blocks common but still
    # resampled when evaluating the nonlinear price-loss difference.
    B=np.load(OUT/'fresh/SPX_blocks.npz');base=frames[1];mask=base.priority_selected.to_numpy();q=base[mask].reset_index(drop=True);parity=base.option_type.eq('put').to_numpy()*(base.forward-base.strike).to_numpy();rng=np.random.default_rng(20261028423);cis=[]
    for c in range(2,7):
        deltas=[]
        for _ in range(512):
            ie=rng.integers(8,size=8);ic=rng.integers(8,size=8);C=B['cheap'][ic].mean(0);p0=B['c1_128'][ie].mean(0)+C-parity;p1=B[f'c{c}_128'][ie].mean(0)+C-parity;deltas.append(loss(q,p1[mask])-loss(q,p0[mask]))
        cis.append(dict(candidate=c,alpha=float(s[s.candidate==c].alpha.iloc[0]),loss_change_CI95=np.quantile(deltas,[.025,.975]).tolist()))
    dump(out/'paired_loss_intervals.json',cis);summary=[]
    for c in range(1,7):
        row=s[s.candidate==c].iloc[0];summary.append(dict(candidate=c,alpha=float(row.alpha),SPX_loss=float(row.loss),SPX_IV_MAE_pp=float(row.IV_MAE_pp),SPX_invalid_IV=int(row.invalid_IV),SPX_median_SE_scales=float(row.median_SE_scales),VIX_training_loss=ledger[c-1]['loss'],VIX_fresh_loss=vs[c]['option_loss'],training_max_forward_gap=ledger[c-1]['max_forward_gap'],checking_max_forward_gap=vs[c]['max_forward_gap']))
    pd.DataFrame(summary).to_csv(out/'comparison.csv',index=False);grid=[]
    for c in [1,3,5,6]:
        fine=frames[c];coarse=pd.read_csv(OUT/f'fresh/SPX_prices_c{c}_64.csv');p=fine.priority_selected;delta=abs(fine.model_price_U-coarse.model_price_U)/fine.price_scale_U;grid.append(dict(candidate=c,median_change_scales=float(np.median(delta[p])),maximum_change_scales=float(np.max(delta[p]))))
    dump(out/'grid_diagnostics.json',grid)
    table='\n'.join(f"| {r['alpha']:g} | {r['SPX_IV_MAE_pp']:.3f} | {r['SPX_loss']:.1f} | {r['VIX_fresh_loss']:.3f} | {r['SPX_invalid_IV']} |" for r in summary)
    report=f'''# Asymmetric whole-response multiplier: results, 28 October 2026

**Tested, not adopted.** The proposed readout gives modest SPX improvements at stronger settings, but damages VIX smiles and raises stock-pricing uncertainty. The VIX-preserving selection retains alpha=0. This is evidence about this fixed functional shape, this existing factor and this bounded seven-level refit—not a model-family impossibility result.

## Equations and exact scope

Keep all exogenous OU dynamics and the previous mixed stock rotation theta=.8. Use the existing feedback projection

$$X_t=(Y_t-E_Q[Y_t])/s_Y,\\qquad s_Y={design_scale():g}.$$

Its instantaneous covariance with the old equity driver is positive; the stock still loads .6 on that driver. Define

$$D_\\alpha(X)=\\exp[\\alpha\\log(1+e^{{-X}})],\\qquad H^\\alpha_t=\\frac{{H^0_tD_\\alpha(X_t)}}{{E_Q[H^0_tD_\\alpha(X_t)]}},\\qquad v_t=f+(\\xi_t-f)H^\\alpha_t.$$

All original M2 parameters, stress rates/amplitudes/weights, floor and nine off-window variance levels are fixed. Alpha=0,.25,.5,1,1.5,2 is shared across every maturity. For each strength, fit only the seven VIX-sensitive variance levels, bounded between the floor and .2, with all four training-estimator forward equalities. Include the original-curve baseline and an alpha=0 curve-refit control. Positive spike shift remains fixed. No endogenous feedback, extra bank, jump or quintic readout is introduced.

## Comparison

The table uses 1024 fresh full-model stock paths on the 128-step/trading-day grid and 4096 fresh full-factor VIX draws per expiry. VIX checking freezes training-fitted levels; it does not force checking forwards back onto target. All 1751 eligible/1294 primary SPX quotes at seven maturities and all 217 VIX quotes/four forwards remain present.

| Alpha | SPX IV MAE, pp | SPX interval loss | VIX interval loss | Invalid primary SPX IVs |
|---:|---:|---:|---:|---:|
{table}

Training forward gaps are below {max(z['max_forward_gap'] for z in ledger):.3g} VIX points for every strength. That is equality of a controlled training estimator, not an exact population-forward guarantee. Fresh forward gaps and both anchor/difference SEs are in `runs/m2_asymmetric_readout_001/fresh_VIX/VIX_forwards_c*.csv`. Signed negative price estimates and unavailable IVs are retained. IV MAE averages available IVs; the missing counts matter. The original unchanged-curve alpha0 baseline on these new stock draws has loss {float(s[s.candidate==0].loss.iloc[0]):.2f} and IV MAE {float(s[s.candidate==0].IV_MAE_pp.iloc[0]):.3f} pp, so comparisons should use this experiment's paired alpha0 controls, not a different experiment's Monte Carlo estimate.

The predeclared selection allowed VIX loss no more than alpha0-refit loss+.1 before minimizing replay SPX IV MAE. No positive tested strength met that VIX screen; candidate1, alpha0, was selected before fresh stocks. The fresh full-grid curves are exploratory response evidence, not independent validation of a newly selected winner. All bounded VIX fits completed with forward feasibility; several curve levels reach their bounds. Wider readout/stress/curve optimization was not attempted.

## Why the normalization matters

Since the stress factors are independent of the old bank,

$$E[H^0D]=(1-w)E[RD]+wE[D].$$

The original R already favours negative feedback states. D favours those states too. Thus E[RD] grows faster than E[D], and whole-product normalization changes the average M2/stress balance despite w staying fixed. At day153 the stress fraction of mean excess variance is 92.3% at alpha0, 85.6% at alpha.5, 72.9% at alpha1, 53.2% at alpha1.5 and 32.1% at alpha2. This mechanism is directly computed under the Gaussian law; attributing the entire VIX loss to it would require an additional controlled experiment.

The volatility reaction is stronger on bad-news states, but the normalization also shifts mean variance away from the stress mechanism that produced the successful VIX fit. Merely matching the four forwards does not recover the VIX distribution or smiles.

## Numerical evidence and limits

The stock sampler rebuilds candidate variances and stock paths from the same factor stream once per batch for every strength. Alpha0 full replay matches previous rotated stocks to {json.loads((OUT/'zero_multiplier_stock_replay.json').read_text())['max_absolute_stock_error']:.3g}. Independent NumPy direct tensor GH versus the C-convolution/YF implementation agrees to {audit['maximum_direct_tensor_quadrature_scaled_error']:.3g} scaled error; a scalar stock recursion agrees to {audit['scalar_stock_recursion_max_error']:.3g}; all seven forward derivatives agree within {audit['forward_derivative_max_absolute_error']:.3g}. Normalizer GH48 versus GH64 discrepancy is at most {max(json.loads((OUT/'normalizer_order_audit.json').read_text())['maximum_relative_GH48_vs64_by_alpha']):.3g} at declared probes. Midpoint time interpolation error is about1.1e-7. These are sampled numerical checks, not rigorous global bounds.

Stock uncertainty uses antithetic pairs and eight independent batches; auxiliary reference uses the previous65536 paths and is shared across screening/checking. Median stock SE rises from about1.66 quote scales at alpha0 to7.88 at alpha1.5; strong-readout numerical precision is not sufficient for a final calibration claim. Coupled64/128 grid changes, paired bootstrap loss intervals, stock means and invalid prices are saved. Eight batches underpin approximate uncertainty, not a certified tail bound.

VIX training uses1024 randomized rows from the earlier65536 state/stress stream and a controlled old-model anchor. Follow-up VIX checking uses genuinely fresh joint38-dimensional Gaussian-mixture draws (eight independent scrambled Sobol batches), with exact likelihood weights and no sample-mass normalization; all strengths share those draws. The old anchor price estimate is still shared, and its uncertainty is retained separately and combined with the independent new correction. Weighted-tail ESS/concentration are saved for every strength and maturity. GH32/time8 checking and matched-state GH48/time12 probes diagnose deterministic VIX error. The same market snapshot is used throughout.

Existing Lean algebraic identities were not changed. The new quadrature/sampler correspondences are numerical audits; neither continuous stock martingality nor a finite-wing-skew response theorem is claimed. No production baseline was overwritten.

## Conclusion to retain

The first whole-response asymmetric multiplier does not solve the joint-fit problem in this bounded trial. Strong settings modestly improve SPX while spoiling VIX and making stock estimation harder. Do not repeat the same strength scan or claim that a lower available-IV MAE alone gives an acceptable double calibration.

A possible next refinement is to keep each component's mean variance allocation fixed when adding asymmetry, using the correct joint normalization inside the M2 component. That is a separate proposed model, not tested or adopted here. A different Gaussian factor/timescale is another hypothesis; this experiment does not establish that all exogenous readouts fail.

Figures: [SPX smiles](visualizations/asymmetric-readout028/SPX-asymmetric-smiles.png), [VIX smiles](visualizations/asymmetric-readout028/VIX-asymmetric-smiles.png), [response and variance allocation](visualizations/asymmetric-readout028/readout-tradeoff.png). Detailed results: `runs/m2_asymmetric_readout_001`.
'''
    (HERE/'ASYMMETRIC_READOUT_RESULTS.md').write_text(report);dump(out/'summary.json',dict(comparison=summary,paired_loss_intervals=cis,grid=grid,adopted=False));print('ASYMMETRIC REPORT COMPLETE',flush=True)

def design_scale():return json.loads((OUT/'design.json').read_text())['scale']

from .vix_only_curve_fit import DAYS
if __name__=='__main__':main()
