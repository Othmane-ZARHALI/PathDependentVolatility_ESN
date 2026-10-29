"""All-maturity smiles, paired uncertainty, and reproducible slow-readout conclusion."""
import hashlib,json,shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .slow_readout import HERE,OUT,arr,dump
from .protected_curve_check import loss
from .plot_joint_curve_fit import implied_vol
WEAK=HERE/'runs/m2_slow_readout_002'
CHECK=WEAK/'fresh_validation'
PLOT=HERE/'visualizations/slow-readout028'
NAMES=['Baseline','Slow signal, fixed curve','Slow signal + curve refit','Curve refit alone','Stronger signal + refit']
COLORS=['#888888','#1d8da4','#d17426','#846bb0','#b84b66']

def main():
    PLOT.mkdir(exist_ok=False,parents=True);sel=json.loads((CHECK/'selection.json').read_text());vix=json.loads((CHECK/'VIX/scores.json').read_text());scores=pd.read_csv(CHECK/'SPX/SPX_scores.csv');fine=scores[scores.steps==256].set_index('candidate');frames={c:pd.read_csv(CHECK/f'SPX/SPX_prices_c{c}_256.csv') for c in range(5)};days=sorted(frames[0].expiry_days.unique());metrics=[]
    def save(fig,name):
        fig.tight_layout()
        for ext in ['png','svg']:fig.savefig(PLOT/f'{name}.{ext}',dpi=170)
        plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(17,8))
    for ax,day in zip(axes.flat,days):
        g=frames[0][(frames[0].expiry_days==day)&frames[0].priority_selected];ax.scatter(g.K_over_F,g.market_IV,s=9,c='black',alpha=.6,label='Market')
        for c in range(5):
            a=frames[c][(frames[c].expiry_days==day)&frames[c].priority_selected];ax.scatter(a.K_over_F,a.model_IV,s=8,c=COLORS[c],alpha=.65,label=NAMES[c]);lo=a.loc[(a.K_over_F-.9).abs().idxmin()];hi=a.loc[(a.K_over_F-1.1).abs().idxmin()];metrics.append(dict(candidate=c,day=int(day),IV_MAE_pp=float(np.nanmean(abs(a.model_IV-a.market_IV))),invalid_IV=int(a.model_IV.isna().sum()),stock_mean=float(a.stock_mean.iloc[0]),stock_mean_SE=float(a.stock_mean_SE.iloc[0]),left_K_over_F=float(lo.K_over_F),right_K_over_F=float(hi.K_over_F),model_right_minus_left_IV_pp=float(hi.model_IV-lo.model_IV),market_right_minus_left_IV_pp=float(hi.market_IV-lo.market_IV)))
        ax.set(title=f'{day} days',xlabel='Strike / forward',ylabel='Implied volatility (%)');ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=6);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,'1024 new stock paths\n256 steps / trading day\nAll seven SPX maturities\n\nWeak signal: 9.48-day half-life, alpha=.125\nStrong contrast: 3.33 days, alpha=.5\nBoth already present in frozen OU bank\n\nOriginal fast M2 retained\nNo endogenous factor feedback\nSigned prices and IV failures retained',va='top',fontsize=10);fig.suptitle('SPX: slower downside signal from the existing bank');save(fig,'SPX-slow-readout');pd.DataFrame(metrics).to_csv(PLOT/'maturity_metrics.csv',index=False)
    vframes={};fig,axes=plt.subplots(1,4,figsize=(17,4.5))
    for c in range(5):
        a=pd.read_csv(CHECK/f'VIX/VIX_prices_c{c}.csv');a['market_IV']=[implied_vol(r.price_mid_U,r.forward,r.strike,r.T_years,r.option_type) for r in a.itertuples()];a['model_IV']=[implied_vol(r.model_price_U,r.forward,r.strike,r.T_years,r.option_type) for r in a.itertuples()];a.to_csv(CHECK/f'VIX/native_IV_c{c}.csv',index=False);vframes[c]=a
    for ax,(day,g) in zip(axes,vframes[0].groupby('expiry_days',sort=True)):
        ax.scatter(g.strike,g.market_IV,c='black',s=12,label='Market')
        for c in range(5):
            a=vframes[c][vframes[c].expiry_days==day];ax.scatter(a.strike,a.model_IV,c=COLORS[c],s=12,label=NAMES[c])
        ax.set(title=f'{day} days',xlabel='VIX strike',ylabel='Implied volatility (%)');ax.grid(alpha=.15)
    axes[0].legend(fontsize=6);fig.suptitle('VIX: preserving conditional variance means does not preserve the smiles');save(fig,'VIX-slow-readout')
    records=[]
    for root,label in [(OUT,'strong'),(WEAK,'weak')]:
        ss=pd.read_csv(root/'SPX/SPX_scores.csv');vv=json.loads((root/'checking_VIX/scores.json').read_text());cc=json.loads((root/'candidates.json').read_text())
        for r in ss[ss.steps==256].itertuples():
            c=int(r.candidate);records.append(dict(stage=label,candidate=c,name=cc['names'][c],SPX_loss=float(r.loss),SPX_IV_MAE_pp=float(r.IV_MAE_pp),VIX_loss=vv[c]['option_loss'],max_forward_gap=vv[c]['max_forward_gap']))
    pd.DataFrame(records).to_csv(PLOT/'pilot_grid.csv',index=False);fig,ax=plt.subplots(figsize=(7,5))
    for stage,color in [('strong','#b84b66'),('weak','#1d8da4')]:
        rr=[r for r in records if r['stage']==stage and r['candidate']>0];ax.scatter([r['VIX_loss'] for r in rr],[r['SPX_IV_MAE_pp'] for r in rr],c=color,label=stage+' signals')
    ax.scatter(records[0]['VIX_loss'],records[0]['SPX_IV_MAE_pp'],c='black',marker='x',s=70,label='Baseline');ax.set(xscale='log',xlabel='VIX interval loss (smaller is better)',ylabel='SPX IV error, percentage points',title='Exploratory grid: stronger SPX changes come with VIX cost');ax.legend();ax.grid(alpha=.15);save(fig,'tradeoff')
    B=np.load(CHECK/'SPX/SPX_blocks.npz');base=frames[0];mask=base.priority_selected.to_numpy();pq=base[mask].reset_index(drop=True);parity=base.option_type.eq('put').to_numpy()*(base.forward-base.strike).to_numpy();rng=np.random.default_rng(20261028614);intervals=[];grid=[]
    for c in range(5):
        co=pd.read_csv(CHECK/f'SPX/SPX_prices_c{c}_128.csv');dv=abs(frames[c].model_price_U-co.model_price_U)/base.price_scale_U;grid.append(dict(candidate=c,median_change_quote_scales=float(np.median(dv[mask])),max_change_quote_scales=float(np.max(dv[mask]))))
        if c==0:continue
        changes=[]
        for _ in range(1024):
            ie=rng.integers(8,size=8);ic=rng.integers(8,size=8);C=B['cheap'][ic].mean(0);p=(B[f'c{c}_256'][ie].mean(0)+C-parity)[mask];p0=(B['c0_256'][ie].mean(0)+C-parity)[mask];changes.append(loss(pq,p)-loss(pq,p0))
        intervals.append(dict(candidate=c,loss_delta=float(fine.loc[c,'loss']-fine.loc[0,'loss']),approximate_batch_bootstrap_CI95=np.quantile(changes,[.025,.975]).tolist()))
    dump(PLOT/'paired_intervals.json',intervals);dump(PLOT/'grid_check.json',grid)
    # Raw target-law payoff concentration; antithetic dependence remains explicit.
    stocks=np.load(CHECK/'stocks.npz')['stocks'];tails=[]
    for c in range(5):
        for j,day in enumerate(days):
            a=base[(base.expiry_days==day)&base.priority_selected];S=stocks[:,j,c,1]
            for side in ['put','call']:
                q=a[a.option_type==side];K=float(q.K_over_F.min() if side=='put' else q.K_over_F.max());pay=np.maximum(K-S,0) if side=='put' else np.maximum(S-K,0);tot=pay.sum();tails.append(dict(candidate=c,day=int(day),side=side,K_over_F=K,raw_payoff_mean=float(pay.mean()),payoff_concentration_ESS=float(tot*tot/(pay@pay)) if tot else 0.,top10_share=float(np.sort(pay)[-10:].sum()/tot) if tot else 0.))
    dump(CHECK/'stock_tail_diagnostics.json',tails);audit=json.loads((OUT/'implementation_audit.json').read_text());quad=json.loads((CHECK/'quadrature_probe.json').read_text());m=pd.DataFrame(metrics);m29=m[m.day==29].set_index('candidate');table='\n'.join(f"| {NAMES[c]} | {fine.loc[c,'loss']:.2f} | {fine.loc[c,'IV_MAE_pp']:.3f} | {vix[c]['option_loss']:.3f} | {vix[c]['max_forward_gap']:.4f} | {int(fine.loc[c,'invalid_IV'])} |" for c in range(5));gain=100*(1-fine.loc[2,'loss']/fine.loc[0,'loss']);train=json.loads((WEAK/'VIX_fit_ledger.json').read_text())[3]
    report=r'''# Slower downside signals from the existing OU bank — 28 October 2026

The slower signals produce a modest SPX improvement, but this particular multiplicative readout does not solve the joint fit. Strong signals improve SPX more while damaging VIX. A weak signal avoids most of that damage but leaves much of the SPX skew missing. Production model and previous results are unchanged.

## Model in plain terms

The existing fast signal continues to react to recent shocks. Add a separate slower signal, so that a negative shock can raise volatility for longer. No OU process is added or given feedback from the simulated price or volatility. Select existing feedback modes nearest in log half-life to generic targets 2, 8 and 32 days. Their actual half-lives are 3.329, 9.481 and 27.004 days; the same construction applies at every option maturity.

Let U_j be the selected existing mode and X_j=U_j-E_Q[U_j]. Each mode has stationary variance one; its dated conditional variance is not assumed to be one. For one mode at a time,

$$dX_j=-\kappa_jX_jdt+\sqrt{2\kappa_j}\,dB_j,\qquad h_j=\log(2)/\kappa_j,$$

$$D_j(X_j)=\exp\{\alpha\log(1+e^{-X_j})\}.$$

D increases as the state becomes more negative. Alpha controls its strength. The current fast M2 component, including its asymmetric positive spike shift, stays exactly in place:

$$\eta=c_CC+q_Y(Y-\beta_Y)^2+q_J(F-\beta_J)^2,\quad h=\ell-k^{-1}\log(1+e^{k(\ell-\eta)}),\quad R=e^{2h}/E_Q[e^{2h}].$$

Keep both lognormal stress factors and their frozen parameters, with A=(1-r)A_s+rA_f and a=1-w. The old multiplier is H0=aR+wA. The new one is

$$m_D(t)=E_Q[D_j],\qquad m_{RD}(t)=E_Q[RD_j],$$

$$H_j=(aR+wA)D_j\frac{a+wA}{a\,m_{RD}(t)+wA\,m_D(t)},\qquad v_t=f+(\xi(t)-f)H_j.$$

R and D are correlated: E[RD] cannot be replaced by E[R]E[D]. The whole old Gaussian bank remains independent of the stress bank. Consequently,

$$E_Q[H_j\mid A_s,A_f]=a+wA,\qquad E_Q[H_j]=1,\qquad E_Q[v_t]=\xi(t).$$

This preserves the conditional mean variance allocated to each stress regime. It does not preserve separate internal component budgets or the entire VIX law. VIX depends on the conditional future variance and a square root, so matching E[v] does not by itself match VIX forwards or options.

Stock loading stays fixed at the previous mixed theta=.8 value:

$$dW^S=.6\,dW^0-.8\{\sqrt{1-r}\,dB_s+\sqrt r\,dB_f\}.$$

The original M2 cap is unchanged. No new cap or maturity-specific exception is introduced. Reusable Gaussian factors are distinguished from rebuilt readouts, normalizers, stock paths and VIX exposures.

## Experiment and controls

First test alpha .5 and 1 for all three memories. Their VIX training losses were too large, so add a reactive exploratory grid .125 and .25. These grids were not jointly preregistered. Both stages share identical 512-path stock draws, eight independent antithetic batches, and identical 512-state-per-expiry training samples. Evaluate the fixed variance curve and refit the same seven active variance levels under four training-forward equalities; nine off-window levels remain frozen. All seven SPX maturities/1751 eligible/1294 primary quotes and all 217 VIX quotes/four forwards remain.

Separate previously examined 4096-state-per-expiry VIX checking is saved for every candidate. A weak candidate is frozen before additional new stock and outer-state draws: require training-forward equality and training VIX loss within .1 of the alpha-zero refit, then minimize the pilot's fine-grid SPX IV MAE among weak refits. This selects the 9.481-day mode, alpha=.125. The .1 threshold is a diagnostic tolerance, not evidence that the resulting VIX smile is satisfactory.

Additional checking uses 1024 newly generated stock paths and 2048 new full-factor importance-sampled VIX states per expiry. The five comparisons are baseline, selected slow signal with fixed curve, selected slow signal with refitted curve, that same refitted curve without the slow signal, and a stronger 3.329-day alpha=.5 refitted contrast. The strong contrast is illustrative, not a second newly selected winner. All four training forwards satisfy numerical equalities; their separate checking gaps below are retained.

## Independent checking results

| Configuration | SPX interval loss | SPX IV error, pp | VIX interval loss | Max VIX forward gap, points | Invalid primary SPX IVs |
|---|---:|---:|---:|---:|---:|
TABLE

IV error equally averages available native-IV errors over seven maturities. Invalid IVs and signed negative price estimates are retained; native-IV metrics alone are not a fit certificate. Forward gaps include the inherited baseline estimator gap. Training equalities are not certified exact population forwards.

For the weak refit, SPX price loss falls GAIN percent in this checking sample. The curve-only control and fixed-curve signal show how much comes from each change. At 29 days, the quoted right-minus-left IV gap near strike/forward .9 and 1.1 is SLOPES percentage points, respectively, for baseline, weak fixed, weak refit, curve alone and stronger refit, versus MARKET in the market. The required downward slope remains missing.

Paired sampling intervals: INTERVALS. These approximate 1024-resample bootstrap intervals use eight antithetic batch means with paired candidate/baseline resampling and independent resampling of the reused cheap reference batches. They exclude deterministic quadrature/grid errors and market uncertainty. New full-model draws do not make the already reused 65536-path auxiliary reference independent of prior research. This is the same market snapshot, not market out-of-sample validation.

## Numerical correspondence and accuracy

Conditional moments use the full Gaussian state law projected onto C,Y,F and the chosen slower mode. Integrate independent C by a 4096-interval cubic-Hermite convolution of the capped readout, Y/F by tensor Gaussian quadrature, and X conditional on Y/F by one-dimensional Gaussian quadrature and an 8192-interval interpolation table. Direct one-dimensional marginal quadrature computes E[D]. This preserves the R/D correlation. Product normalizers use GH48 and adaptive logarithmic-time cubic splines; conditional VIX checking uses GH32 and time order8 with positive stress quadrature32.

Alpha-zero full-horizon stock replay error is REPLAY. Independent scalar conditional tensor quadrature relative errors for the three memories are ORACLE. GH48/64 product-normalizer checks are below ORDER. An independent 97-step scalar stock recurrence for five candidates, both grids and both antithetic partners has error RECURRENCE. These are numerical implementation checks, not Lean certification or uniform error bounds.

Coupled stock grids use 128 and 256 steps per trading day; their median absolute price changes in quote-scale units are GRID. A matched 64-state-per-expiry conditional VIX refinement, GH32/time8 to GH48/time12, gives maximum VIX change QUAD points and maximum absolute weighted mean change QMEAN points. These probes do not bound unobserved tails. Raw stock payoff concentration, importance-weighted VIX tail concentration, sample mass, stock means and SEs, quote membership, forward SEs and negative/invalid estimates are preserved in the run. No clipping, smoothing, sample-weight normalization or sample-stock-forward renormalization.

Use the previously validated stress-only auxiliary stock reference with its conditional lognormal call prices and fixed deterministic delta controls. Adding the new readout changes the expensive target payoff correction, not the auxiliary model. Forward and option checking SEs combine the independent fresh difference and old anchor SEs; separate components remain available.

## Durable conclusion

The slow-mode hypothesis has now been tested for three generic existing memories and four strengths, with fixed-curve, constrained-curve and curve-only controls. This simple one-sided multiplier supplies some downside response, but stronger changes also alter the VIX distribution; protecting the first moment does not protect option smiles. Do not interpret the weak aggregate improvement as a solved SPX skew or repeat the same strength scan as the next default step. The result is specific to this readout and tested range, not a no-go theorem for slower projections or exogenous Gaussian banks.

No candidate is adopted. The production baseline and original fast/Zumbach-motivated readout remain intact. No endogenous return/volatility feedback, new OU bank, new constant-theta scan or quintic/jump extension has been introduced.

Figures: [all seven SPX smiles](visualizations/slow-readout028/SPX-slow-readout.png), [all four VIX smiles](visualizations/slow-readout028/VIX-slow-readout.png), [full exploratory trade-off](visualizations/slow-readout028/tradeoff.png). Data: `runs/m2_slow_readout_001` and `runs/m2_slow_readout_002`, with independent checking in the latter's `fresh_validation` folder.
'''
    replacements={'TABLE':table,'GAIN':f'{gain:.2f}','SLOPES':', '.join(f"{m29.loc[c,'model_right_minus_left_IV_pp']:.2f}" for c in range(5)),'MARKET':f"{m29.loc[0,'market_right_minus_left_IV_pp']:.2f}",'INTERVALS':json.dumps(intervals),'REPLAY':f"{audit['alpha_zero_existing_stock_max_error']:.3g}",'ORACLE':json.dumps(audit['scalar_conditional_tensor_GH24_relative_errors']),'ORDER':f"{max(json.loads((r/'moment_order_audit.json').read_text())['max_GH48_64_relative_error'] for r in [OUT,WEAK]):.3g}",'RECURRENCE':str(json.loads((WEAK/'recurrence_audit.json').read_text())['maximum_scalar_stock_error']),'GRID':', '.join(f"{g['median_change_quote_scales']:.3f}" for g in grid),'QUAD':f"{max(q['max_VIX_change'] for q in quad):.4g}",'QMEAN':f"{max(abs(q['weighted_mean_VIX_change']) for q in quad):.4g}"}
    for a,b in replacements.items():report=report.replace(a,b)
    (HERE/'SLOW_READOUT_RESULTS.md').write_text(report);dump(PLOT/'summary.json',dict(scores=fine.reset_index().to_dict('records'),VIX_scores=vix,price_loss_reduction_percent=gain,paired_intervals=intervals,grid=grid,adopted=False));print('SLOW REPORT COMPLETE',gain,flush=True)
if __name__=='__main__':main()
