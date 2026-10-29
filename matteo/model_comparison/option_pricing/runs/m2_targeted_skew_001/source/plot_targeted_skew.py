"""Paired SPX/VIX diagnostics for the two isolated mechanisms."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .targeted_skew import OUT,HERE,ALPHA
from .full_readout_fit import dump
from .protected_curve_check import loss
from .plot_joint_curve_fit import implied_vol
from .vix_only_curve_fit import DAYS

def main():
    out=HERE/'visualizations/targeted-skew028';out.mkdir(exist_ok=False,parents=True);config=json.loads((OUT/'candidates.json').read_text());scores=pd.read_csv(OUT/'fresh/SPX_scores.csv');scores=scores[scores.steps==128];vs={r['candidate']:r for r in json.loads((OUT/'checking_VIX/scores.json').read_text())};vs[1]=vs[0];selected=[0,1,4,9];labels=['Original','Adaptive stock loading','Regime α=1.5, fixed curve','Regime α=1.5, refitted curve'];colors=['#999999','#147d96','#b6732a','#8756a1'];mainframes={c:pd.read_csv(OUT/f'fresh/SPX_prices_c{c}_128.csv') for c in selected};refined={c:pd.read_csv(OUT/f'grid_refinement/prices/SPX_prices_c{j}_256.csv') for j,c in enumerate(selected)};days=sorted(mainframes[0].expiry_days.unique());metrics=[]
    def save(fig,name):
        fig.tight_layout()
        for ext in ['png','svg']:fig.savefig(out/f'{name}.{ext}',dpi=170)
        plt.close(fig)
    def smiles(frames,name,title,note):
        fig,axes=plt.subplots(2,4,figsize=(16,8))
        for ax,day in zip(axes.flat,days):
            base=frames[0][(frames[0].expiry_days==day)&frames[0].priority_selected];ax.scatter(base.K_over_F,base.market_IV,color='black',s=8,alpha=.55,label='Market')
            for c,label,color in zip(selected,labels,colors):
                a=frames[c][(frames[c].expiry_days==day)&frames[c].priority_selected];ax.scatter(a.K_over_F,a.model_IV,color=color,s=7,alpha=.75,label=label)
            ax.set_title(f'{day} days');ax.set_xlabel('Strike / forward');ax.set_ylabel('Implied volatility (%)');ax.grid(alpha=.15)
        axes.flat[0].legend(fontsize=6.4);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,note,va='top');fig.suptitle(title);save(fig,name)
    smiles(mainframes,'SPX-targeted-128','Two isolated mechanisms: main paired SPX check','1024 fresh stock paths\n128 steps / trading day\nAll seven maturities\nSame auxiliary reference sample\n\nVIX unchanged for adaptive stock\nRegime normalization reprices VIX\nSigned prices and invalid IVs retained\nSampling error remains material')
    smiles(refined,'SPX-targeted-256','Finer-grid check: adaptive stock loading and conditional normalization','2048 fresh stock paths\n256 steps / trading day\nPaired128/256 grids\nFrozen illustrative α=1.5\n\nSame market snapshot\nNew factor construction on finer grid\nNo candidate reselection\nSigned estimates retained')
    for stage,frames in [('main128',mainframes),('refined256',refined)]:
        for c,frame in frames.items():
            for day,g in frame[frame.priority_selected].groupby('expiry_days'):
                lo=g.loc[(g.K_over_F-.9).abs().idxmin()];hi=g.loc[(g.K_over_F-1.1).abs().idxmin()];metrics.append(dict(stage=stage,candidate=c,name=config['names'][c],day=int(day),IV_MAE_pp=float(np.nanmean(abs(g.model_IV-g.market_IV))),invalid_IV=int(g.model_IV.isna().sum()),stock_mean=float(g.stock_mean.iloc[0]),stock_mean_SE=float(g.stock_mean_SE.iloc[0]),left_K_over_F=float(lo.K_over_F),right_K_over_F=float(hi.K_over_F),model_right_minus_left_IV_pp=float(hi.model_IV-lo.model_IV),market_right_minus_left_IV_pp=float(hi.market_IV-lo.market_IV)))
    pd.DataFrame(metrics).to_csv(out/'SPX_maturity_metrics.csv',index=False)
    fig,axes=plt.subplots(1,4,figsize=(16,4.3));vixframes=[]
    for ax,day in zip(axes,DAYS):
        for c,label,color in [(0,'Original = adaptive stock','#999999'),(4,'Regime1.5 fixed','#b6732a'),(9,'Regime1.5 refit','#8756a1')]:
            g=pd.read_csv(OUT/f'checking_VIX/VIX_prices_c{c}.csv');g=g[g.expiry_days==day].copy();g['model_IV']=[implied_vol(r.model_price_U,r.forward,r.strike,r.T_years,r.option_type) for r in g.itertuples()];g['market_IV']=[implied_vol(r.price_mid_U,r.forward,r.strike,r.T_years,r.option_type) for r in g.itertuples()];ax.plot(g.strike,g.model_IV,'.-',color=color,lw=1,label=label);g['candidate']=c;vixframes.append(g)
            if c==0:ax.scatter(g.strike,g.market_IV,s=10,color='black',label='Market')
        ax.set_title(f'{day} days');ax.set_xlabel('VIX strike');ax.set_ylabel('Implied volatility (%)');ax.grid(alpha=.15)
    axes[0].legend(fontsize=7);fig.suptitle('VIX: stock loading preserves the law; conditional normalization preserves less\nSupplied market forward used for IV; checking forward gaps saved separately');save(fig,'VIX-targeted');pd.concat(vixframes,ignore_index=True).to_csv(out/'VIX_native_IV.csv',index=False)
    # Entire readout response grid: fixed curves and forward-constrained refits.
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for candidates,label,color in [(range(2,6),'Original curve','#b6732a'),(range(6,11),'Forward-constrained refit','#8756a1')]:
        a=scores[scores.candidate.isin(candidates)];axes[0].plot(a.alpha,a.IV_MAE_pp,'.-',color=color,label=label);axes[1].plot(a.alpha,[vs[int(c)]['option_loss'] for c in a.candidate],'.-',color=color,label=label)
    axes[0].axhline(float(scores[scores.candidate==0].IV_MAE_pp.iloc[0]),color='#999999',ls='--');axes[1].axhline(vs[0]['option_loss'],color='#999999',ls='--')
    for ax in axes:ax.set_xlabel('α');ax.grid(alpha=.15);ax.legend(fontsize=8)
    axes[0].set_ylabel('SPX IV MAE, percentage points');axes[1].set_ylabel('VIX price interval loss');fig.suptitle('Stress-regime normalization: finite-grid response, not an adopted fit');save(fig,'regime-response')
    grids=[];intervals=[];rng=np.random.default_rng(20261028673)
    for stage,folder,steps,frames,cmap in [('main',OUT/'fresh',[64,128],mainframes,{c:c for c in selected}),('refined',OUT/'grid_refinement/prices',[128,256],refined,{c:j for j,c in enumerate(selected)})]:
        B=np.load(folder/'SPX_blocks.npz');base=frames[0];mask=base.priority_selected.to_numpy();q=base[mask].reset_index(drop=True);parity=base.option_type.eq('put').to_numpy()*(base.forward-base.strike).to_numpy();E0=B[f'c0_{steps[1]}']
        for c in selected:
            j=cmap[c];coarse=pd.read_csv(folder/f'SPX_prices_c{j}_{steps[0]}.csv');fine=frames[c];delta=abs(fine.model_price_U-coarse.model_price_U)/fine.price_scale_U;grids.append(dict(stage=stage,candidate=c,coarse_steps=steps[0],fine_steps=steps[1],median_grid_change_scales=float(np.median(delta[mask])),maximum_grid_change_scales=float(np.max(delta[mask]))))
            if c==0:continue
            E=B[f'c{j}_{steps[1]}'];diffs=[]
            for _ in range(512):
                ie=rng.integers(8,size=8);ic=rng.integers(8,size=8);C=B['cheap'][ic].mean(0);p0=E0[ie].mean(0)+C-parity;p1=E[ie].mean(0)+C-parity;diffs.append(loss(q,p1[mask])-loss(q,p0[mask]))
            intervals.append(dict(stage=stage,candidate=c,loss_delta_CI95=np.quantile(diffs,[.025,.975]).tolist()))
    dump(out/'grid_diagnostics.json',grids);dump(out/'paired_loss_intervals.json',intervals)
    summary=[]
    for c in range(11):
        z=scores[scores.candidate==c].iloc[0];summary.append(dict(candidate=c,name=config['names'][c],alpha=float(z.alpha),SPX_loss=float(z.loss),SPX_IV_MAE_pp=float(z.IV_MAE_pp),invalid_IV=int(z.invalid_IV),SPX_median_SE_scales=float(z.median_SE_scales),VIX_loss=vs[c]['option_loss'],VIX_max_forward_gap=vs[c]['max_forward_gap']))
    pd.DataFrame(summary).to_csv(out/'comparison.csv',index=False);fine_scores=pd.read_csv(OUT/'grid_refinement/prices/SPX_scores.csv');fine_scores=fine_scores[fine_scores.steps==256];dump(out/'summary.json',dict(main=summary,refined=fine_scores.to_dict('records'),grid=grids,paired=intervals,adopted=False));print('TARGETED PLOTS COMPLETE',flush=True)

if __name__=='__main__':main()
