"""Protected-fit native SPX/VIX smiles and unchanged-window variance curves."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .full_readout_fit import HERE,dump
from .plot_joint_curve_fit import implied_vol
from .two_speed_joint_core import physical
from .two_speed_precision import FREE
from .vix_only_curve_fit import PRIOR


def plot(run,out):
    out.mkdir(parents=True,exist_ok=False);fresh=run/'fresh_check';baseline=pd.read_csv(fresh/'prices_baseline_128.csv');selected=pd.read_csv(fresh/'prices_selected_128.csv');days=sorted(baseline.expiry_days.unique());frames=[];metrics=[]
    for name,q in [('baseline',baseline),('selected',selected)]:
        q=q.copy();q['market_native_IV']=[implied_vol(r.price_mid_U,r.forward,r.strike,r.T_years,r.option_type) for r in q.itertuples()];q['model_native_IV']=[implied_vol(r.model_price_U,r.forward,r.strike,r.T_years,r.option_type) for r in q.itertuples()];q['model_IV_lower']=[implied_vol(r.model_price_U-1.96*r.mc_se_U,r.forward,r.strike,r.T_years,r.option_type) for r in q.itertuples()];q['model_IV_upper']=[implied_vol(r.model_price_U+1.96*r.mc_se_U,r.forward,r.strike,r.T_years,r.option_type) for r in q.itertuples()];q['invalid_model_IV']=~np.isfinite(q.model_native_IV);q.to_csv(out/f'native_IV_{name}_all.csv',index=False);frames.append(q)
        for day,g in q[q.priority_selected].groupby('expiry_days'):
            metrics.append(dict(model=name,day=int(day),quotes=len(g),loss=float(np.mean(g.interval_residual**2)),IV_MAE_pp=float(np.nanmean(abs(g.model_native_IV-g.market_native_IV))),inside=int(g.interval_residual.eq(0).sum()),negative=int((g.model_price_U<0).sum()),invalid_IV=int(g.invalid_model_IV.sum()),median_SE_scales=float(np.median(g.mc_se_U/g.price_scale_U))))
    pd.DataFrame(metrics).to_csv(out/'SPX_metrics.csv',index=False)
    fig,axes=plt.subplots(2,4,figsize=(16,8))
    for ax,day in zip(axes.flat,days):
        a=frames[0][(frames[0].expiry_days==day)&frames[0].priority_selected];b=frames[1][(frames[1].expiry_days==day)&frames[1].priority_selected]
        ax.scatter(a.K_over_F,a.market_native_IV,color='black',s=7,alpha=.55,label='Market')
        ax.scatter(a.K_over_F,a.model_native_IV,color='#a0a0a0',s=7,alpha=.65,label='VIX baseline')
        ax.scatter(b.K_over_F,b.model_native_IV,color='#267b96',s=7,alpha=.75,label='Protected curve fit')
        for r in b.iloc[::max(1,len(b)//20)].itertuples():
            if np.isfinite(r.model_IV_lower) and np.isfinite(r.model_IV_upper):ax.plot([r.K_over_F]*2,[r.model_IV_lower,r.model_IV_upper],color='#267b96',alpha=.4,lw=.8)
        ax.set_title(f'{day} days');ax.set_xlabel('Strike / forward');ax.set_ylabel('Native implied vol (%)');ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=8);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,'Fresh 512 full-model paths\n256 antithetic pairs\n128 steps / trading day\n\nBlue error bars: ±1.96 MC SE\nTimestep / normalizer error excluded\n\nOnly nine variance levels changed\nVIX model unchanged',va='top',fontsize=11);fig.suptitle('Protected SPX fit: improvement without sacrificing VIX',fontsize=15);fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(out/f'SPX-protected-smiles.{ext}',dpi=170)
    plt.close(fig)
    q=pd.read_csv(HERE/'runs/m2_two_speed_001/fresh_check/native_iv_two_speed.csv');fig,axes=plt.subplots(1,4,figsize=(15,3.9))
    for ax,(day,g) in zip(axes,q.groupby('expiry_days',sort=True)):
        ax.scatter(g.strike,g.market_native_IV,color='black',s=12,label='Market');ax.plot(g.strike,g.model_native_IV,color='#267b96',marker='.',lw=1,label='Baseline = protected fit');ax.set_title(f'{day} days');ax.set_xlabel('VIX strike');ax.set_ylabel('Native implied vol (%)');ax.grid(alpha=.15)
    axes[0].legend(fontsize=8);fig.suptitle('VIX smiles are inherited unchanged: identical VIX model inputs');fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(out/f'VIX-protected-smiles.{ext}',dpi=170)
    plt.close(fig)
    base=np.array(json.loads((run/'design.json').read_text())['baseline']);chosen=np.array(json.loads((run/'selection.json').read_text())['parameters']);ends=365*np.array(json.loads((PRIOR/'run.json').read_text())['ends_years']);a,f=physical(base);b,_=physical(chosen);table=pd.DataFrame(dict(index=np.arange(16),start_days=np.r_[0,ends[:-1]],end_days=ends,free=np.isin(np.arange(16),FREE),baseline_xi=a,candidate_xi=b,baseline_vol_percent=100*np.sqrt(a),candidate_vol_percent=100*np.sqrt(b)));table.to_csv(out/'curve_levels.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    for ax,horizon in zip(axes,[211,1030]):
        for levels,color,label in [(a,'#a0a0a0','VIX baseline'),(b,'#267b96','Protected curve fit')]:ax.stairs(100*np.sqrt(levels),np.r_[0,ends],color=color,label=label,lw=1.7)
        for left,right in [(27,57),(55,85),(118,148),(153,183)]:ax.axvspan(left,right,color='#e9a23b',alpha=.1)
        ax.axhline(100*np.sqrt(f),color='#ad2929',linestyle='--',lw=.8,label='Frozen floor');ax.set_xlim(0,horizon);ax.set_xlabel('Days');ax.set_ylabel('sqrt(forward variance), %');ax.grid(alpha=.15)
    axes[0].legend(fontsize=8);fig.suptitle('Only off-window levels move; shaded VIX windows remain fixed');fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(out/f'protected-variance-curve.{ext}',dpi=170)
    plt.close(fig)
    gridrows=[]
    for name in ['baseline','selected']:
        g=pd.read_csv(fresh/f'paired_grid_{name}.csv');p=g[g.priority_selected]
        gridrows.append(dict(model=name,median_absolute_grid_change_scales=float(np.median(abs(p.fine_minus_coarse_U)/p.price_scale_U)),max_absolute_grid_change_scales=float(np.max(abs(p.fine_minus_coarse_U)/p.price_scale_U)),median_grid_SE_scales=float(np.median(p.paired_SE_U/p.price_scale_U))))
    dump(out/'grid_summary.json',gridrows)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();plot(args.run.resolve(),args.output.resolve())
