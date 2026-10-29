"""Frozen-model precision comparison and native SPX smiles with MC error bars."""
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


def run(out):
    out.mkdir(parents=True,exist_ok=False);parent=HERE/'runs/m2_two_speed_joint_001';dynamic=HERE/'runs/m2_two_speed_precision_002';conditional=HERE/'runs/m2_two_speed_precision_conditional_002'
    old=pd.read_csv(parent/'checking/prices_spx_vix_only.csv');mid=pd.read_csv(dynamic/'check/prices_dynamic.csv');new=pd.read_csv(conditional/'check/prices_conditional.csv');days=sorted(new.expiry_days.unique());rows=[]
    for day in days:
        a=old[old.expiry_days==day];b=mid[(mid.expiry_days==day)&mid.priority_selected];c=new[(new.expiry_days==day)&new.priority_selected]
        rows.append(dict(day=int(day),old=float(np.median(a.mc_se_U/a.price_scale_U)),dynamic=float(np.median(b.mc_se_U/b.price_scale_U)),conditional=float(np.median(c.mc_se_U/c.price_scale_U)),residual=float(np.median(c.residual_SE_U/c.price_scale_U)),cheap=float(np.median(c.cheap_SE_U/c.price_scale_U))))
    frame=pd.DataFrame(rows);frame.to_csv(out/'precision.csv',index=False);fig,ax=plt.subplots(figsize=(9,4.6));x=np.arange(len(days))
    for shift,label,color,key in [(-.25,'Previous estimator','#9a9a9a','old'),(0,'Daily hedge control','#e9a23b','dynamic'),(.25,'One-step conditional control','#267b96','conditional')]:ax.bar(x+shift,frame[key],width=.24,label=label,color=color)
    ax.axhline(.125,color='#ad2929',linestyle='--',label='Required SE ≤ 0.125 price-scale units');ax.set_yscale('log');ax.set_xticks(x,[str(d) for d in days]);ax.set_xlabel('SPX maturity (days)');ax.set_ylabel('Median estimated SE / price scale');ax.set_title('Same frozen VIX model — SPX pricing noise falls, gate still fails');ax.legend(fontsize=8);fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(out/f'precision-comparison.{ext}',dpi=170)
    plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(16,8));iv_frames=[]
    for ax,day in zip(axes.flat,days):
        g=new[(new.expiry_days==day)&new.priority_selected].copy()
        def iv(price,row):return implied_vol(price,row.forward,row.strike,row.T_years,row.option_type)
        for name,column in [('market','price_mid_U'),('model','model_price_U')]:g[f'{name}_IV']=[iv(getattr(r,column),r) for r in g.itertuples()]
        g['lower_IV']=[iv(r.model_price_U-1.96*r.mc_se_U,r) for r in g.itertuples()];g['upper_IV']=[iv(r.model_price_U+1.96*r.mc_se_U,r) for r in g.itertuples()];g['invalid_model_IV']=~np.isfinite(g.model_IV)
        # No smoothing or interpolation: every native quote remains in the CSV.
        ax.scatter(g.K_over_F,g.market_IV,s=7,color='black',label='Market native quotes',alpha=.55)
        ax.scatter(g.K_over_F,g.model_IV,s=7,color='#267b96',label='Frozen VIX model',alpha=.7)
        for r in g.iloc[::max(1,len(g)//24)].itertuples():
            if np.isfinite(r.lower_IV) and np.isfinite(r.upper_IV):ax.plot([r.K_over_F]*2,[r.lower_IV,r.upper_IV],color='#267b96',alpha=.45,lw=.8)
        ax.set_title(f'{day} days');ax.set_xlabel('Strike / forward');ax.set_ylabel('Native implied vol (%)');ax.grid(alpha=.15);iv_frames.append(g)
    axes.flat[-1].axis('off');axes.flat[-1].text(0,1,'No parameter refit\n\nError bars: ±1.96 estimated SE\nSampling error only; excludes timestep\nand normalizer uncertainty.\n\nAll seven SPX maturities retained.\nAll VIX parameters unchanged.',va='top',fontsize=11);axes.flat[0].legend(fontsize=8);fig.suptitle('Protected two-speed VIX model: more precise SPX diagnostics',fontsize=15);fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(out/f'frozen-model-spx-smiles.{ext}',dpi=170)
    plt.close(fig);pd.concat(iv_frames).to_csv(out/'native-IV-diagnostics.csv',index=False)
    # Allocation projection uses independent-batch sqrt(N) scaling only. QMC
    # gains from enlarging a net are not assumed: increase independent scrambles.
    primary=new[new.priority_selected].copy();component_target=.125/np.sqrt(2)
    primary['projected_expensive_paths']=128*np.ceil(512*(primary.residual_SE_U/primary.price_scale_U/component_target)**2/128)
    primary['projected_reference_scrambles']=np.maximum(8,np.ceil(8*(primary.cheap_SE_U/primary.price_scale_U/component_target)**2))
    primary['reference_points_per_scramble']=4096;primary.to_csv(out/'sampling-budget-projection.csv',index=False)
    budget=[]
    for day,g in primary.groupby('expiry_days',sort=True):budget.append(dict(day=int(day),median_expensive_paths=float(np.median(g.projected_expensive_paths)),max_expensive_paths=int(g.projected_expensive_paths.max()),median_reference_scrambles=float(np.median(g.projected_reference_scrambles)),max_reference_scrambles=int(g.projected_reference_scrambles.max())))
    dump(out/'sampling-budget-projection.json',dict(per_maturity=budget,meaning='Pilot-based allocation projection, not a guaranteed budget or achieved precision. Equal variance budget for expensive correction and cheap reference; no bias allowance.',full_model_stored_path_limit=512))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();run(args.output.resolve())
