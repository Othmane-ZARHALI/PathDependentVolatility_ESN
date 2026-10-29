"""Native SPX/VIX IV panels for the joint two-speed experiment."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .plot_joint_curve_fit import implied_vol


def render(run,out):
    out.mkdir(parents=True,exist_ok=True);plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    specs=[('vix_only','#236eb0','--','VIX-only parameters'),('curve_floor','#939393',':','Curve/floor joint fit'),('joint','#d04a35','-','All23 joint fit')]
    def save(fig,name):
        for ext in ['png','svg']:fig.savefig(out/f'{name}.{ext}',dpi=155,bbox_inches='tight')
        plt.close(fig)
    for instrument,days in [('spx',[29,57,119,211,302,666,1030]),('vix',[27,55,118,153])]:
        fig,axes=plt.subplots(3 if instrument=='spx' else 2,3 if instrument=='spx' else 2,figsize=(15,12) if instrument=='spx' else (13,9),layout='constrained')
        for ax,day in zip(axes.flat,days):
            q=pd.read_csv(run/'checking'/f'prices_{instrument}_vix_only.csv').query('expiry_days == @day')
            def iv(col):return np.array([implied_vol(getattr(a,col),a.forward,a.strike,a.T_years,a.option_type) for a in q.itertuples()],float)
            ax.fill_between(q.K_over_F,iv('price_bid_U'),iv('price_ask_U'),color='.3',alpha=.17,label='Market bid–ask');ax.plot(q.K_over_F,q.market_native_IV,'k-',lw=1.3,label='Market mid')
            for name,color,ls,label in specs:
                g=pd.read_csv(run/'checking'/f'prices_{instrument}_{name}.csv').query('expiry_days == @day');ax.plot(g.K_over_F,g.model_native_IV,color=color,ls=ls,label=label,lw=1.7);bad=~np.isfinite(g.model_native_IV)
                if bad.any():ax.plot(g.K_over_F[bad],np.full(bad.sum(),.02),'|',color=color,transform=ax.get_xaxis_transform())
            ax.set(xlabel='Strike / market forward',ylabel='Native Black IV (%)',title=f'{instrument.upper()} {day} days');ax.grid(alpha=.18)
            if instrument=='vix':ax.set_xscale('log')
        for ax in list(axes.flat)[len(days):]:ax.set_visible(False)
        handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=3,frameon=False);fig.suptitle('Two-speed model: joint SPX/VIX trade-off\n64 steps/trading day · checking partly reuses streams · sampling and timestep gates assessed separately',fontsize=12);save(fig,f'{instrument}-joint-smiles')
    fig,ax=plt.subplots(figsize=(10,4.5),layout='constrained');base=json.loads((run/'design.json').read_text())['baseline']
    from .two_speed_joint_core import physical
    ends=np.array(json.loads((run/'experimental_model_joint.json').read_text())['ends_years'])*365;records=json.loads((run/'optimizer_ledger.json').read_text());curves=[('vix_only',np.array(base)),('curve_floor',np.array(records[0]['parameters'])),('joint',np.array(records[1]['parameters']))]
    for (name,x),(_,color,ls,label) in zip(curves,specs):
        levels,floor=physical(x);ax.stairs(100*np.sqrt(levels),np.r_[0,ends],color=color,ls=ls,label=label)
    ax.set(xlabel='Days from pricing date',ylabel='Forward variance level, shown as volatility (%)',title='All16 curve levels are free in the joint fit');ax.legend(frameon=False);ax.grid(alpha=.2);save(fig,'joint-variance-curves')
    fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    for name,color,ls,label in specs:
        f=pd.read_csv(run/'checking'/f'forwards_{name}.csv');ax.errorbar(f.expiry_days,f.forward_error,yerr=2*f.mc_se,fmt='o',color=color,capsize=4,label=label)
    ax.axhline(0,color='k',lw=1);ax.set(xlabel='VIX expiry days',ylabel='Model minus target forward',title='Checking forwards: ±2 block SE');ax.legend(frameon=False);ax.grid(alpha=.2);save(fig,'joint-forward-errors')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();render(a.run,a.output)
