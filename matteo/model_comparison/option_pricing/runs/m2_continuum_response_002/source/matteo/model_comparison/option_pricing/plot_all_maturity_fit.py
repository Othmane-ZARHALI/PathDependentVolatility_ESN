"""All seven SPX and four VIX fit panels, on the fixed market-forward IV scale."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from .plot_joint_curve_fit import implied_vol


def render(run,out,role='check'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    out.mkdir(parents=True,exist_ok=True)
    record_paths=[run/'run.json',run/'score_audit.json'] if role=='check' else [run/'design.json',run/'selected.json']
    if role=='check':r=json.loads(record_paths[0].read_text())
    else:r=json.loads(record_paths[0].read_text())|{'selected_levels':json.loads(record_paths[1].read_text())['levels']}
    p=[pd.read_csv(run/f'prices_{role}_{name}.csv') for name in ['baseline','selected']]
    f=[pd.read_csv(run/f'forwards_{role}_{name}.csv') for name in ['baseline','selected']]
    note=('Same market quotes used for calibration; separate simulation draws check numerical stability only.' if role=='check'
          else 'SEARCH-SAMPLE PREVIEW: parameters selected on these simulations; numerical checking pending.')
    labels=['Previous working parameters','Refitted on all retained maturities'];colors=['#bc7628','#2269ad'];styles=['--','-']
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    outputs=[]
    def save(fig,name):
        for ext in ['png','svg']:
            path=out/f'{name}.{ext}';fig.savefig(path,dpi=145,bbox_inches='tight');outputs.append(path)
        plt.close(fig)
    def ivs(frame,column):
        return np.array([implied_vol(float(getattr(t,column)),float(t.forward),float(t.strike),float(t.T_years),t.option_type) for t in frame.itertuples()],dtype=float)
    for instrument,shape in [('SPX',(4,2)),('VIX',(2,2))]:
        fig,axes=plt.subplots(*shape,figsize=(13,4.05*shape[0]),layout='constrained')
        dates=sorted(p[0][p[0].instrument==instrument].expiry_days.unique())
        for ax,days in zip(axes.flat,dates):
            groups=[t[(t.instrument==instrument)&(t.expiry_days==days)].sort_values('strike') for t in p]
            unavailable=0
            for side in ['put','call']:
                market=groups[0][groups[0].option_type==side]
                ax.fill_between(market.K_over_F,ivs(market,'price_bid_U'),ivs(market,'price_ask_U'),color='.2',alpha=.18,label='Market bid–ask' if side=='put' else None)
                ax.plot(market.K_over_F,ivs(market,'price_mid_U'),color='.15',lw=1.2)
                for k,(g,label,color,style) in enumerate(zip(groups,labels,colors,styles)):
                    t=g[g.option_type==side];v=ivs(t,'model_price_U');unavailable+=int(np.isnan(v).sum()) if k==1 else 0
                    ax.plot(t.K_over_F,v,lw=1.8,color=color,ls=style,label=label if side=='put' else None)
                    if k==1:
                        lower=t.copy();upper=t.copy()
                        lower['model_price_U']=t.model_price_U-2*t.mc_se_U;upper['model_price_U']=t.model_price_U+2*t.mc_se_U
                        ax.fill_between(t.K_over_F,ivs(lower,'model_price_U'),ivs(upper,'model_price_U'),color=color,alpha=.13,label='Refit ±2 sampling SE' if side=='put' else None)
            ax.set_xscale('log');xs=groups[0].K_over_F.to_numpy();ticks=np.geomspace(xs.min(),xs.max(),5)
            ax.set_xticks(ticks,[f'{t:.2f}' for t in ticks]);ax.set_xticks([],minor=True)
            ax.set(xlabel='Strike / market forward (log scale)',ylabel='Black implied volatility (%)',
                title=f'{instrument} {days} days · {len(groups[0])} options · IN-SAMPLE FIT')
            ax.axvline(1,color='.6',lw=.8);ax.grid(alpha=.18);ax.yaxis.set_major_locator(MaxNLocator(5))
            if unavailable:ax.text(.02,.03,f'{unavailable} model IVs unavailable; prices retained',transform=ax.transAxes,fontsize=9)
        if len(dates)<axes.size:
            ax=axes.flat[-1];ax.axis('off')
            ax.text(.04,.94,'All seven SPX maturities enter this fit.\n\nCutoff: |log(K/F)| / (ATM proxy × √T) ≤ 3\n\nExcluded only: 1, 1401, 1765, 2129 days.\n\nOriginal split labels are historical metadata.\nThere is no market out-of-sample panel.\n\nAll IVs use the observed market forwards.\nSE shading excludes numerical bias.',va='top',fontsize=11)
        handles,text=axes.flat[0].get_legend_handles_labels();fig.legend(handles,text,loc='outside lower center',ncols=2,frameon=False)
        fig.suptitle(f'{instrument}: fit on every retained maturity\n'+note,fontsize=14)
        save(fig,instrument.lower()+'-all-maturity-fit')
    fig,axes=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
    axes[0].plot(f[0].expiry_days,f[0].forward_vix_points,'D-',color='.1',label='Observed market forwards')
    for t,label,color,style in zip(f,labels,colors,styles):
        axes[0].errorbar(t.expiry_days,t.model_forward,yerr=2*t.mc_se,marker='o',color=color,ls=style,capsize=3,label=label)
    axes[0].set(xlabel='VIX expiry (days)',ylabel='VIX forward (points)',title='All four VIX forwards are fitting targets')
    for xi,label,color,style in zip([r['baseline_levels'],r['selected_levels']],labels,colors,styles):
        end=np.array(r['ends_years'])*365;x=np.r_[0,end[end<1030],1030];vol=100*np.sqrt(xi)
        axes[1].step(x,np.r_[vol,vol[-1]],where='post',color=color,ls=style,lw=1.8,label=label)
    axes[1].set(xlabel='Future calendar day',ylabel='Square root of annual variance (%)',title='Variance input across the full SPX fitting horizon')
    for ax in axes:ax.grid(alpha=.2)
    handles,text=axes[0].get_legend_handles_labels();fig.legend(handles,text,loc='outside lower center',ncols=1,frameon=False)
    fig.suptitle('In-sample joint SPX/VIX fitting\n'+note,fontsize=12)
    save(fig,'forwards-and-variance')
    sources=[Path(__file__),*record_paths,
        *[run/f'{kind}_{role}_{name}.csv' for kind in ['prices','forwards'] for name in ['baseline','selected']]]
    if role=='check':
        summary=pd.read_csv(run/'expiry_scores.csv');summary.to_csv(out/'expiry-fit-scores.csv',index=False)
        sources.append(run/'expiry_scores.csv')
    manifest=dict(source_sha256={str(t):hashlib.sha256(t.read_bytes()).hexdigest() for t in sources},
        output_sha256={t.name:hashlib.sha256(t.read_bytes()).hexdigest() for t in outputs},
        role=role,market_out_of_sample=False,iv_forward='Observed market forwards, fixed for all model and market IVs.',
        uncertainty='Sampling only; excludes timestep, normalizer/GH/time bias and pilot uncertainty.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--role',choices=['search','check'],default='check')
    a=p.parse_args();render(a.run,a.output,a.role)
