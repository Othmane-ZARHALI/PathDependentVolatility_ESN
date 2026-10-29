"""Scientific plots of audited pooled state and five-piece checking results."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from .plot_joint_curve_fit import implied_vol


def render(state_run,curve_run,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    out.mkdir(parents=True,exist_ok=True)
    states=json.loads((state_run/'run.json').read_text());curves=json.loads((curve_run/'run.json').read_text())
    labels=['Zero state · 3 pieces','Selected state · 3 pieces','Selected state · 5 pieces']
    colors=['#737373','#c17621','#246cb1'];styles=[':', '--','-']
    datasets=[(state_run,'baseline'),(curve_run,'baseline'),(curve_run,'selected')]
    qs=[pd.read_csv(p/f'prices_fresh_pooled_{name}.csv') for p,name in datasets]
    fs=[pd.read_csv(p/f'forwards_fresh_pooled_{name}.csv') for p,name in datasets]
    for q in qs[1:]:pd.testing.assert_frame_equal(q[['instrument','expiry','strike','option_type']],qs[0][['instrument','expiry','strike','option_type']])
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    outputs=[]
    def save(fig,stem):
        for suffix in ['png','svg']:
            p=out/f'{stem}.{suffix}';fig.savefig(p,dpi=150,bbox_inches='tight');outputs.append(p)
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
    axes[0].plot(fs[0].expiry_days,fs[0].forward_vix_points,'D-',color='0.15',label='Market')
    for f,label,color,style in zip(fs,labels,colors,styles):
        axes[0].errorbar(f.expiry_days,f.model_forward,yerr=2*f.mc_se,ls=style,marker='o',color=color,capsize=3,label=label)
    axes[0].set(xlabel='VIX expiry (days)',ylabel='VIX forward (points)',title='All four forward targets · ±2 sampling SE')
    e3=[.25*365,.75*365,365];e5=np.array(curves['ends_years'])*365
    for levels,ends,label,color,style in zip([states['baseline_levels'],curves['three_levels'],curves['selected_levels']],
            [e3,e3,e5],labels,colors,styles):
        x=np.r_[0,np.array(ends)[np.array(ends)<365],365]
        v=100*np.sqrt(levels);y=[v[min(i,len(v)-1)] for i in range(len(x))]
        axes[1].step(x,y,where='post',color=color,ls=style,lw=1.8,label=label)
    axes[1].set(xlabel='Time (calendar days)',ylabel='Square root of annual variance (%)',title='Variance input · first year')
    for ax in axes:ax.grid(alpha=.2)
    handles,legend=axes[0].get_legend_handles_labels();fig.legend(handles,legend,loc='outside lower center',ncols=2,frameon=False)
    fig.suptitle('Joint-fit progression · pooled checking ensembles\nForward and option improvements remain conditional numerical evidence, not a qualified market fit.',fontsize=13)
    save(fig,'forwards-and-curve-progression')
    for instrument in ['SPX','VIX']:
        fig,axes=plt.subplots(2,2,figsize=(12,8.5),sharey=True,layout='constrained')
        for ax,days in zip(axes.flat,sorted(qs[0][qs[0].instrument==instrument].expiry_days.unique())):
            groups=[q[(q.instrument==instrument)&(q.expiry_days==days)].sort_values('strike') for q in qs]
            market=groups[0]
            for side in ['put','call']:
                m=market[market.option_type==side]
                def invs(frame,column):
                    return np.array([implied_vol(float(getattr(r,column)),float(r.forward),float(r.strike),float(r.T_years),r.option_type) for r in frame.itertuples()],dtype=float)
                x=m.K_over_F.to_numpy();bid=invs(m,'price_bid_U');ask=invs(m,'price_ask_U');mid=invs(m,'price_mid_U')
                ax.fill_between(x,bid,ask,color='0.2',alpha=.16,label='Market bid–ask' if side=='put' else None)
                ax.plot(x,mid,color='0.15',lw=1)
                for k,(group,label,color,style) in enumerate(zip(groups,labels,colors,styles)):
                    g=group[group.option_type==side];iv=invs(g,'model_price_U')
                    ax.plot(g.K_over_F,iv,color=color,ls=style,lw=1.5,label=label if side=='put' else None)
                    if k==2:
                        low=[];high=[]
                        for r in g.itertuples():
                            inv=lambda v:implied_vol(v,float(r.forward),float(r.strike),float(r.T_years),r.option_type)
                            low.append(inv(r.model_price_U-2*r.mc_se_U));high.append(inv(r.model_price_U+2*r.mc_se_U))
                        ax.fill_between(g.K_over_F,np.array(low,dtype=float),np.array(high,dtype=float),color=color,alpha=.12)
            ax.set_xscale('log');x=market.K_over_F.to_numpy();ticks=np.linspace(x.min(),x.max(),4)
            ax.set_xticks(ticks,[f'{v:.2f}' for v in ticks]);ax.set_xticks([],minor=True)
            ax.set_xlabel('Strike / market forward (log scale)');ax.set_ylabel('Annualized Black IV (%)')
            ax.set_title(f'{instrument} · {days} days · {len(market)} options',loc='left')
            ax.axvline(1,color='0.6',lw=.7);ax.grid(alpha=.2);ax.yaxis.set_major_locator(MaxNLocator(5))
        handles,legend=axes[0,0].get_legend_handles_labels();fig.legend(handles,legend,loc='outside lower center',ncols=2,frameon=False)
        fig.suptitle(f'{instrument}: successive joint candidates on pooled checking draws\nAll IVs use fixed market forwards; selected shading is ±2 sampling SE only, excluding numerical bias.',fontsize=13)
        save(fig,instrument.lower()+'-fit-progression')
    sources=[Path(__file__),state_run/'run.json',curve_run/'run.json',
             *[p/f'{kind}_fresh_pooled_{name}.csv' for p,name in datasets for kind in ['prices','forwards']]]
    manifest=dict(source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs},
        caveat='The state and refinement comparisons reuse the same checking primitives. Sampling bands exclude all numerical bias and pilot variability; IV gaps are retained.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--state-run',type=Path,required=True)
    parser.add_argument('--curve-run',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();render(args.state_run,args.curve_run,args.output)
