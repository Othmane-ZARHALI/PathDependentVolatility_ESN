"""Full-strike VIX smiles and all sixteen independently represented curve levels."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .plot_maturity_curve import iv


def render(run,out):
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    specs=[('baseline','#a67b3a','--','Previous curve, no new factor'),('curve_only','#2679ad','-','Free curve only'),('fixed_factor_curve','#af6db0',':','Free curve + previous strong factor'),('factor_curve','#d14c36','-','Free curve + fitted factor')]
    tables={name:pd.read_csv(run/'fresh_check'/f'prices_canonical_{name}.csv') for name,*_ in specs}
    def save(fig,name):
        for ext in ['png','svg']:fig.savefig(out/f'{name}.{ext}',dpi=155,bbox_inches='tight')
        plt.close(fig)
    for wing in [False,True]:
        fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
        for ax,day in zip(axes.flat,[27,55,118,153]):
            q=tables['baseline'].query('expiry_days == @day')
            if wing:q=q[q.K_over_F>=2]
            ax.fill_between(q.K_over_F,iv(q,q.price_bid_U),iv(q,q.price_ask_U),alpha=.2,color='.25',label='Market bid–ask')
            ax.plot(q.K_over_F,iv(q,q.price_mid_U),color='.15',label='Market mid',lw=1.4)
            missing=[]
            for name,color,style,label in specs:
                t=tables[name].query('expiry_days == @day')
                if wing:t=t[t.K_over_F>=2]
                y=iv(t,t.model_price_U);missing.append(int(np.sum(~np.isfinite(y))))
                ax.plot(t.K_over_F,y,color=color,ls=style,label=label,lw=1.6)
                if name=='factor_curve':ax.fill_between(t.K_over_F,iv(t,t.model_price_U-2*t.mc_se_U),iv(t,t.model_price_U+2*t.mc_se_U),color=color,alpha=.12,label='Fitted factor ±2 sampling SE')
            ax.set(xlabel='Strike / observed forward',ylabel='Black implied volatility (%)',title=f'VIX expiry {day} days · {len(q)} strikes shown')
            ax.set_xscale('log');ticks=np.geomspace(q.K_over_F.min(),q.K_over_F.max(),5);ax.set_xticks(ticks,[f'{x:.2f}' for x in ticks]);ax.set_xticks([],minor=True);ax.grid(alpha=.18)
            ax.text(.02,.02,'Unavailable IVs by curve: '+', '.join(map(str,missing)),transform=ax.transAxes,fontsize=8)
        handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=3,frameon=False)
        fig.suptitle(('Right wing' if wing else 'All retained strikes')+' · VIX-only fit with the original 16-piece curve\nFresh 8192-point simulations per expiry; same calibration market data; no SPX fit',fontsize=13)
        save(fig,'vix-right-wing' if wing else 'vix-all-strikes')
    table=pd.read_csv(run/'curve_levels.csv');fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for ax,limit in zip(axes,[211,1030]):
        for name,color,style,label in specs:
            ax.stairs(100*np.sqrt(table[name]),np.r_[table.start_day.iloc[0],table.end_day],color=color,ls=style,label=label,lw=1.5)
        for a,b in [(27,57),(55,85),(118,148),(153,183)]:ax.axvspan(a,b,color='.5',alpha=.1)
        ax.set(xlim=(0,limit),xlabel='Calendar days from valuation',ylabel='Square root of forward variance (%)',title='16 pieces retained · shaded VIX averaging periods');ax.grid(alpha=.18)
    axes[0].legend(fontsize=8,frameon=False)
    fig.suptitle('Seven pieces affect these VIX options; nine retain their old values\nA narrow piece hits the 0.20 variance bound (44.7% volatility): exploratory fit, not an adopted curve',fontsize=12)
    save(fig,'variance-curve')
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for name,color,style,label in specs:
        f=pd.read_csv(run/'fresh_check'/f'forwards_{name}.csv');axes[0].errorbar(f.expiry_days,f.model_forward,yerr=2*f.mc_se,color=color,ls=style,marker='o',capsize=2,label=label)
        s=json.loads((run/'fresh_check'/f'score_{name}.json').read_text());axes[1].plot([r['day'] for r in s['moments']],[r['exact_second_moment'] for r in s['moments']],color=color,ls=style,marker='o',label=label)
    axes[0].plot(f.expiry_days,f.forward_vix_points,'kD-',label='Observed market forward')
    axes[1].plot([27,55,118,153],[468.99932893,556.67298547,648.62569161,695.29544408],'kD--',label='Necessary bound for all quotes + exact forward')
    axes[0].set(title='VIX forwards',xlabel='Expiry (days)',ylabel='VIX points');axes[1].set(title='Exact curve-implied second moments',xlabel='Expiry (days)',ylabel='E[VIX²] (VIX points squared)')
    for ax in axes:ax.grid(alpha=.18);ax.legend(fontsize=8,frameon=False)
    fig.suptitle('Average variance remains a fitting parameter; a lower overall loss does not fit every quote',fontsize=12)
    save(fig,'forwards-and-moments')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();render(a.run,a.output)
