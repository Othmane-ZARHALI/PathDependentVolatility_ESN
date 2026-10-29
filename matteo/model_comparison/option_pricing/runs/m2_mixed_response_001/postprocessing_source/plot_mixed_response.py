"""VIX mixture fits, full strikes, tails and the full sixteen-piece curve."""
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
    out.mkdir(parents=True,exist_ok=True);plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    specs=[('baseline','#9b7b47','--','Old curve, original M2'),('previous_single_factor','#7f7f7f',':','Previous fitted single factor'),('product_price','#256daf','-','Two sensitivities × M2, price fit'),('additive_price','#d04a35','-','M2 + stress component, price fit')]
    shapes=[('previous_single_factor','#7f7f7f',':','Previous fitted single factor'),('product_shape','#256daf','-','Two sensitivities × M2, shape selection'),('additive_shape','#d04a35','-','M2 + stress component, shape selection')]
    allnames=set(n for n,*_ in specs+shapes);tables={n:pd.read_csv(run/'fresh_check'/f'prices_canonical_{n}.csv') for n in allnames}
    def save(fig,name):
        for ext in ['png','svg']:fig.savefig(out/f'{name}.{ext}',dpi=155,bbox_inches='tight')
        plt.close(fig)
    for selection,wing,filename in [(specs,False,'vix-all-strikes'),(specs,True,'vix-right-wing'),(shapes,False,'vix-shape-selections')]:
        fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
        for ax,day in zip(axes.flat,[27,55,118,153]):
            q=tables[selection[0][0]].query('expiry_days == @day')
            if wing:q=q[q.K_over_F>=2]
            ax.fill_between(q.K_over_F,iv(q,q.price_bid_U),iv(q,q.price_ask_U),alpha=.2,color='.25',label='Market bid–ask');ax.plot(q.K_over_F,iv(q,q.price_mid_U),color='.1',label='Market mid',lw=1.4)
            missing=[]
            for name,color,style,label in selection:
                t=tables[name].query('expiry_days == @day')
                if wing:t=t[t.K_over_F>=2]
                y=iv(t,t.model_price_U);missing.append(int(np.sum(~np.isfinite(y))));ax.plot(t.K_over_F,y,color=color,ls=style,label=label,lw=1.7)
                if name.startswith('additive'):ax.fill_between(t.K_over_F,iv(t,t.model_price_U-2*t.mc_se_U),iv(t,t.model_price_U+2*t.mc_se_U),color=color,alpha=.12,label='Additive fit ±2 sampling SE')
                bad=~np.isfinite(y)
                if bad.any():ax.plot(t.K_over_F.to_numpy()[bad],np.full(bad.sum(),.02),'|',transform=ax.get_xaxis_transform(),color=color)
            ax.set(title=f'VIX {day} days · {len(q)} strikes shown',xlabel='Strike / observed forward',ylabel='Black implied volatility (%)');ax.set_xscale('log');ticks=np.geomspace(q.K_over_F.min(),q.K_over_F.max(),5);ax.set_xticks(ticks,[f'{x:.2f}' for x in ticks]);ax.set_xticks([],minor=True);ax.grid(alpha=.18);ax.text(.02,.02,'Unavailable IVs by curve: '+', '.join(map(str,missing)),transform=ax.transAxes,fontsize=8)
        handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=2,frameon=False)
        fig.suptitle('Two variance responses · forward variance fitted in both models\nFresh 16384-point simulations per expiry; all original VIX strikes retained; SPX paused',fontsize=13);save(fig,filename)
    levels=pd.read_csv(run/'curve_levels.csv');fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for ax,limit in zip(axes,[211,1030]):
        for name,color,style,label in specs:ax.stairs(100*np.sqrt(levels[name]),np.r_[levels.start_day.iloc[0],levels.end_day],baseline=None,color=color,ls=style,label=label,lw=1.6)
        for a,b in [(27,57),(55,85),(118,148),(153,183)]:ax.axvspan(a,b,color='.5',alpha=.1)
        ax.set(xlim=(0,limit),xlabel='Calendar days from valuation',ylabel='Square root of forward variance (%)');ax.grid(alpha=.18)
    axes[0].legend(fontsize=8,frameon=False);fig.suptitle('All 16 original intervals retained · seven VIX-sensitive levels fitted independently\nShading shows the four 30-day VIX averaging periods',fontsize=12);save(fig,'variance-curve')
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for name,color,style,label in specs:
        f=pd.read_csv(run/'fresh_check'/f'forwards_{name}.csv');axes[0].errorbar(f.expiry_days,f.model_forward,yerr=2*f.mc_se,color=color,ls=style,marker='o',capsize=2,label=label)
        score=json.loads((run/'fresh_check'/f'score_{name}.json').read_text());axes[1].plot([r['day'] for r in score['moments']],[r['exact_second_moment'] for r in score['moments']],color=color,ls=style,marker='o',label=label)
    axes[0].plot(f.expiry_days,f.forward_vix_points,'kD-',label='Observed market forward');axes[1].plot([27,55,118,153],[468.99932893,556.67298547,648.62569161,695.29544408],'kD--',label='Necessary bound: all quotes + exact forward')
    axes[0].set(title='VIX forwards',xlabel='Expiry (days)',ylabel='VIX points');axes[1].set(title='Exact curve-implied second moments',xlabel='Expiry (days)',ylabel='E[VIX²] (VIX points squared)')
    for ax in axes:ax.grid(alpha=.18);ax.legend(fontsize=8,frameon=False)
    fig.suptitle('Mean variance fitted together with each mixture; same calibration market data',fontsize=12);save(fig,'forwards-and-moments')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();render(a.run,a.output)
