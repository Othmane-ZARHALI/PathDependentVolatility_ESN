"""Include the training-selected genuinely spread continuous measure."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import betainc
from .plot_maturity_curve import iv


def render(run,out):
    out.mkdir(parents=True,exist_ok=True);plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    specs=[('previous','#999999','--','Previous fitted-floor model'),('boundary','#256daf',':','Two-atom refit'),('continuous','#d04a35','-','Continuous: training winner'),('interior','#168265','-','Continuous: broader measure')]
    def save(fig,name):
        for ext in ['png','svg']:fig.savefig(out/f'{name}.{ext}',dpi=155,bbox_inches='tight')
        plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for ax,day in zip(axes.flat,[27,55,118,153]):
        market=pd.read_csv(run/'fresh_check/prices_canonical_continuous.csv').query('expiry_days == @day');ax.fill_between(market.K_over_F,iv(market,market.price_bid_U),iv(market,market.price_ask_U),color='.25',alpha=.18,label='Market bid–ask');ax.plot(market.K_over_F,iv(market,market.price_mid_U),'k-',lw=1.5,label='Market mid')
        for name,color,style,label in specs:
            q=pd.read_csv(run/'fresh_check'/f'prices_canonical_{name}.csv').query('expiry_days == @day');y=iv(q,q.model_price_U);ax.plot(q.K_over_F,y,color=color,ls=style,lw=1.7,label=label)
            if np.any(~np.isfinite(y)):ax.plot(q.K_over_F.to_numpy()[~np.isfinite(y)],np.full(np.sum(~np.isfinite(y)),.02),'|',color=color,transform=ax.get_xaxis_transform())
        ax.set_xscale('log');ticks=np.geomspace(market.K_over_F.min(),market.K_over_F.max(),5);ax.set_xticks(ticks,[f'{t:.2f}' for t in ticks]);ax.set_xticks([],minor=True);ax.set(title=f'VIX {day} days',xlabel='Strike / observed forward',ylabel='Black implied volatility (%)');ax.grid(alpha=.18)
    h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle('Continuous response spectrum: endpoint-like and broader measures\nCanonical calls at market forward · reference65536 + difference16384/expiry · unavailable IVs marked at foot',fontsize=12);save(fig,'vix-continuous-comparison')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained');bins=np.linspace(0,1,11);u=np.linspace(0,1,1001)
    for name,color in [('continuous','#d04a35'),('interior','#168265')]:
        z=json.loads((run/'fresh_check'/f'score_{name}.json').read_text())['model'];m=z['measure_mean'];k=z['kappa'];cdf=betainc(m*k,(1-m)*k,u);axes[0].plot(u,cdf,color=color,lw=2,label=f'{name}: m={m:.3f}, kappa={k:.3g}');mass=np.diff(betainc(m*k,(1-m)*k,bins));axes[1].plot((bins[:-1]+bins[1:])/2,mass,'o-',color=color,label=name)
    axes[0].set(xlabel='Component position u: M2 (0) to stress (1)',ylabel='Cumulative measure mass',title='Probability measure over variance components');axes[1].set(xlabel='Component position u',ylabel='Mass in each width0.1 interval',title='Expected non-floor variance contributions')
    for ax in axes:ax.grid(alpha=.18);ax.legend(frameon=False,fontsize=9)
    fig.suptitle('Both measures are continuous; smaller kappa concentrates mass near the endpoints',fontsize=12);save(fig,'continuous-measures-comparison')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();render(a.run,a.output)
