"""Four-maturity comparisons without filtering native prices or missing IVs."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .plot_maturity_curve import iv


def render(run,out):
    out.mkdir(parents=True,exist_ok=True);plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    specs=[('two_speed','#236eb0','--','Previous fitted two-speed model'),('bank','#d04a35','-','Fixed geometric four-speed bank')]
    for kind in ['smiles','right-wing','prices']:
        fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
        for ax,day in zip(axes.flat,[27,55,118,153]):
            market=pd.read_csv(run/'fresh_check/prices_canonical_two_speed.csv').query('expiry_days == @day')
            if kind=='right-wing':market=market[market.K_over_F>2]
            x=market.K_over_F
            def yy(q,col):return q[col] if kind=='prices' else iv(q,q[col])
            ax.fill_between(x,yy(market,'price_bid_U'),yy(market,'price_ask_U'),color='.3',alpha=.18,label='Market bid–ask')
            ax.plot(x,yy(market,'price_mid_U'),'k-',lw=1.4,label='Market mid')
            for name,color,ls,label in specs:
                q=pd.read_csv(run/'fresh_check'/f'prices_canonical_{name}.csv').query('expiry_days == @day')
                if kind=='right-wing':q=q[q.K_over_F>2]
                y=yy(q,'model_price_U');ax.plot(q.K_over_F,y,color=color,ls=ls,lw=1.8,label=label)
                if np.any(~np.isfinite(y)):ax.plot(q.K_over_F.to_numpy()[~np.isfinite(y)],np.full(np.sum(~np.isfinite(y)),.02),'|',color=color,transform=ax.get_xaxis_transform())
            ax.set_xscale('log');ax.set(title=f'VIX {day} days',xlabel='Strike / market forward',ylabel='Undiscounted call price' if kind=='prices' else 'Black implied volatility (%)');ax.grid(alpha=.18)
            if kind=='prices':ax.set_yscale('symlog',linthresh=.005)
        handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=4,frameon=False)
        fig.suptitle('One shared model across all four maturities\nFixed half-lives 7, 28, 112, 448 days · checking partly reuses states · forward equalities apply to training',fontsize=12)
        for ext in ['png','svg']:fig.savefig(out/f'vix-stress-bank-{kind}.{ext}',dpi=155,bbox_inches='tight')
        plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    for name,color,ls,label in specs:
        q=pd.read_csv(run/'fresh_check'/f'forwards_{name}.csv');ax.errorbar(q.expiry_days,q.forward_error,yerr=2*q.mc_se,fmt='o',color=color,label=label,capsize=4)
    ax.axhline(0,color='k',lw=1);ax.set(xlabel='Expiry days',ylabel='Model − target VIX forward',title='Checking forward errors with ±2 block SE');ax.legend();ax.grid(alpha=.2)
    for ext in ['png','svg']:fig.savefig(out/f'vix-stress-bank-forwards.{ext}',dpi=155,bbox_inches='tight')
    plt.close(fig)

    model=json.loads((run/'selection.json').read_text())['model'];h=np.array(model['half_lives_days']);t=np.linspace(0,500,501)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for half in h:axes[0].plot(t,2**(-t/half),label=f'{half:g} days')
    axes[0].axvline(30,color='.5',ls=':',label='30-day VIX window');axes[0].set(xlabel='Days since a shock',ylabel='Fraction of state shock remaining',title='Fixed geometric persistence grid');axes[0].legend(frameon=False);axes[0].grid(alpha=.2)
    labels=['M2']+[f'{v:g}d' for v in h];ww=[model['M2_weight']]+model['stress_weights']
    axes[1].bar(labels,ww,color=['.5','#b8d9ec','#77b6d9','#367fa8','#164c70']);axes[1].set(ylabel='Weight in normalized variance response',title='Jointly fitted mixture weights');axes[1].grid(axis='y',alpha=.2)
    for j,v in enumerate(ww):axes[1].text(j,v+.01,f'{100*v:.1f}%',ha='center',fontsize=9)
    axes[1].set_ylim(0,max(ww)+.12)
    for ext in ['png','svg']:fig.savefig(out/f'stress-bank-persistence-weights.{ext}',dpi=155,bbox_inches='tight')
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();render(a.run,a.output)
