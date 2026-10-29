"""Continuous-measure fit against current and refitted two-atom references."""
import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .plot_maturity_curve import iv


def render(run,out):
    out.mkdir(parents=True,exist_ok=True);plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    specs=[('previous','#888888','--','Previous fitted floor'),('boundary','#256daf','-','Two-atom refit'),('continuous','#d04a35','-','Continuous Beta measure')];table={n:pd.read_csv(run/'fresh_check'/f'prices_canonical_{n}.csv') for n,*_ in specs};native={n:pd.read_csv(run/'fresh_check'/f'native_iv_{n}.csv') for n,*_ in specs}
    def save(fig,name):
        for ext in ['png','svg']:fig.savefig(out/f'{name}.{ext}',dpi=155,bbox_inches='tight')
        plt.close(fig)
    for mode in ['canonical','native','wing']:
        fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
        for ax,day in zip(axes.flat,[27,55,118,153]):
            q=table['previous'].query('expiry_days == @day')
            if mode=='wing':q=q[q.K_over_F>=2]
            if mode!='native':ax.fill_between(q.K_over_F,iv(q,q.price_bid_U),iv(q,q.price_ask_U),color='.25',alpha=.2,label='Market bid–ask')
            ax.plot(q.K_over_F,iv(q,q.price_mid_U),color='.1',lw=1.5,label='Market mid');missing=[]
            for n,color,style,label in specs:
                t=(native[n] if mode=='native' else table[n]).query('expiry_days == @day')
                if mode=='wing':t=t[t.K_over_F>=2]
                y=t.model_native_IV.to_numpy() if mode=='native' else iv(t,t.model_price_U);missing.append(int(np.sum(~np.isfinite(y))));ax.plot(t.K_over_F,y,color=color,ls=style,lw=1.7,label=label)
                if mode!='native' and n=='continuous':ax.fill_between(t.K_over_F,iv(t,t.model_price_U-2*t.mc_se_U),iv(t,t.model_price_U+2*t.mc_se_U),color=color,alpha=.12,label='Continuous fit ±2 price sampling SE')
                bad=~np.isfinite(y)
                if bad.any():ax.plot(t.K_over_F.to_numpy()[bad],np.full(bad.sum(),.02),'|',transform=ax.get_xaxis_transform(),color=color)
            ax.set(title=f'VIX {day} days',xlabel='Strike / observed forward',ylabel='Black implied volatility (%)');ax.set_xscale('log');ticks=np.geomspace(q.K_over_F.min(),q.K_over_F.max(),5);ax.set_xticks(ticks,[f'{x:.2f}' for x in ticks]);ax.set_xticks([],minor=True);ax.grid(alpha=.18);ax.text(.02,.02,'Unavailable IVs by curve: '+', '.join(map(str,missing)),transform=ax.transAxes,fontsize=8)
        handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=2,frameon=False)
        convention='Native OTM prices; each model uses its own forward' if mode=='native' else 'Canonical calls; market forward used for every IV inversion'
        fig.suptitle('Continuous M2/stress response spectrum\n'+convention+' · fresh16384points/expiry; SPX paused',fontsize=12);save(fig,{'canonical':'vix-all-strikes','native':'vix-native-IV','wing':'vix-right-wing'}[mode])
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,day in zip(axes.flat,[27,55,118,153]):
        q=native['previous'].query('expiry_days == @day and option_type == "put"');ax.fill_between(q.K_over_F,q.price_bid_U,q.price_ask_U,color='.25',alpha=.2,label='Market bid–ask');ax.plot(q.K_over_F,q.price_mid_U,'k.-',label='Market mid')
        for n,color,style,label in specs:
            t=native[n].query('expiry_days == @day and option_type == "put"');ax.plot(t.K_over_F,t.model_price_U,color=color,ls=style,lw=1.7,label=label)
        ax.set(title=f'VIX {day} days · original OTM puts',xlabel='Strike / observed forward',ylabel='Undiscounted native put price');ax.grid(alpha=.18)
    handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=2,frameon=False);fig.suptitle('Actual low-strike put prices, with no call-parity conversion',fontsize=13);save(fig,'native-put-prices')
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained');levels=pd.read_csv(run/'curve_levels.csv')
    for n,color,style,label in specs:
        f=pd.read_csv(run/'fresh_check'/f'forwards_{n}.csv');axes[0].errorbar(f.expiry_days,f.forward_error,yerr=2*f.mc_se,color=color,ls=style,marker='o',capsize=2,label=label)
        s=json.loads((run/'fresh_check'/f'score_{n}.json').read_text());axes[1].stairs(100*np.sqrt(levels[n]),np.r_[levels.start_day.iloc[0],levels.end_day],baseline=None,color=color,ls=style,label=label,lw=1.5);axes[1].axhline(100*s['model']['floor_volatility'],color=color,ls=':',lw=1)
    axes[0].axhline(0,color='.1',lw=1);axes[0].set(title='Fresh forward errors ±2 sampling SE',xlabel='Expiry (days)',ylabel='Model minus observed forward (VIX points)');axes[1].set(title='16-piece curve, first211days; dotted floors',xlim=(0,211),xlabel='Calendar days from valuation',ylabel='Square root of forward variance (%)')
    for ax in axes:ax.grid(alpha=.18);ax.legend(fontsize=8,frameon=False)
    fig.suptitle('Forward equality constraints apply to training estimates; fresh residuals remain visible',fontsize=12);save(fig,'forwards-and-curve')


    from scipy.special import betainc
    selection=json.loads((run/'selection.json').read_text());m=selection['model']['measure_mean'];k=selection['model']['kappa'];u=np.linspace(0,1,501)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained');axes[0].plot(u,betainc(m*k,(1-m)*k,u),color='#d04a35',lw=2)
    axes[0].set(xlabel='Component position u: M2 (0) to stress (1)',ylabel='Cumulative expected non-floor variance share',title=f'Fitted Beta measure: m={m:.4g}, kappa={k:.4g}')
    bins=np.linspace(0,1,11);mass=np.diff(betainc(m*k,(1-m)*k,bins));axes[1].bar((bins[:-1]+bins[1:])/2,mass,width=.09,color='#d04a35');axes[1].set(xlabel='Component position u',ylabel='Measure mass in each interval',title='Contributions spread across the continuous family')
    for ax in axes:ax.grid(alpha=.18)
    fig.suptitle('Measure weights describe expected variance contributions; they are not regime probabilities',fontsize=11);save(fig,'fitted-measure')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();render(a.run,a.output)
