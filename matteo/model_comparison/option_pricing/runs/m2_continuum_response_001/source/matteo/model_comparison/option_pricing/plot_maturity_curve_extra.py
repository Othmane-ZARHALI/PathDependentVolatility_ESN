"""Curve-level and separate-stream VIX overlays; true prices only."""
import hashlib
import json
import numpy as np
import pandas as pd
from .plot_parity_parameter_fit import iv


def extra(run,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    r=json.loads((run/'run.json').read_text());outputs=[]
    def save(fig,name):
        for ext in ['png','svg']:
            p=out/f'{name}.{ext}';fig.savefig(p,dpi=145,bbox_inches='tight');outputs.append(p)
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(14,5),layout='constrained');ends=np.array(r['ends_days']);edges=np.r_[0,ends]
    for ax,horizon in zip(axes,[320,1030]):
        for name,color,style,label in [('baseline','#bc7628','--','Previous curve'),('selected','#2269ad','-','New maturity-aligned curve')]:
            p=np.array(r[name]);ax.stairs(100*np.sqrt(p[7:-1]),edges,baseline=None,color=color,ls=style,lw=1.8,label=label)
        ax.axhline(100*np.sqrt(r['selected'][-1]),color='.2',ls=':',label='Selected variance-floor volatility')
        for day in [27,55,118,153]:ax.axvspan(day,day+30,color='#7f9e84',alpha=.09)
        ax.set(xlim=(0,horizon),xlabel='Calendar days',ylabel='Square root of expected variance (%)',title=f'Forward variance curve: first {horizon} days');ax.grid(alpha=.2)
    axes[0].legend(frameon=False,fontsize=9);fig.suptitle('16 curve pieces · green bands are VIX 30-day windows');save(fig,'forward-variance-curve')
    tables={n:pd.read_csv(run/'checking_vix'/f'prices_canonical_{n}.csv') for n in ['baseline','selected']}
    fig,axes=plt.subplots(2,2,figsize=(13,8.5),layout='constrained')
    for ax,day in zip(axes.flat,[27,55,118,153]):
        q=tables['selected'].query('expiry_days == @day').sort_values('strike');x=q.K_over_F
        ax.fill_between(x,iv(q,q.price_bid_U),iv(q,q.price_ask_U),color='.2',alpha=.18,label='Market bid–ask')
        ax.plot(x,iv(q,q.price_mid_U),color='.1',label='Market mid')
        for name,color,style,label in [('baseline','#bc7628','--','Previous curve/readout'),('selected','#2269ad','-','Selected 24-parameter calibration')]:
            t=tables[name].query('expiry_days == @day').sort_values('strike');ax.plot(t.K_over_F,iv(t,t.model_price_U),color=color,ls=style,label=label)
        ax.fill_between(x,iv(q,q.model_price_U-2*q.mc_se_U),iv(q,q.model_price_U+2*q.mc_se_U),color='#2269ad',alpha=.14,label='Selected ±2 pointwise sampling SE')
        missing=~np.isfinite(iv(q,q.model_price_U))
        if missing.any():ax.plot(x[missing],np.full(missing.sum(),.015),'|',transform=ax.get_xaxis_transform(),color='red')
        ax.set_xscale('log');ticks=np.geomspace(x.min(),x.max(),5);ax.set_xticks(ticks,[f'{t:.2f}' for t in ticks]);ax.set_xticks([],minor=True)
        ax.set(title=f'VIX {day} days · {int(missing.sum())} selected IVs unavailable',xlabel='Strike / observed market forward',ylabel='Black call implied volatility (%)');ax.grid(alpha=.18)
    handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=2,frameon=False)
    fig.suptitle('VIX: separate reused numerical streams, parameters frozen before checking\nAll 4096 points · observed-forward call IV convention · no smoothing',fontsize=13)
    save(fig,'vix-separate-checking')
    path=out/'manifest.json';m=json.loads(path.read_text());m['output_sha256'].update({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs})
    m['source_sha256'][str(__file__)]=hashlib.sha256(open(__file__,'rb').read()).hexdigest()
    for name in ['baseline','selected']:
        source=run/'checking_vix'/f'prices_canonical_{name}.csv'
        m['source_sha256'][str(source)]=hashlib.sha256(source.read_bytes()).hexdigest()
    path.write_text(json.dumps(m,indent=2)+'\n')
    with (out/'README.md').open('a') as f:
        f.write('\n- [Forward variance curve](forward-variance-curve.png)\n- [VIX on separate reused numerical streams](vix-separate-checking.png)\n\nThe six original overlays use training draws. The separate VIX figure uses previously used checking streams; this is not a pristine holdout or market OOS test.\n')
