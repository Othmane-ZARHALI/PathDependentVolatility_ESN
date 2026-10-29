"""Raw price-derived smiles for the post-cap amplitude experiment."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .plot_maturity_curve import iv
from .cap_controlled_core import shape_metrics
from .full_readout_fit import dump


def render(run,out):
    out.mkdir(parents=True,exist_ok=True);selection=json.loads((run/'selection.json').read_text());names=selection['full_candidates'];profiles=json.loads((run/'profiles.json').read_text());outputs=[];sources=[Path(__file__),run/'selection.json',run/'profiles.json',run/'run.json']
    colors=dict(zip(names,['#2865ad','#de8a2b','#168958']));labels={r['name']:('Original model' if r['name']=='baseline' else f"gamma={r['spec'][1]:g}, half-life={r['spec'][3]}d") for r in profiles}
    def table(folder,name):
        path=folder/f'prices_canonical_{name}.csv';sources.append(path);return pd.read_csv(path)
    def save(fig,name):
        for ext in ['png','svg']:
            path=out/f'{name}.{ext}';fig.savefig(path,dpi=140,bbox_inches='tight');outputs.append(path)
        plt.close(fig)
    def finish(fig,axes,title):
        h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle(title)
    train={n:table(run/'full',n) for n in names};fresh={n:table(run/'fresh_check',n) for n in names}
    for scale in ['log','linear']:
        fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
        for ax,day in zip(axes.flat,[27,55,118,153]):
            q=fresh['baseline'];q=q[q.expiry_days.eq(day)].sort_values('strike');ax.fill_between(q.K_over_F,iv(q,q.price_bid_U),iv(q,q.price_ask_U),color='.25',alpha=.15,label='Market bid–ask');ax.plot(q.K_over_F,iv(q,q.price_mid_U),'k-',label='Market mid')
            for name in names:
                g=fresh[name];g=g[g.expiry_days.eq(day)].sort_values('strike');y=iv(g,g.model_price_U);ax.plot(g.K_over_F,y,color=colors[name],ls='--' if name=='baseline' else '-',label=labels[name]);ax.fill_between(g.K_over_F,iv(g,g.model_price_U-2*g.mc_se_U),iv(g,g.model_price_U+2*g.mc_se_U),color=colors[name],alpha=.08)
                bad=~np.isfinite(y)
                if bad.any():ax.plot(g.K_over_F[bad],np.full(bad.sum(),.01),'|',color=colors[name],transform=ax.get_xaxis_transform())
            ax.set(xscale=scale,title=f'VIX {day} days',xlabel='Strike / market forward',ylabel='Black call IV (%)');ax.grid(alpha=.18)
        finish(fig,axes,'Independent multiplier AFTER original cap · new VIX simulations\nFixed variance curve and floor; ±2 pointwise eight-scramble sampling SE');save(fig,'vix-fresh-'+scale)
    for kind in ['iv','signed-errors']:
        fig,axes=plt.subplots(4,2,figsize=(13,15),layout='constrained')
        for ax,day in zip(axes.flat,[29,57,119,211,302,666,1030]):
            q=train['baseline'];q=q[q.instrument.eq('SPX')&q.expiry_days.eq(day)].sort_values('strike')
            if kind=='iv':ax.fill_between(q.K_over_F,iv(q,q.price_bid_U),iv(q,q.price_ask_U),color='.3',alpha=.15,label='Market bid–ask');ax.plot(q.K_over_F,iv(q,q.price_mid_U),'k-',label='Market mid')
            else:ax.axhline(0,color='.3',lw=.7)
            for name in names:
                g=train[name];g=g[g.instrument.eq('SPX')&g.expiry_days.eq(day)].sort_values('strike');y=iv(g,g.model_price_U) if kind=='iv' else (g.model_price_U-g.price_mid_U)/g.price_scale_U
                ax.plot(g.K_over_F,y,color=colors[name],ls='--' if name=='baseline' else '-',label=labels[name]);bad=~np.isfinite(y)
                if bad.any():ax.plot(g.K_over_F[bad],np.full(bad.sum(),.01),'|',color=colors[name],transform=ax.get_xaxis_transform())
            ax.set(title=f'SPX {day} days',xlabel='Strike / market forward',ylabel='Black call IV (%)' if kind=='iv' else 'Signed error / quote scale');ax.grid(alpha=.18)
        axes.flat[-1].axis('off');axes.flat[-1].text(.04,.95,'Original cap, all 16 variance levels and floor fixed.\nIndependent normalized exponential OU amplitude.\n4096 paths at every strike.\nStock volatility and paths recomputed.\n\nNo clipping or smoothing; missing IVs marked.\nSame market quotes and old stock random primitives.\nNo adopted joint fit.',va='top',fontsize=11)
        finish(fig,axes,'SPX: full repricing after changing the amplitude channel\nTwo VIX selections frozen before stock pricing');save(fig,'spx-all-maturities-'+kind)
    fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained');q=fresh['baseline'];market=shape_metrics(q,q.price_mid_U.to_numpy());days=[r['days'] for r in market['rows']]
    for ax,key,title in [(axes[0],'wing_slope','Right-wing slope'),(axes[1],'wing_curvature','Right-wing curvature')]:
        ax.plot(days,[r[key] for r in market['rows']],'kD-',label='Market')
        for name in names:
            g=fresh[name];metric=shape_metrics(g,g.model_price_U.to_numpy());ax.plot(days,[r[key] for r in metric['rows']],marker='o',color=colors[name],ls='--' if name=='baseline' else '-',label=labels[name])
        ax.set(title=title,xlabel='VIX maturity (days)',ylabel='IV points / log-strike' if key.endswith('slope') else 'IV points / log-strike²');ax.grid(alpha=.18)
    finish(fig,axes,'Smile shape on fresh simulations · K / market forward > 2\nDescriptive regressions; curvature depends on strike coordinate');save(fig,'vix-wing-shape')
    fig,axes=plt.subplots(3,2,figsize=(13,13),layout='constrained');grid={r['name']:table(run,r['name']) for r in profiles}
    for row,gamma in enumerate([.25,.5,1.]):
        for col,scale in enumerate(['log','linear']):
            ax=axes[row,col];q=grid['baseline'];q=q[q.expiry_days.eq(153)].sort_values('strike');ax.plot(q.K_over_F,iv(q,q.price_mid_U),'k-',label='Market');ax.plot(q.K_over_F,iv(q,q.model_price_U),'--',color='#2865ad',label='Original')
            for half,color in zip([7,30,90],['#168958','#de8a2b','#9b54b8']):
                g=grid[f'g{gamma:g}_h{half}'];g=g[g.expiry_days.eq(153)].sort_values('strike');ax.plot(g.K_over_F,iv(g,g.model_price_U),color=color,label=f'{half}d half-life')
            ax.set(xscale=scale,title=f'153 days; gamma={gamma:g}; {scale} strike axis',xlabel='Strike / market forward',ylabel='Black call IV (%)');ax.grid(alpha=.18);ax.legend(frameon=False,fontsize=9)
    fig.suptitle('Complete strength/persistence comparison at the longest VIX expiry\nTraining simulations; all variants retained');save(fig,'vix-153-grid')
    fig,ax=plt.subplots(figsize=(9,6),layout='constrained')
    for name in names:
        path=run/'fresh_check'/f'forwards_{name}.csv';sources.append(path);f=pd.read_csv(path);ax.errorbar(f.expiry_days,f.model_forward,yerr=2*f.mc_se,marker='o',color=colors[name],label=labels[name],capsize=2)
    ax.plot(f.expiry_days,f.forward_vix_points,'kD-',label='Market');ax.set(xlabel='Expiry (days)',ylabel='VIX forward (points)',title='VIX forwards on fresh simulations · ±2 sampling SE');ax.grid(alpha=.18);ax.legend(frameon=False);save(fig,'vix-forwards')
    dump(out/'manifest.json',dict(source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in set(sources)},output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs},prices_clipped=False,market_out_of_sample=False))
    (out/'README.md').write_text('''# Independent amplitude after the original cap

The original capped normalized variance component is multiplied by an independent
mean-one OU exponential. The original curve, floor, reservoir and cap stay fixed.
Gamma is log-variance sensitivity; it is not a bounded multiplier range.

`vix-fresh-log` and `vix-fresh-linear` show new joint numerical draws, on both
strike axes. `vix-wing-shape` compares slope and curvature over K/marketF>2.
`spx-all-maturities-*` shows actual stock repricing. `vix-153-grid` retains all
nine training variants at the longest VIX expiry. `vix-forwards` checks forwards.

Lines connect raw price-derived Black call IVs using market forwards, without
smoothing/clipping. Missing IVs are marked at the foot of panels. Pointwise bands
are sampling SE only, excluding deterministic, pilot and selection uncertainty.
New numerical draws use the same market data; there is no market OOS claim.
See POSTCAP_VOV_RESULTS.md and the run's audit for results and limitations.
''')
    print('EXPORTED',len(outputs),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();render(a.run.resolve(),a.output.resolve())
