"""Raw all-maturity comparisons for the dedicated sensitivity experiment."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .plot_maturity_curve import iv
from .full_readout_fit import dump
from .cap_controlled_core import shape_metrics


def render(run,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    out.mkdir(parents=True,exist_ok=True)
    profiles=json.loads((run/'profiles.json').read_text());selection=json.loads((run/'selection.json').read_text());names=selection['full_candidates'];sources=[run/'run.json',run/'profiles.json',run/'selection.json',Path(__file__)];outputs=[]
    colors=dict(zip(names,['#2865ad','#de8a2b','#168958','#9b54b8','#da5369','#4899a6']))
    def label(name):
        row=next(x for x in profiles if x['name']==name);cap='Current cap' if row['cap']=='cap0' else 'Slow growth'
        return cap+'; no new state' if not row['layer'] else f"{cap}; {row['spec'][4]:g}×, {row['spec'][3]}d half-life"
    def table(folder,name):
        p=folder/f'prices_canonical_{name}.csv';sources.append(p);return pd.read_csv(p)
    train={n:table(run/'full',n) for n in names};fresh={n:table(run/'fresh_check',n) for n in names}
    def save(fig,name):
        for ext in ['png','svg']:
            p=out/f'{name}.{ext}';fig.savefig(p,dpi=140,bbox_inches='tight');outputs.append(p)
        plt.close(fig)
    def finish(fig,axes,title):
        h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle(title)
    for ins,layout,days in [('VIX',(2,2),[27,55,118,153]),('SPX',(4,2),[29,57,119,211,302,666,1030])]:
        for kind in ['iv','signed-errors']:
            fig,axes=plt.subplots(*layout,figsize=(13,4*layout[0]),layout='constrained')
            for ax,day in zip(axes.flat,days):
                q=train['cap0_ref'];q=q[q.instrument.eq(ins)&q.expiry_days.eq(day)].sort_values('strike')
                if kind=='iv':
                    ax.fill_between(q.K_over_F,iv(q,q.price_bid_U),iv(q,q.price_ask_U),color='.3',alpha=.15,label='Market bid–ask');ax.plot(q.K_over_F,iv(q,q.price_mid_U),'k-',label='Market mid')
                else:ax.axhline(0,color='.2',lw=.7)
                for name in names:
                    g=train[name];g=g[g.instrument.eq(ins)&g.expiry_days.eq(day)].sort_values('strike')
                    y=iv(g,g.model_price_U) if kind=='iv' else (g.model_price_U-g.price_mid_U)/g.price_scale_U
                    ax.plot(g.K_over_F,y,color=colors[name],ls='--' if name.endswith('_ref') else '-',label=label(name),lw=1.5)
                    if kind=='iv' and (~np.isfinite(y)).any():
                        bad=~np.isfinite(y);ax.plot(g.K_over_F[bad],np.full(bad.sum(),.01),'|',color=colors[name],transform=ax.get_xaxis_transform())
                ax.set(title=f'{ins} {day} days',xlabel='Strike / market forward',ylabel='Black call IV (%)' if kind=='iv' else 'Signed call error / quote scale')
                if ins=='VIX':ax.set_xscale('log')
                ax.grid(alpha=.18)
            for ax in axes.flat[len(days):]:
                ax.axis('off');ax.text(.02,.96,'FIXED EXPECTED VARIANCE CURVE\n\nOne new independent OU sensitivity state.\nAll 16 variance levels and floor fixed.\nOriginal reservoir and correlations fixed.\n4096 paths at every strike.\nStock volatility and paths recomputed.\n\nTraining simulations; no accepted joint fit.\nMissing IVs marked; signed prices retained.',va='top',fontsize=11)
            finish(fig,axes,f'{ins}: dedicated stochastic sensitivity · all agreed maturities\nMatched cap comparisons, no curve or floor refit');save(fig,ins.lower()+'-all-maturities-'+kind)
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for ax,day in zip(axes.flat,[27,55,118,153]):
        q=fresh['cap0_ref'];q=q[q.expiry_days.eq(day)].sort_values('strike');ax.plot(q.K_over_F,iv(q,q.price_mid_U),'k-',label='Market mid')
        for name in names:
            g=fresh[name];g=g[g.expiry_days.eq(day)].sort_values('strike');y=iv(g,g.model_price_U)
            ax.plot(g.K_over_F,y,color=colors[name],ls='--' if name.endswith('_ref') else '-',label=label(name))
            ax.fill_between(g.K_over_F,iv(g,g.model_price_U-2*g.mc_se_U),iv(g,g.model_price_U+2*g.mc_se_U),color=colors[name],alpha=.08)
            if (~np.isfinite(y)).any():
                bad=~np.isfinite(y);ax.plot(g.K_over_F[bad],np.full(bad.sum(),.01),'|',color=colors[name],transform=ax.get_xaxis_transform())
        ax.set(xscale='log',title=f'VIX {day} days',xlabel='Strike / market forward',ylabel='Black call IV (%)');ax.grid(alpha=.18)
    finish(fig,axes,'NEW joint Gaussian VIX draws generated after selection\n±2 pointwise eight-scramble SE · same in-sample market quotes');save(fig,'vix-new-checking')
    # Show all six strength/persistence cases, separately for each cap.
    for cap in ['cap0','cap05']:
        rows=[r for r in profiles if r['cap']==cap];tables={r['name']:table(run,r['name']) for r in rows}
        fig,axes=plt.subplots(2,2,figsize=(13,10),layout='constrained')
        palette=plt.get_cmap('viridis')(np.linspace(.1,.9,6));gridcolors={r['name']:palette[i] for i,r in enumerate(rows[1:])}
        for ax,day in zip(axes.flat,[27,55,118,153]):
            q=tables[cap+'_ref'];q=q[q.expiry_days.eq(day)].sort_values('strike');ax.plot(q.K_over_F,iv(q,q.price_mid_U),'k-',label='Market mid');ax.plot(q.K_over_F,iv(q,q.model_price_U),'--',color='#2865ad',label='Same cap, no new state')
            for row in rows[1:]:
                name=row['name'];g=tables[name];g=g[g.expiry_days.eq(day)].sort_values('strike')
                ax.plot(g.K_over_F,iv(g,g.model_price_U),color=gridcolors[name],label=f"{row['spec'][4]:g}× sensitivity, {row['spec'][3]}d half-life",lw=1.25)
            ax.set(xscale='log',title=f'VIX {day} days',xlabel='Strike / market forward',ylabel='Black call IV (%)');ax.grid(alpha=.18)
        finish(fig,axes,('Current cap' if cap=='cap0' else 'Logarithmic slow growth')+': complete fixed-curve VIX grid\nTwo strengths × three persistence half-lives; 4096 points per strike');save(fig,'vix-grid-'+cap)
    fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
    for cap,ax in zip(['cap0','cap05'],axes):
        rows=[r for r in profiles if r['cap']==cap and r['layer']]
        for factor in [1.5,2.]:
            rr=sorted([r for r in rows if r['spec'][4]==factor],key=lambda x:x['spec'][3]);ax.plot([r['spec'][3] for r in rr],[r['shape']['wing_MAE'] for r in rr],'o-',label=f'{factor:g}× sensitivity')
        base=next(r for r in profiles if r['name']==cap+'_ref');ax.axhline(base['shape']['wing_MAE'],ls='--',color='.4',label='Same cap, no state')
        ax.set(title='Current cap' if cap=='cap0' else 'Logarithmic slow growth',xlabel='Sensitivity-state half-life (days)',ylabel='Right-wing IV MAE (percentage points)',xticks=[7,30,90]);ax.grid(alpha=.18);ax.legend(frameon=False)
    fig.suptitle('Does persistence help the wing? Lower error is better\nFixed entire variance curve and floor; averages give each VIX maturity equal weight');save(fig,'strength-persistence-grid')
    fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
    q=fresh['cap0_ref'];market=shape_metrics(q,q.price_mid_U.to_numpy());days=[x['days'] for x in market['rows']]
    for ax,key,title in [(axes[0],'wing_slope','Right-wing slope'),(axes[1],'wing_curvature','Right-wing curvature')]:
        ax.plot(days,[x[key] for x in market['rows']],'kD-',label='Market')
        for name in names:
            g=fresh[name];metrics=shape_metrics(g,g.model_price_U.to_numpy())
            ax.plot(days,[x[key] for x in metrics['rows']],marker='o',color=colors[name],ls='--' if name.endswith('_ref') else '-',label=label(name))
        ax.set(title=title,xlabel='VIX expiry (days)',ylabel='IV percentage points / log-strike' if key=='wing_slope' else 'IV percentage points / log-strike²');ax.grid(alpha=.18)
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle('Closer wing IV levels do not necessarily mean a steeper wing\nDescriptive fits over K / market forward > 2 · NEW numerical draws');save(fig,'vix-wing-slope-curvature')
    fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
    for ax,folder,title in [(axes[0],run/'full','Training'),(axes[1],run/'fresh_check','New numerical checking')]:
        for name in names:
            path=folder/f'forwards_{name}.csv';sources.append(path);f=pd.read_csv(path);ax.errorbar(f.expiry_days,f.model_forward,yerr=2*f.mc_se,marker='o',color=colors[name],ls='--' if name.endswith('_ref') else '-',label=label(name),capsize=2)
        ax.plot(f.expiry_days,f.forward_vix_points,'kD-',label='Market forward');ax.set(title=title,xlabel='Expiry (days)',ylabel='VIX forward (points)');ax.grid(alpha=.18)
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle('VIX forwards · ±2 pointwise sampling SE');save(fig,'vix-forwards')
    fig,axes=plt.subplots(2,2,figsize=(13,10),layout='constrained')
    for ax,day in zip(axes.flat,[27,55,118,153]):
        for name in names:
            path=run/'fresh_check'/f'paths_{name}.npz';sources.append(path)
            with np.load(path) as a:
                v=a[f'vix_{day}'];w=a[f'weights_vix_{day}'];K=np.geomspace(20,250,90)
                b=(w[:,None]*np.maximum(v[:,None]-K,0)).reshape(8,512,len(K)).mean(1);m=b.mean(0);se=b.std(0,ddof=1)/np.sqrt(8)
                ax.plot(K,np.where(m>0,m,np.nan),color=colors[name],ls='--' if name.endswith('_ref') else '-',label=label(name));ax.fill_between(K,np.where(m-2*se>0,m-2*se,np.nan),m+2*se,color=colors[name],alpha=.08)
        q=fresh['cap0_ref'];q=q[q.expiry_days.eq(day)&q.K_over_F.gt(1)];ax.plot(q.strike,q.price_mid_U,'k.',label='Market calls',ms=3)
        ax.set(xscale='log',yscale='log',title=f'VIX {day} days',xlabel='VIX strike',ylabel='Raw likelihood-weighted call');ax.grid(alpha=.18)
    finish(fig,axes,'New-stream weighted tail contributions · ±2 scramble SE\nRaw payoffs exclude signed controls; complete controlled estimates saved separately');save(fig,'vix-weighted-tails')
    manifest=dict(source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in set(sources)},output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs},market_out_of_sample=False,uncertainty='Pointwise sampling only; excludes deterministic, pilot and parameter uncertainty',prices_smoothed_or_clipped=False)
    dump(out/'manifest.json',manifest)
    (out/'README.md').write_text('''# Dedicated stochastic sensitivity: reading the figures

The new state randomly changes spike sensitivity. Its half-life controls how long
an episode lasts. The entire expected variance curve, floor and original reservoir
are fixed. The two outer transformations are the current cap and a logarithmic
continuation with strength .5. This is a bounded diagnostic, not a new accepted fit.

- `vix-grid-cap0` / `vix-grid-cap05`: every proposed strength/persistence case.
- `strength-persistence-grid`: wing errors; lower is closer to market.
- `vix-all-maturities-iv`: matched training selections and cap-only controls.
- `spx-all-maturities-iv`: actual full stock repricing at all seven dates.
- `*-signed-errors`: signed call errors in quote scales, including failures.
- `vix-new-checking`: new joint scrambles generated after selection, with bands.
- `vix-forwards` / `vix-weighted-tails`: level and rare-tail diagnostics.

Solid lines add the new state; dashed lines are no-layer references. Market inputs
remain original OTM quotes converted to common calls. IVs use market forwards.
All 4096 observations enter every strike. Curves connect raw price-derived IVs;
there is no smoothing or price clipping. Missing IVs are marked. Sampling bands
exclude integration, stock timestep, pilot and parameter uncertainty. New numerical
streams use the same in-sample market targets, with no market out-of-sample claim.
See STOCHASTIC_VOV_RESULTS.md and the run's audit for evidence and limitations.
''')
    print(len(outputs),'exports',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();render(a.run.resolve(),a.output.resolve())
