"""Unsmoothed all-maturity plots of the controlled cap experiment."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from .full_readout_fit import HERE
from .plot_maturity_curve import iv


def render(run,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    out.mkdir(parents=True,exist_ok=True);d=json.loads((run/'design.json').read_text());selection=json.loads((run/'selection.json').read_text())['candidates'];prior=Path(d['parent']);outputs=[];sources=[Path(__file__),run/'design.json',run/'selection.json',run/'run.json',run/'profiles.json'];names=['baseline']+[r['name'] for r in selection];colors=dict(zip(names,['#2269ad','#168c50','#ef8630']));labels={'baseline':'Selected baseline','shape00':'Current cap, four window adjustments'}
    for row in selection:labels[row['name']]='Higher cap, limited curve changes' if row['family']=='capped' else 'Slow growth, limited curve changes'
    def save(fig,name):
        for ext in ['png','svg']:
            p=out/f'{name}.{ext}';fig.savefig(p,dpi=140,bbox_inches='tight');outputs.append(p)
        plt.close(fig)
    def read(folder,name):
        p=folder/f'prices_canonical_{name}.csv';sources.append(p);return pd.read_csv(p)
    baseline=read(prior,'selected');train={'baseline':baseline};fresh={n:read(run/'fresh_check',n) for n in names}
    for n in names[1:]:train[n]=read(run/'full',n)
    control=read(run,'shape00')
    for ins,layout in [('VIX',(2,2)),('SPX',(4,2))]:
        dates=[27,55,118,153] if ins=='VIX' else [29,57,119,211,302,666,1030]
        for kind in ['iv','signed-errors']:
            fig,axes=plt.subplots(*layout,figsize=(13,4*layout[0]),layout='constrained')
            for ax,day in zip(axes.flat,dates):
                q=baseline[baseline.instrument.eq(ins)&baseline.expiry_days.eq(day)].sort_values('strike');x=q.K_over_F.to_numpy()
                if kind=='iv':ax.fill_between(x,iv(q,q.price_bid_U),iv(q,q.price_ask_U),color='.3',alpha=.15,label='Market bid–ask');ax.plot(x,iv(q,q.price_mid_U),'k-',label='Market mid')
                else:ax.axhline(0,color='.2',lw=.7)
                for n in names:
                    g=train[n];g=g[g.instrument.eq(ins)&g.expiry_days.eq(day)].sort_values('strike');y=iv(g,g.model_price_U) if kind=='iv' else (g.model_price_U-g.price_mid_U)/g.price_scale_U
                    ax.plot(g.K_over_F,y,color=colors[n],lw=1.6,label=labels[n]);bad=~np.isfinite(y)
                    if kind=='iv' and bad.any():ax.plot(g.K_over_F[bad],np.full(bad.sum(),.01),'|',color='#bd3030',transform=ax.get_xaxis_transform())
                if ins=='VIX':
                    g=control[control.expiry_days.eq(day)].sort_values('strike');y=iv(g,g.model_price_U) if kind=='iv' else (g.model_price_U-g.price_mid_U)/g.price_scale_U;ax.plot(g.K_over_F,y,'--',color='.55',lw=1.1,label=labels['shape00']);ax.set_xscale('log')
                ax.set(title=f'{ins} {day} days',xlabel='Strike / observed market forward',ylabel='Black call IV (%)' if kind=='iv' else 'Signed call error / quote scale');ax.grid(alpha=.18)
            for ax in axes.flat[len(dates):]:
                ax.axis('off');ax.text(.02,.95,'FLOOR AND OFF-WINDOW CURVE FIXED\n\nOnly four broad VIX-window adjustments.\nAll 4096 stock paths recomputed per expiry.\nOU reservoir and equity shocks reused.\nOriginal OTM market sources, cutoff 3.\nNo price clipping or IV smoothing.\nNo new accepted joint calibration.\nSampling and deterministic errors remain.',va='top',fontsize=11)
            h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle(f'{ins}: a cleaner cap comparison · all agreed maturities\nTraining simulations; exact floor and off-window curve preserved');save(fig,ins.lower()+'-all-maturities-'+kind)
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for ax,day in zip(axes.flat,[27,55,118,153]):
        q=fresh['baseline'];q=q[q.expiry_days.eq(day)].sort_values('strike');ax.plot(q.K_over_F,iv(q,q.price_mid_U),'k-',label='Market mid')
        for n in names:
            g=fresh[n];g=g[g.expiry_days.eq(day)].sort_values('strike');y=iv(g,g.model_price_U);ax.plot(g.K_over_F,y,color=colors[n],label=labels[n]);ax.fill_between(g.K_over_F,iv(g,g.model_price_U-2*g.mc_se_U),iv(g,g.model_price_U+2*g.mc_se_U),color=colors[n],alpha=.12);bad=~np.isfinite(y)
            if bad.any():ax.plot(g.K_over_F[bad],np.full(bad.sum(),.01),'|',color='#bd3030',transform=ax.get_xaxis_transform())
        ax.set(xscale='log',title=f'VIX {day} days',xlabel='Strike / market forward',ylabel='Black call IV (%)');ax.grid(alpha=.18)
    h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle('Genuinely NEW VIX draws, after selection froze\n4096 points per strike · ±2 pointwise scramble SE · same in-sample market quotes');save(fig,'vix-fresh-checking')
    # Show that no unrelated variance curve pieces moved.
    ends=np.array(json.loads((prior/'run.json').read_text())['ends_days']);base=np.array(d['baseline']);fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for ax,xlim in zip(axes,[(0,211),(0,1030)]):
        for n,p in [('baseline',base)]+[(z['name'],np.array(z['parameters'])) for z in selection]:
            xx=np.r_[0,ends];yy=100*np.sqrt(np.r_[p[7:-1],p[-2]]);ax.step(xx,yy,where='post',color=colors[n],label=labels[n])
        for l,r in [(27,57),(55,85),(118,148),(153,183)]:ax.axvspan(l,r,color='.6',alpha=.1)
        ax.set(xlim=xlim,xlabel='Calendar day',ylabel='Square root of expected variance (%)',title='Only shaded VIX windows may change');ax.grid(alpha=.2)
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle('Floor remains 5.138% · off-window levels are exactly unchanged');save(fig,'bounded-variance-curve')
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for ax,folder,title in [(axes[0],run,'Training'),(axes[1],run/'fresh_check','New numerical checking')]:
        for n in names:
            p=(prior/'forwards_selected.csv') if folder==run and n=='baseline' else ((run/'full'/f'forwards_{n}.csv') if folder==run else folder/f'forwards_{n}.csv');sources.append(p);f=pd.read_csv(p);ax.errorbar(f.expiry_days,f.model_forward,yerr=2*f.mc_se,marker='o',color=colors[n],capsize=3,label=labels[n])
        ax.plot(f.expiry_days,f.forward_vix_points,'kD--',label='Market');ax.set(title=title,xlabel='Expiry (days)',ylabel='Actual model VIX forward');ax.grid(alpha=.2)
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle('VIX forwards · ±2 pointwise sampling SE');save(fig,'vix-forwards')
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for ax,day in zip(axes.flat,[27,55,118,153]):
        for n in names:
            path=run/'fresh_check'/f'paths_{n}.npz';sources.append(path)
            with np.load(path) as a:
                v=a[f'vix_{day}'];w=a[f'weights_vix_{day}'];K=np.geomspace(20,250,100);blocks=(w[:,None]*np.maximum(v[:,None]-K,0)).reshape(8,512,len(K)).mean(1);m=blocks.mean(0);se=blocks.std(0,ddof=1)/np.sqrt(8);ax.plot(K,m,color=colors[n],label=labels[n]);ax.fill_between(K,np.where(m-2*se>0,m-2*se,np.nan),m+2*se,color=colors[n],alpha=.12)
        g=fresh['baseline'];g=g[g.expiry_days.eq(day)&g.K_over_F.gt(1)];ax.plot(g.strike,g.price_mid_U,'k.',ms=3,label='Market calls');ax.set(xscale='log',yscale='log',title=f'VIX {day} days',xlabel='VIX strike',ylabel='Raw likelihood-weighted call');ax.grid(alpha=.2)
    h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle('New-stream tail payoff contributions · ±2 scramble SE\nRaw estimates exclude signed control corrections; controlled prices are saved separately');save(fig,'fresh-weighted-tail-calls')
    profiles=json.loads((run/'profiles.json').read_text());fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for row in profiles:
        color='#168c50' if row['continuation']==0 else '#ef8630';label=f"L+{row['cap_shift']:g}, a={row['continuation']:g}";marker='o' if row['preferred'] else 'x'
        axes[0].scatter(row['shape']['wing_MAE'],row['fit']['native_loss'],color=color,marker=marker);axes[0].annotate(label,(row['shape']['wing_MAE'],row['fit']['native_loss']),xytext=((-40,-25) if row['name']=='shape02' else ((6,10) if row['name']=='shape03' else (3,3))),textcoords='offset points',fontsize=7)
        axes[1].scatter(row['shape']['wing_MAE'],row['shape']['ATM_MAE'],color=color,marker=marker);axes[1].annotate(label,(row['shape']['wing_MAE'],row['shape']['ATM_MAE']),xytext=(3,3),textcoords='offset points',fontsize=7)
    b=json.loads((run/'baseline_replay.json').read_text());axes[0].scatter(b['shape']['wing_MAE'],b['native_vix_plus_forward'],marker='D',color=colors['baseline'],label='Unchanged baseline');axes[1].scatter(b['shape']['wing_MAE'],b['shape']['ATM_MAE'],marker='D',color=colors['baseline'])
    axes[0].set(xlabel='Right-wing IV MAE (percentage points)',ylabel='Native VIX + forward loss',title='Better wing shape can cost price fit');axes[1].set(xlabel='Right-wing IV MAE (percentage points)',ylabel='Near-ATM IV MAE (percentage points)',title='Compare both parts of the smile')
    for ax in axes:ax.grid(alpha=.2);ax.margins(.18)
    axes[0].legend(frameon=False);fig.suptitle('Nine bounded shapes · four modest curve adjustments each\nAll nine fits reach the common-call loss safeguard; circles meet preference screens');save(fig,'nine-shape-tradeoffs')
    manifest=dict(source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in set(sources)},output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs},uncertainty='Pointwise scramble sampling bands; excludes pilot, parameter and deterministic errors',all_market_data_in_sample=True,fresh_numerical_streams=True)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(len(outputs),'exports')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();render(a.run.resolve(),a.output.resolve())
