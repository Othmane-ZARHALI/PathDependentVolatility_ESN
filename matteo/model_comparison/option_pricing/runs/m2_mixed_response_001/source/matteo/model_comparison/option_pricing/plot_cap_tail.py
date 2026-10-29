"""All-maturity, unsmoothed cap profiles, prices, weighted tails and tradeoffs."""
import json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from .full_readout_fit import HERE
from .plot_maturity_curve import iv
from .cap_tail_study import VARIANTS


def render():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    run=HERE/'runs/m2_cap_tail_001';out=HERE/'visualizations/cap-tail-014';out.mkdir(parents=True,exist_ok=True);outputs=[];sources=[Path(__file__)]
    def save(fig,name):
        for ext in ['png','svg']:
            path=out/f'{name}.{ext}';fig.savefig(path,dpi=140,bbox_inches='tight');outputs.append(path)
        plt.close(fig)
    def table(folder,name):
        path=folder/f'prices_canonical_{name}.csv';sources.append(path);return pd.read_csv(path)
    base=table(run,'baseline');alltables={n:table(run,n) for n,*_ in VARIANTS};labels={'baseline':'Selected baseline','lower_L':'L − 0.75','higher_L':'L + 0.75','soft_k1':'κ = 1','sharp_k4':'κ = 4','algebraic':'Algebraic saturation','logtail025':'Log continuation a=.25','logtail050':'Log continuation a=.50','separate':'Separate channel caps','cluster050':'Cluster spike boost b=.5','cluster100':'Cluster spike boost b=1','higher_L_fit':'L+.75, VIX-only curve fit','logtail025_fit':'Log a=.25, VIX-only curve fit'}
    colors={n:plt.get_cmap('tab20')(i) for i,n in enumerate(labels)}
    colors.update(baseline='#2269ad',higher_L='#168c50',lower_L='#db8c20',soft_k1='#c33131',sharp_k4='#9467bd',algebraic='#a77b3a',logtail025='#ef8630',logtail050='#c33131',separate='#9467bd',cluster050='#ef8630',cluster100='#168c50',higher_L_fit='#168c50',logtail025_fit='#ef8630')
    groups=[('vix-cap-location',['baseline','lower_L','higher_L','soft_k1','sharp_k4']),('vix-cap-shapes',['baseline','algebraic','logtail025','logtail050','separate'])]
    if (run/'cluster_vov/run.json').exists():
        for n in ['cluster050','cluster100']:alltables[n]=table(run/'cluster_vov',n)
        groups.append(('vix-cluster-vol-of-vol',['baseline','cluster050','cluster100']))
    full=(run/'followup/run.json').exists()
    if full:
        for n in ['higher_L_fit','logtail025_fit']:alltables[n]=table(run/'followup',n)
        groups.append(('vix-focused-fits',['baseline','higher_L_fit','logtail025_fit']))
    smile=[]
    for name,names in groups:
        fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
        for ax,day in zip(axes.flat,[27,55,118,153]):
            q=base[base.expiry_days.eq(day)].sort_values('strike');x=q.K_over_F.to_numpy();ax.fill_between(x,iv(q,q.price_bid_U),iv(q,q.price_ask_U),color='.4',alpha=.18);ax.plot(x,iv(q,q.price_mid_U),color='black',label='Market (OTM converted to calls)')
            for n in names:
                t=alltables[n];t=t[t.instrument.eq('VIX')&t.expiry_days.eq(day)].sort_values('strike');y=iv(t,t.model_price_U);ax.plot(t.K_over_F,y,label=labels[n],lw=1.4,color=colors[n]);bad=~np.isfinite(y)
                if bad.any():ax.plot(t.K_over_F[bad],np.full(bad.sum(),.01),'|',transform=ax.get_xaxis_transform(),color='#b52e2e')
                for region,mask in [('ATM',abs(np.log(t.K_over_F))<=.2),('right',t.K_over_F>2)]:
                    valid=mask&np.isfinite(y);xm=np.log(t.K_over_F[valid]);yv=y[valid];market=iv(t,t.price_mid_U);slope=curv=np.nan
                    if len(xm)>=3:
                        pol=np.polyfit(xm,yv,2);slope=float(np.polyfit(xm,yv,1)[0]);curv=float(2*pol[0])
                    smile.append(dict(candidate=n,days=day,region=region,IV_mean=float(np.mean(yv)) if len(yv) else None,IV_MAE=float(np.nanmean(abs(y[mask]-market[mask]))),slope_per_logstrike=slope,quadratic_curvature_per_logstrike_squared=curv,available=int(valid.sum()),quotes=int(mask.sum())))
            ax.set(title=f'VIX {day} days',xlabel='Strike / observed market forward',ylabel='Black call IV (%)',xscale='log');ax.grid(alpha=.2)
        handles,lab=axes.flat[0].get_legend_handles_labels();fig.legend(handles,lab,loc='outside lower center',ncols=3,frameon=False);fig.suptitle('4096 points per strike · fixed OU reservoir · training simulations\n'+('VIX-focused fits; SPX tradeoffs shown separately' if 'fits' in name else 'Curve/floor fixed · every candidate renormalized under Q'));save(fig,name)
    pd.DataFrame(smile).drop_duplicates(['candidate','days','region']).to_csv(out/'smile-shape.csv',index=False)
    # Forwards for every profile; exact four market targets.
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for j,names in enumerate([['baseline','lower_L','higher_L','soft_k1','sharp_k4'],['baseline','algebraic','logtail025','logtail050','separate']]):
        for n in names:
            f=pd.read_csv(run/f'forwards_{n}.csv');sources.append(run/f'forwards_{n}.csv');axes[j].errorbar(f.expiry_days,f.model_forward,yerr=2*f.mc_se,marker='o',capsize=2,label=labels[n],color=colors[n])
        axes[j].plot(f.expiry_days,f.forward_vix_points,'kD--',label='Market');axes[j].set(xlabel='Expiry (days)',ylabel='VIX forward',title='Actual E[VIX], with ±2 scramble SE');axes[j].grid(alpha=.2);axes[j].legend(fontsize=8,ncols=2)
    save(fig,'vix-forwards')
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,day in zip(axes.flat,[27,55,118,153]):
        for n in ['baseline','higher_L','logtail025','logtail050']:
            with np.load(run/f'paths_{n}.npz') as a:
                v=a[f'vix_{day}'];w=a[f'weights_vix_{day}'];K=np.geomspace(20,250,100);pay=w[:,None]*np.maximum(v[:,None]-K,0);blocks=pay.reshape(8,512,len(K)).mean(1);mean=blocks.mean(0);se=blocks.std(0,ddof=1)/np.sqrt(8)
                ax.plot(K,mean,label=labels[n],color=colors[n]);ax.fill_between(K,np.where(mean-2*se>0,mean-2*se,np.nan),mean+2*se,alpha=.08,color=colors[n])
        q=base[base.expiry_days.eq(day)];wing=q.K_over_F>1;ax.plot(q.strike[wing],q.price_mid_U[wing],'k.',ms=3,label='Market calls');ax.set(xscale='log',yscale='log',xlabel='VIX strike',ylabel='Raw weighted call expectation',title=f'{day} days · ±2 eight-net SE');ax.grid(alpha=.2)
    h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=3,frameon=False);fig.suptitle('Likelihood-weighted tail payoffs, before signed control correction');save(fig,'weighted-tail-calls')
    # Signed actual controlled prices remain visible even when IV is undefined.
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,day in zip(axes.flat,[27,55,118,153]):
        for n in ['baseline','higher_L','logtail025','logtail050','separate']:
            t=alltables[n];t=t[t.expiry_days.eq(day)&t.instrument.eq('VIX')].sort_values('strike');ax.plot(t.K_over_F,(t.model_price_U-t.price_mid_U)/t.price_scale_U,label=labels[n],color=colors[n])
        ax.axhline(0,color='k',lw=.6);ax.grid(alpha=.2);ax.set(title=f'VIX {day} days',xlabel='Strike / market forward',ylabel='Signed call error / quote scale')
    h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=3,frameon=False);save(fig,'vix-signed-price-errors')
    norms=pd.read_csv(run/'normalizers.csv');fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    eta=np.linspace(-1,15,600);p=np.array(json.loads((run/'design.json').read_text())['baseline'])
    for n,mode,a,dl,kap in VARIANTS[:-1]:
        L=p[5]+dl;k=kap or p[6];h=L-np.logaddexp(0,k*(L-eta))/k
        if mode==1:h=L-.5*((L-eta)+np.hypot(L-eta,2*np.log(2)/k))
        if mode==2:h+=a*np.log1p(np.logaddexp(0,k*(eta-L))/k)
        M=float(norms[(norms.candidate==n)&norms.days.eq(55)].M.iloc[0]);axes[0].plot(eta,h,label=labels[n],color=colors[n]);axes[1].plot(eta,np.exp(2*h)/M,label=labels[n],color=colors[n])
    axes[0].set(xlabel='Score η',ylabel='Tempered score h',title='Transition and growth beyond it');axes[1].set(xlabel='Score η',ylabel='G / M_Q(day55)',yscale='log',title='Candidate-specific normalization included')
    for ax in axes:ax.grid(alpha=.2)
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=4,frameon=False);save(fig,'cap-and-normalization')
    if full:
        parent=HERE/'runs/m2_maturity_curve_fit_001';old=table(parent,'selected');names=['higher_L_fit','logtail025_fit']
        for kind in ['iv','price-errors']:
            fig,axes=plt.subplots(4,2,figsize=(12,13),layout='constrained')
            for ax,day in zip(axes.flat,[29,57,119,211,302,666,1030]):
                q=old[old.instrument.eq('SPX')&old.expiry_days.eq(day)].sort_values('strike');x=q.K_over_F.to_numpy()
                if kind=='iv':ax.plot(x,iv(q,q.price_mid_U),color='black',label='Market')
                else:ax.axhline(0,color='black',lw=.6)
                for n,t in [('baseline',old)]+[(n,alltables[n]) for n in names]:
                    g=t[t.instrument.eq('SPX')&t.expiry_days.eq(day)].sort_values('strike');y=iv(g,g.model_price_U) if kind=='iv' else (g.model_price_U-g.price_mid_U)/g.price_scale_U;ax.plot(g.K_over_F,y,label=labels[n],color=colors[n]);
                    if kind=='iv':
                        bad=~np.isfinite(y)
                        if bad.any():ax.plot(g.K_over_F[bad],np.full(bad.sum(),.01),'|',transform=ax.get_xaxis_transform(),color='#b52e2e')
                ax.set(title=f'SPX {day} days',xlabel='Strike / market forward',ylabel='Black call IV (%)' if kind=='iv' else 'Signed price error / quote scale');ax.grid(alpha=.2)
            axes.flat[-1].axis('off');axes.flat[-1].text(.02,.95,'VIX-focused curve fits, NOT accepted joint fits.\nAll 4096 stock paths regenerated from cached OU\nfeatures and equity increments at each expiry.\nMissing IVs remain unavailable; signed prices retained.\nSampling/timestep/normalizer errors unresolved.',va='top');h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle('Actual SPX tradeoffs after VIX-focused fitting · all retained maturities');save(fig,'spx-all-maturities-'+kind)
        # Numerical checking overlays after parameter freeze.
        fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained');old=table(parent/'checking_vix','selected')
        for ax,day in zip(axes.flat,[27,55,118,153]):
            q=old[old.expiry_days.eq(day)].sort_values('strike');ax.plot(q.K_over_F,iv(q,q.price_mid_U),'k-',label='Market')
            for n,t in [('baseline',old)]+[(n,table(run/'followup/checking_vix',n)) for n in names]:
                g=t[t.expiry_days.eq(day)].sort_values('strike');y=iv(g,g.model_price_U);ax.plot(g.K_over_F,y,label=labels[n],color=colors[n]);lo=iv(g,g.model_price_U-2*g.mc_se_U);hi=iv(g,g.model_price_U+2*g.mc_se_U);ax.fill_between(g.K_over_F,lo,hi,alpha=.10,color=colors[n])
            ax.set(xscale='log',title=f'VIX {day} days',xlabel='Strike / market forward',ylabel='Black call IV (%)');ax.grid(alpha=.2)
        h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle('Previously used separate VIX streams · ±2 pointwise scramble SE\nSame in-sample market quotes; not untouched holdout or market out-of-sample');save(fig,'vix-separate-checking')
    manifest=dict(source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in set(sources)},output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs},SPX_complete=full,notes='No IV smoothing/clipping. Profile charts omit bands for legibility; all SEs saved. Slope/curvature are descriptive least-squares summaries of raw IVs, not smoothed prices. All market data in sample.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(len(outputs),'exports',flush=True)

if __name__=='__main__':render()
