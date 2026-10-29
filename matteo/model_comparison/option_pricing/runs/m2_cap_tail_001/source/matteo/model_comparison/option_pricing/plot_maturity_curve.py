"""Inspect all retained maturities on one canonical-call IV convention."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .plot_joint_curve_fit import implied_vol


def iv(q,price):
    return np.array([implied_vol(float(x),float(r.forward),float(r.strike),float(r.T_years),'call')
                     if x>0 else np.nan for r,x in zip(q.itertuples(),price)],float)


def render(run,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    out.mkdir(parents=True,exist_ok=True)
    tables={n:pd.read_csv(run/f'prices_canonical_{n}.csv') for n in ['baseline','selected']}
    f={n:pd.read_csv(run/f'forwards_{n}.csv') for n in tables}
    outputs=[]
    def save(fig,name):
        for ext in ['png','svg']:
            p=out/f'{name}.{ext}';fig.savefig(p,dpi=145,bbox_inches='tight');outputs.append(p)
        plt.close(fig)
    for instrument,layout in [('SPX',(4,2)),('VIX',(2,2))]:
        fig,axes=plt.subplots(*layout,figsize=(13,4.05*layout[0]),layout='constrained')
        dates=sorted(tables['selected'].query('instrument == @instrument').expiry_days.unique())
        for ax,days in zip(axes.flat,dates):
            groups={n:q[(q.instrument==instrument)&(q.expiry_days==days)].sort_values('strike') for n,q in tables.items()}
            q=groups['selected'];x=q.K_over_F.to_numpy()
            ax.fill_between(x,iv(q,q.price_bid_U),iv(q,q.price_ask_U),color='.2',alpha=.18,label='Market bid–ask (original OTM sources)')
            ax.plot(x,iv(q,q.price_mid_U),color='.1',lw=1.3,label='Market mid, parity converted')
            for name,color,style,label in [('baseline','#bc7628','--','Previous fitted parameters'),('selected','#2269ad','-','New 24-parameter calibration')]:
                t=groups[name];ax.plot(t.K_over_F,iv(t,t.model_price_U),color=color,ls=style,lw=1.6,label=label)
            low=iv(q,q.model_price_U-2*q.mc_se_U);high=iv(q,q.model_price_U+2*q.mc_se_U)
            ax.fill_between(x,low,high,color='#2269ad',alpha=.14,label='Selected ±2 pointwise sampling SE')
            unavailable=~np.isfinite(iv(q,q.model_price_U))
            if unavailable.any():
                ax.plot(x[unavailable],np.full(unavailable.sum(),.015),'|',transform=ax.get_xaxis_transform(),color='#ba3030',markersize=8)
            ax.text(.02,.03,f'{int(unavailable.sum())} selected IVs unavailable · effective paths {q.effective_paths.min():.0f}–{q.effective_paths.max():.0f}',transform=ax.transAxes,fontsize=8)
            ax.axvline(1,color='.6',lw=.8);ax.set_xscale('log');ticks=np.geomspace(x.min(),x.max(),5)
            ax.set_xticks(ticks,[f'{t:.2f}' for t in ticks]);ax.set_xticks([],minor=True);ax.grid(alpha=.18)
            ax.set(title=f'{instrument} {days} days · {len(q)} original OTM quotes',xlabel='Strike / observed market forward',ylabel='Black call implied volatility (%)')
        for ax in axes.flat[len(dates):]:
            ax.axis('off');ax.text(.02,.97,'ONE CALL PAYOFF AT EVERY STRIKE\n\nMarket OTM puts converted using market F.\nNo additional ITM observations.\n\n4096 paths at EVERY strike.\nEight independent VIX scrambles.\n\nSixteen variance pieces; 24 parameters.\nCalibration simulations.\nNo IV smoothing.\n\nMissing IVs marked in red; prices retained.\nSampling bands exclude deterministic bias\nand independent-pilot uncertainty.',va='top',fontsize=11)
        handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=2,frameon=False)
        note='Model and market VIX forwards can differ; see forward diagnostics.' if instrument=='VIX' else 'Same OTM cutoff-3 targets; separate put/call IV branches removed.'
        fig.suptitle(f'{instrument}: maturity-curve calibration with 4096 paths everywhere\n'+note,fontsize=14)
        save(fig,instrument.lower()+'-all-maturity-fit')
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    axes[0].plot(f['selected'].expiry_days,f['selected'].forward_vix_points,'D-',color='.1',label='Observed market forward')
    for n,color,label in [('baseline','#bc7628','Previous fit'),('selected','#2269ad','New calibration')]:
        t=f[n];axes[0].errorbar(t.expiry_days,t.model_forward,yerr=2*t.mc_se,marker='o',color=color,capsize=3,label=label)
    axes[0].set(title='Actual model VIX expectation remains a fitting target',xlabel='Expiry (days)',ylabel='VIX forward (points)');axes[0].legend(frameon=False)
    d=np.linspace(0,4,401);eff=np.full_like(d,4096.)
    axes[1].plot(d,eff,color='#2269ad',lw=2,label='Paths / VIX points used per strike')
    axes[1].set_ylim(0,4600)
    axes[1].set(title='All simulated paths used at every strike',xlabel='|log(K/F)| / (ATM proxy × √T)',ylabel='Number of paths');axes[1].legend(frameon=False)
    for ax in axes:ax.grid(alpha=.2)
    save(fig,'forwards-and-path-budgets')
    for instrument,layout in [('SPX',(4,2)),('VIX',(2,2))]:
        fig,axes=plt.subplots(*layout,figsize=(13,3.2*layout[0]),layout='constrained');q=tables['selected'];dates=sorted(q[q.instrument==instrument].expiry_days.unique())
        for ax,days in zip(axes.flat,dates):
            for n,color,label in [('baseline','#bc7628','Previous fit'),('selected','#2269ad','New calibration')]:
                t=tables[n];t=t[(t.instrument==instrument)&(t.expiry_days==days)].sort_values('strike')
                ax.plot(t.K_over_F,(t.model_price_U-t.price_mid_U)/t.price_scale_U,color=color,label=label)
            ax.axhline(0,color='.3',lw=.8);ax.grid(alpha=.2);ax.set(title=f'{instrument} {days} days',xlabel='Strike / market forward',ylabel='Canonical-call error / quote scale')
        for ax in axes.flat[len(dates):]:ax.axis('off')
        handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=2,frameon=False)
        fig.suptitle('Signed price errors: all estimates remain visible',fontsize=13);save(fig,instrument.lower()+'-price-errors')
    u=pd.read_csv(run/'uniform_budget_prices_selected.csv');rows=[]
    for (ins,day,n),g in u.groupby(['instrument','expiry_days','uniform_paths']):
        rows.append(dict(instrument=ins,days=day,paths=n,median_scaled_se=float(np.median(g.mc_se_U/g.price_scale_U)),
                         maximum_scaled_se=float(np.max(g.mc_se_U/g.price_scale_U))))
    b=pd.DataFrame(rows);b.to_csv(out/'sampling-precision.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for ax,ins in zip(axes,['SPX','VIX']):
        for day,g in b[b.instrument==ins].groupby('days'):ax.plot(g.paths,g.median_scaled_se,'o-',label=f'{day} days')
        ax.set(xscale='log',yscale='log',xlabel='Uniform nested path budget',ylabel='Median SE / quote scale',title=f'{ins}: observed precision at the same parameters')
        ax.set_xticks([512,1024,2048,4096],['512','1024','2048','4096']);ax.set_xticks([],minor=True);ax.grid(alpha=.2);ax.legend(frameon=False,ncols=2)
    save(fig,'sampling-precision')
    sources=[Path(__file__),run/'run.json',*[run/f'prices_canonical_{n}.csv' for n in tables],*[run/f'forwards_{n}.csv' for n in tables],run/'uniform_budget_prices_selected.csv']
    manifest=dict(source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs},
                  convention='All IVs invert calls with observed market forwards. Sources remain OTM only.',fitted_in_this_run=True,uncertainty='Pointwise sampling only; no bias/pilot/parameter uncertainty.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (out/'README.md').write_text('''# Maturity-curve calibration with 4096 paths

Twenty-four-parameter calibration; full-price-guarded local search on the calibration simulations. The original OTM
observations are converted to calls by parity. Blue uses the new selected
parameters; orange uses the previous fitted vector on the SAME fixed training samples.

- [All seven SPX maturities](spx-all-maturity-fit.png)
- [All four VIX maturities](vix-all-maturity-fit.png)
- [Actual VIX forwards and path allocation](forwards-and-path-budgets.png)
- [SPX signed price errors](spx-price-errors.png)
- [VIX signed price errors](vix-price-errors.png)
- [Sampling precision by path budget](sampling-precision.png)

All figures also have SVG exports. Bands show ±2 pointwise sampling SE, not
parameter uncertainty or a simultaneous confidence band. Invalid IVs are marked;
signed prices remain in the tables and loss. Both market and model forward
conventions are preserved in the run; a common call plot does not establish
that the model VIX forward is correct. Original-OTM fit errors are also retained.
''')

from .plot_maturity_curve_extra import extra

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--run',type=Path,required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();render(args.run.resolve(),args.output.resolve());extra(args.run.resolve(),args.output.resolve())
