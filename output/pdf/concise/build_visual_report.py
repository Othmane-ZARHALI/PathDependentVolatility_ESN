"""Replot archived prices only: no simulation, fitting, smoothing or clipping."""
from pathlib import Path
import hashlib, json, os, re, zipfile
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/m2-report-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import ndtr
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
RESEARCH=ROOT/'matteo/model_comparison/option_pricing'
RUN=RESEARCH/'runs'
FIG=HERE/'figures'; FIG.mkdir(parents=True,exist_ok=True)
INPUTS={}; POINTS=[]; DIAGNOSTICS=[]
COLORS=['#3a6fa5','#d66a2e','#32936f','#9a59a6','#c69a26']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':12,
 'axes.labelsize':11,'legend.fontsize':10.5,'axes.spines.top':False,'axes.spines.right':False,
 'pdf.fonttype':42,'axes.grid':True,'grid.alpha':.16})

def read(p):
    p=Path(p); INPUTS[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    return pd.read_csv(p)

def iv(p,f,k,t,side):
    if not np.isfinite(p) or min(f,k,t)<=0:return np.nan
    if side=='put':p=p+f-k
    intrinsic=max(f-k,0.)
    if p<intrinsic or p>=f:return np.nan
    if p==intrinsic:return 0.
    def err(s):
        z=s*np.sqrt(t); d1=np.log(f/k)/z+z/2
        return f*ndtr(d1)-k*ndtr(d1-z)-p
    try:return 100*brentq(err,1e-9,50,xtol=1e-11)
    except ValueError:return np.nan

def load(prices,forwards=None):
    q=read(prices).copy()
    if 'original_option_type' in q and (q.option_type=='call').all():
        # Prefer explicit native files, but support stored canonical tables.
        side=q.original_option_type
        if 'native_otm_model_price_U' in q:q['model_price_U']=q.native_otm_model_price_U
        else:raise ValueError('Native price needed')
        for c in ['price_bid_U','price_ask_U','price_mid_U']:
            q[c]=q[c]-(side=='put')*(q.forward-q.strike)
        q['option_type']=side
    q['mf']=q.forward
    if 'model_forward' in q:q['mf']=q.model_forward
    if forwards:
        f=read(forwards).set_index('expiry_days').model_forward
        mask=q.instrument.eq('VIX');q.loc[mask,'mf']=q.loc[mask,'expiry_days'].map(f)
    for c,out in [('price_mid_U','market'),('price_bid_U','bid'),('price_ask_U','ask')]:
        q[out]=[iv(p,f,k,t,s) for p,f,k,t,s in zip(q[c],q.forward,q.strike,q.T_years,q.option_type)]
    q['model']=[iv(p,f,k,t,s) for p,f,k,t,s in zip(q.model_price_U,q.mf,q.strike,q.T_years,q.option_type)]
    return q

def native(run,folder,name):
    p=RUN/run/folder
    return load(p/f'prices_otm_{name}.csv',p/f'forwards_{name}.csv')

def modern(run,sp,vi,vf):
    return pd.concat([load(RUN/run/sp),load(RUN/run/vi,RUN/run/vf)],ignore_index=True)

def plate(name,variants,note,vix_only=False,delta=False):
    """Native IVs use each model's forward; all x axes use market forwards."""
    rows=1 if vix_only else 2
    fig,axs=plt.subplots(rows,2,figsize=(10.7,3.7 if vix_only else 7.0),squeeze=False)
    for row,ins in enumerate(['VIX'] if vix_only else ['SPX','VIX']):
        days=[29,211] if ins=='SPX' else [27,153]
        for ax,day in zip(axs[row],days):
            base=variants[0][1];b=base[(base.instrument==ins)&(base.expiry_days==day)].copy()
            if ins=='SPX' and 'priority_selected' in b:b=b[b.priority_selected==True]
            b=b.sort_values(['strike','option_type'])
            # Distinct source calls/puts can share a strike: do not average them.
            if delta:
                ax.axhline(0,color='.4',ls='--',lw=.8)
            else:
                for side,g in b.groupby('option_type',sort=True):
                    ax.fill_between(g.K_over_F,g.bid,g.ask,color='.25',alpha=.13,lw=0)
                    ax.plot(g.K_over_F,g.market,'o',ms=2.6,color='#20252a',alpha=.7,
                            label='Market mid / bid-ask' if side==b.option_type.iloc[0] else None)
            for j,(label,q) in enumerate(variants):
                g=q[(q.instrument==ins)&(q.expiry_days==day)].copy()
                if ins=='SPX' and 'priority_selected' in g:g=g[g.priority_selected==True]
                g=g.sort_values(['strike','option_type'])
                assert len(g)>0,(name,label,ins,day)
                assert np.all(g.mf>0)
                color=COLORS[j];bad=g.model.isna()
                if delta:
                    base_iv=b.set_index(['strike','option_type']).model
                    g['display']=g.model-np.array([base_iv.loc[(k,s)] for k,s in zip(g.strike,g.option_type)])
                else:g['display']=g.model
                for k,(side,z) in enumerate(g.groupby('option_type',sort=True)):
                    ax.plot(z.K_over_F,z.display,color=color,lw=1.45,marker='.',ms=2.3,
                            label=label if k==0 else None)
                if bad.any():ax.plot(g.loc[bad,'K_over_F'],np.full(bad.sum(),.025+j*.025),
                     '|',color=color,ms=6,transform=ax.get_xaxis_transform())
                if 'mc_se_U' in g:se=g.mc_se_U
                elif 'combined_SE_U' in g:se=g.combined_SE_U
                else:se=pd.Series(np.nan,index=g.index)
                err=np.where(g.model_price_U<g.price_bid_U,g.model_price_U-g.price_bid_U,
                    np.where(g.model_price_U>g.price_ask_U,g.model_price_U-g.price_ask_U,0))/g.price_scale_U
                DIAGNOSTICS.append(dict(figure=name,variant=label,instrument=ins,days=day,
                    quotes=len(g),invalid_iv=int(bad.sum()),interval_loss=float(np.mean(err**2)),
                    native_iv_mae_pp=float(np.nanmean(abs(g.model-g.market))),
                    max_forward_gap=float(np.max(abs(g.mf-g.forward))),
                    median_sampling_SE_quote_scales=float(np.nanmedian(se/g.price_scale_U))
                      if np.isfinite(se).any() else None))
                z=g[['instrument','expiry_days','strike','option_type','K_over_F','forward','mf',
                    'model_price_U','price_mid_U','model','market','bid','ask']].copy()
                z['figure']=name;z['variant']=label;POINTS.append(z)
            ax.set_title(f'{ins}: {day} days',loc='left',fontweight='bold')
            ax.set_xlabel('Strike / market forward')
            ax.set_ylabel('IV change from baseline (pp)' if delta else 'Implied volatility (%)')
            if ins=='VIX':
                from matplotlib.ticker import NullFormatter
                ax.set_xscale('log')
                lo,hi=ax.get_xlim(); ticks=[x for x in [1,2,4,8] if lo<=x<=hi]
                ax.set_xticks(ticks);ax.set_xticklabels([str(x) for x in ticks])
                ax.xaxis.set_minor_formatter(NullFormatter());ax.set_xlim(lo,hi)
            ax.margins(x=.025)
    h,l=axs[0,0].get_legend_handles_labels()
    # Explicit handles avoid call/put ordering affecting legend membership.
    from matplotlib.lines import Line2D
    h=[Line2D([],[],color='#20252a',marker='o',lw=0,ms=4)]+[
        Line2D([],[],color=COLORS[j],lw=1.8) for j in range(len(variants))]
    l=['Market mid / bid-ask']+[v[0] for v in variants]
    if delta:h=h[1:];l=l[1:]
    fig.legend(h,l,loc='lower center',ncols=2 if len(h)>3 else 3,frameon=False,
               bbox_to_anchor=(.5,.005))
    fig.subplots_adjust(left=.075,right=.985,top=.96,bottom=.23 if rows==2 else .3,
                        hspace=.42,wspace=.25)
    fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(FIG/f'{name}.png',dpi=160,bbox_inches='tight')
    plt.close(fig)

def main():
    old=native('m2_maturity_curve_fit_001','','baseline')
    baseline=native('m2_maturity_curve_fit_001','','selected')
    cap=native('m2_cap_controlled_001','full','shape03')
    tail=native('m2_cap_controlled_001','full','shape02')
    plate('01-curve-cap',[('Before aligned curve',old),('16-piece curve',baseline),
          ('Higher cap + window refit',cap),('Slow-growth cap + window refit',tail)],'Training')
    inside=native('m2_stochastic_vov_001','full','cap0_z2_h90')
    after=native('m2_postcap_vov_001','full','g1_h7')
    plate('02-stochastic-response',[('Original M2',baseline),('Inside-cap sensitivity: 2x / 90d',inside),
          ('After-cap multiplier: gamma=1 / 7d',after)],'Training')
    single=native('m2_two_speed_001','fresh_check','single_speed')
    two=native('m2_two_speed_001','fresh_check','two_speed')
    plate('03-two-speeds',[('One stress speed',single),('Shared fast + slow stress',two)],'Checking',True)
    def joint(n):
        return modern('m2_two_speed_joint_001',f'checking/prices_spx_{n}.csv',
           f'checking/prices_vix_{n}.csv',f'checking/forwards_{n}.csv')
    plate('04-joint-tradeoff',[('VIX-only parameters',joint('vix_only')),
          ('Joint curve / floor fit',joint('curve_floor')),('Joint 23-parameter fit',joint('joint'))],'Checking')
    def protected(n):
        return modern('m2_protected_curve_001',f'fresh_check/prices_{n}_128.csv',
             'VIX_protection/prices_otm_two_speed.csv','VIX_protection/forwards_two_speed.csv')
    plate('05-protected-curve',[('Before off-window refit',protected('baseline')),
          ('Nine off-window levels refitted',protected('selected'))],'Fresh SPX; VIX inherited')
    def lev(t):
        s=load(RUN/f'm2_leverage_scan_001/fresh/prices_mixed_{t}_128.csv')
        return pd.concat([s,two],ignore_index=True)
    plate('06-negative-correlation',[('Independent stress / stock',lev('0')),
          ('Negative stress leverage: theta=0.8',lev('0.8'))],'Fresh SPX; VIX inherited')
    # Show global and conditional normalization in their common VIX states;
    # SPX examples below share the global experiment's original factor draws.
    def asym(c):
        return modern('m2_asymmetric_readout_001',f'fresh/SPX_prices_c{c}_128.csv',
             f'fresh_VIX/VIX_prices_c{c}.csv',f'fresh_VIX/VIX_forwards_c{c}.csv')
    plate('07-asymmetric-readout',[('Alpha=0: curve-refit control',asym(1)),
          ('Alpha=0.5 + curve refit',asym(3)),('Alpha=1.5 + curve refit',asym(5))],'Fresh draws, shared anchor')
    def slow(c):
        return modern('m2_slow_readout_002',f'fresh_validation/SPX/SPX_prices_c{c}_256.csv',
             f'fresh_validation/VIX/VIX_prices_c{c}.csv',f'fresh_validation/VIX/VIX_forwards_c{c}.csv')
    plate('08-slower-signal',[('Baseline',slow(0)),('Weak 9.48d signal + refit',slow(2)),
          ('Stronger 3.33d signal + refit',slow(4))],'New factor draws, reused anchor',delta=True)
    pd.concat(POINTS,ignore_index=True).to_csv(HERE/'plotted_points.csv',index=False)
    pd.DataFrame(DIAGNOSTICS).to_csv(HERE/'panel_diagnostics.csv',index=False)
    (HERE/'figure_provenance.json').write_text(json.dumps(dict(
        inputs_sha256=INPUTS,maturities=dict(SPX=[29,211],VIX=[27,153]),
        iv_convention='Native option IV: quoted forward for market, model forward for model; SPX theoretical forward equals quoted forward',
        difference_figure='08-slower-signal: native IV minus paired baseline native IV',
        prices_clipped=False,prices_smoothed=False,simulations_rerun=False,
        missing_iv_policy='NaN line gaps and colored ticks at panel foot'),indent=2)+'\n')
    print(f'Built 8 plates from {len(INPUTS)} archived price/forward files.')

if __name__=='__main__':main()
