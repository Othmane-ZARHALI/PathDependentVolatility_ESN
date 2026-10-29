"""Read-only diagnostics for the continuous-model reduction and 27-day wing.

No calibration, model change, clipping or alteration of archived run artifacts.
"""
import argparse,hashlib,importlib.util,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import betainc
from .full_readout_fit import HERE,dump
from .plot_maturity_curve import iv
from .plot_joint_curve_fit import implied_vol

RUN=HERE/'runs/m2_continuum_response_003'
CERT=HERE/'runs/m2_vix_design_research_001/moment_certificates.json'


def run(out):
    out.mkdir(parents=True,exist_ok=True);scores=json.loads((RUN/'fresh_check/scores.json').read_text());summaries=[];selected_rows=[];ou=[];all_prices=[]
    sources=[Path(__file__),CERT,HERE/'runs/m2_vix_design_research_001/moment_certificate.py',RUN/'selection.json',RUN/'curve_levels.csv',RUN/'fresh_check/scores.json',RUN/'audit/audit.json',RUN/'audit/quadrature_refinement.json',RUN/'training/forwards_continuous.csv']
    for z in scores:
        name=z['name'];path=RUN/'fresh_check'/f'prices_canonical_{name}.csv';sources.append(path);q=pd.read_csv(path);q['market_IV']=iv(q,q.price_mid_U);q['model_IV_market_F']=iv(q,q.model_price_U)
        mf=dict(zip([27,55,118,153],pd.read_csv(RUN/'fresh_check'/f'forwards_{name}.csv').model_forward));sources.append(RUN/'fresh_check'/f'forwards_{name}.csv')
        q['model_IV_model_F']=[implied_vol(float(r.model_price_U),mf[int(r.expiry_days)],float(r.strike),float(r.T_years),'call') if r.model_price_U>0 else np.nan for r in q.itertuples()]
        q['candidate']=name;q['price_error']=q.model_price_U-q.price_mid_U;q['outside_bidask_error_in_SE']=np.where(q.model_price_U<q.price_bid_U,q.model_price_U-q.price_bid_U,np.where(q.model_price_U>q.price_ask_U,q.model_price_U-q.price_ask_U,0.))/q.mc_se_U;q['wing']=q.K_over_F.gt(2);all_prices.append(q)
        for day,g in q.groupby('expiry_days'):
            w=g[g.wing];x=np.log(w.K_over_F);mk=w.market_IV.to_numpy();md=w.model_IV_market_F.to_numpy();shape=np.polyfit(x,mk,2);modelshape=np.polyfit(x,md,2)
            sources.append(RUN/'fresh_check'/f'prices_otm_{name}.csv');native=pd.read_csv(RUN/'fresh_check'/f'prices_otm_{name}.csv').query('expiry_days == @day');outside=np.where(native.model_price_U<native.price_bid_U,native.model_price_U-native.price_bid_U,np.where(native.model_price_U>native.price_ask_U,native.model_price_U-native.price_ask_U,0.))/native.price_scale_U
            summaries.append(dict(candidate=name,day=int(day),wing_quotes=len(w),market_slope=float(np.polyfit(x,mk,1)[0]),model_slope=float(np.polyfit(x,md,1)[0]),market_curvature=float(2*shape[0]),model_curvature=float(2*modelshape[0]),wing_IV_MAE_pp=float(np.mean(abs(md-mk))),max_wing_IV_forward_convention_difference_pp=float(np.max(abs(w.model_IV_model_F-w.model_IV_market_F))),far_call_price_gap=float(w.iloc[-1].price_error),far_call_gap_in_SE=float(w.iloc[-1].price_error/w.iloc[-1].mc_se_U),native_option_loss_contribution=float(np.mean(outside**2)/4)))
            if day==27:
                ids=sorted(set([int(np.argmin(abs(g.strike-20))),int(np.argmin(abs(g.strike-40))),int(np.argmin(abs(g.strike-60))),len(g)-1]))
                for j in ids:selected_rows.append(g.iloc[j].to_dict())
        if 'kappa' in z['model']:
            m=z['model'];gamma=m['b'];a=np.log(2)/m['half_life_days']
            for d in [27,55,118,153]:
                T=d/365;qT=-np.expm1(-2*a*d);load=gamma*np.sqrt(qT);mid=load*np.exp(-a*15)
                ou.append(dict(candidate=name,day=d,q_T=float(qT),standard_normal_stress_loading_at_expiry=float(load),standard_normal_stress_loading_at_mid_window=float(mid),mid_window_pure_stress_logVIX_loading=float(mid/2),mid_window_single_time_pure_stress_Black_IV_percent=float(100*mid/2/np.sqrt(T))))
    table=pd.concat(all_prices,ignore_index=True);table.to_csv(out/'quote_diagnostics.csv',index=False);pd.DataFrame(summaries).to_csv(out/'maturity_shape_metrics.csv',index=False);pd.DataFrame(selected_rows).to_csv(out/'selected_27day_quotes.csv',index=False);pd.DataFrame(ou).to_csv(out/'OU_clock_diagnostics.csv',index=False)
    common=[]
    for (name,day),g in table.groupby(['candidate','expiry_days']):
        w=g[(g.K_over_F>2)&(g.K_over_F<=100/20.1)]
        common.append(dict(candidate=name,day=int(day),quotes=len(w),market_slope=float(np.polyfit(np.log(w.K_over_F),w.market_IV,1)[0]),model_slope=float(np.polyfit(np.log(w.K_over_F),w.model_IV_market_F,1)[0])))
    pd.DataFrame(common).to_csv(out/'common_range_wing_slopes.csv',index=False)
    # Payoff-weighted decomposition of the matching two-atom numerical reference.
    # This is an exact decomposition for that reference, not an attribution of
    # the continuous model's small paired-difference estimator.
    model=json.loads((RUN/'experimental_model_continuous.json').read_text());pp=np.array(model['parameters']);ends=np.array(model['ends_years']);levels=pp[7:-1];floor=pp[-1];m=model['readout'];b=m['b'];w=m['measure_mean'];rate=np.log(2)*365/m['half_life_days'];tail=[]
    sources.append(RUN/'experimental_model_continuous.json')
    for day in [27,55,118,153]:
        path=HERE/'runs/m2_continuum_response_002/fresh_power_tables'/f'd{day}.npz';drawpath=HERE/'runs/m2_continuum_response_002/fresh_cache/production'/f'vix_d{day}.npz';sources.extend([path,drawpath])
        with np.load(path) as z:times=z['times'];tw=z['weights'];old=z['ratios']
        with np.load(drawpath) as z:eps=z['z_normals'];iw=z['weights']
        T=day/365;load=b*np.exp(-rate*(times-T))*np.sqrt(-np.expm1(-2*rate*T));L=np.exp(load[:,None]*eps[None,:]-.5*load[:,None]**2);seg=np.minimum(np.searchsorted(ends,times,side='right'),len(ends)-1);cw=(levels[seg]-floor)*tw
        stress=w*(cw@L);base=floor*(30/365)+(1-w)*(cw@old);H=base+stress;V=100*np.sqrt(H/(30/365));share=stress/H;elastic=.5*w*((cw*load)@L)/H
        maxK=table.query('candidate == "continuous" and expiry_days == @day').strike.max()
        for K in sorted(set([20.,40.,60.,float(maxK)])):
            pay=iw*np.maximum(V-K,0);total=pay.sum();assert total>0
            tail.append(dict(day=day,strike=K,reference_raw_call=float(pay.mean()),payoff_weighted_stress_variance_share=float(pay@share/total),payoff_weighted_logVIX_elasticity_to_standard_normal=float(pay@elastic/total),payoff_fraction_from_states_stress_share_above_90pct=float(pay[share>.9].sum()/total),payoff_ESS=float(total*total/(pay@pay)),top10_payoff_fraction=float(np.sort(pay)[-10:].sum()/total)))
    pd.DataFrame(tail).to_csv(out/'payoff_weighted_reference_components.csv',index=False)

    # Reverify archived certificates in exact rational arithmetic; do not rerun LP.
    spec=importlib.util.spec_from_file_location('cert_check',HERE/'runs/m2_vix_design_research_001/moment_certificate.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);cert=json.loads(CERT.read_text());cc=[mod.verify(z) for z in cert['records']];assert hashlib.sha256(Path(cert['source']).read_bytes()).hexdigest()==cert['source_sha256']
    c27=next(z for z in cert['records'] if z['day']==27 and z['forward_tolerance']==0);selected=next(z for z in scores if z['name']=='continuous');E2=next(z for z in selected['moments'] if z['day']==27)['exact_second_moment'];F=c27['forward'];needed=c27['certified_lower_bound_float'];forward_slope=float(mod.Q(c27['native_forward_coefficient']));fcheck=pd.read_csv(RUN/'fresh_check/forwards_continuous.csv').query('expiry_days == 27').iloc[0]
    budget=dict(model_second_moment=E2,necessary_second_moment_at_exact_market_forward=needed,minimum_shortfall=needed-E2,relative_second_moment_shortfall=(needed-E2)/E2,model_variance_at_exact_market_forward=E2-F**2,necessary_variance_at_exact_market_forward=needed-F**2,required_minimum_average_xi_increase=(needed-E2)/10000,certificate_forward_coefficient=forward_slope,necessary_bound_if_mean_equal_checking_forward=needed+forward_slope*(fcheck.model_forward-F),qualification='A necessary condition for all native bid/ask intervals and the stated mean, not proof of a unique wing mechanism or global family infeasibility.')
    # Curve allocation within the first VIX window; numerical zero is displayed,
    # never clipped in prices or used to alter the fitted model.
    curve=pd.read_csv(RUN/'curve_levels.csv');curve['continuous']=np.array(model['parameters'])[7:-1];floor=selected['model']['floor_variance'];assert np.all(curve.continuous>=floor);length=np.maximum(0,np.minimum(curve.end_day,57)-np.maximum(curve.start_day,27));c=curve[length>0].copy();c['overlap_days']=length[length>0];c['expected_nonfloor_variance_integral']=c.overlap_days*(c.continuous-floor);c['expected_nonfloor_variance_share']=c.expected_nonfloor_variance_integral/c.expected_nonfloor_variance_integral.sum();c.to_csv(out/'first_window_curve_allocation.csv',index=False)
    dates=pd.to_datetime(table.expiry_date)-pd.to_timedelta(table.expiry_days,unit='D');assert dates.nunique()==1;date=dates.iloc[0].date().isoformat()
    measures=[]
    for z in scores:
        m=z['model']
        if 'kappa' in m:measures.append(dict(candidate=z['name'],m=m['measure_mean'],kappa=m['kappa'],interior_mass_01_to_09=float(betainc(m['alpha'],m['beta'],.9)-betainc(m['alpha'],m['beta'],.1))))
    # Exact single-future-time loading calculation, not a calibration or full
    # smile prediction. Check its stationary point against a numerical grid.
    mm=selected['model'];T0=27.;lag=15.;astar=np.log1p(T0/lag)/(2*T0);hstar=np.log(2)/astar
    load=lambda days,half:mm['b']*np.sqrt(-np.expm1(-2*np.log(2)*days/half))*np.exp(-np.log(2)*lag/half)
    grid=np.geomspace(1,365,10001);gg=load(T0,grid);assert abs(np.log(grid[np.argmax(gg)]/hstar))<.001
    proxy=dict(qualification='Single future-time stress loading at lag15days, fixed b. Not a full VIX or implied-volatility recalibration.',half_life_current_days=mm['half_life_days'],half_life_maximizing_27day_loading_days=float(hstar),loading_27_current=float(load(27,mm['half_life_days'])),loading_27_at_proxy_optimum=float(load(27,hstar)),loading_153_current=float(load(153,mm['half_life_days'])),loading_153_at_27day_proxy_optimum=float(load(153,hstar)))

    summary=dict(status='completed diagnostics only',market_snapshot_date_inferred_from_expiry_dates=date,experiment_date='2026-10-27',scope='All217 quotes retained; wing means K/F>2. Common-range slopes also retained. No calibration or parameter changes.',second_moment_budget=budget,midwindow_loading_proxy=proxy,shape_metrics=summaries,common_range_slopes=common,payoff_weighted_reference_components=tail,OU_clock=ou,measures=measures,exact_certificate_rechecks=cc,limits=['Black IV slope is measured in percentage points per log(K/F), curvature per log(K/F)^2.','Smaller q(T) alone does not imply a flatter annualized smile; clock table includes annualization.','Current first27day window is[27,57], with active level[27,29] at its floor; this is not a single variance level.','Conditional quadrature refinements cover sampled states, not uniform error bounds.'])
    dump(out/'diagnostics.json',summary);dump(out/'input_and_source_hashes.json',{str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(13,4.8),layout='constrained');q=table.query('candidate == "continuous" and expiry_days == 27')
    axes[0].plot(q.K_over_F,q.market_IV,'k-',label='Market');axes[0].plot(q.K_over_F,q.model_IV_market_F,color='#d04a35',label='Continuous fit');axes[0].set(xlabel='Strike / market forward',ylabel='Black IV (%)',title='27-day smile crosses the market');axes[0].set_xscale('log')
    calls=q[q.wing];axes[1].plot(calls.strike,calls.price_mid_U,'k.-',label='Market call');axes[1].errorbar(calls.strike,calls.model_price_U,yerr=2*calls.mc_se_U,color='#d04a35',marker='.',capsize=2,label='Model call ±2 sampling SE');axes[1].fill_between(calls.strike,calls.price_bid_U,calls.price_ask_U,color='.2',alpha=.15,label='Market bid–ask');axes[1].set(xlabel='Strike (VIX points)',ylabel='Undiscounted call price',yscale='log',title='Far-right price shortfall persists before IV inversion')
    for ax in axes:ax.grid(alpha=.2);ax.legend(frameon=False)
    for ext in ['png','svg']:fig.savefig(out/f'27day_shape_and_prices.{ext}',dpi=160,bbox_inches='tight')
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,4.6),layout='constrained')
    for day in [27,55,118,153]:ax.plot(grid,load(day,grid),label=f'{day}-day expiry')
    ax.axvline(mm['half_life_days'],ls='--',color='.25',label='Fitted shared half-life');ax.axvline(hstar,ls=':',color='#d04a35',label='27-day midpoint-loading optimum');ax.set(xscale='log',xlabel='Stress half-life (calendar days)',ylabel='Stress loading on a standard normal',title='One shared persistence creates a maturity trade-off\nSingle future date at expiry +15days; fixed b; not a smile calibration');ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=9)
    for ext in ['png','svg']:fig.savefig(out/f'OU_midwindow_tradeoff.{ext}',dpi=160,bbox_inches='tight')
    plt.close(fig);print(json.dumps(dict(date=date,budget=budget,continuous_metrics=[z for z in summaries if z['candidate']=='continuous']),indent=2),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);run(ap.parse_args().output.resolve())
