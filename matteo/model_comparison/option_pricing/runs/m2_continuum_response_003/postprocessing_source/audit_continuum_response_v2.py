"""Independent powered-readout, conditional and integration diagnostics."""
import argparse
import ctypes as ct
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.special import roots_jacobi
from .continuum_response_study import context,tables,subset,PARENT,metadata
from .continuum_response_core import ContinuumFit,measure
from .forward_floor_core import ForwardFloorFit,unpack,LOG_HALF_MAX
from .vix_only_curve_fit import PRIOR,ARCHIVE,DAYS,ACTIVE,load_draws,DELTA
from .full_readout_fit import HERE,arr,dump
from .m2_spx import ROOT
from .postcap_vov_core import conditional_multiplier


def numpy_power(law,p,u):
    m,o,w=law;x=m[:,None,:]+o[None,:,:]
    eta=p[0]*x[:,:,0]+p[1]*(x[:,:,1]-p[2])**2+p[3]*(x[:,:,2]-p[4])**2
    h=p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6]
    return np.exp(2*(1-u)*h)@w


def run(run):
    out=run/'audit';out.mkdir(exist_ok=True);prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);initial=np.array(prior['initial_state']);sel=json.loads((run/'selection.json').read_text());x=np.array(sel['parameters']);bx=np.array(sel['boundary_parameters']);parent=np.array(sel['previous_parameters'])
    market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv');forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv');controls=json.loads((PRIOR/'vix_controls.json').read_text());draws=load_draws(run/'fresh_cache');ctx,k,dated=context(out,p,ends)
    coarse={'u_grid':np.linspace(0,1,33)};probes={};ids_by_day={};oracles=[];interpolation=[]
    for day,a in draws.items():
        with np.load(run/'fresh_power_tables'/f'd{day}.npz') as z:
            for key,v in z.items():coarse[f'd{day}_{key}']=v
        V=np.load(run/'fresh_check/paths_continuous.npz')[f'd{day}_V'];w=a['weights'];K=market[market.expiry_days.eq(day)].strike.max()
        ids=np.unique(np.r_[0,1,np.argmax(abs(a['states']),axis=0),np.argmin(a['z_normals']),np.argmax(a['z_normals']),np.argsort(w*np.maximum(V-K,0))[-4:],np.argsort(V)[:3]])
        ids_by_day[day]=ids;probes[day]={key:v[ids] for key,v in a.items()}
        sp=CubicSpline(coarse['u_grid'],coarse[f'd{day}_log_ratios'][:,:,ids],axis=0)
        times=coarse[f'd{day}_times'];stateids=np.arange(min(4,len(ids)))
        for it in [0,len(times)//2,len(times)-1]:
            t=times[it];law=ctx.law(float(t-day/365),a['states'][ids],32);nlaw=ctx.law(float(t),initial,48)
            for u in [.017,.1234,.555,.879,.998]:
                k.lib.configure(0,1-u);actual=k.expectation(*law,p)[:,0]/k.expectation(*nlaw,p)[0,0];pred=np.exp(sp(u)[it]);err=float(np.max(abs(pred/actual-1)));interpolation.append(dict(day=day,time_index=it,u=u,relative_error=err))
                # NumPy independently implements the original soft cap and power.
                if it==0:
                    small=(arr(law[0][stateids]),law[1],law[2]);val=numpy_power(small,p,u);cpp=k.expectation(*small,p)[:,0];oracles.append(dict(day=day,u=u,relative_error=float(max(abs(cpp/val-1)))))
    dump(out/'powered_numpy_oracles.json',oracles);dump(out/'off_grid_power_checks.json',interpolation)
    assert max(r['relative_error'] for r in oracles)<2e-11
    assert max(r['relative_error'] for r in interpolation)<2e-5
    # Rebuild every powered conditional expectation at doubled u/time resolution,
    # using original-state GH48, on boundary and largest payoff-contributing states.
    # Full combined refinement is costly under direct GH48, so use three
    # deliberately stressed states/expiry; the larger state set above checks u
    # interpolation at several off-grid powers separately.
    for day,a in draws.items():
        paths=np.load(run/'fresh_check/paths_continuous.npz');V=paths[f'd{day}_V'];K=market[market.expiry_days.eq(day)].strike.max()
        ids=np.unique([np.argmax(a['weights']*np.maximum(V-K,0)),np.argmin(V),np.argmax(np.linalg.norm(a['states'],axis=1))]);ids_by_day[day]=ids;probes[day]={key:v[ids] for key,v in a.items()}
    fine=tables(ctx,k,p,initial,probes,out/'refined_power_tables',grid_size=65,time_order=16,order=48)
    cn={'u_grid':coarse['u_grid']}
    for day,ids in ids_by_day.items():
        for key in ['times','weights']:cn[f'd{day}_{key}']=coarse[f'd{day}_{key}']
        cn[f'd{day}_ratios']=coarse[f'd{day}_ratios'][:,ids];cn[f'd{day}_log_ratios']=coarse[f'd{day}_log_ratios'][:,:,ids]
    ce=ContinuumFit(cn,probes,market,forwards,controls,p,ends,measure_order=24);fe=ContinuumFit(fine,probes,market,forwards,controls,p,ends,measure_order=48)
    c=ce.predict(x,True);f=fe.predict(x,True);refine=[]
    for day in DAYS:
        delta=f['raw'][day]['V']-c['raw'][day]['V'];refine.append(dict(day=day,probe_count=len(ids_by_day[day]),max_VIX_difference=float(max(abs(delta))),ids=ids_by_day[day].tolist(),coarse='u33/measure24/GH32/time8',fine='u65/measure48/GH48/time16'))
    dump(out/'quadrature_refinement.json',refine)
    # Separate independent physical-parameter integration: positive Gauss-Jacobi
    # sum, log-power splines, OU conditional formula and variance time integration.
    m,kap=measure(x);zz,ww=roots_jacobi(24,(1-m)*kap-1,m*kap-1);uu=(zz+1)/2;ww/=ww.sum();levels,floor,theta=unpack(x,p);b=theta[0];rate=np.log(2)*365/np.exp(theta[2]);replays=[];rawputs=[]
    for day,a in draws.items():
        times=coarse[f'd{day}_times'];tw=coarse[f'd{day}_weights'];lr=coarse[f'd{day}_log_ratios'];sp=CubicSpline(coarse['u_grid'],lr,axis=0);T=day/365;z=np.sqrt(-np.expm1(-2*rate*T))*a['z_normals'];integ=np.full(len(z),floor*DELTA)
        for u,wt in zip(uu,ww):
            power=np.exp(sp(u));L=np.array([conditional_multiplier(z,T,float(t-T),b*u,rate) for t in times]);seg=np.minimum(np.searchsorted(ends,times,side='right'),len(ends)-1);integ+=wt*((levels[seg]-floor)*tw)@(power*L)
        V=100*np.sqrt(integ/DELTA);saved=np.load(run/'fresh_check/paths_continuous.npz');err=float(max(abs(V-saved[f'd{day}_V'])));rel=float(max(abs(V-saved[f'd{day}_V'])/np.maximum(1,abs(saved[f'd{day}_V']))));assert err<1e-9 or rel<1e-11;replays.append(dict(day=day,max_VIX_difference=err,max_relative_difference=rel,weighted_mean_difference=float(abs(np.mean(a['weights']*(V-saved[f'd{day}_V']))))))
        for q in market[market.expiry_days.eq(day)&market.option_type.eq('put')].itertuples():
            pay=a['weights']*np.maximum(q.strike-V,0);block=np.array([pay[a['nets']==j].mean() for j in range(8)]);rawputs.append(dict(day=day,strike=float(q.strike),raw_price=float(pay.mean()),net_SE=float(block.std(ddof=1)/np.sqrt(8)),positive_payoffs=int(np.sum(pay>0))))
    pd.DataFrame(rawputs).to_csv(out/'raw_native_puts.csv',index=False);dump(out/'independent_replay.json',replays)
    # Regenerate every fresh importance scramble from the frozen design.
    from .mixed_response_study import ZSHIFTS,ZPROB,seed_for
    from esn_eq.importance_sampling import DefensiveGaussianMixture
    from esn_eq.outer_sampling import GaussianOuterSampler
    des=json.loads((run/'fresh_cache/design.json').read_text());stream_checks=[]
    for day,a in draws.items():
        dec=des['old_proposals'][str(day)];shifts=np.array([np.r_[old,z] for old in dec['shifts'] for z in ZSHIFTS]);prob=np.array([p0*q for p0 in dec['probabilities'] for q in ZPROB]);mix=DefensiveGaussianMixture(shifts,prob);sampler=GaussianOuterSampler(dated,day/365,initial)
        for net in range(8):
            ids=np.flatnonzero(a['nets']==net);zz,ww,cc=mix.draw(11,seed_for(des['seed_base'],day,net));state=sampler.mean+zz[:,:-1]@sampler.root.T
            assert np.array_equal(zz[:,:-1],a['normals'][ids]);assert np.array_equal(zz[:,-1],a['z_normals'][ids]);assert np.array_equal(ww,a['weights'][ids]);assert np.array_equal(cc,a['components'][ids]);assert np.array_equal(state,a['states'][ids]);stream_checks.append(dict(day=day,net=net,pass_exact_redraw=True))
    dump(out/'fresh_scramble_redraw.json',stream_checks)

    # Source/input integrity, including unrelated captured sources unchanged.
    design=json.loads((run/'design.json').read_text());source=[]
    for path,h in design['source_sha256'].items():
        target=ROOT/path;source.append(dict(path=path,pass_hash=target.exists() and hashlib.sha256(target.read_bytes()).hexdigest()==h))
    inputs=[]
    for path,h in design['input_sha256'].items():inputs.append(dict(path=path,pass_hash=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h))
    assert all(t['pass_hash'] for t in source+inputs)
    levels,floor,_=unpack(x,p);assert np.all(levels>=floor);assert np.array_equal(levels[np.setdiff1d(np.arange(16),ACTIVE)],p[7:-1][np.setdiff1d(np.arange(16),ACTIVE)])
    quotes=pd.read_csv(run/'fresh_check/prices_otm_continuous.csv');negative=quotes[quotes.model_price_U<0].copy();negative['negative_in_SE']=-negative.model_price_U/negative.mc_se_U;negative.to_csv(out/'negative_native_prices.csv',index=False)
    dump(out/'audit.json',dict(status='PASS',source_checks=source,input_checks=inputs,original_readout_unchanged=True,beta_J=float(p[4]),floor_support=True,inactive_nine_unchanged=True,negative_price_count=len(negative),powered_numpy_max_error=max(r['relative_error'] for r in oracles),power_interpolation_max_relative_error=max(r['relative_error'] for r in interpolation),refinements=refine,limitations=['Sampled deterministic refinements, no uniform certificate.','No new Lean theorem or exact population forward certificate.']))
    print('AUDIT PASS',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);run(ap.parse_args().run.resolve())
