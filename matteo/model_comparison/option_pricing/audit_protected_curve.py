"""Independent fresh-kernel replay, Gaussian law checks and VIX protection audit."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from .full_readout_fit import HERE,arr,dump
from .vix_only_curve_fit import PRIOR,ACTIVE
from .two_speed_joint_core import physical,vix_x
from .two_speed_precision import FREE
from .protected_curve_check import FreshKernel,context
from .audit_two_speed_precision import normal


def audit(run):
    fresh=run/'fresh_check';design=json.loads((run/'design.json').read_text());base=arr(design['baseline']);chosen=arr(json.loads((run/'selection.json').read_text())['parameters']);models=arr([base,chosen]);checked=0
    for folder in [run,fresh]:
        hashes=json.loads((folder/'design.json').read_text())['source_hashes']
        for name,h in hashes.items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==h,name;checked+=1
    np.testing.assert_array_equal(vix_x(base),vix_x(chosen));fixed=np.setdiff1d(np.arange(23),FREE);np.testing.assert_array_equal(base[fixed],chosen[fixed]);assert np.all((chosen[FREE]>=0)&(chosen[FREE]<=1));prior=json.loads((PRIOR/'run.json').read_text());p=arr(prior['selected']);assert p[4]>0 and p[3]>0
    derivative=json.loads((run/'derivative_checks.json').read_text());assert len(derivative)==9 and {a['index'] for a in derivative}==set(FREE) and max(a['relative_error'] for a in derivative)<3e-5
    full=dict(np.load(fresh/'grid_and_transitions.npz'));k=FreshKernel(fresh/'fresh.dylib');d={name:v for name,v in full.items() if name not in ['reference_grid','reference_segments']};d['expiries']=int(d['expiries']);nt=17;pairs=3;probe=d.copy()
    for name in ['transition','time','dt','seg','coarse','reference','mean','lognorm']:probe[name]=np.ascontiguousarray(d[name][:nt])
    probe['expiries']=2;probe['expiry']=np.full(nt,-1,np.int32);probe['expiry'][7]=0;probe['expiry'][-1]=1;probe['reference']=((np.arange(nt)%8)==0).astype('int32');seed=2026102859;actual=k.simulate(probe,pairs,p,models,seed);expected=np.zeros_like(actual);nd=d['projection'].shape[1];state=np.zeros((pairs,nd));zs=np.zeros((pairs,2));logs=np.zeros((2,pairs,2,2));logref=np.zeros((2,pairs,2));iv=np.zeros_like(logref);hold=np.zeros_like(logref);qr=np.ones((2,pairs));hf=365**base[21];hs=hf*(365/hf)**base[18];rates=np.log(2)*365/np.array([hs,hf]);b=5*base[[16,20]];levels,f=physical(base)
    for t in range(nt):
        tr=probe['transition'][t];dt=probe['dt'][t];eps=np.array([[normal(seed,(t*pairs+i)*(nd+3)+j) for j in range(nd+3)] for i in range(pairs)]);dw=np.sqrt(dt)*eps[:,0];projected=state@d['projection'].T
        for anti,sign in enumerate([1,-1]):
            X=probe['mean'][t]+sign*projected;eta=p[0]*X[:,0]+p[1]*(X[:,1]-p[2])**2+p[3]*(X[:,2]-p[4])**2;old=np.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6])-probe['lognorm'][t]);A=np.exp(sign*zs*b-.5*b*b*(-np.expm1(-2*rates*probe['time'][t])));stress=(1-base[22])*A[:,0]+base[22]*A[:,1]
            if probe['reference'][t]:qr[anti]=1-base[17]+base[17]*stress
            Q=(1-base[17])*old+base[17]*stress
            for m,x in enumerate(models):
                v=f+(.2-f)*x[probe['seg'][t]]*Q;vr=f+(.2-f)*x[probe['seg'][t]]*qr[anti]
                if probe['coarse'][t]:hold[anti,:,m]=v
                vv=np.stack([hold[anti,:,m],v],axis=1);logs[anti,:,m,:]+=-.5*vv*dt+np.sqrt(vv)*sign*dw[:,None];logref[anti,:,m]+=-.5*vr*dt+np.sqrt(vr)*sign*dw;iv[anti,:,m]+=vr*dt
                if probe['expiry'][t]>=0:
                    col=probe['expiry'][t];expected[anti*pairs:(anti+1)*pairs,col,m,:2]=np.exp(logs[anti,:,m,:]);expected[anti*pairs:(anti+1)*pairs,col,m,2]=np.exp(logref[anti,:,m]);expected[anti*pairs:(anti+1)*pairs,col,m,3]=iv[anti,:,m]
        state=state*d['decay'][tr]+eps[:,1:nd+1]@d['root'][tr].T+dw[:,None]*d['cross'][tr];zs=zs*np.exp(-rates*dt)+np.sqrt(-np.expm1(-2*rates*dt))*eps[:,nd+1:]
    np.testing.assert_allclose(actual,expected,rtol=3e-12,atol=3e-14)
    # A full deterministic mean/covariance recurrence checks risk-neutral drift,
    # retained-state marginalization and all transition sizes over the horizon.
    prior,dated,_,_=context();keep=np.array(json.loads((fresh/'transition_audit.json').read_text())['retained_modes']);rates0=dated.layout.rates[keep];stationary=dated.stationary[np.ix_(keep,keep)];initial=arr(prior['initial_state']);cov=np.zeros_like(stationary);innov=np.array([r@r.T+np.outer(c,c)*dt for r,c,dt in zip(d['root'],d['cross'],np.unique(d['dt']))]);probes=set(np.linspace(0,len(d['dt'])-1,23,dtype=int));cov_errors=[];mean_errors=[]
    for t,dt in enumerate(d['dt']):
        tr=d['transition'][t];dd=d['decay'][tr];cov=dd[:,None]*cov*dd[None,:]+innov[tr]
        if t in probes:
            elapsed=d['time'][t]+dt;oracle=stationary*(-np.expm1(-(rates0[:,None]+rates0[None,:])*elapsed));error=float(np.max(abs(cov-oracle)));assert error<2e-9;cov_errors.append(error);mu,_=dated.gaussian_projection_law(float(d['time'][t]),initial);error=float(np.max(abs(mu-d['mean'][t])));assert error<2e-12;mean_errors.append(error)
    stored=np.load(fresh/'stocks.npz');seeds=stored['seeds'];stocks=stored['stock'];N=len(stocks);pairbatch=N//16;replays=[]
    for batch,seed in enumerate(seeds):
        replay=k.simulate(d,pairbatch,p,models,int(seed));oracle=np.concatenate([stocks[batch*pairbatch:(batch+1)*pairbatch],stocks[N//2+batch*pairbatch:N//2+(batch+1)*pairbatch]]);np.testing.assert_array_equal(replay,oracle);replays.append(int(seed))
    moments=[]
    for day,j in zip([29,57,119,211,302,666,1030],range(7)):
        for m,name in enumerate(['baseline','selected']):
            for level,steps in enumerate([64,128]):
                S=stocks[:,j,m,level];paired=(S[:N//2]+S[N//2:])/2;moments.append(dict(day=day,model=name,steps=steps,stock_mean=float(S.mean()),stock_pair_SE=float(paired.std(ddof=1)/np.sqrt(N//2)),minimum_stock=float(S.min()),maximum_stock=float(S.max())))
    pd.DataFrame(moments).to_csv(fresh/'stock_moments.csv',index=False)
    protection=run/'VIX_protection';protection.mkdir(exist_ok=True);selected=HERE/'runs/m2_two_speed_001';records={}
    for name in ['selection.json','experimental_model_two_speed.json','fresh_check/prices_otm_two_speed.csv','fresh_check/forwards_two_speed.csv','fresh_check/score_two_speed.json']:
        src=selected/name;dst=protection/Path(name).name;shutil.copy2(src,dst);records[str(src)]=hashlib.sha256(src.read_bytes()).hexdigest()
    dump(protection/'invariance.json',dict(baseline_VIX_vector=vix_x(base).tolist(),candidate_VIX_vector=vix_x(chosen).tolist(),exact_array_equality=True,source_hashes=records,interpretation='Identical VIX parameter inputs and frozen exogenous state law imply identical VIX distributions. Copied existing prices/forwards are inherited, not newly estimated or newly certified population values.'))
    result=dict(status='implementation_checks_PASS; fit_accuracy_separate',verified_source_hashes=checked,all9_price_derivative_checks=True,independent_fresh_kernel_numpy_replay_error=float(np.max(abs(actual-expected))),retained_OU_covariance_semigroup_max_error=max(cov_errors),Q_projection_mean_max_error=max(mean_errors),all8_fresh_batch_seed_replays=replays,exact_VIX_map_and_fixed_parameter_invariance=True,positive_spike_shift=float(p[4]),free_indices=FREE.tolist(),limitations=['Independent fresh simulation checking at one market snapshot, not market OOS.','Eight-batch bootstrap is approximate; fitting variability and numerical biases excluded.','Two stock grids do not prove time-grid convergence.','No final bid-ask calibration, model-family infeasibility or Lean theorem.'])
    dump(run/'audit.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);args=ap.parse_args();audit(args.run.resolve())
