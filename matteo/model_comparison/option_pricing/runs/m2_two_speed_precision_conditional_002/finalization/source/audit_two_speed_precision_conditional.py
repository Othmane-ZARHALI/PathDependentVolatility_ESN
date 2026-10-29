"""Independent conditional-value, OU-bridge covariance and curve-gradient audits."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtri
from scipy.stats import qmc

from .full_readout_fit import HERE,arr,dump
from .two_speed_precision import PrecisionKernel,reference_grid
from .two_speed_precision_conditional import ConditionalKernel,bridge
from .two_speed_joint_study import load_spx
from .two_speed_joint_core import physical
from esn_eq.conditional_spx import lognormal_payoffs


def value(logS,I,K):
    return lognormal_payoffs(logS[:,None]-.5*I[:,None],I[:,None],K[None,:],np.full((1,len(K)),'call'))


def audit(run):
    design=json.loads((run/'design.json').read_text());parent=Path(design['parent']);x=arr(design['baseline']);k=ConditionalKernel(run/'conditional.dylib');p=PrecisionKernel(parent/'precision.dylib');checked=0
    for path in [run,parent]:
        de=json.loads((path/'design.json').read_text());manifest=json.loads((path/'artifact_manifest.json').read_text())
        for name,h in de['source_hashes'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==h,name;checked+=1
        for name,h in manifest.items():assert hashlib.sha256((path/name).read_bytes()).hexdigest()==h,name;checked+=1
    rng=np.random.default_rng(2026102823);nt=17;n=5;dt=arr(rng.uniform(.0005,.003,nt));t=arr(np.r_[0,np.cumsum(dt)]);d=dict(old=arr(np.exp(rng.normal(size=(nt,n))*.3-.045)),equity=arr(rng.normal(size=(nt,n))*np.sqrt(dt)[:,None]),eps=arr(rng.normal(size=(nt,n,2))),time=arr(t[:-1]),dt=dt,seg=np.minimum(np.arange(nt),15).astype('int32'),grid=t)
    target,_,grid=p.stock(d,x);K=arr([.8,1.,1.2]);actual=k.correction(d,x,target[:,2],K);oracle=np.zeros((n,len(K)));zs=np.zeros((n,2));logs=np.zeros(n);prefix=np.zeros(n);_,mask,*_=reference_grid(d)
    hf=365**x[21];hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf]);b=5*x[[16,20]];levels,f=physical(x)
    for j,dt0 in enumerate(dt):
        A=np.exp(zs*b-.5*b*b*(-np.expm1(-2*rates*t[j])));qs=(1-x[17])+x[17]*((1-x[22])*A[:,0]+x[22]*A[:,1])
        if mask[j]:qr=qs.copy()
        q=qs+(1-x[17])*(d['old'][j]-1);v=f+(levels[d['seg'][j]]-f)*q;vr=f+(levels[d['seg'][j]]-f)*qr;remaining=target[:,2]-prefix
        oracle+=value(logs,remaining+(v-vr)*dt0,K)-value(logs,remaining,K);prefix+=vr*dt0;logs+=-.5*v*dt0+np.sqrt(v)*d['equity'][j];zs=zs*np.exp(-rates*dt0)+np.sqrt(-np.expm1(-2*rates*dt0))*d['eps'][j]
    np.testing.assert_allclose(actual[:,:-1],oracle,rtol=2e-10,atol=3e-14);np.testing.assert_allclose(actual[:,-1],target[:,0],rtol=2e-12)
    # OU-bridge coefficients reproduce the exact zero-initialized Gaussian law.
    grid=(t,dt,d['seg']);nodes,left,right,c=bridge(grid,x);bridge_errors=[]
    for factor,a in enumerate(rates):
        M=np.zeros((len(dt),len(nodes)))
        for j,node in enumerate(nodes):M[node]=c[j,3*factor]*M[left[j]]+c[j,3*factor+1]*M[right[j]];M[node,j]+=c[j,3*factor+2]
        ts=grid[0][:-1];cov=np.exp(-a*abs(ts[:,None]-ts[None,:]))-np.exp(-a*(ts[:,None]+ts[None,:]));error=float(np.max(abs(M@M.T-cov)));assert error<2e-14;bridge_errors.append(error)
    H=k.integrals(grid,x,4,207);normals=ndtri(qmc.Sobol(d=2*len(nodes),scramble=True,seed=207).random_base2(4));z=np.zeros((len(dt),16,2));expected=np.zeros_like(H)
    for j,node in enumerate(nodes):
        for factor in range(2):z[node,:,factor]=c[j,3*factor]*z[left[j],:,factor]+c[j,3*factor+1]*z[right[j],:,factor]+c[j,3*factor+2]*normals[:,2*j+factor]
    for j in range(len(dt)):
        A=np.exp(z[j]*b-.5*b*b*(-np.expm1(-2*rates*t[j])));Q=1-x[17]+x[17]*((1-x[22])*A[:,0]+x[22]*A[:,1]);expected[:,d['seg'][j]]+=Q*dt[j]
    np.testing.assert_allclose(H,expected,rtol=2e-13,atol=2e-15)
    # Full-horizon probe covers every curve interval, not just early coordinates.
    data,_=load_spx('search');d0=data[1030][0];ids=np.array([0,7,19,47]);probe={name:arr(d0[name][:,ids]) for name in ['old','equity','eps']};probe.update({name:np.ascontiguousarray(d0[name]) for name in ['time','dt','seg','grid']});xx=x.copy();xx[:16]=.3;aa=p.stock(probe,xx)[0];derivatives=[]
    for j in range(16):
        dx=np.zeros(23);dx[j]=1e-6;fd=(p.stock(probe,xx+dx)[0][:,:3]-p.stock(probe,xx-dx)[0][:,:3])/2e-6;analytic=aa[:,[3+j,19+j,35+j]];error=float(np.max(abs(fd-analytic)/(1+abs(fd))));assert error<3e-6 and np.max(abs(analytic))>1e-8;derivatives.append(error)
    # Reconstruct one stored scramble at every expiry and role, including the
    # full long-horizon bridge (rather than only a synthetic short-grid test).
    replay=[]
    for role in ['search','check']:
        dd,_=load_spx(role)
        for day in sorted(dd):
            stored=np.load(run/role/f'd{day}_reference.npz');grid=(stored['grid'],np.diff(stored['grid']),stored['segments']);h=k.integrals(grid,x,int(np.log2(len(stored['H'][0]))),int(stored['seeds'][0]));np.testing.assert_array_equal(h,stored['H'][0]);replay.append(dict(role=role,day=day))
    # Nontrivial importance-weighted paired expectation test of the new control.
    n=65536;nt=8;dt=arr(np.full(nt,.03));times=arr(np.r_[0,np.cumsum(dt)]);z=rng.normal(size=(nt,n));tilt=.3;dw=arr(np.sqrt(dt)[:,None]*z+tilt*dt[:,None]);weight=np.exp(-tilt*dw.sum(0)+.5*tilt*tilt*times[-1]);state=np.r_[np.zeros((1,n)),np.cumsum(dw,axis=0)[:-1]];old=arr(np.exp(.9*state-.5*.9**2*times[:-1,None]));toy=dict(old=old,equity=dw,eps=arr(rng.normal(size=(nt,n,2))),time=arr(times[:-1]),dt=dt,seg=np.zeros(nt,np.int32),grid=times);target,_,_=p.stock(toy,x);cor=k.correction(toy,x,target[:,2],K);raw=np.maximum(target[:,0,None]-K,0);B=value(np.zeros(n),target[:,2],K);difference=weight[:,None]*(raw-cor[:,:-1])-B;mean=difference.mean(0);se=difference.std(0,ddof=1)/np.sqrt(n);assert np.max(abs(mean/se))<5
    result=dict(status='implementation_checks_PASS; pricing_precision_separate',verified_hashes=checked,independent_Python_conditional_correction_error=float(np.max(abs(actual[:,:-1]-oracle))),OU_bridge_covariance_errors=bridge_errors,OU_bridge_integral_replay_error=float(np.max(abs(H-expected))),full_horizon_curve_derivative_errors=derivatives,stored_scramble_replays=replay,finite_law_weighted_control_test=dict(mean=mean.tolist(),SE=se.tolist(),z_scores=(mean/se).tolist()),limitations=['Finite-law checks are not universal formal proofs.','Existing512 expensive paths per expiry remain the limiting expensive sample.','No optimizer, correlation change or model adoption.'])
    dump(run/'audit.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);args=ap.parse_args();audit(args.run.resolve())
