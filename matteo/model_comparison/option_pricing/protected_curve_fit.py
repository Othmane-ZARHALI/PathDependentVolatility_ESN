"""Nine-curve-level exploratory SPX fit with the entire VIX law protected."""
import argparse
import hashlib
import json
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.special import ndtr

from .full_readout_fit import HERE,arr,dump
from .two_speed_precision import PrecisionKernel,FREE
from .two_speed_joint_core import physical,vix_x,metadata
from .two_speed_joint_study import load_spx
from .vix_only_curve_fit import PRIOR

PRECISION=HERE/'runs/m2_two_speed_precision_002'
CONDITIONAL=HERE/'runs/m2_two_speed_precision_conditional_002'
PARENT=HERE/'runs/m2_two_speed_joint_001'
SELECTED=HERE/'runs/m2_two_speed_001'


def black_parts(iv,K):
    sd=np.sqrt(iv[:,None]);d1=-np.log(K)[None,:]/sd+.5*sd;delta=ndtr(d1);phi=np.exp(-.5*d1*d1)/np.sqrt(2*np.pi)
    return delta-K[None,:]*ndtr(d1-sd),phi/(2*sd)


class ProtectedFit:
    def __init__(self,data,quotes,reference,points=512):
        self.data=data;self.market=quotes;self.k=PrecisionKernel(PRECISION/'precision.dylib');self.reference={};self.points=points
        for day in data:
            a=np.load(reference/f'd{day}_reference.npz');H=arr(a['H'][:,:points].reshape(-1,16));self.reference[day]=(H,float(a['grid'][-1]))
    def predict(self,x,detail=False):
        levels,f=physical(x);prices=[];jac=[];exp_blocks=[];cheap_blocks=[];raw={}
        for day,g in self.market.groupby('expiry_days',sort=True):
            a=np.concatenate([self.k.stock(d,x)[0] for d in self.data[day]]);w=np.concatenate([d['weights'] for d in self.data[day]]);S,R,iv=a[:,:3].T;ds=a[:,3:19];dr=a[:,19:35];di=a[:,35:51];K=g.K_over_F.to_numpy();F=g.forward.to_numpy();sd=np.sqrt(iv[:,None]);d1=-np.log(K)[None,:]/sd+.5*sd;delta=ndtr(d1);phi=np.exp(-.5*d1*d1)/np.sqrt(2*np.pi);ddelta=phi*(np.log(K)[None,:]/(2*iv[:,None]**1.5)+.25/sd)
            y=w[:,None]*(np.maximum(S[:,None]-K,0)-np.maximum(R[:,None]-K,0)-delta*(S-R)[:,None]);dy=w[:,None,None]*(((S[:,None]>K)-delta)[:,:,None]*ds[:,None,:]-((R[:,None]>K)-delta)[:,:,None]*dr[:,None,:]-(S-R)[:,None,None]*ddelta[:,:,None]*di[:,None,:])
            H,T=self.reference[day];I=f*T+(.2-f)*(H@x[:16]);B,vega=black_parts(I,K);cheap=B.reshape(8,-1,len(K)).mean(1);dB=(.2-f)*(vega.T@H)/len(H)
            price=(y.mean(0)+B.mean(0))*F-g.option_type.eq('put').to_numpy()*(F-g.strike.to_numpy());prices.extend(price);jac.extend((dy.mean(0)+dB)*F[:,None]);exp_blocks.append(y.reshape(8,-1,len(K)).mean(1)*F);cheap_blocks.append(cheap*F)
            if detail:raw[day]=dict(out=a,weights=w,residual=y,cheap=cheap)
        prices=np.array(prices);jac=np.array(jac);q=self.market;lo=q.price_bid_U.to_numpy();hi=q.price_ask_U.to_numpy();scale=q.price_scale_U.to_numpy();dist=np.where(prices<lo,prices-lo,np.where(prices>hi,prices-hi,0));wt=np.concatenate([np.full(len(g),1/(7*len(g))) for _,g in q.groupby('expiry_days',sort=True)]);r=np.sqrt(wt)*dist/scale;J=np.sqrt(wt)[:,None]*((prices<lo)|(prices>hi))[:,None]*jac/scale[:,None]
        return dict(prices=prices,price_jacobian=jac,residuals=r,jacobian=J,loss=float(r@r),expensive_blocks=np.concatenate(exp_blocks,axis=1),cheap_blocks=np.concatenate(cheap_blocks,axis=1),raw=raw)


def record(engine,x,out,name):
    r=engine.predict(x,True);q=engine.market.copy();q['model_price_U']=r['prices'];q['mc_se_U']=np.sqrt(r['expensive_blocks'].var(0,ddof=1)/8+r['cheap_blocks'].var(0,ddof=1)/8);q['interval_residual']=np.where(q.model_price_U<q.price_bid_U,q.model_price_U-q.price_bid_U,np.where(q.model_price_U>q.price_ask_U,q.model_price_U-q.price_ask_U,0))/q.price_scale_U;q.to_csv(out/f'prices_{name}.csv',index=False);np.savez_compressed(out/f'blocks_{name}.npz',expensive=r['expensive_blocks'],cheap=r['cheap_blocks']);dump(out/f'score_{name}.json',dict(loss=r['loss'],parameters=x.tolist(),inside=int(q.interval_residual.eq(0).sum()),negative=int((q.model_price_U<0).sum())));return r


def run(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();baseline=arr(json.loads((PARENT/'design.json').read_text())['baseline']);quotes=pd.read_csv(PARENT/'spx_contracts.csv');assert len(quotes)==1294
    sources=[Path(__file__),HERE/'PROTECTED_CURVE_PLAN.md',HERE/'two_speed_precision.py',HERE/'two_speed_precision_kernel.cpp',HERE/'two_speed_joint_core.py',SELECTED/'selection.json',SELECTED/'experimental_model_two_speed.json',PARENT/'design.json',CONDITIONAL/'design.json'];hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources};(out/'source').mkdir()
    for f in sources:shutil.copy2(f,out/'source'/f.name)
    data,_=load_spx('search');e=ProtectedFit(data,quotes,CONDITIONAL/'search',args.points);tick=time.perf_counter();r=e.predict(baseline);print('BASELINE loss',r['loss'],'seconds',time.perf_counter()-tick,flush=True)
    # Check native-price derivatives on every free direction using a four-path,
    # full-horizon kernel probe; include nonconstant cheap Black/delta terms.
    probe=baseline.copy();probe[FREE]=.2;check=[]
    small={day:[{**d,**{key:arr(d[key][:,:4]) for key in ['old','equity','eps']},'weights':arr(d['weights'][:4])}] for day,ds in data.items() for d in ds[:1]}
    # Eight uncertainty blocks require a multiple of8. Derivative test only:
    # repeat four selected paths once, leaving finite-difference identity intact.
    for ds in small.values():
        for d in ds:
            for key in ['old','equity','eps']:d[key]=arr(np.concatenate([d[key],d[key]],axis=1))
            d['weights']=arr(np.tile(d['weights'],2))
    pe=ProtectedFit(small,quotes,CONDITIONAL/'search',64);a=pe.predict(probe)
    for j in FREE:
        dx=np.zeros(23);dx[j]=1e-6;fd=(pe.predict(probe+dx)['prices']-pe.predict(probe-dx)['prices'])/2e-6;error=float(np.max(abs(fd-a['price_jacobian'][:,j])/(1+abs(fd))));assert error<3e-5,(j,error);check.append(dict(index=int(j),relative_error=error))
    dump(out/'derivative_checks.json',check);print('ALL9 PRICE DERIVATIVES PASS',flush=True)
    dump(out/'design.json',dict(baseline=baseline.tolist(),source_hashes=hashes,free_indices=FREE.tolist(),training_expensive_paths=512,training_reference_points=8*args.points,quote_count=1294,objective='unchanged native bid-ask interval loss, equal maturity weights',estimator='terminal-delta independent-stress auxiliary; analytic curve derivatives',exploratory=True,all_quote_accuracy_gate='final acceptance only, superseded as a fitting prerequisite by explicit user authorization'))
    train=out/'training';train.mkdir();baseline_record=record(e,baseline,train,'baseline');ledger=[];best=baseline.copy();bestloss=r['loss']
    for label,u in [('baseline_start',baseline[FREE]),('floor_start',np.zeros(9))]:
        memo={};count=[0]
        def predict(v):
            key=tuple(v)
            if key not in memo:
                x=baseline.copy();x[FREE]=v;np.testing.assert_array_equal(vix_x(x),vix_x(baseline));memo.clear();memo[key]=e.predict(x);count[0]+=1
                if count[0]%5==0:print('FIT',label,'evaluation',count[0],'loss',memo[key]['loss'],'seconds',round(time.perf_counter()-begin,1),flush=True)
            return memo[key]
        opt=least_squares(lambda u:predict(u)['residuals'],u,jac=lambda u:predict(u)['jacobian'][:,FREE],bounds=(np.zeros(9),np.ones(9)),max_nfev=args.evaluations,ftol=1e-6,xtol=1e-6,gtol=1e-6)
        candidate=baseline.copy();candidate[FREE]=opt.x;rr=predict(opt.x);ledger.append(dict(label=label,loss=rr['loss'],parameters=candidate.tolist(),success=bool(opt.success),message=opt.message,nfev=opt.nfev,optimality=opt.optimality));dump(out/'optimizer_ledger.json',ledger)
        if rr['loss']<bestloss:best=candidate.copy();bestloss=rr['loss']
        print('END',label,rr['loss'],opt.message,flush=True)
    np.testing.assert_array_equal(vix_x(best),vix_x(baseline));dump(out/'selection.json',dict(parameters=best.tolist(),training_loss=bestloss,selected_before_checking=True,adopted=False,VIX_map_identical=True));record(e,best,train,'selected')
    prior=json.loads((PRIOR/'run.json').read_text());p=arr(prior['selected']);levels,f=physical(best);p[7:-1]=levels;p[-1]=f;old=json.loads((SELECTED/'experimental_model_two_speed.json').read_text());old.update(parameters=p.tolist(),readout=metadata(best,p),adopted=False);dump(out/'experimental_model_protected.json',old)
    for f,h in hashes.items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h
    dump(out/'run.json',dict(status='protected training completed; fresh validation pending',baseline_loss=r['loss'],selected_loss=bestloss,wall_seconds=time.perf_counter()-begin,production_model_changed=False,VIX_parameters_changed=False))
    print('PROTECTED FIT COMPLETE',bestloss,'seconds',time.perf_counter()-begin,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--points',type=int,default=512);ap.add_argument('--evaluations',type=int,default=40);args=ap.parse_args();run(args)
