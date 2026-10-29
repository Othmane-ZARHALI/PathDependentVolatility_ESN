"""Fixed-factor, paired leverage scan; prices remain signed and quotes retained."""
import ctypes as ct
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import ndtr
from .full_readout_fit import HERE,arr,dump
from .protected_curve_check import loss
from .protected_curve_fit import black_parts
from .plot_joint_curve_fit import implied_vol
from .vix_only_curve_fit import PRIOR

PARENT=HERE/'runs/m2_protected_curve_001'
OUT=HERE/'runs/m2_leverage_scan_001'
THETAS=[0,.2,.4,.6,.8,.95]

class Kernel:
    def __init__(self,lib):
        self.k=ct.CDLL(str(lib));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.k.leverage_stats.argtypes=[ct.c_int]*4+[I]+[P]*7+[I]*4+[P]*3+[ct.c_uint64,P]
        self.k.leverage_reference.argtypes=[ct.c_int]*3+[P]*2+[I]*2+[P,ct.c_uint64,P]
    def stats(self,d,pairs,p,x,seed):
        out=np.zeros((2*pairs,int(d['expiries']),3,4));self.k.leverage_stats(len(d['dt']),pairs,d['projection'].shape[1],int(d['expiries']),d['transition'],d['decay'],d['root'],d['cross'],d['mean'],d['projection'],d['time'],d['dt'],d['seg'],d['coarse'],d['reference'],d['expiry'],d['lognorm'],arr(p),arr(x),int(seed),out);return out
    def reference(self,d,pairs,x,seed):
        idx=np.flatnonzero(d['reference']);grid=np.r_[d['time'][idx],d['time'][-1]+d['dt'][-1]];expiry=np.full(len(idx),-1,np.int32)
        for t,e in enumerate(d['expiry']):
            if e>=0:
                end=d['time'][t]+d['dt'][t];j=np.searchsorted(grid,end);assert abs(grid[j]-end)<1e-12;expiry[j-1]=e
        out=np.zeros((2*pairs,int(d['expiries']),4));self.k.leverage_reference(len(idx),pairs,int(d['expiries']),arr(grid[:-1]),arr(np.diff(grid)),d['seg'][idx].copy(),expiry,arr(x),int(seed),out);return out

def combine(batches):
    n=len(batches[0])//2
    return np.concatenate([np.concatenate([b[:n] for b in batches]),np.concatenate([b[n:] for b in batches])])

def stock(s,theta,cs,cf):
    return np.exp(-.5*s[...,0]+np.sqrt(1-theta*theta)*s[...,1]-theta*(cs*s[...,2]+cf*s[...,3]))

def blocks(y):
    n=len(y)//2
    return ((y[:n]+y[n:])/2).reshape(8,-1,y.shape[-1]).mean(1)

def evaluate(stats,ref,q,x,label,out):
    days=sorted(q.expiry_days.unique());directions={'mixed':(np.sqrt(1-x[22]),np.sqrt(x[22])),'slow':(1.,0.),'fast':(0.,1.)};summaries=[];frames={};allblocks={}
    floor=(.05+.07*x[19])**2;ends=np.array(json.loads((PRIOR/'run.json').read_text())['ends_years']);levels=floor+(.2-floor)*x[:16]
    for direction,(cs,cf) in directions.items():
        for theta in THETAS:
            cheap=[];res=[[],[]];frames_by_grid=[[],[]]
            for j,day in enumerate(days):
                g=q[q.expiry_days==day].copy();K=g.K_over_F.to_numpy();F=g.forward.to_numpy();a=ref[:,j];I=a[:,0];v=(1-theta**2)*I;m=-.5*I-theta*(cs*a[:,2]+cf*a[:,3]);sd=np.sqrt(v)[:,None];d1=(m[:,None]+v[:,None]-np.log(K)[None,:])/sd;mu=np.exp(m+v/2);B=mu[:,None]*ndtr(d1)-K*ndtr(d1-sd);cb=blocks(B)*F;cheap.append(cb)
                T=day/365;lengths=np.maximum(0,np.minimum(ends,T)-np.minimum(np.r_[0,ends[:-1]],T));detI=float(lengths@levels);delta=ndtr(-np.log(K)/np.sqrt(detI)+.5*np.sqrt(detI));R=stock(stats[:,j,2],theta,cs,cf)
                for grid in range(2):
                    S=stock(stats[:,j,grid],theta,cs,cf);y=np.maximum(S[:,None]-K,0)-np.maximum(R[:,None]-K,0)-delta*(S-R)[:,None];eb=blocks(y)*F;res[grid].append(eb);h=g.copy();h['model_price_U']=eb.mean(0)+cb.mean(0)-h.option_type.eq('put').to_numpy()*(F-h.strike.to_numpy());h['mc_se_U']=np.sqrt(eb.var(0,ddof=1)/8+cb.var(0,ddof=1)/8);h['stock_mean']=S.mean();h['stock_mean_SE']=blocks(S[:,None]).std(ddof=1)/np.sqrt(8);frames_by_grid[grid].append(h)
            C=np.concatenate(cheap,axis=1)
            for grid,steps in enumerate([64,128]):
                E=np.concatenate(res[grid],axis=1);frame=pd.concat(frames_by_grid[grid],ignore_index=True);primary=frame.priority_selected.to_numpy();p=frame[primary].reset_index(drop=True);iv=np.array([np.nan if (z:=implied_vol(r.model_price_U,r.forward,r.strike,r.T_years,r.option_type)) is None else z for r in frame.itertuples()]);market=np.array([np.nan if (z:=implied_vol(r.price_mid_U,r.forward,r.strike,r.T_years,r.option_type)) is None else z for r in frame.itertuples()]);frame['model_IV']=iv;frame['market_IV']=market;frame['direction']=direction;frame['theta']=theta;frame['steps']=steps;key=f'{direction}_{theta:g}_{steps}';frames[key]=frame;allblocks[key]=(E,C);frame.to_csv(out/f'prices_{key}.csv',index=False)
                iverr=np.abs(iv-market);ivmae=float(np.mean([np.nanmean(iverr[primary & frame.expiry_days.eq(day)]) for day in days]));inside=(p.model_price_U>=p.price_bid_U)&(p.model_price_U<=p.price_ask_U);summaries.append(dict(direction=direction,theta=theta,steps=steps,loss=loss(p,p.model_price_U.to_numpy()),inside=int(inside.sum()),IV_MAE_pp=ivmae,invalid_IV=int((~np.isfinite(iv[primary])).sum()),negative=int((p.model_price_U<0).sum()),median_SE_scales=float(np.median(p.mc_se_U/p.price_scale_U))))
            print(label,direction,theta,'loss',summaries[-1]['loss'],'IVMAE',summaries[-1]['IV_MAE_pp'],flush=True)
    base=frames['mixed_0_128'];E0,C0=allblocks['mixed_0_128'];rng=np.random.default_rng(20261028091)
    for s in summaries:
        key=f"{s['direction']}_{s['theta']:g}_{s['steps']}";frame=frames[key];E,C=allblocks[key];d=frame[['expiry_days','strike','option_type','priority_selected','price_scale_U']].copy();d['difference_from_theta0_U']=frame.model_price_U-base.model_price_U;d['paired_SE_U']=np.sqrt((E-E0).var(0,ddof=1)/8+(C-C0).var(0,ddof=1)/8);d.to_csv(out/f'paired_{key}.csv',index=False)
        if s['steps']==128:
            diffs=[];mask=base.priority_selected.to_numpy();pq=base[mask].reset_index(drop=True);parity=base.option_type.eq('put').to_numpy()*(base.forward-base.strike).to_numpy()
            for _ in range(256):
                ie=rng.integers(8,size=8);ic=rng.integers(8,size=8);diffs.append(loss(pq,(E[ie].mean(0)+C[ic].mean(0)-parity)[mask])-loss(pq,(E0[ie].mean(0)+C0[ic].mean(0)-parity)[mask]))
            s['paired_loss_delta_CI95']=np.quantile(diffs,[.025,.975]).tolist()
    pd.DataFrame(summaries).to_csv(out/'scores.csv',index=False);dump(out/'scores.json',summaries);np.savez_compressed(out/'blocks.npz',**{k+'_E':v[0] for k,v in allblocks.items()},**{k+'_C':v[1] for k,v in allblocks.items()});return summaries

def run():
    start=time.perf_counter();OUT.mkdir(exist_ok=False);d=dict(np.load(PARENT/'fresh_check/grid_and_transitions.npz'));x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);p=arr(json.loads((PRIOR/'run.json').read_text())['selected']);q=pd.read_csv(PARENT/'spx_membership.csv');q.to_csv(OUT/'spx_membership.csv',index=False);(OUT/'source').mkdir()
    sources=[Path(__file__),HERE/'leverage_scan.cpp',HERE/'two_speed_precision_kernel.cpp',PARENT/'selection.json',PARENT/'fresh_check/grid_and_transitions.npz',PRIOR/'run.json',HERE/'LEVERAGE_SCAN_PLAN.md'];hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources}
    for f in sources:
        if f.suffix!='.npz':shutil.copy2(f,OUT/'source'/f.name)
    dump(OUT/'design.json',dict(parameters=x.tolist(),theta=THETAS,directions=['mixed','slow','fast'],factor_law_changed=False,VIX_law_changed=False,primary_quotes=int(q.priority_selected.sum()),eligible_quotes=len(q),source_hashes=hashes,estimator='conditional reference call + paired target/reference payoff difference - fixed deterministic Black delta times stock difference; antithetic pairs and 8 independent batches',selection_metric='minimum fine-grid equal-maturity price interval loss on replay scan; freeze before fresh checking',martingale_scope='finite predictable log-Euler stocks have unit expectation; no continuous-time martingale theorem claimed'))
    lib=OUT/'leverage.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'leverage_scan.cpp'),'-o',str(lib)],check=True);k=Kernel(lib);old=np.load(PARENT/'fresh_check/stocks.npz');replay=[]
    for seed in old['seeds']:
        t=time.perf_counter();replay.append(k.stats(d,32,p,x,int(seed)));print('REPLAY',len(replay),'seconds',time.perf_counter()-t,flush=True)
    stats=combine(replay);np.savez_compressed(OUT/'replay_stats.npz',stats=stats,seeds=old['seeds']);actual=np.exp(-.5*stats[:,:,:,0]+stats[:,:,:,1]);expected=old['stock'][:,:,1,:3];err=float(np.max(abs(actual-expected)));assert err<1e-10
    ref=[];seeds=[]
    for b in range(8):
        seed=int(np.random.SeedSequence([2026102819,100,b]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed);t=time.perf_counter();ref.append(k.reference(d,4096,x,seed));print('REFERENCE',b,'seconds',time.perf_counter()-t,flush=True)
    ref=combine(ref);np.savez_compressed(OUT/'reference.npz',stats=ref,seeds=np.array(seeds,dtype=np.uint64));scan=OUT/'replay';scan.mkdir();scores=evaluate(stats,ref,q,x,'REPLAY',scan);best=min([s for s in scores if s['steps']==128],key=lambda s:s['loss']);dump(OUT/'selection.json',dict(selected=best,selection_frozen_before_fresh=True,not_adopted=True));print('FROZEN SELECTION',best,flush=True)
    fresh=[];seeds=[]
    for b in range(8):
        seed=int(np.random.SeedSequence([2026102819,200,b]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed);t=time.perf_counter();fresh.append(k.stats(d,64,p,x,seed));print('FRESH',b,'seconds',time.perf_counter()-t,flush=True)
    stats=combine(fresh);np.savez_compressed(OUT/'fresh_stats.npz',stats=stats,seeds=np.array(seeds,dtype=np.uint64));check=OUT/'fresh';check.mkdir();freshscores=evaluate(stats,ref,q,x,'FRESH',check)
    dump(OUT/'run.json',dict(status='completed fixed-factor exploratory scan',wall_seconds=time.perf_counter()-start,replay_theta0_max_absolute_stock_error=err,fresh_expensive_paths=len(stats),cheap_reference_paths=len(ref),checking_reference_reused_from_scan=True,selection=best,scores=freshscores,no_parameters_adopted=True,VIX_note='factor law and readout exactly preserved algebraically; inherited forward accuracy, no new population forward certificate'))
    print('COMPLETE',time.perf_counter()-start,flush=True)

if __name__=='__main__':run()
