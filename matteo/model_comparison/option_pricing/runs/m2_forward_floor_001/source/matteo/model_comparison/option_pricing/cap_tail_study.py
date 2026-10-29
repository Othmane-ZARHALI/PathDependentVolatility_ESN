"""Bounded cap profiles, exact same cached primitives; no production mutation."""
import argparse,ctypes as ct,json,time,subprocess,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from .full_readout_fit import HERE,arr,dump
from .m2_spx import ROOT,M2,frozen_stage3
from .maturity_curve import CurveKernel,CurveEngine,vix_only_loss
from .parity_adaptive import canonical_market
from .fit_parity_parameters_v3 import save_result
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from esn_eq.pricing import DatedM2Pricer

VARIANTS=[('baseline',0,0.,0.,None),('lower_L',0,0.,-.75,None),('higher_L',0,0.,.75,None),('soft_k1',0,0.,0.,1.),('sharp_k4',0,0.,0.,4.),('algebraic',1,0.,0.,None),('logtail025',2,.25,0.,None),('logtail050',2,.5,0.,None),('separate',3,0.,0.,None)]

def tails(folder,name,quotes,p,ends):
    rows=[];moments=[]
    with np.load(folder/f'paths_{name}.npz') as a:
        for days,q in quotes.groupby('expiry_days',sort=True):
            v=a[f'vix_{days}'];w=a[f'weights_vix_{days}'];second=(w*v*v).reshape(8,512).mean(1)
            length=np.maximum(0,np.minimum(ends,days+30)-np.maximum(np.r_[0,ends[:-1]],days))
            moments.append(dict(candidate=name,days=int(days),maximum_sample=float(v.max()),mass=float(w.mean()),second=float(second.mean()),second_SE=float(second.std(ddof=1)/np.sqrt(8)),target=float(10000*length@p[7:-1]/30)))
            for K in sorted(set([float(q.forward.iloc[0]),40.,60.,float(q.strike.max())])):
                hit=v>K;prob=w*hit;pay=w*np.maximum(v-K,0);pn=pay.reshape(8,512).mean(1);pr=prob.reshape(8,512).mean(1);total=pay.sum();order=np.argsort(pay)[::-1]
                rows.append(dict(candidate=name,days=int(days),strike=K,hits=int(hit.sum()),probability=float(prob.mean()),prob_SE=float(pr.std(ddof=1)/np.sqrt(8)),raw_call=float(pay.mean()),call_SE=float(pn.std(ddof=1)/np.sqrt(8)),severity=float(pay.mean()/prob.mean()) if prob.sum()>0 else None,payoff_ESS=float(total**2/(pay@pay)) if total>0 else 0,largest_share=float(pay.max()/total) if total>0 else None,top10_share=float(pay[order[:10]].sum()/total) if total>0 else None,net_calls=pn.tolist()))
    return rows,moments

def engine(run,cache,k,dated,initial,spx,vix,f,vc,controls,design,vix_only=True):
    e=CurveEngine(cache,k,dated,initial,spx,vix,f,controls=vc,frozen_controls=controls,ends=np.array(design['ends_years']),reference_levels=design['reference_levels'])
    if vix_only:
        e.spx=spx.iloc[:0];e.original_quotes=vix.reset_index(drop=True);e.quotes=canonical_market(e.original_quotes);e.loss=vix_only_loss(e.quotes,f);e.otm_loss=vix_only_loss(e.original_quotes,f)
    return e

def run(out):
    out.mkdir(parents=True,exist_ok=False);start=time.perf_counter();prior=HERE/'runs/m2_maturity_curve_fit_001';design=json.loads((prior/'run.json').read_text());p=np.array(design['selected']);ends=np.array(design['ends_days'])
    sources=[Path(__file__),HERE/'cap_tail_kernel.cpp',HERE/'CAP_TAIL_PLAN.md',prior/'run.json',prior/'selected_model.json',prior/'vix_controls.json',HERE/'active_scope.json',*sorted(HERE.glob('*.py')),*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(sources);save_sources(out,captured,ROOT)
    lib=out/'cap.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'cap_tail_kernel.cpp'),'-o',str(lib)],check=True)
    k=CurveKernel(lib,lib,np.array(design['ends_years']));k.lib.configure.argtypes=[ct.c_int,ct.c_double];k.lib.direct.argtypes=k.lib.expectation.argtypes
    cache=Path(design['cache']);ctrl=Path(design['controls']);spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv');vc=json.loads((prior/'vix_controls.json').read_text())
    a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial=np.array(design['initial_state'])
    config=dict(status='running cap profiles',baseline=p.tolist(),variants=VARIANTS,source_sha256=source_hashes(captured,ROOT),parent=str(prior),cache=str(cache),checking_cache=design['check_cache'],randomness='Exact saved OU states, shocks, likelihood weights, eight512-point VIX scrambles; no new random streams',paths_per_strike=4096,market_out_of_sample=False)
    dump(out/'design.json',config);e=engine(out,cache,k,dated,initial,spx,vix,f,vc,ctrl,design);scores=[];tailrows=[];moments=[];normrows=[];checks=[]
    try:
        for name,mode,strength,dl,kap in VARIANTS:
            tick=time.perf_counter();k.lib.configure(mode,strength);x=p.copy();x[5]+=dl
            if kap is not None:x[6]=kap
            e.last_readout=None;e.exposure_readout_key=None
            print('VIX',name,flush=True);r=e.predict(x,True,False);score=save_result(e,r,x,name,out);score['seconds']=time.perf_counter()-tick;scores.append(score)
            if name=='baseline':
                saved=pd.read_csv(prior/'prices_canonical_selected.csv');saved=saved[saved.instrument.eq('VIX')]
                discrepancy=float(np.max(abs(saved.model_price_U.to_numpy()-r['prices'])));assert discrepancy<2e-6,discrepancy;dump(out/'baseline_replay.json',dict(maximum_price_difference=discrepancy))
            norm,_=e.normalizer(x)
            for day in ends:
                nn=norm(day/365)[0];normrows.append(dict(candidate=name,days=float(day),M=float(np.exp(nn)),G_sup=float(np.exp(2*x[5])) if mode in [0,1] else None))
            # Representative exact same finite-GH conditional checks, including largest VIX states.
            for day in [27,153]:
                aa,laws=e.vixlaws[(day,'production')];v=r['raw'][f'vix_{day}'];ids=np.unique(np.r_[0,1,np.argmax(v),np.argsort(v)[-2]])
                for index in [0,len(laws)//2,len(laws)-1]:
                    mean,offset,w=laws[index];outd=np.empty((len(ids),8));k.lib.direct(len(ids),len(w),arr(mean[ids]),offset,w,arr(x),outd);pred=k.expectation(*laws[index],x)[ids,0];err=float(np.max(abs(pred/outd[:,0]-1)));assert err<2e-7,(name,day,index,err);checks.append(dict(candidate=name,day=day,node=index,maximum_relative_error=err))
            t,m=tails(out,name,vix,x,ends);tailrows+=t;moments+=m
            pd.DataFrame(tailrows).to_csv(out/'tails.csv',index=False);pd.DataFrame(moments).to_csv(out/'moments.csv',index=False);pd.DataFrame(normrows).to_csv(out/'normalizers.csv',index=False);dump(out/'direct_checks.json',checks);dump(out/'scores.json',scores)
            print('DONE',name,score['original_otm_components'],score['seconds'],flush=True)
    finally:e.close()
    check_sources(captured);dump(out/'run.json',config|dict(status='cap VIX profiles complete; follow-up selection pending',scores=scores,wall_seconds=time.perf_counter()-start))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();run(args.output.resolve())
