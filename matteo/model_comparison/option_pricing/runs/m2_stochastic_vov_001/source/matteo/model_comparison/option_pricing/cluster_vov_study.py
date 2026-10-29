"""After cap profiles, fixed-curve cluster-state vol-of-vol diagnostics."""
import ctypes as ct,json,subprocess,time
import numpy as np,pandas as pd
from .full_readout_fit import HERE,arr,dump
from .m2_spx import ROOT,frozen_stage3
from .maturity_curve import CurveKernel
from .cap_tail_study import engine,tails
from .fit_parity_parameters_v3 import save_result
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from esn_eq.pricing import DatedM2Pricer

def run():
    parent=HERE/'runs/m2_cap_tail_001';assert (parent/'run.json').exists()
    out=parent/'cluster_vov';out.mkdir(exist_ok=False);prior=HERE/'runs/m2_maturity_curve_fit_001';d=json.loads((prior/'run.json').read_text());p=np.array(d['selected']);cache=HERE/'runs/m2_parity_adaptive_cache_001';spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv');vc=json.loads((prior/'vix_controls.json').read_text())
    files=[__file__,HERE/'cluster_vov_kernel.cpp',HERE/'cap_tail_study.py',prior/'selected_model.json',parent/'run.json'];captured=capture_sources([__import__('pathlib').Path(s).resolve() for s in files]);save_sources(out,captured,ROOT)
    dump(out/'design.json',dict(parameters=p.tolist(),strengths=[.5,1.],equation='G=(1-w(C))*exp(2h(eta))+w(C)*exp(2h(eta+b*qJ*(F-betaJ)^2)); w(C)=0.5+0.45*tanh(0.5*C)',status='Fixed-curve diagnostic, not calibrated or adopted. Original OU unchanged; no new randomness.',source_sha256=source_hashes(captured,ROOT)))
    lib=out/'cluster.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'cluster_vov_kernel.cpp'),'-o',str(lib)],check=True);k=CurveKernel(lib,lib,np.array(d['ends_years']));k.lib.configure.argtypes=[ct.c_int,ct.c_double];k.lib.direct.argtypes=k.lib.expectation.argtypes;a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial=np.array(d['initial_state']);scores=[];checks=[];ts=[];ms=[]
    e=engine(out,cache,k,dated,initial,spx,vix,f,vc, __import__('pathlib').Path(d['controls']),d)
    try:
        for strength in [.5,1.]:
            name=f'cluster{int(strength*100):03}';k.lib.configure(0,strength);e.last_readout=None;e.exposure_readout_key=None;tick=time.perf_counter();print('CLUSTER',name,flush=True);r=e.predict(p,True,False);scores.append(save_result(e,r,p,name,out));print('SCORE',name,r['original_otm_components'],time.perf_counter()-tick,flush=True)
            for day in [27,55,118,153]:
                aa,laws=e.vixlaws[(day,'production')];v=r['raw'][f'vix_{day}'];ids=np.unique(np.r_[0,1,np.argmax(v),np.argsort(v)[-2]])
                for index in [0,len(laws)//2,len(laws)-1]:
                    mean,o,w=laws[index];direct=np.empty((len(ids),8));k.lib.direct(len(ids),len(w),arr(mean[ids]),o,w,arr(p),direct);pred=k.expectation(*laws[index],p)[ids,0];err=float(np.max(abs(pred/direct[:,0]-1)));assert err<1e-6,(day,index,err);checks.append(dict(candidate=name,day=day,node=index,error=err))
            t,m=tails(out,name,vix,p,np.array(d['ends_days']));ts+=t;ms+=m;dump(out/'scores.json',scores);dump(out/'direct_checks.json',checks);pd.DataFrame(ts).to_csv(out/'tails.csv',index=False);pd.DataFrame(ms).to_csv(out/'moments.csv',index=False)
    finally:e.close()
    check_sources(captured);dump(out/'run.json',dict(status='completed two cluster-modulated VIX diagnostics',scores=scores,stock_repriced=False,not_joint_calibration=True))
if __name__=='__main__':run()
