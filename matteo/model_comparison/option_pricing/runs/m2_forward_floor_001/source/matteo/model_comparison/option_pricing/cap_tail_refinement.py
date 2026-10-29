"""Representative extreme-state GH/time refinement, not a surface certificate."""
import ctypes as ct,json
from pathlib import Path
import numpy as np,pandas as pd
from .full_readout_fit import HERE,arr,dump
from .m2_spx import frozen_stage3
from .maturity_curve import CurveKernel,window_nodes
from .cap_tail_study import engine
from esn_eq.pricing import DatedM2Pricer,_normal_rule
from esn_eq.states import symmetric_square_root


def run():
    out=HERE/'runs/m2_cap_tail_001';prior=HERE/'runs/m2_maturity_curve_fit_001';d=json.loads((prior/'run.json').read_text());p=np.array(d['selected']);cache=Path(d['cache']);spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv');vc=json.loads((prior/'vix_controls.json').read_text());a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial=np.array(d['initial_state']);k=CurveKernel(out/'cap.dylib',out/'cap.dylib',np.array(d['ends_years']));k.lib.configure.argtypes=[ct.c_int,ct.c_double];records=[];nr=[]
    for name,mode,strength,dl in [('baseline',0,0.,0.),('higher_L',0,0.,.75),('logtail025',2,.25,0.)]:
        print('REFINEMENT',name,flush=True);k.lib.configure(mode,strength);x=p.copy();x[5]+=dl;e=engine(out,cache,k,dated,initial,spx,vix,f,vc,Path(d['controls']),d);norm,_=e.normalizer(x)
        try:
            for day in [27,55,118,153]:
                aa=np.load(cache/'production'/f'vix_d{day}.npz')
                with np.load(out/f'paths_{name}.npz') as saved:v=saved[f'vix_{day}'];weights=saved[f'weights_vix_{day}']
                K=float(vix[vix.expiry_days.eq(day)].strike.max());ids=np.unique(np.r_[0,np.argmax(v),np.argmax(weights*np.maximum(v-K,0))]);states=aa['states'][ids];vv={}
                for gh,order in [(32,8),(48,8),(32,16)]:
                    times,tw=window_nodes(day/365,np.array(d['ends_years']),e.tau,order);integ=np.zeros(len(ids))
                    for t,wt in zip(times,tw):
                        law=e.law(t-day/365,states,gh);ev=k.expectation(*law,x)[:,0];xi=x[7:-1][min(np.searchsorted(d['ends_years'],t,side='right'),15)];integ+=wt*(x[-1]+(xi-x[-1])*ev/np.exp(norm(t)[0]))
                    vv[(gh,order)]=100*np.sqrt(integ/(30/365))
                for j,ix in enumerate(ids):records.append(dict(candidate=name,day=day,path_index=int(ix),saved_value=float(v[ix]),GH32_time8=float(vv[(32,8)][j]),GH48_time8=float(vv[(48,8)][j]),GH32_time16=float(vv[(32,16)][j]),GH_change=float(vv[(48,8)][j]-vv[(32,8)][j]),time_change=float(vv[(32,16)][j]-vv[(32,8)][j])))
                mean,c=dated.gaussian_projection_law((day+15)/365,initial);root=symmetric_square_root(c);ms=[]
                for order in [48,64,96]:
                    z,w=_normal_rule(order);ms.append(float(k.expectation(arr(mean),arr(z@root.T),arr(w),x)[0,0]))
                nr.append(dict(candidate=name,day=day+15,M48=ms[0],M64=ms[1],M96=ms[2],relative48to96=ms[0]/ms[2]-1))
                pd.DataFrame(records).to_csv(out/'representative_refinement.csv',index=False);pd.DataFrame(nr).to_csv(out/'normalizer_refinement.csv',index=False)
        finally:e.close()
    dump(out/'refinement.json',dict(status='completed',scope='Three representative outer states per expiry (first, largestVIX, largest weighted highest-strike payoff), all four dates, three fixed-curve profiles; no full-surface/global convergence claim',maximum_conditional_GH_change=float(max(abs(r['GH_change']) for r in records)),maximum_time_change=float(max(abs(r['time_change']) for r in records)),source_sha256=__import__('hashlib').sha256(Path(__file__).read_bytes()).hexdigest()))
if __name__=='__main__':run()
