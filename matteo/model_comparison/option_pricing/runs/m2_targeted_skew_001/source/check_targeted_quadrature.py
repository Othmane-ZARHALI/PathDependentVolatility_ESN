"""Matched-state refinement of both Gaussian quadratures and VIX time rule."""
import json
import numpy as np
from .targeted_skew import OUT,Kernel,ComponentMeans,ALPHA,exposures
from .asymmetric_readout import OUT as OLD,PRIOR,arr,dump
from .protected_curve_check import context
from .vix_only_curve_fit import DAYS,DELTA
from .two_speed_joint_core import physical

def main():
    prior,dated,norm,_=context();initial=arr(prior['initial_state']);p=arr(prior['selected']);design=json.loads((OUT/'design.json').read_text());scale=design['scale'];x=arr(design['baseline']);k=Kernel(OUT/'targeted.dylib');M=ComponentMeans(k,dated,initial,p,norm,ALPHA,scale).build();indices=[0,3,4];alpha=ALPHA[indices];draws={};selection={}
    for day in DAYS:
        a=dict(np.load(OLD/f'fresh_VIX_draws/d{day}.npz'));idx=np.concatenate([np.flatnonzero(a['nets']==j)[:16] for j in range(8)]);selection[day]=idx;draws[day]={key:value[idx] for key,value in a.items()}
    H=exposures(k,dated,initial,p,norm,lambda t:M(t)[indices],x,draws,alpha,scale,48,48,12,OUT/'quadrature_probe');rows=[];config=json.loads((OUT/'candidates.json').read_text());models=arr(config['models']);_,floor=physical(x)
    for day in DAYS:
        old=np.load(OUT/f'checking_exposures/d{day}.npz')['H'][indices][:,selection[day]];fine=H[day];levels,_=physical(models[9]);coeff=levels-floor;I0=floor*DELTA+old@coeff;I1=floor*DELTA+fine@coeff;V0=100*np.sqrt(I0/DELTA);V1=100*np.sqrt(I1/DELTA);rows.append(dict(day=day,alpha=alpha.tolist(),maximum_VIX_change_points=np.max(abs(V1-V0),axis=1).tolist(),scope='128matched states/expiry, representative regime-alpha1.5-refit curve; GH32/stress32/time8 versus GH48/stress48/time12, not a uniform bound'))
    dump(OUT/'quadrature_probe.json',rows);print('TARGETED QUADRATURE PROBE',rows,flush=True)

if __name__=='__main__':main()
