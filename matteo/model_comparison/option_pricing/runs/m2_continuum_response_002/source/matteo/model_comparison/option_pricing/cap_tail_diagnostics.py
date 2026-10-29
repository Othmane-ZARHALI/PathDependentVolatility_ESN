"""Targeted cap response/moment diagnostics, distinct from option convergence."""
import json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from scipy.special import ndtri,expit,logsumexp
from scipy.stats import qmc
from .full_readout_fit import HERE,Kernel,arr,dump
from .m2_spx import frozen_stage3
from .maturity_curve import CurveEngine,window_nodes
from .cap_tail_study import engine
from esn_eq.pricing import DatedM2Pricer,_normal_rule
from esn_eq.states import symmetric_square_root

def run(out):
    prior=HERE/'runs/m2_maturity_curve_fit_001';d=json.loads((prior/'run.json').read_text());p=np.array(d['selected']);cache=Path(d['cache']);spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv');vc=json.loads((prior/'vix_controls.json').read_text())
    a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial=np.array(d['initial_state']);k=Kernel(prior/'convolution.dylib');e=engine(out,cache,k,dated,initial,spx,vix,f,vc,Path(d['controls']),d);norm,_=e.normalizer(p)
    marginal=[];conditional=[];bounds=[];eigen=[]
    try:
        with np.load(prior/'paths_selected.npz') as paths:
            for day in [27,55,118,153]:
                t=(day+15)/365;m,c=dated.gaussian_projection_law(t,initial);root=symmetric_square_root(c);A=np.diag([0,p[1],p[3]]);lam=float(np.linalg.eigvalsh(4*root.T@A@root).max());eigen.append(dict(day=day+15,uncapped_quadratic_eigenvalue=lam,integrable_uncapped=lam<1))
                gh=[]
                for order in [32,48,64,96]:
                    z,w=_normal_rule(order);gh.append(float(k.expectation(arr(m),arr(z@root.T),arr(w),p)[0,0]))
                for net in range(8):
                    seed=202610250+day*10+net;u=qmc.Sobol(4,scramble=True,seed=seed).random_base2(15);z=ndtri(u[:,:3]);scale=np.where(u[:,3]<.5,1.,2.5);y=z*scale[:,None];x=m+y@root.T
                    logbase=-.5*np.sum(y*y,axis=1);loginflate=-3*np.log(2.5)-.5*np.sum((y/2.5)**2,axis=1);w=np.exp(logbase-np.logaddexp(logbase,loginflate)+np.log(2))
                    eta=p[0]*x[:,0]+p[1]*(x[:,1]-p[2])**2+p[3]*(x[:,2]-p[4])**2;g=k.multiplier(x,p)[:,0];s=expit(p[6]*(p[5]-eta));wg=w*g
                    marginal.append(dict(days=day+15,net=net,seed=seed,points=len(w),mass=float(w.mean()),M=float(wg.mean()),M_GH32=gh[0],M_GH48=gh[1],M_GH64=gh[2],M_GH96=gh[3],prob_over_transition=float(np.mean(w*(eta>p[5]))),mean_share_over_transition=float(np.sum(wg*(eta>p[5]))/wg.sum()),prob_slope_below_01=float(np.mean(w*(s<.1))),mean_share_slope_below_01=float(np.sum(wg*(s<.1))/wg.sum()),G_weighted_response=float(np.sum(wg*s)/wg.sum())))
                times,wt=window_nodes(day/365,np.array(d['ends_years']),e.tau);xi=p[7:-1][np.minimum(np.searchsorted(d['ends_years'],times,side='right'),15)]
                bound=100*np.sqrt(p[-1]+np.sum(wt*(xi-p[-1])*np.exp(2*p[5]-norm(times)[:,0]))/(30/365))
                bounds.append(dict(days=day,deterministic_VIX_upper_bound=float(bound),max_strike=float(vix[vix.expiry_days.eq(day)].strike.max()),maximum_saved=float(paths[f'vix_{day}'].max())))
                aa=np.load(cache/'production'/f'vix_d{day}.npz');law=e.law(15/365,aa['states'],32);ev=k.expectation(*law,p);response=1-ev[:,6]/(2*ev[:,0]);v=paths[f'vix_{day}'];weights=aa['weights']
                for threshold in [float(vix[vix.expiry_days.eq(day)].forward.iloc[0]),60.,float(vix[vix.expiry_days.eq(day)].strike.max())]:
                    pay=weights*np.maximum(v-threshold,0);conditional.append(dict(days=day,strike=threshold,payoff_weighted_G_response_midwindow=float(pay@response/pay.sum()),minimum_state_response=float(response.min()),maximum_state_response=float(response.max())))
        pd.DataFrame(marginal).to_csv(out/'marginal_response.csv',index=False);pd.DataFrame(conditional).to_csv(out/'conditional_response.csv',index=False);dump(out/'cap_bounds_and_integrability.json',dict(bounds=bounds,uncapped=eigen,randomness='8 independent new scrambled Sobol nets per time,32768 each; mixture half unit half2.5-standard-deviation Gaussian. Weights p/q, no self normalization.',limitation='Direct marginal moments and representative midpoint conditional responses; not VIX-price convergence or global error bounds.'))
    finally:e.close()
    print(json.dumps(dict(bounds=bounds,uncapped=eigen),indent=2))

if __name__=='__main__':run(HERE/'runs/m2_cap_tail_001')
