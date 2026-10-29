"""Time-interpolated conditional means and new forward-fit Jacobian checks."""
import json
import numpy as np
import pandas as pd
from .targeted_skew import OUT,Kernel,ALPHA,VIX,load_outer,arr,dump
from .asymmetric_readout import PRIOR
from .protected_curve_check import context

def main():
    prior,dated,norm,_=context();p=arr(prior['selected']);initial=arr(prior['initial_state']);design=json.loads((OUT/'design.json').read_text());scale=design['scale'];x=arr(design['baseline']);k=Kernel(OUT/'targeted.dylib');data=np.load(OUT/'component_means.npz');indices=[0,500,10000,len(data['times'])-1];rows=[]
    for i in indices:
        t=float(data['times'][i]);mean,cov=dated.gaussian_projection_law(t,initial);actual=k.moments(mean,cov,p,mean[1],scale,ALPHA,64)[0];actual[:,1]/=np.exp(norm(t)[0]);actual[0]=1.;interpolated=data['moments'][i];errs=[]
        for A in [1e-6,.1,1,10,1e4]:
            value=((1-x[17])*actual[:,1]+x[17]*A*actual[:,0])*((1-x[17])+x[17]*A)/((1-x[17])*interpolated[:,1]+x[17]*A*interpolated[:,0]);target=(1-x[17])+x[17]*A;errs.append(float(np.max(abs(value/target-1))))
        assert max(errs)<2e-6,errs;rows.append(dict(time=t,maximum_relative_conditional_mean_error=max(errs)))
    H={day:np.load(OUT/f'training_exposures/d{day}.npz')['H'] for day in [27,55,118,153]};draws=load_outer(128,20261028491);q=pd.read_csv(OUT/'VIX_membership.csv');f=pd.read_csv(OUT/'forward_targets.csv');v=VIX(H,draws,x,q,f);y=np.array([.3,.2,.8,.3,.8,.4,.4]);r=v.predict(3,y);errors=[]
    for j in range(7):
        delta=np.zeros(7);delta[j]=1e-6;fd=(v.predict(3,y+delta)['forwards']-v.predict(3,y-delta)['forwards'])/(2e-6);err=float(np.max(abs(fd-r['fjac'][:,j])));assert err<1e-7;errors.append(err)
    dump(OUT/'extended_audit.json',dict(all_checks_passed=True,conditional_mean_interpolation_checks=rows,new_regime_forward_derivative_max_absolute_error=max(errors),semantic_note='Adaptive direction uses normalized sensitivity of H0. For xi>floor this equals normalized sensitivity of variance. At xi=floor it uses that continuous limiting direction; all unit directions tie at zero instantaneous covariance. If the H0 gradient itself vanishes, the old fixed rotation is used. This refines the initial plan fallback convention without changing factor dynamics.'))
    print('EXTENDED AUDIT PASS',max(r['maximum_relative_conditional_mean_error'] for r in rows),max(errors),flush=True)

if __name__=='__main__':main()
