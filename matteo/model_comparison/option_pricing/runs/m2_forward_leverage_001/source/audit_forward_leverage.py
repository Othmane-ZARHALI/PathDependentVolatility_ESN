"""Independent finite-difference future-variance sensitivity and rule checks."""
import json
import numpy as np
from .forward_leverage import OUT,TARGET,L,setup,basis,signals,arr,dump
from .asymmetric_readout import Kernel as MomentKernel
from esn_eq.outer_sampling import GaussianOuterSampler

def field(k,states,stress,signal,b,p,x,t):return k.probe(states,stress,{key:signal[key][0] for key in ['futuremean','inverse_norm','hat']},b,p,x,t)

def oracle_value(states,stress,signal,b,p,x,t,dated,k):
    values=np.zeros(len(states));w=x[17];r=x[22];bs=5*x[16];bf=5*x[20];L0=np.log(365);hf=np.exp(L0*x[21]);hs=hf*(365/hf)**x[18];a=np.log(2)*365/np.array([hs,hf]);q=-np.expm1(-2*a*t)
    for j,lag in enumerate(b['lags']):
        means=signal['futuremean'][0,j]+states@b['projection_lag'][j].T;_,cov=dated.gaussian_projection_law(float(lag),np.zeros(len(dated.layout.rates)));EG=k.moments(means,cov,p,0,1,arr([0]),32)[:,0,1];dec=np.exp(-a*lag);EA=(1-r)*np.exp(bs*dec[0]*stress[:,0]-.5*bs*bs*dec[0]**2*q[0])+r*np.exp(bf*dec[1]*stress[:,1]-.5*bf*bf*dec[1]**2*q[1]);values+=signal['hat'][0,j]*((1-w)*EG*signal['inverse_norm'][0,j]+w*EA)
    return values

def main():
    prior,dated,norm,_,p,x,d,k=setup(OUT);b12=dict(np.load(OUT/'basis.npz'));b24=basis(k,dated,d,p,22,24);b32=basis(k,dated,d,p,22,32);b44=basis(k,dated,d,p,44,24);rng=np.random.default_rng(20261028577);rows=[];largest=0.;laglargest=0.;probe_fields=[]
    for day in [0,29,119,666,1030]:
        t=day/365;sig=signals(d,dated,prior,norm,x,b12,[t]);sig44=signals(d,dated,prior,norm,x,b44,[t]);sampler=GaussianOuterSampler(dated,t,arr(prior['initial_state']));draw=arr(rng.standard_normal((12,len(sampler.mean))));states=(draw@sampler.root.T)[:,:24];stress=rng.standard_normal((12,2));L0=np.log(365);hf=np.exp(L0*x[21]);hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf]);stress*=np.sqrt(-np.expm1(-2*rates*t));states[8:]*=2;stress[8:]*=2
        K12=field(k,states,stress,sig,b12,p,x,t);K24=field(k,states,stress,sig,b24,p,x,t);K32=field(k,states,stress,sig,b32,p,x,t);K44=field(k,states,stress,sig44,b44,p,x,t);directions=[-z/np.linalg.norm(z,axis=1)[:,None] for z in [K12,K24,K32,K44]];e12=float(np.max(np.linalg.norm(directions[0]-directions[2],axis=1)));e24=float(np.max(np.linalg.norm(directions[1]-directions[2],axis=1)));elag=float(np.max(np.linalg.norm(directions[1]-directions[3],axis=1)));largest=max(largest,e12);laglargest=max(laglargest,elag);rows.append(dict(day=day,max_direction_distance_GH12_vs32=e12,max_direction_distance_GH24_vs32=e24,max_direction_distance_22_vs44_lag_knots=elag));probe_fields.append((t,states,stress,K32,sig))
    # Full primitive finite-difference of independently repriced F_t, GH32.
    t,states,stress,K,sig=probe_fields[2];state=states[0];zs=stress[0];ob=d['diffusion'].shape[1];n=ob+2;perturbed=np.repeat(state[None,:],2*n,axis=0);zstress=np.repeat(zs[None,:],2*n,axis=0);rates=dated.layout.rates[:24];load=dated.layout.driver_loadings[:24];used=np.flatnonzero(np.any(load!=0,axis=0));B=np.sqrt(2*rates)[:,None]*load[:,used];L0=np.log(365);hf=np.exp(L0*x[21]);hs=hf*(365/hf)**x[18];a=np.log(2)*365/np.array([hs,hf]);h=1e-5
    for j in range(n):
        if j<ob:perturbed[2*j]+=h*B[:,j];perturbed[2*j+1]-=h*B[:,j]
        else:zstress[2*j,j-ob]+=h*np.sqrt(2*a[j-ob]);zstress[2*j+1,j-ob]-=h*np.sqrt(2*a[j-ob])
    mk=MomentKernel(OUT/'forward.dylib');V=oracle_value(perturbed,zstress,sig,b32,p,x,t,dated,mk);fd=(V[::2]-V[1::2])/(2*h);fd_error=float(np.max(abs(fd-K[0])/(1+abs(fd))));assert fd_error<2e-6,fd_error
    # Exact integration of the curve under the piecewise-linear hat partition.
    ends=arr(prior['ends_years']);f=(.05+.07*x[19])**2;coeff=(.2-f)*x[:16];sum_errors=[]
    for day in [0,26,27,28,54,117,150,1015,1030]:
        t=day/365;s=signals(d,dated,prior,norm,x,b12,[t]);left=np.r_[0,ends[:-1]];exact=np.maximum(0,np.minimum(ends,t+L)-np.maximum(left,t))@coeff+max(0,t+L-max(t,ends[-1]))*coeff[-1];sum_errors.append(abs(s['hat'].sum()-exact))
    assert max(sum_errors)<1e-13
    dump(OUT/'integration_audit.json',dict(rows=rows,finite_difference_primitive_sensitivity_max_scaled_error=fd_error,hat_curve_integral_max_error=max(sum_errors),GH12_vs32_max_direction_distance=largest,lag22_vs44_max_direction_distance=laglargest,production_rule_approved=bool(largest<.03 and laglargest<.03),scope='60representative Q/stress states including amplified tails; not uniform integration bounds. Stock direction is valid unit predictable loading even when numerically approximated.'))
    np.savez_compressed(OUT/'probe_fields.npz',fd=fd,K=K[0]);print('FORWARD INTEGRATION AUDIT',largest,laglargest,fd_error,'approved',largest<.03 and laglargest<.03,flush=True)

if __name__=='__main__':main()
