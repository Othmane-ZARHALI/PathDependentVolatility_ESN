"""Independent RNG replay and finite-law audits of the stress-reference control."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr

from .full_readout_fit import HERE, dump
from .two_speed_precision import PrecisionKernel, residuals, PARENT, SELECTED, FREE
from .two_speed_joint_core import physical, vix_x


def normal(seed,i):
    mask=(1<<64)-1
    def mix(v):
        v&=mask;v=((v^(v>>30))*0xbf58476d1ce4e5b9)&mask;v=((v^(v>>27))*0x94d049bb133111eb)&mask;return v^(v>>31)
    a=mix(seed+2*i*0x9e3779b97f4a7c15);b=mix(seed+(2*i+1)*0x9e3779b97f4a7c15)
    return math.sqrt(-2*math.log(((a>>11)+.5)*2**-53))*math.cos(2*math.pi*((b>>11)+.5)*2**-53)


def audit(run):
    design=json.loads((run/'design.json').read_text());x=np.array(design['baseline']);k=PrecisionKernel(run/'precision.dylib');checked=0
    for record in ['design.json','input_hashes.json','artifact_manifest.json']:
        d=json.loads((run/record).read_text());hashes=d['source_hashes'] if record=='design.json' else d
        for name,h in hashes.items():
            f=run/name if record=='artifact_manifest.json' else Path(name)
            assert hashlib.sha256(f.read_bytes()).hexdigest()==h,str(f);checked+=1
    np.testing.assert_allclose(vix_x(x),json.loads((SELECTED/'selection.json').read_text())['parameters'],rtol=0,atol=2e-16)
    old=json.loads((SELECTED/'experimental_model_two_speed.json').read_text());levels,floor=physical(x);np.testing.assert_allclose(levels,old['parameters'][7:-1],rtol=0,atol=3e-17);assert abs(floor-old['parameters'][-1])<1e-17 and old['parameters'][4]>0
    # Scalar Python replay of both signs of the counter-Gaussian cheap OU paths.
    dt=np.array([.001,.007,.013,.009,.02]);times=np.r_[0,np.cumsum(dt)];seg=np.array([0,1,1,4,7],np.int32);pairs=5;seed=2026102891
    actual=k.integrals((times,dt,seg),x,pairs,seed);expected=np.zeros_like(actual)
    hf=365**x[21];hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf]);b=5*x[[16,20]];z=np.zeros((pairs,2))
    for t,d in enumerate(dt):
        for sign_index,sign in enumerate([1,-1]):
            A=np.exp(sign*z*b-.5*b*b*(-np.expm1(-2*rates*times[t])));Q=1-x[17]+x[17]*((1-x[22])*A[:,0]+x[22]*A[:,1]);expected[:,sign_index,seg[t]]+=Q*d
        eps=np.array([[normal(seed,2*(t*pairs+i)+j) for j in range(2)] for i in range(pairs)])
        z=z*np.exp(-rates*d)+np.sqrt(-np.expm1(-2*rates*d))*eps
    np.testing.assert_allclose(actual,expected,rtol=3e-13,atol=2e-15)
    # Nontrivial adapted target with importance-tilted equity and independent
    # future stress: each claimed zero-mean hedge is checked under the full law.
    rng=np.random.default_rng(2026102893);n=131072;nt=8;dt=np.full(nt,.03);z=rng.normal(size=(nt,n));stress=rng.normal(size=(nt,n));tilt=.3;dw=np.sqrt(dt)[:,None]*(z+tilt*np.sqrt(dt)[:,None]);weight=np.exp(-tilt*dw.sum(0)+.5*tilt*tilt*dt.sum())
    old_state=np.r_[np.zeros((1,n)),np.cumsum(dw,axis=0)[:-1]];R=np.exp(.9*old_state-.5*.9**2*np.r_[0,np.cumsum(dt)[:-1]][:,None]);H=np.exp(.6*stress-.5*.6**2);vr=.02+.05*(.1+.9*H);v=.02+.05*(.1*R+.9*H)
    Sa=np.exp(np.r_[np.zeros((1,n)),np.cumsum(-.5*v*dt[:,None]+np.sqrt(v)*dw,axis=0)]);Sr=np.exp(np.r_[np.zeros((1,n)),np.cumsum(-.5*vr*dt[:,None]+np.sqrt(vr)*dw,axis=0)]);iv=np.r_[np.zeros((1,n)),np.cumsum(vr*dt[:,None],axis=0)];path=np.stack([Sa,Sr,iv],axis=2);out=np.c_[Sa[-1],Sr[-1],iv[-1]];strikes=np.array([.85,1.,1.15]);rr=residuals(out,path,strikes,weight)
    records=[]
    for name in ['terminal','dynamic']:
        control=rr['unhedged']-rr[name];mean=control.mean(0);se=control.std(0,ddof=1)/np.sqrt(n);zscore=mean/se;assert np.max(abs(zscore))<5,(name,zscore)
        records.append(dict(method=name,weighted_control_mean=mean.tolist(),SE=se.tolist(),z_scores=zscore.tolist()))
    # Record component precision and paired differences to the earlier estimator.
    summary=pd.read_csv(run/'summary.csv');g=summary[(summary.role=='check')&(summary.method=='dynamic')];comparisons=[]
    for day in sorted(g.day):
        q=pd.read_csv(run/'check/prices_dynamic.csv');q=q[q.expiry_days==day];q=q[q.priority_selected]
        oldq=pd.read_csv(PARENT/'checking/prices_spx_vix_only.csv');oldq=oldq[oldq.expiry_days==day]
        merged=q.merge(oldq,on=['expiry_days','strike','option_type'],suffixes=('_new','_old'),validate='one_to_one');assert len(merged)==len(q)
        comparisons.append(dict(day=int(day),old_median_SE_scales=float(np.median(merged.mc_se_U_old/merged.price_scale_U_old)),new_median_SE_scales=float(np.median(merged.mc_se_U_new/merged.price_scale_U_new)),median_SE_reduction_factor=float(np.median(merged.mc_se_U_old/merged.mc_se_U_new)),maximum_price_difference_scales=float(np.max(abs(merged.model_price_U_new-merged.model_price_U_old)/merged.price_scale_U_old))))
    result=dict(status='implementation_checks_PASS; numerical_precision_separate',verified_hashes=checked,cheap_RNG_and_antithetic_replay_error=float(np.max(abs(actual-expected))),finite_adapted_law_control_checks=records,protected_VIX_parameter_and_model_replay=True,positive_spike_shift=old['parameters'][4],free_indices=FREE.tolist(),checking_comparisons=comparisons,limitations=['Finite-law statistical control checks are not a universal proof.','Existing expensive streams reused; no market OOS or untouched checking claim.','Sampling screen does not certify timestep, normalization, forward or calibration accuracy.'])
    dump(run/'audit.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);args=ap.parse_args();audit(args.run.resolve())
