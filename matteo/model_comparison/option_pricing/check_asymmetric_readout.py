"""Fresh full-factor VIX draws for every frozen multiplier strength."""
import json,time
import numpy as np
import pandas as pd
from .asymmetric_readout import OUT,PRIOR,ALPHA,Kernel,Normalizer,exposures,VIX,save_vix,arr,dump
from .protected_curve_check import context
from .vix_only_curve_fit import ARCHIVE,DAYS,ACTIVE
from .mixed_response_study import ZSHIFTS,ZPROB
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture

def main():
    begin=time.perf_counter();prior,dated,norm,_=context();initial=arr(prior['initial_state']);p=arr(prior['selected']);design=json.loads((OUT/'design.json').read_text());scale=design['scale'];x=arr(design['baseline']);k=Kernel(OUT/'asymmetric.dylib');Z=Normalizer(k,dated,initial,p,norm,ALPHA,scale).build(x[17]);props=json.loads((ARCHIVE/'m2_vix_targeted_check_001/proposals.json').read_text());draws={};folder=OUT/'fresh_VIX_draws';folder.mkdir(exist_ok=False);seeds=[]
    for day in DAYS:
        dec=props[str(day)];shifts=arr([np.r_[old,s,f] for old in dec['shifts'] for s in ZSHIFTS for f in ZSHIFTS]);prob=arr([po*ps*pf for po in dec['probabilities'] for ps in ZPROB for pf in ZPROB]);mix=DefensiveGaussianMixture(shifts,prob);sampler=GaussianOuterSampler(dated,day/365,initial);zz=[];ww=[]
        for batch in range(8):
            seed=int(np.random.SeedSequence([20261028343,day,batch]).generate_state(1,dtype=np.uint64)[0]);z,w,_=mix.draw(9,seed);zz.append(z);ww.append(w);seeds.append(dict(day=day,batch=batch,seed=seed))
        z=np.concatenate(zz);a=dict(states=sampler.mean+z[:,:-2]@sampler.root.T,z_normals=z[:,-2],fast_normals=z[:,-1],weights=np.concatenate(ww),nets=np.repeat(np.arange(8),512),source_indices=np.arange(len(z)));draws[day]=a;np.savez_compressed(folder/f'd{day}.npz',**a);assert abs(np.linalg.norm(sampler.root@sampler.root.T-sampler.covariance))<1e-9
    dump(folder/'design.json',dict(seeds=seeds,points_per_expiry=4096,proposal='full joint defensive Gaussian mixture on old36 + slow + fast coordinates, exact target/proposal weights; no sample mass normalization',new_factor_draws=True,all_strengths_reuse_same_states=True,selection_already_frozen=True))
    H=exposures(k,dated,initial,p,norm,Z,x,draws,ALPHA,scale,32,8,OUT/'fresh_VIX_exposures');q=pd.read_csv(OUT/'VIX_membership.csv');f=pd.read_csv(OUT/'forward_targets.csv');v=VIX(H,draws,x,q,f);models=arr(json.loads((OUT/'candidates.json').read_text())['models']);folder=OUT/'fresh_VIX';folder.mkdir(exist_ok=False);scores=[];anchor=np.load(OUT.parent/'m2_two_speed_001/fresh_check/blocks_two_speed.npz')
    for candidate in range(len(models)):
        idx=max(0,candidate-1);y=models[candidate,ACTIVE];r=v.predict(idx,y);score=save_vix(v,idx,y,f'c{candidate}',folder);score['candidate']=candidate;score['alpha']=0. if candidate==0 else float(ALPHA[idx]);scores.append(score);g=pd.read_csv(folder/f'VIX_forwards_c{candidate}.csv');g['combined_SE']=np.sqrt(g.difference_SE**2+g.anchor_SE**2);g.to_csv(folder/f'VIX_forwards_c{candidate}.csv',index=False);g=pd.read_csv(folder/f'VIX_prices_c{candidate}.csv');g['combined_SE_U']=np.sqrt(g.difference_SE_U**2+g.anchor_SE_U**2);g.to_csv(folder/f'VIX_prices_c{candidate}.csv',index=False)
        # Weighted-tail diagnostics for the new actual VIX variable.
        tails=[]
        for day,a in r['raw'].items():
            V=a['V'];w=a['weights'];K=float(q[q.expiry_days==day].strike.max());tail=w*np.maximum(V-K,0);total=tail.sum();tails.append(dict(day=int(day),mass=float(w.mean()),maximum_VIX=float(V.max()),weighted_tail_mean=float(tail.mean()),tail_ESS=float(total**2/np.sum(tail**2)) if total else 0.,top10_payoff_share=float(np.sort(tail)[-10:].sum()/total) if total else 0.))
        score['tails']=tails;print('FRESH VIX',candidate,score['alpha'],score['option_loss'],score['max_forward_gap'],flush=True)
    dump(folder/'scores.json',scores);pd.DataFrame([{k:v for k,v in row.items() if k!='tails'} for row in scores]).to_csv(folder/'scores.csv',index=False)
    # Same fresh states, small deterministic quadrature refinement diagnostic.
    subset={day:{key:value[np.concatenate([np.flatnonzero(a['nets']==j)[:16] for j in range(8)])] for key,value in a.items()} for day,a in draws.items()};aa=arr([0,1.5]);indices=[0,4];fine=exposures(k,dated,initial,p,norm,lambda t:Z(t)[indices],x,subset,aa,scale,48,12,OUT/'quadrature_probe');errors=[]
    for day in DAYS:
        ids=np.concatenate([np.flatnonzero(draws[day]['nets']==j)[:16] for j in range(8)]);old=H[day][indices][:,ids];new=fine[day];levels=(.05+.07*x[19])**2+(.2-(.05+.07*x[19])**2)*models[5,:16];coeff=levels-(.05+.07*x[19])**2;i0=(.05+.07*x[19])**2*(30/365)+old@coeff;i1=(.05+.07*x[19])**2*(30/365)+new@coeff;dv=100*(np.sqrt(i1/(30/365))-np.sqrt(i0/(30/365)));errors.append(dict(day=int(day),max_VIX_change_points_by_alpha=np.max(abs(dv),axis=1).tolist(),note='128 matched states; GH32/time8 versus GH48/time12; representative alpha1.5 curve, not uniform error bound'))
    dump(OUT/'quadrature_probe.json',errors);dump(folder/'run.json',dict(status='completed fresh-factor VIX checking for all strengths',wall_seconds=time.perf_counter()-begin,scores=scores,independent_of_training_factor_states=True,old_anchor_price_estimator_shared=True,forward_equality_not_refitted=True));print('FRESH VIX COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':main()
