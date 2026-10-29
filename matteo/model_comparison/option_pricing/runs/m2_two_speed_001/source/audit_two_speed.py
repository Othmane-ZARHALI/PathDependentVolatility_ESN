"""Independent scalar checks and weighted-tail diagnostics for the two-speed pilot."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.special import roots_hermitenorm,logsumexp
from .full_readout_fit import HERE,dump
from .two_speed_core import metadata
from .forward_floor_core import unpack
from .vix_only_curve_fit import PRIOR,DAYS,DELTA
from .mixed_response_study import ZSHIFTS,ZPROB
from esn_eq.importance_sampling import DefensiveGaussianMixture
from .two_speed_study import PREVIOUS,augmented


def audit(run):
    design=json.loads((run/'design.json').read_text())
    for key in ['source_hashes','input_hashes']:
        for fn,sha in design[key].items():assert hashlib.sha256(Path(fn).read_bytes()).hexdigest()==sha,fn
    prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);sel=json.loads((run/'selection.json').read_text());x=np.array(sel['parameters']);m=metadata(x,p);levels,floor,_=unpack(x,p);assert p[4]>0 and np.all(levels>=floor)
    eps,ghw=roots_hermitenorm(80);ghw/=np.sqrt(2*np.pi);rows=[];tails=[]
    nodes=dict(np.load(PREVIOUS/'fresh_base_nodes.npz'));fast_design=json.loads((run/'checking_fast/design.json').read_text());mixture=DefensiveGaussianMixture(ZSHIFTS[:,None],ZPROB)
    for day in DAYS:
        T=day/365;times=nodes[f'd{day}_times'];tw=nodes[f'd{day}_weights'];ratios=nodes[f'd{day}_ratios'];seg=np.minimum(np.searchsorted(ends,times,side='right'),len(ends)-1)
        base=dict(np.load(PREVIOUS/'fresh_base_cache/production'/f'vix_d{day}.npz'));fast=dict(np.load(run/'checking_fast'/f'd{day}.npz'))
        zf=fast['fast_normals'];zs=base['z_normals']
        for rec in fast_design['seeds']:
            if rec['day']!=day:continue
            ids_net=np.flatnonzero(base['nets']==rec['net']);zz,ww,_=mixture.draw(int(np.log2(rec['count'])),rec['seed']);perm=np.random.default_rng(rec['permutation_seed']).permutation(rec['count'])
            assert np.array_equal(zf[ids_net],zz[perm,0]);assert np.max(abs(fast['fast_weight'][ids_net]-ww[perm]))<1e-14
        weight=base['weights']*fast['fast_weight']
        # Independent likelihood evaluation in log space, no sample renormalization.
        wf=np.exp(-logsumexp(zf[:,None]*ZSHIFTS[None,:]-.5*ZSHIFTS[None,:]**2+np.log(ZPROB)[None,:],axis=1));assert np.max(abs(wf-fast['fast_weight']))<1e-14
        with np.load(run/'fresh_check/paths_two_speed.npz') as a:V=a[f'd{day}_V'];savedw=a[f'd{day}_weights']
        assert np.max(abs(savedw-weight))<1e-14
        q=pd.read_csv(run/'fresh_check/prices_canonical_two_speed.csv').query('expiry_days == @day');K=q.strike.max();pay=weight*np.maximum(V-K,0);ids=np.unique(np.r_[0,np.argsort(pay)[-24:],np.argmax(zs),np.argmax(zf),np.argmin(zs),np.argmin(zf)])
        coeff=levels[seg]-floor;old=floor*DELTA+(1-x[8])*(coeff*tw)@ratios[:,ids];integral=old.copy();condmax=0.
        for b,h,z,mix in [(m['b_slow'],m['slow_half_days'],zs,x[8]*(1-x[13])),(m['b_fast'],m['fast_half_days'],zf,x[8]*x[13])]:
            a=np.log(2)*365/h;qT=-np.expm1(-2*a*T)
            def cond(s,zvalue):
                load=b*np.exp(-a*(s-T))*np.sqrt(qT)
                return np.exp(load*zvalue-.5*load*load)
            for lag in [0.,15/365,30/365]:
                s=T+lag;qs=-np.expm1(-2*a*s);qinc=-np.expm1(-2*a*lag)
                for z0 in [-1.,0.,2.]:
                    direct=ghw@np.exp(b*(np.exp(-a*lag)*np.sqrt(qT)*z0+np.sqrt(qinc)*eps)-.5*b*b*qs)
                    oracle=cond(s,z0);condmax=max(condmax,abs(direct-oracle)/max(1,abs(oracle)))
            for j,idx in enumerate(ids):
                total=0.
                for k,end in enumerate(ends):
                    lo=max(T,0 if k==0 else ends[k-1]);hi=min(T+DELTA,end)
                    if hi>lo:total+=(levels[k]-floor)*quad(lambda s:cond(s,z[idx]),lo,hi,epsabs=1e-11,epsrel=1e-11)[0]
                integral[j]+=mix*total
        vr=100*np.sqrt(integral/DELTA);err=abs(vr-V[ids]);assert np.max(err/np.maximum(1,abs(vr)))<1e-8,(day,err.max())
        assert condmax<1e-10
        rows.append(dict(day=day,probe_count=len(ids),stress_time_quad_max_VIX_error=float(err.max()),stress_time_quad_max_relative_error=float(np.max(err/np.maximum(1,abs(vr)))),conditional_GH80_max_error=float(condmax)))
        total=pay.sum();top=np.sort(pay)[::-1];mass=np.array([weight[base['nets']==j].mean() for j in range(8)])
        tails.append(dict(day=day,mass=float(mass.mean()),mass_se=float(mass.std(ddof=1)/np.sqrt(8)),last_strike=float(K),weighted_call=float(pay.mean()),payoff_ESS=float(total**2/(pay@pay)) if total else 0,top10_fraction=float(top[:10].sum()/total) if total else 0,minimum_sampled_VIX=float(V.min()),maximum_sampled_VIX=float(V.max())))
        assert V.min()>=100*np.sqrt(floor)-1e-9
    result=dict(status='PASS',source_and_input_hashes=True,independent_likelihood_reconstruction=True,all_32_checking_fast_seed_replays=True,conditional_and_stress_quadrature=rows,weighted_tail_checks=tails,limitations=['M2 conditional tables reused; their previous audits are not repeated.','Scalar stress checks are sampled numerical checks, not uniform error bounds.','No Lean proof or population-forward certificate.'])
    dump(run/'audit.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();audit(a.run)
