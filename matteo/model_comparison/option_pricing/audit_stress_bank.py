"""Independent scalar and likelihood checks for the geometric stress bank."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.special import logsumexp,roots_hermitenorm
from .full_readout_fit import dump
from .vix_only_curve_fit import PRIOR,DAYS,DELTA
from .stress_bank_core import HALF_LIVES,unpack,metadata,weights
from .stress_bank_study import PREVIOUS,OLD_CHECK,EXTRA_SHIFTS,EXTRA_PROB
from esn_eq.importance_sampling import DefensiveGaussianMixture


def audit(run):
    design=json.loads((run/'design.json').read_text())
    for key in ['source_hashes','input_hashes']:
        for name,sha in design[key].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha,name
    prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);x=np.array(json.loads((run/'selection.json').read_text())['parameters']);levels,floor,_=unpack(x,p);alpha,m2,_=weights(x)
    assert p[4]>0 and np.all(levels>=floor) and np.all(alpha>=0) and m2>=0 and abs(m2+alpha.sum()-1)<1e-14
    nodes=dict(np.load(OLD_CHECK/'fresh_base_nodes.npz'));fd=json.loads((run/'checking_bank/design.json').read_text());mixture=DefensiveGaussianMixture(EXTRA_SHIFTS,EXTRA_PROB)
    gh,gw=roots_hermitenorm(80);gw/=np.sqrt(2*np.pi);rows=[];tails=[]
    for day in DAYS:
        T=day/365;bn=dict(np.load(run/'checking_bank'/f'd{day}.npz'));base=dict(np.load(OLD_CHECK/'fresh_base_cache/production'/f'vix_d{day}.npz'));fast=dict(np.load(PREVIOUS/'checking_fast'/f'd{day}.npz'));z=bn['bank_normals']
        assert np.array_equal(z[:,2],fast['fast_normals']) and np.array_equal(z[:,3],base['z_normals'])
        for rec in fd['seeds']:
            if rec['day']!=day:continue
            ids=np.flatnonzero(base['nets']==rec['net']);zz,ww,_=mixture.draw(int(np.log2(rec['count'])),rec['seed']);perm=np.random.default_rng(rec['permutation_seed']).permutation(rec['count'])
            assert np.array_equal(z[ids,:2],zz[perm]);assert np.array_equal(bn['extra_weights'][ids],ww[perm])
        lograt=z[:,:2]@EXTRA_SHIFTS.T-.5*np.sum(EXTRA_SHIFTS**2,axis=1);ew=np.exp(-logsumexp(lograt+np.log(EXTRA_PROB),axis=1));w=base['weights']*fast['fast_weight']*ew
        assert np.max(abs(w-bn['combined_weights']))<1e-13
        with np.load(run/'fresh_check/paths_bank.npz') as a:V=a[f'd{day}_V'];savedw=a[f'd{day}_weights']
        assert np.max(abs(w-savedw))<1e-13
        q=pd.read_csv(run/'fresh_check/prices_canonical_bank.csv').query('expiry_days == @day');K=q.strike.max();pay=w*np.maximum(V-K,0);ids=np.unique(np.r_[0,np.argsort(pay)[-24:],np.argmax(z,axis=0),np.argmin(z,axis=0),np.argmax(V)])
        times=nodes[f'd{day}_times'];tw=nodes[f'd{day}_weights'];ratios=nodes[f'd{day}_ratios'];seg=np.minimum(np.searchsorted(ends,times,side='right'),len(ends)-1)
        integral=floor*DELTA+m2*((levels[seg]-floor)*tw)@ratios[:,ids];conditional_error=0.;normalization_error=0.
        for j,h in enumerate(HALF_LIVES):
            a=np.log(2)*365/h;b=5*x[8+j];qT=-np.expm1(-2*a*T)
            normalization_error=max(normalization_error,abs(gw@np.exp(b*np.sqrt(qT)*gh-.5*b*b*qT)-1))
            def response(s,zvalue):
                load=b*np.exp(-a*(s-T))*np.sqrt(qT)
                return np.exp(load*zvalue-.5*load*load)
            for lag in [0.,15/365,30/365]:
                s=T+lag;qs=-np.expm1(-2*a*s);qi=-np.expm1(-2*a*lag)
                for z0 in [-1.,0.,2.]:
                    numeric=gw@np.exp(b*(np.exp(-a*lag)*np.sqrt(qT)*z0+np.sqrt(qi)*gh)-.5*b*b*qs);exact=response(s,z0)
                    conditional_error=max(conditional_error,abs(numeric-exact)/max(1,abs(exact)))
            for k,idx in enumerate(ids):
                total=0.
                for segj,end in enumerate(ends):
                    lo=max(T,0 if segj==0 else ends[segj-1]);hi=min(T+DELTA,end)
                    if hi>lo:total+=(levels[segj]-floor)*quad(lambda s:response(s,z[idx,j]),lo,hi,epsabs=1e-11,epsrel=1e-11)[0]
                integral[k]+=alpha[j]*total
        vr=100*np.sqrt(integral/DELTA);err=abs(vr-V[ids]);rel=float(max(err/np.maximum(1,abs(vr))));assert rel<1e-8,(day,rel);assert normalization_error<1e-10 and conditional_error<1e-10
        rows.append(dict(day=day,probe_count=len(ids),scalar_replay_max_abs_VIX_error=float(err.max()),scalar_replay_max_relative_error=rel,conditional_GH80_error=float(conditional_error),normalization_GH80_error=float(normalization_error)))
        total=pay.sum();top=np.sort(pay)[::-1];mass=np.array([w[base['nets']==k].mean() for k in range(8)]);moment=np.array([(w*V*V)[base['nets']==k].mean() for k in range(8)])
        lengths=np.maximum(0,np.minimum(ends,T+DELTA)-np.maximum(np.r_[0,ends[:-1]],T));exactmoment=10000/DELTA*(lengths@levels)
        tails.append(dict(day=day,mass=float(mass.mean()),mass_se=float(mass.std(ddof=1)/np.sqrt(8)),last_strike=float(K),raw_call=float(pay.mean()),payoff_ESS=float(total**2/(pay@pay)) if total else 0.,top10_fraction=float(top[:10].sum()/total) if total else 0.,sample_second_moment=float(moment.mean()),exact_second_moment=float(exactmoment),second_moment_se=float(moment.std(ddof=1)/np.sqrt(8)),minimum_VIX=float(V.min()),maximum_VIX=float(V.max())))
        assert V.min()>=100*np.sqrt(floor)-1e-9
    result=dict(status='PASS',source_and_input_hashes=True,all_32_extra_seed_replays=True,independent_likelihood_reconstruction=True,positive_spike_shift=True,mixture_normalization=float(m2+alpha.sum()),scalar_checks=rows,weighted_tails=tails,limitations=['M2 conditional tables reused with their previous audits.','Checking reuses old/node2/node3 coordinates.','Scalar probes and eight-block SEs do not certify population-forward equality or uniform error.','No new Lean theorem.'])
    dump(run/'audit.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);args=ap.parse_args();audit(args.run)
