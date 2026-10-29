"""Independent nodal replay, exact proposal weights, nesting and curve preservation."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .mixed_response_study import training_inputs,ZSHIFTS,ZPROB,CHECK_SEED,seed_for
from .vix_only_curve_fit import PRIOR,ARCHIVE,ACTIVE,DAYS,DELTA,load_draws
from .postcap_vov_core import conditional_multiplier
from .parity_adaptive import gaussian_controls,canonical_market
from .cap_controlled_core import shape_metrics
from .full_readout_fit import dump
from .m2_spx import ROOT,frozen_stage3
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler


def audit(run):
    assert json.loads((run/'run.json').read_text())['status']=='completed'
    design=json.loads((run/'design.json').read_text());checks=[]
    for group in ['source_sha256','input_sha256']:
        for path,h in design[group].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h,path
            if group=='source_sha256':assert hashlib.sha256((run/'source'/path).read_bytes()).hexdigest()==h,path
    d=json.loads((PRIOR/'run.json').read_text());p=np.array(d['selected']);ends=np.array(d['ends_years']);left=np.r_[0,ends[:-1]]
    selected=json.loads((run/'selection.json').read_text());controls=json.loads((PRIOR/'vix_controls.json').read_text());market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv')
    assert len(market)==217 and set(market.expiry_days)==set(DAYS) and p[4]>0
    assert selected['frozen_before_check'] and not design['SPX_priced'] and not design['baseline_replaced']
    # Recreate training in a separate audit directory; never overwrite the run's design.
    folder=run/'audit_replay';folder.mkdir(exist_ok=True);tn,td=training_inputs(folder)
    assert json.loads((folder/'training_design.json').read_text())==json.loads((run/'training_design.json').read_text())
    for role,nodes,draws in [('training',tn,td),('fresh_check',dict(np.load(run/'fresh_conditional_nodes.npz')),load_draws(run/'fresh_cache'))]:
        for name,candidate in selected['candidates'].items():
            model=json.loads((run/f'experimental_model_{name}.json').read_text());pp=np.array(model['parameters']);m=model['readout'];family=model['family']
            np.testing.assert_array_equal(pp[:7],p[:7]);assert pp[-1]==p[-1];assert len(pp[7:-1])==16
            off=np.setdiff1d(np.arange(16),ACTIVE);np.testing.assert_array_equal(pp[7+off],p[7+off]);assert np.all(pp[7:-1]>=p[-1]) and max(pp[7:-1])<=.2
            assert 0<=m['weight']<=1 and 0<=m['b1']<=m['b2']<=5 and 1-1e-12<=m['half_life_days']<=365+1e-9
            stored=pd.read_csv(run/role/f'prices_otm_{name}.csv');fw=pd.read_csv(run/role/f'forwards_{name}.csv');score=json.loads((run/role/f'score_{name}.json').read_text());errors=[]
            with np.load(run/role/f'paths_{name}.npz') as paths:
                for day,g in market.groupby('expiry_days',sort=True):
                    a=draws[day];rate=np.log(2)*365/m['half_life_days'];T=day/365;z=np.sqrt(-np.expm1(-2*rate*T))*a['z_normals'];I=np.full(len(z),p[-1]*DELTA)
                    for t,w,old in zip(nodes[f'd{day}_times'],nodes[f'd{day}_weights'],nodes[f'd{day}_ratios']):
                        j=min(np.searchsorted(ends,t,side='right'),15)
                        L2=conditional_multiplier(z,T,t-T,m['b2'],rate)
                        if family=='product':
                            L1=conditional_multiplier(z,T,t-T,m['b1'],rate);Q=old*((1-m['weight'])*L1+m['weight']*L2)
                        else:Q=(1-m['weight'])*old+m['weight']*L2
                        I+=w*(pp[7+j]-p[-1])*Q
                    V=100*np.sqrt(I/DELTA);np.testing.assert_allclose(V,paths[f'd{day}_V'],rtol=1e-13,atol=2e-9)
                    c=gaussian_controls(a['normals'],a['weights']);fixed=controls[str(day)];b=np.array([np.interp(np.log(g.strike),np.log(fixed['anchors']),row) for row in fixed['beta']])
                    forward=a['weights']*V-c@np.array(fixed['forward_beta']);call=a['weights'][:,None]*np.maximum(V[:,None]-g.strike.to_numpy(),0)-c@b
                    native=call-g.option_type.eq('put').to_numpy()[None,:]*(forward[:,None]-g.strike.to_numpy())
                    errors.append(float(np.max(abs(stored.loc[g.index,'model_price_U'].to_numpy()-native.mean(0)))))
                    assert abs(fw.loc[fw.expiry_days.eq(day),'model_forward'].iloc[0]-forward.mean())<1e-9
                    lengths=np.maximum(0,np.minimum(ends,T+DELTA)-np.maximum(left,T));np.testing.assert_allclose(lengths.sum(),DELTA,rtol=1e-14)
                    exact=10000/DELTA*(lengths@pp[7:-1]);row=next(r for r in score['moments'] if r['day']==day);assert abs(exact-row['exact_second_moment'])<1e-9
            assert max(errors)<1e-9;checks.append(dict(role=role,candidate=name,maximum_price_error=max(errors)))
    # Proposal factorization and independent complete-scramble regeneration.
    fresh=load_draws(run/'fresh_cache');props=json.loads((run/'fresh_cache/design.json').read_text())['old_proposals'];zm=DefensiveGaussianMixture(ZSHIFTS[:,None],ZPROB);seeds=[]
    architecture,_,premium,_=frozen_stage3();dated=DatedM2Pricer(architecture,premium);initial=np.array(d['initial_state'])
    for day in DAYS:
        a=fresh[day];dec=props[str(day)];old=DefensiveGaussianMixture(dec['shifts'],dec['probabilities']);weights=old.weights(a['normals'])*zm.weights(a['z_normals'][:,None])
        np.testing.assert_allclose(weights,a['weights'],rtol=3e-13,atol=1e-14)
        joint=DefensiveGaussianMixture(np.array([np.r_[x,z] for x in dec['shifts'] for z in ZSHIFTS]),np.array([p*q for p in dec['probabilities'] for q in ZPROB]));sampler=GaussianOuterSampler(dated,day/365,initial)
        for net in range(8):
            z,w,c=joint.draw(11,seed_for(CHECK_SEED,day,net));ids=a['nets']==net
            np.testing.assert_array_equal(z[:,:-1],a['normals'][ids]);np.testing.assert_array_equal(z[:,-1],a['z_normals'][ids]);np.testing.assert_array_equal(w,a['weights'][ids]);np.testing.assert_array_equal(c,a['components'][ids]);np.testing.assert_allclose(sampler.mean+z[:,:-1]@sampler.root.T,a['states'][ids],atol=1e-13,rtol=1e-13);seeds.append(dict(day=day,net=net))
    # Paired wing sensitivity across eight independent complete scrambles.
    quotes=canonical_market(market);shape_full={};shape_jk={}
    for name in selected['candidates']:
        a=dict(np.load(run/'fresh_check'/f'blocks_{name}.npz'));shape_full[name]=shape_metrics(quotes,a['call'].mean(0))
        shape_jk[name]=[shape_metrics(quotes,a['call'][np.arange(8)!=j].mean(0)) for j in range(8)]
    pairs=[]
    for name in selected['candidates']:
        if name in ['baseline','curve_only','previous_single_factor']:continue
        for ref in ['curve_only','previous_single_factor']:
            for k,day in enumerate(DAYS):
                for metric in ['wing_slope','wing_curvature']:
                    change=shape_full[name]['rows'][k][metric]-shape_full[ref]['rows'][k][metric]
                    values=[shape_jk[name][j]['rows'][k][metric]-shape_jk[ref][j]['rows'][k][metric] for j in range(8)]
                    pairs.append(dict(candidate=name,reference=ref,day=day,metric=metric,change=change,paired_jackknife_se=float(np.std(values,ddof=0)*np.sqrt(7))))
    dump(run/'fresh_check/paired_wing_shape.json',pairs)
    result=dict(status='PASS',source_files=len(design['source_sha256']),input_files=len(design['input_sha256']),price_replays=checks,fresh_scramble_regenerations=len(seeds),proposal_weights_factorization_checked=True,all_16_curve_pieces_preserved=True,off_window_levels_unchanged=True,positive_beta_J=True,SPX_priced=False,baseline_replaced=False,limitations='Independent implementation/numerical checks; not a population-convergence, global-optimum or Lean certificate.')
    dump(run/'audit.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();audit(a.run)
