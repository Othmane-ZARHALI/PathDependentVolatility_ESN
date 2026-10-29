"""Independent physical-parameter replay of constrained-forward/floor experiments."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .forward_floor_study import PARENT,CHECK_SEED
from .mixed_response_study import ZSHIFTS,ZPROB,seed_for
from .vix_only_curve_fit import PRIOR,ARCHIVE,ACTIVE,DAYS,DELTA,load_draws
from .postcap_vov_core import conditional_multiplier
from .parity_adaptive import gaussian_controls
from .full_readout_fit import dump
from .m2_spx import ROOT,frozen_stage3
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler


def audit(run):
    assert json.loads((run/'run.json').read_text())['status']=='completed';design=json.loads((run/'design.json').read_text())
    for group in ['source_sha256','input_sha256']:
        for path,h in design[group].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h,path
            if group=='source_sha256':assert hashlib.sha256((run/'source'/path).read_bytes()).hexdigest()==h,path
    old=json.loads((PRIOR/'run.json').read_text());p=np.array(old['selected']);ends=np.array(old['ends_years']);left=np.r_[0,ends[:-1]];controls=json.loads((PRIOR/'vix_controls.json').read_text());market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv');forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv');selection=json.loads((run/'selection.json').read_text());checks=[];negative_checks=[]
    assert len(market)==217 and set(market.expiry_days)==set(DAYS) and p[4]>0 and selection['frozen_before_fresh_check'];assert not design['SPX_priced'] and not design['baseline_replaced']
    for role,nodes,draws in [('training',dict(np.load(PARENT/'fresh_conditional_nodes.npz')),load_draws(PARENT/'fresh_cache')),('fresh_check',dict(np.load(run/'fresh_conditional_nodes.npz')),load_draws(run/'fresh_cache'))]:
        for name in selection['parameters']:
            model=json.loads((run/f'experimental_model_{name}.json').read_text());pp=np.array(model['parameters']);m=model['readout'];f=pp[-1];levels=pp[7:-1]
            np.testing.assert_array_equal(pp[:7],p[:7]);off=np.setdiff1d(np.arange(16),ACTIVE);np.testing.assert_array_equal(levels[off],p[7+off]);assert len(levels)==16 and np.all(levels>=f-1e-14) and max(levels)<=.2+1e-13;assert .05-1e-13<=np.sqrt(f)<=.12+1e-13
            if name=='tight_forwards':assert abs(f-p[-1])<1e-15
            stored=pd.read_csv(run/role/f'prices_otm_{name}.csv');common=pd.read_csv(run/role/f'prices_canonical_{name}.csv');fw=pd.read_csv(run/role/f'forwards_{name}.csv');score=json.loads((run/role/f'score_{name}.json').read_text());errors=[];put_rows=[]
            with np.load(run/role/f'paths_{name}.npz') as paths:
                for day,g in market.groupby('expiry_days',sort=True):
                    a=draws[day];rate=np.log(2)*365/m['half_life_days'];T=day/365;z=np.sqrt(-np.expm1(-2*rate*T))*a['z_normals'];I=np.full(len(z),f*DELTA)
                    for t,w,oldratio in zip(nodes[f'd{day}_times'],nodes[f'd{day}_weights'],nodes[f'd{day}_ratios']):
                        j=min(np.searchsorted(ends,t,side='right'),15);factor=conditional_multiplier(z,T,t-T,m['b'],rate);Q=(1-m['weight'])*oldratio+m['weight']*factor;I+=w*(levels[j]-f)*Q
                    V=100*np.sqrt(I/DELTA);np.testing.assert_allclose(V,paths[f'd{day}_V'],rtol=1e-13,atol=2e-9);assert min(V)>=100*np.sqrt(f)-1e-10
                    c=gaussian_controls(a['normals'],a['weights']);fixed=controls[str(day)];b=np.array([np.interp(np.log(g.strike),np.log(fixed['anchors']),row) for row in fixed['beta']]);F=a['weights']*V-c@np.array(fixed['forward_beta']);C=a['weights'][:,None]*np.maximum(V[:,None]-g.strike.to_numpy(),0)-c@b;P=C-g.option_type.eq('put').to_numpy()[None,:]*(F[:,None]-g.strike.to_numpy())
                    errors.append(float(np.max(abs(stored.loc[g.index,'model_price_U'].to_numpy()-P.mean(0)))));assert np.max(abs(common.loc[g.index,'model_price_U'].to_numpy()-C.mean(0)))<1e-9;assert abs(fw[fw.expiry_days.eq(day)].model_forward.iloc[0]-F.mean())<1e-9
                    lengths=np.maximum(0,np.minimum(ends,T+DELTA)-np.maximum(left,T));exact=10000/DELTA*(lengths@levels);mr=next(x for x in score['moments'] if x['day']==day);assert abs(exact-mr['exact_second_moment'])<1e-9
                    puts=g.option_type.eq('put').to_numpy();ce=common.loc[g.index,'model_price_U'].to_numpy()-common.loc[g.index,'price_mid_U'].to_numpy();pe=stored.loc[g.index,'model_price_U'].to_numpy()-stored.loc[g.index,'price_mid_U'].to_numpy();fg=fw[fw.expiry_days.eq(day)].iloc[0];np.testing.assert_allclose(ce[puts],pe[puts]+fg.model_forward-fg.forward_vix_points,atol=1e-10,rtol=1e-10)
                    # Positive raw put payoffs are an independent diagnostic; retain signed
                    # controlled estimates and the original calibration objective unchanged.
                    raw_put=a['weights'][:,None]*np.maximum(g.strike.to_numpy()[None,:]-V[:,None],0)
                    for j,(index,quote) in enumerate(g.iterrows()):
                        if quote.option_type!='put':continue
                        blocks=raw_put[:,j].reshape(8,-1).mean(1);correction=P[:,j]-raw_put[:,j];cb=correction.reshape(8,-1).mean(1)
                        np.testing.assert_allclose(blocks.mean()+cb.mean(),stored.loc[index,'model_price_U'],atol=1e-10)
                        put_rows.append(dict(day=int(day),strike=float(quote.strike),market_put=float(quote.price_mid_U),controlled_put=float(stored.loc[index,'model_price_U']),controlled_SE=float(stored.loc[index,'mc_se_U']),raw_put=float(blocks.mean()),raw_put_SE=float(blocks.std(ddof=1)/np.sqrt(8)),control_correction=float(cb.mean()),control_correction_SE=float(cb.std(ddof=1)/np.sqrt(8))))
                        if stored.loc[index,'model_price_U']<0:negative_checks.append(dict(role=role,candidate=name,day=int(day),strike=float(quote.strike),controlled_price=float(stored.loc[index,'model_price_U']),SE=float(stored.loc[index,'mc_se_U']),negative_magnitude_in_SE=float(-stored.loc[index,'model_price_U']/stored.loc[index,'mc_se_U']),raw_put=float(blocks.mean())))
            assert max(errors)<1e-9
            if role=='training' and name!='previous_mixture':assert score['maximum_forward_error']<=1e-5
            checks.append(dict(role=role,candidate=name,maximum_native_price_error=max(errors),maximum_forward_error=score['maximum_forward_error']))
            pd.DataFrame(put_rows).to_csv(run/role/f'raw_native_puts_{name}.csv',index=False)
    fresh=load_draws(run/'fresh_cache');props=json.loads((run/'fresh_cache/design.json').read_text())['old_proposals'];zm=DefensiveGaussianMixture(ZSHIFTS[:,None],ZPROB);architecture,_,premium,_=frozen_stage3();dated=DatedM2Pricer(architecture,premium);initial=np.array(old['initial_state']);regenerated=0
    for day in DAYS:
        a=fresh[day];dec=props[str(day)];om=DefensiveGaussianMixture(dec['shifts'],dec['probabilities']);np.testing.assert_allclose(om.weights(a['normals'])*zm.weights(a['z_normals'][:,None]),a['weights'],atol=1e-14,rtol=3e-13);joint=DefensiveGaussianMixture(np.array([np.r_[x,z] for x in dec['shifts'] for z in ZSHIFTS]),np.array([p*q for p in dec['probabilities'] for q in ZPROB]));sampler=GaussianOuterSampler(dated,day/365,initial)
        for net in range(8):
            z,w,c=joint.draw(12,seed_for(CHECK_SEED,day,net));ids=a['nets']==net;np.testing.assert_array_equal(z[:,:-1],a['normals'][ids]);np.testing.assert_array_equal(z[:,-1],a['z_normals'][ids]);np.testing.assert_array_equal(w,a['weights'][ids]);np.testing.assert_array_equal(c,a['components'][ids]);np.testing.assert_allclose(sampler.mean+z[:,:-1]@sampler.root.T,a['states'][ids],rtol=1e-13,atol=1e-13);regenerated+=1
    dump(run/'negative_price_diagnostics.json',negative_checks)
    result=dict(status='PASS',source_files=len(design['source_sha256']),input_files=len(design['input_sha256']),price_replays=checks,fresh_scramble_regenerations=regenerated,training_constraints_checked=True,floor_support_checked=True,all_16_curve_levels_preserved=True,off_window_levels_unchanged=True,positive_beta_J=True,SPX_priced=False,baseline_replaced=False,negative_estimates_retained=len(negative_checks),maximum_negative_magnitude_in_SE=max((x['negative_magnitude_in_SE'] for x in negative_checks),default=0.),limitations='Training constraints and finite numerical replay; no exact population-forward, convergence, global optimum or Lean certificate.')
    dump(run/'audit.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();audit(a.run)
