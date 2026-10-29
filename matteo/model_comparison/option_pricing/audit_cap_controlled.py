"""Saved-draw, four-window fit and genuinely new sampler correspondence audit."""
import json,hashlib,ctypes as ct
from pathlib import Path
import numpy as np,pandas as pd
from .full_readout_fit import HERE,dump,Kernel,arr
from .m2_spx import ROOT,frozen_stage3
from .cap_controlled_core import window_basis,adjusted,shape_metrics
from .audit_maturity_curve import replay
from .audit_cap_tail import paired
from .parity_adaptive import gaussian_controls,canonical_market
from .maturity_curve import vix_only_loss
from esn_eq.pricing import DatedM2Pricer,_normal_rule
from esn_eq.states import symmetric_square_root
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture


def exposure_prices(folder,name,p,cache,vix,vc):
    calls=[];native=[];fw=[]
    with np.load(folder/f'exposures_{name}.npz') as expo:
        for day,g in vix.groupby('expiry_days',sort=True):
            a=np.load(cache/'production'/f'vix_d{day}.npz');v=100*np.sqrt(p[-1]+expo[f'd{day}']@(p[7:-1]-p[-1])/(30/365));w=a['weights'];c=gaussian_controls(a['normals'],w);record=vc[str(day)]
            beta=np.stack([np.interp(np.log(g.strike),np.log(record['anchors']),b) for b in record['beta']]);y=w[:,None]*np.maximum(v[:,None]-g.strike.to_numpy(),0)-c@beta
            fwd=float(np.mean(w*v-c@np.array(record['forward_beta'])));call=y.mean(0);calls.extend(call);native.extend(call-g.option_type.eq('put').to_numpy()*(fwd-g.strike.to_numpy()));fw.append(fwd)
    return np.array(calls),np.array(native),np.array(fw)

def run(out):
    r=json.loads((out/'run.json').read_text());d=json.loads((out/'design.json').read_text());base=np.array(d['baseline']);prior=Path(d['parent']);cache=Path(d['cache']);pd0=json.loads((prior/'run.json').read_text());B=window_basis(pd0['ends_days']);np.testing.assert_array_equal(B,d['basis']);spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv');q=pd.concat([spx,vix],ignore_index=True);vc=json.loads((prior/'vix_controls.json').read_text());canonical=canonical_market(vix);native_loss=vix_only_loss(vix,f);common_loss=vix_only_loss(canonical,f);records=[];histories=0
    replay(out,'baseline',base[7:-1],vix)
    profiles=json.loads((out/'profiles.json').read_text());baseline=json.loads((out/'baseline_replay.json').read_text())
    for row in profiles:
        x=np.array(row['parameters']);p=base.copy();p[5]+=row['cap_shift'];xx=adjusted(p,row['fit']['adjustments'],B);np.testing.assert_array_equal(x,xx);np.testing.assert_array_equal(x[7:-1][B.sum(1)==0],base[7:-1][B.sum(1)==0]);assert x[-1]==base[-1] and x[4]>0 and np.all(x[7:-1]>=x[-1])
        rr=replay(out,row['name'],x[7:-1],vix);records.append(dict(name=row['name'],role='training profile',native=rr['original_otm_total']))
        saved=pd.read_csv(out/f"prices_canonical_{row['name']}.csv");m=shape_metrics(saved,saved.model_price_U.to_numpy());np.testing.assert_allclose([m['wing_MAE'],m['ATM_MAE']],[row['shape']['wing_MAE'],row['shape']['ATM_MAE']],rtol=1e-11)
        for h in row['fit']['history']:
            xx=adjusted(p,h['adjustments'],B);call,native,forward=exposure_prices(out,row['name'],xx,cache,vix,vc);a=native_loss.evaluate(native,forward,xx[7:-1]);b=common_loss.evaluate(call,forward,xx[7:-1]);value=lambda z:z['components']['VIX']+z['components']['forwards'];np.testing.assert_allclose([value(a),value(b)],[h['native_loss'],h['common_loss']],rtol=1e-10,atol=1e-9);np.testing.assert_allclose(forward,h['forwards'],rtol=1e-12);histories+=1
        preferred=row['fit']['feasible'] and row['fit']['native_loss']<baseline['native_vix_plus_forward'] and m['ATM_MAE']<=baseline['shape']['ATM_MAE']+5 and m['wing_complete'];assert preferred==row['preferred']
    selected=json.loads((out/'selection.json').read_text())['candidates'];assert selected==r['selection']
    for fam in ['capped','continued']:
        pool=[z for z in profiles if (z['continuation']==0)==(fam=='capped') and z['fit']['feasible']];pref=[z for z in pool if z['preferred']]
        if pool:assert next(s['name'] for s in selected if s['family']==fam)==min(pref or pool,key=lambda z:z['shape']['wing_MAE'])['name']
    for row in selected:
        x=np.array(row['parameters']);rr=replay(out/'full',row['name'],x[7:-1],q);records.append(dict(name=row['name'],role='actual full SPX/VIX',native=rr['original_otm_total']));dump(out/'full'/f"paired_{row['name']}.json",paired(out/'full',row['name'],prior,'selected',q,f,x[7:-1],base[7:-1]))
        full=pd.read_csv(out/'full'/f"prices_canonical_{row['name']}.csv");train=pd.read_csv(out/f"prices_canonical_{row['name']}.csv");np.testing.assert_allclose(full[full.instrument.eq('VIX')].model_price_U,train.model_price_U,atol=2e-12,rtol=1e-12)
    for row in [dict(name='baseline',parameters=base.tolist()),*selected]:
        x=np.array(row['parameters']);rr=replay(out/'fresh_check',row['name'],x[7:-1],vix);records.append(dict(name=row['name'],role='genuinely fresh VIX checking',native=rr['original_otm_total']))
        if row['name']!='baseline':dump(out/'fresh_check'/f"paired_{row['name']}.json",paired(out/'fresh_check',row['name'],out/'fresh_check','baseline',vix,f,x[7:-1],base[7:-1]))
    # Regenerate each entire NEW Sobol mixture draw and exact terminal-state map.
    fresh=out/'fresh_cache';fd=json.loads((fresh/'design.json').read_text());old=set()
    for path in (HERE/'runs').rglob('vix_d*_seeds.json'):
        if out in path.parents:continue
        for z in json.loads(path.read_text()):old.add(int(z['seed']))
    old.update(pd.read_csv(HERE/'runs/m2_cap_tail_001/marginal_response.csv').seed.astype(int));allseeds=[];samplers=[];a,_,prem,_=frozen_stage3();dated=DatedM2Pricer(a,prem);initial=np.array(pd0['initial_state'])
    for day in [27,55,118,153]:
        dec=fd['proposal'][str(day)];mixture=DefensiveGaussianMixture(dec['shifts'],dec['probabilities']);sampler=GaussianOuterSampler(dated,day/365,initial);path=fresh/'production'/f'vix_d{day}.npz';saved=np.load(path)
        for net,record in enumerate(json.loads((fresh/'production'/f'vix_d{day}_seeds.json').read_text())):
            seed=int(np.random.SeedSequence([fd['seed_base'],2,day,net]).generate_state(1)[0]);assert record['seed']==seed and seed not in old and seed not in allseeds;allseeds.append(seed);z,w,c=mixture.draw(9,seed);ix=slice(net*512,(net+1)*512);np.testing.assert_array_equal(z,saved['normals'][ix]);np.testing.assert_array_equal(w,saved['weights'][ix]);np.testing.assert_array_equal(c,saved['components'][ix]);np.testing.assert_array_equal(sampler.mean+z@sampler.root.T,saved['states'][ix])
        entry=next(z for z in fd['files'] if z['days']==day);assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256'];samplers.append(dict(day=day,exact_complete_scramble_replays=8))
    direct_checks=[];kernel=Kernel(out/'cap.dylib');kernel.lib.configure.argtypes=[ct.c_int,ct.c_double]
    for row in selected:
        p=np.array(row['parameters']);kernel.lib.configure(0 if row['continuation']==0 else 2,row['continuation'])
        for day in [27,153]:
            saved=np.load(fresh/'production'/f'vix_d{day}.npz');paths=np.load(out/'fresh_check'/f"paths_{row['name']}.npz");v=paths[f'vix_{day}'];K=float(vix[vix.expiry_days.eq(day)].strike.max());ids=np.unique([0,1,int(np.argmax(v)),int(np.argmax(saved['weights']*np.maximum(v-K,0)))])
            for lag in [1/max(dated.layout.rates),15/365,30/365]:
                mean,c=dated.gaussian_projection_law(lag,saved['states']);z,w=_normal_rule(32);offset=arr(z@symmetric_square_root(c).T);pred=kernel.expectation(mean,offset,arr(w),p)[ids,0];direct=kernel.expectation(mean[ids],offset,arr(w),p)[:,0];error=float(np.max(abs(pred/direct-1)));assert error<2e-7;direct_checks.append(dict(candidate=row['name'],day=day,lag_days=float(lag*365),relative_error=error))
    for rel,sha in d['source_sha256'].items():assert hashlib.sha256((out/'source'/rel).read_bytes()).hexdigest()==sha
    resume_checks=0
    if (out/'resume/design.json').exists():
        rd=json.loads((out/'resume/design.json').read_text())
        for rel,sha in rd['source_sha256'].items():assert hashlib.sha256((out/'resume/source'/rel).read_bytes()).hexdigest()==sha;resume_checks+=1
    frozen=out/'source'/str((prior/'selected_model.json').relative_to(ROOT));assert frozen.read_bytes()==(prior/'selected_model.json').read_bytes()
    result=dict(status='PASS',saved_score_records=len(records)+1,records=records,four_window_fit_trials_replayed=histories,frozen_source_checks=len(d['source_sha256'])+resume_checks,exact_fresh_scramble_replays=len(allseeds),fresh_seeds_disjoint_from_recorded_history=True,historical_seeds_checked=len(old),samplers=samplers,direct_conditional_GH32_checks=direct_checks,parent_selection_unchanged=True,limitation='Numerical correspondence and saved evidence replay; not continuous-price convergence, market OOS, a new Lean certificate or physical compatibility.')
    dump(out/'audit.json',result);print('PASS',len(records)+1,'scores',histories,'fit trials',len(allseeds),'exact fresh scramble replays')

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);run(p.parse_args().run.resolve())
