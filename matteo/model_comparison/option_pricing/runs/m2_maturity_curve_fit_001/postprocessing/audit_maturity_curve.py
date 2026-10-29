"""Full saved-price replay for arbitrary curve size and VIX-only checking."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT
from .full_readout_fit import dump
from .all_maturity_core import InSampleLoss
from .maturity_curve import validate,vix_only_loss
from .maturity_curve_local import CurvePriceLocal
from .parity_parameter_core import accept_proposal
from .parity_checkpoint import load_result


def loss_factory(q,f):
    return vix_only_loss(q,f) if q.instrument.eq('VIX').all() else InSampleLoss(q,f)


def replay(folder,name,levels,original):
    q=pd.read_csv(folder/f'prices_canonical_{name}.csv');o=pd.read_csv(folder/f'prices_otm_{name}.csv');f=pd.read_csv(folder/f'forwards_{name}.csv')
    a=np.load(folder/f'paths_{name}.npz');prices=[];errors=[];tails=[]
    assert len(q)==len(original) and (q.actual_paths==4096).all()
    for (ins,day),g in q.groupby(['instrument','expiry_days'],sort=False):
        ix=g.index.to_numpy();raw=a[f'spx_{day}' if ins=='SPX' else f'vix_values_{day}']
        if ins=='SPX':
            price=raw.mean(0);se=raw.std(0,ddof=1)/64;forward=float(g.forward.iloc[0]);underlying=a[f'stock_{day}']*forward
        else:
            b=raw.reshape(8,512,len(g)).mean(1);price=b.mean(0);se=b.std(0,ddof=1)/np.sqrt(8)
            fb=a[f'vix_forward_values_{day}'].reshape(8,512).mean(1);forward=float(fb.mean());row=f[f.expiry_days==day].iloc[0]
            np.testing.assert_allclose([row.model_forward,row.mc_se],[forward,fb.std(ddof=1)/np.sqrt(8)],rtol=1e-12,atol=1e-12)
            np.testing.assert_allclose(o.mc_se_U.iloc[ix],a[f'vix_native_blocks_{day}'].std(0,ddof=1)/np.sqrt(8),rtol=1e-11,atol=1e-12);underlying=a[f'vix_{day}']
        prices.extend(price);errors.extend(se)
        original_side=original.option_type.iloc[ix];np.testing.assert_allclose(q.model_price_U.iloc[ix]-o.model_price_U.iloc[ix],np.where(original_side.eq('put'),forward-g.strike,0),atol=3e-12)
        native_pay=np.maximum(np.where(original_side.eq('call'),1.,-1.)*(underlying[:,None]-g.strike.to_numpy()),0)
        hits=np.sum(native_pay>0,axis=0);tails.append(dict(instrument=ins,days=int(day),zero_hit_contracts=int(np.sum(hits==0)),minimum_raw_hits=int(hits.min())))
    np.testing.assert_allclose(q.model_price_U,prices,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(q.mc_se_U,errors,rtol=1e-12,atol=1e-12)
    c=loss_factory(q,f).evaluate(np.array(prices),f.model_forward.to_numpy(),levels);n=loss_factory(o,f).evaluate(o.model_price_U.to_numpy(),f.model_forward.to_numpy(),levels)
    for result in [c,n]:result['total']+=.1/3;result['components']['state_penalty']=.1/3
    saved=json.loads((folder/f'score_{name}.json').read_text())
    np.testing.assert_allclose([c['total'],n['total']],[saved['canonical_loss'],saved['original_otm_loss']],rtol=1e-11)
    np.testing.assert_allclose(sum(saved['canonical_components'].values()),saved['canonical_loss'],rtol=1e-11)
    np.testing.assert_allclose(sum(saved['original_otm_components'].values()),saved['original_otm_loss'],rtol=1e-11)
    return dict(total=c['total'],components=c['components'],original_otm_total=n['total'],original_otm_components=n['components'],tails=tails)


def audit(run):
    r=json.loads((run/'run.json').read_text());cache=Path(r['cache']);cc=Path(r['check_cache'])
    original=pd.concat([pd.read_csv(cache/'spx_contracts.csv'),pd.read_csv(cache/'vix_contracts.csv')],ignore_index=True)
    vix=original[original.instrument.eq('VIX')].reset_index(drop=True)
    base=validate(r['baseline']);selected=validate(r['selected']);anchor=replay(run,'baseline',base[7:-1],original);p=base.copy();accepted=0
    for row in r['outer_proposals']:
        x=validate(row['parameters']);trial=replay(run,f"proposal{row['cycle']}",x[7:-1],original)
        d=accept_proposal(anchor,trial,row['decision']['predicted_reduction']);assert d['accepted']==row['decision']['accepted']
        screen=next(z for z in json.loads((run/'vix_screens.json').read_text()) if z['cycle']==row['cycle'])
        actual=[trial['components']['VIX']+trial['components']['forwards'],trial['original_otm_components']['VIX']+trial['original_otm_components']['forwards']]
        np.testing.assert_allclose(screen['true_VIX_criteria'],actual,rtol=1e-10,atol=1e-10)
        if d['accepted']:anchor=trial;p=x.copy();accepted+=1
    np.testing.assert_array_equal(p,selected);train=replay(run,'selected',selected[7:-1],original)
    checks={name:replay(run/'checking_vix',name,x[7:-1],vix) for name,x in [('baseline',base),('selected',selected)]}
    q=pd.read_csv(run/'prices_canonical_baseline.csv');f=pd.read_csv(run/'forwards_baseline.csv');anchor=load_result(run,'baseline');p=base.copy();count=0
    for step in json.loads((run/'local_searches.json').read_text()):
        np.testing.assert_array_equal(step['anchor_parameters'],p);model=CurvePriceLocal(p,anchor,q,original,f)
        for h in step['history']:
            np.testing.assert_allclose(model.evaluate(h['parameters'])['criteria'],h['criteria'],rtol=1e-12,atol=1e-10);count+=1
        row=next((z for z in r['outer_proposals'] if z['cycle']==step['cycle']),None)
        if row is not None and row['decision']['accepted']:anchor=load_result(run,f"proposal{step['cycle']}");p=np.array(row['parameters'])
    for path,digest in r['source_sha256'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,('Source changed',path)
    assert (Path(r['controls'])/'coefficients.json').read_bytes()==(Path(r['check_controls'])/'coefficients.json').read_bytes()
    seeds=[]
    for root in [cache,cc]:
        z=[]
        for day in [27,55,118,153]:z.extend(x['seed'] for x in json.loads((root/'production'/f'vix_d{day}_seeds.json').read_text()))
        seeds.append(set(z))
    assert len(seeds[0])==len(seeds[1])==32 and not seeds[0]&seeds[1]
    assert (run/'selected.json').stat().st_mtime_ns<(run/'checking_vix/score_baseline.json').stat().st_mtime_ns
    config=json.loads((run/'selected_model.json').read_text())
    assert config['variance_floor']==selected[-1] and config['forward_variance']['levels']==selected[7:-1].tolist()
    assert config['forward_variance']['ends_years']==r['ends_years'] and config['initial_state']==r['initial_state']
    result=dict(status='maturity-curve price/SE/parity/loss/local/acceptance/source/policy/checking replay passed',accepted_steps=accepted,local_trials_replayed=count,
                sources_verified=len(r['source_sha256']),VIX_seed_overlap=0,controls_refitted=False,selected_before_checking=True,training=train,checking=checks)
    dump(run/'audit.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['training','checking']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();audit(a.run.resolve())

