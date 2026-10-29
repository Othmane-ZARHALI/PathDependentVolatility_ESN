"""Independent saved-draw, objective, acceptance and sampling-split replay."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT
from .all_maturity_core import InSampleLoss
from .full_readout_fit import dump
from .parity_parameter_core import accept_proposal


def replay(folder,name,levels,original):
    q=pd.read_csv(folder/f'prices_canonical_{name}.csv');o=pd.read_csv(folder/f'prices_otm_{name}.csv');f=pd.read_csv(folder/f'forwards_{name}.csv')
    a=np.load(folder/f'paths_{name}.npz');prices=[];errors=[];tails=[]
    assert len(q)==1511 and (q.actual_paths==4096).all()
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
    c=InSampleLoss(q,f).evaluate(np.array(prices),f.model_forward.to_numpy(),levels);n=InSampleLoss(o,f).evaluate(o.model_price_U.to_numpy(),f.model_forward.to_numpy(),levels)
    for result in [c,n]:result['total']+=.1/3;result['components']['state_penalty']=.1/3
    saved=json.loads((folder/f'score_{name}.json').read_text())
    np.testing.assert_allclose([c['total'],n['total']],[saved['canonical_loss'],saved['original_otm_loss']],rtol=1e-11)
    np.testing.assert_allclose(sum(saved['canonical_components'].values()),saved['canonical_loss'],rtol=1e-11)
    np.testing.assert_allclose(sum(saved['original_otm_components'].values()),saved['original_otm_loss'],rtol=1e-11)
    return dict(total=c['total'],components=c['components'],original_otm_total=n['total'],original_otm_components=n['components'],tails=tails)


def audit(run):
    r=json.loads((run/'run.json').read_text());cache=Path(r['cache']);check=Path(r['check_cache']);control=Path(r['controls']);cc=Path(r['check_controls']);source_count=set()
    for root in [run,check,cc]:
        d=json.loads((root/'design.json').read_text())
        for p,digest in d['source_sha256'].items():
            actual=ROOT/p;assert hashlib.sha256(actual.read_bytes()).hexdigest()==digest,('Changed source',actual);source_count.add(str(actual))
    original=pd.concat([pd.read_csv(cache/'spx_contracts.csv'),pd.read_csv(cache/'vix_contracts.csv')],ignore_index=True)
    assert not ((original.option_type.eq('put')&(original.strike>original.forward))|(original.option_type.eq('call')&(original.strike<original.forward))).any()
    assert (control/'coefficients.json').read_bytes()==(cc/'coefficients.json').read_bytes()
    training_seeds=[];checking_seeds=[]
    for root,output in [(cache,training_seeds),(check,checking_seeds)]:
        for p in sorted((root/'production').glob('d*/batch*/completed.json')):output.append(json.loads(p.read_text())['seed'])
    assert len(training_seeds)==len(checking_seeds)==224
    assert not set(training_seeds)&set(checking_seeds),'Training/checking seed overlap'
    base=np.array(r['baseline']);selected=np.array(r['selected']);lo=np.array(r['bounds_lower']);hi=np.array(r['bounds_upper'])
    assert np.all(selected>=lo) and np.all(selected<=hi) and selected[4]>0
    train_base=replay(run,'baseline',base[7:10],original);anchor=train_base;anchor_p=base.copy();acceptances=[]
    for row in r['outer_proposals']:
        trial=replay(run,f"proposal{row['cycle']}",np.array(row['parameters'])[7:10],original)
        decision=accept_proposal(anchor,trial,row['decision']['predicted_reduction']);assert decision['accepted']==row['decision']['accepted']
        acceptances.append(decision)
        if decision['accepted']:anchor=trial;anchor_p=np.array(row['parameters'])
    np.testing.assert_array_equal(selected,anchor_p)
    train_final=replay(run,'selected',selected[7:10],original)
    check_base=replay(run/'checking','baseline',base[7:10],original);check_final=replay(run/'checking','selected',selected[7:10],original)
    assert (run/'selected.json').stat().st_mtime_ns<(run/'checking/score_baseline.json').stat().st_mtime_ns
    result=dict(status='fit/independent-check replay passed',sources_verified=len(source_count),distinct_training_seeds=len(set(training_seeds)),distinct_checking_seeds=len(set(checking_seeds)),
                seed_overlap=0,controls_fitted_on_checking=False,all_accepted_steps_replayed=True,accepted_steps=sum(x['accepted'] for x in acceptances),
                training=dict(baseline=train_base,selected=train_final),checking=dict(baseline=check_base,selected=check_final),
                limitation='Finite signed Monte Carlo prices and pointwise SEs; no deterministic error/market-OOS/global optimality certificate.')
    dump(run/'audit.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ['training','checking']},indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--run',type=Path,required=True);args=a.parse_args();audit(args.run.resolve())
