"""Independent replay of saved canonical-call prices, SEs and forward conventions."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .parity_adaptive import canonical_market,iid_statistics,scrambled_statistics
from .parity_adaptive_engine import distances
from .all_maturity_core import InSampleLoss
from .full_readout_fit import dump
from .plot_joint_curve_fit import implied_vol


def audit(run):
    design=json.loads((run/'run.json').read_text());cache=Path(design['cache']);controls=Path(design['controls'])
    sources=[]
    roots=[run,cache,controls]
    if (cache/'resume_001/design.json').exists():roots.append(cache/'resume_001')
    for root in roots:
        d=json.loads((root/'design.json').read_text());key='source_sha256'
        for name,digest in d[key].items():
            p=Path(name)
            if not p.is_absolute():p=Path(__file__).resolve().parents[3]/p
            assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,('Source changed',p)
            sources.append(str(p))
    verified=json.loads((run/'immutable_inputs_verified.json').read_text());assert verified['status']=='all unchanged'
    original=pd.concat([pd.read_csv(cache/'spx_contracts.csv'),pd.read_csv(cache/'vix_contracts.csv')],ignore_index=True)
    assert len(original)==1511
    assert not ((original.option_type.eq('put')&(original.strike>original.forward))|(original.option_type.eq('call')&(original.strike<original.forward))).any()
    canonical=canonical_market(original);forwards=pd.read_csv(cache/'forward_targets.csv');records=[];tails=[]
    for name in ['selected','baseline']:
        q=pd.read_csv(run/f'prices_canonical_{name}.csv');o=pd.read_csv(run/f'prices_otm_{name}.csv');f=pd.read_csv(run/f'forwards_{name}.csv')
        saved=np.load(run/f'paths_{name}.npz');p=[];se=[];maxdiff=0.
        for (ins,day),g in original.groupby(['instrument','expiry_days'],sort=False):
            y=saved[f'spx_{day}' if ins=='SPX' else f'vix_values_{day}'];d=np.full(len(g),3.)
            stat=iid_statistics(y,d) if ins=='SPX' else scrambled_statistics(y.reshape(8,512,len(g)),d)
            p.extend(stat['price']);se.extend(stat['se'])
            indices=g.index.to_numpy();np.testing.assert_allclose(q.actual_paths.iloc[indices],stat['actual_paths'],atol=0,rtol=0)
            np.testing.assert_allclose(q.effective_paths.iloc[indices],stat['effective_paths'],rtol=1e-12)
            if ins=='SPX':
                underlying=saved[f'stock_{day}']*float(g.forward.iloc[0]);forward=float(g.forward.iloc[0])
                np.testing.assert_allclose(q.model_price_U.iloc[indices]-o.model_price_U.iloc[indices],np.where(g.option_type.eq('put'),forward-g.strike,0),atol=3e-12)
                # Per-strike active prefixes: tails must be monitored using the
                # actually used draws, not all 4096 when allocation is smaller.
                native_pay=np.maximum(np.where(g.option_type.eq('call'),1.,-1.)*(underlying[:,None]-g.strike.to_numpy()),0)
                hits=np.array([np.count_nonzero(native_pay[:n,j]) for j,n in enumerate(stat['actual_paths'])])
            else:
                underlying=saved[f'vix_{day}'];yf=saved[f'vix_forward_values_{day}'];bf=yf.reshape(8,512).mean(1)
                row=f[f.expiry_days==day].iloc[0];np.testing.assert_allclose(row.model_forward,bf.mean(),rtol=1e-13)
                np.testing.assert_allclose(row.mc_se,bf.std(ddof=1)/np.sqrt(8),rtol=1e-12)
                forward=float(row.model_forward)
                np.testing.assert_allclose(q.model_price_U.iloc[indices]-o.model_price_U.iloc[indices],np.where(g.option_type.eq('put'),forward-g.strike,0),atol=1e-12)
                native_pay=np.maximum(np.where(g.option_type.eq('call'),1.,-1.)*(underlying[:,None]-g.strike.to_numpy()),0).reshape(8,512,len(g))
                hits=np.array([np.count_nonzero(native_pay[:,:n//8,j]) for j,n in enumerate(stat['actual_paths'])])
                blocks=saved[f'vix_native_blocks_{day}'];np.testing.assert_allclose(o.mc_se_U.iloc[indices],blocks.std(0,ddof=1)/np.sqrt(8),rtol=1e-11)
            tails.append(dict(candidate=name,instrument=ins,days=int(day),zero_hit_contracts=int(np.sum(hits==0)),minimum_hits=int(hits.min())))
        np.testing.assert_allclose(q.model_price_U,p,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(q.mc_se_U,se,rtol=1e-12,atol=1e-12)
        for column in ['price_bid_U','price_mid_U','price_ask_U']:
            np.testing.assert_allclose(q[column],canonical[column],rtol=1e-13,atol=2e-12)
        for frame in [q,o]:
            loss=InSampleLoss(frame,forwards).evaluate(frame.model_price_U,f.model_forward,design[name][7:10]);loss['total']+=.1/3
            summary=json.loads((run/f'score_{name}.json').read_text())
            expected=summary['canonical_loss' if frame is q else 'original_otm_loss']
            np.testing.assert_allclose(loss['total'],expected,rtol=1e-11)
        uniform=pd.read_csv(run/f'uniform_budget_prices_{name}.csv');curves=[]
        for (ins,day,n),g in uniform.groupby(['instrument','expiry_days','uniform_paths']):
            curves.append(dict(instrument=ins,days=int(day),paths=int(n),median_se_in_quote_scales=float(np.median(g.mc_se_U/g.price_scale_U)),
                               negative_calls=int(np.sum(g.model_price_U<0)),iv_unavailable=int(np.sum(~np.isfinite(g.model_iv_percent)))))
        records.append(dict(candidate=name,prices_and_SE_replay=True,parity_and_forward_replay=True,market_conversion_replay=True,
                            both_objectives_replay=True,uniform_precision=curves))
    result=dict(status='saved numerical replay passed',distinct_sources_verified=len(set(sources)),records=records,tails=tails,
                limitations=['Sampling SE excludes independent-pilot uncertainty and numerical bias.','One fresh ensemble; no parameter re-estimation or population-convergence proof.',
                             'Common canonical-call objective changes VIX put residuals when model and market forwards differ.'])
    dump(run/'audit.json',result);print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--run',type=Path,required=True);args=a.parse_args();audit(args.run.resolve())
