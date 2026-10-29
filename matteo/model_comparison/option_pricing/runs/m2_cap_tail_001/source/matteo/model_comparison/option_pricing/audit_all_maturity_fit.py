"""Replay expanded in-sample membership, estimates and fixed-scale scores."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .all_maturity_core import populations,nearest_same_side,SPX_DAYS,VIX_DAYS,EXCLUDED_DAYS
from .audit_joint_curve_fit import score,paired_loss
from .m2_spx import ROOT


def audit(run):
    record=json.loads((run/'run.json').read_text())
    for name,h in record['implementation_sha256'].items():
        assert hashlib.sha256((run/'source'/name).read_bytes()).hexdigest()==h,name
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    raw=ROOT/'data/processed/options/2026-02-19/quotes_calibration.csv'
    source=pd.read_csv(raw)
    assert set(source.loc[source.instrument=='SPX','expiry_days'])-set(EXCLUDED_DAYS)==set(SPX_DAYS)
    spx,vix,m,a=populations(source)
    for name,table in [('spx_contracts',spx),('vix_contracts',vix),('spx_membership',m),('atm_proxies',a)]:
        pd.testing.assert_frame_equal(pd.read_csv(run/(name+'.csv')),table,check_dtype=False,rtol=1e-12,atol=1e-12)
    q=pd.concat([spx,vix],ignore_index=True);targets=pd.read_csv(run/'forward_targets.csv')
    assert q.split.eq('in_sample_fit').all() and len(q)==1511
    assert set(spx.expiry_days)==set(SPX_DAYS) and set(vix.expiry_days)==set(VIX_DAYS)
    assert json.loads((run/'selected.json').read_text())['checking_not_used']
    levels={'baseline':record['baseline_levels'],'selected':record['selected_levels']}
    archive=np.load(run/'comparison_paths.npz');summary=[];expiry=[];byrole={}
    for row in record['scores']:
        role,name=row['role'],row['candidate'];prefix=f'{role}_{name}'
        p=pd.read_csv(run/f'prices_{role}_{name}.csv');f=pd.read_csv(run/f'forwards_{role}_{name}.csv')
        pd.testing.assert_frame_equal(p[['instrument','expiry','strike','option_type']],q[['instrument','expiry','strike','option_type']])
        total,components,residual=score(q,targets,p.model_price_U.to_numpy(),f.model_forward.to_numpy(),levels[name])
        total+=record['state_penalty'];components['state_penalty']=record['state_penalty']
        np.testing.assert_allclose(total,row['total'],rtol=1e-12)
        for key,value in components.items():np.testing.assert_allclose(value,row['components'][key],rtol=1e-12)
        blocks=archive[prefix+'_blocks'];fb=archive[prefix+'_forward_blocks']
        np.testing.assert_allclose(blocks.mean(axis=0),p.model_price_U,rtol=1e-12,atol=1e-10)
        np.testing.assert_allclose(fb.mean(axis=0),f.model_forward,rtol=1e-12,atol=1e-10)
        for (instrument,days),g in p.groupby(['instrument','expiry_days']):
            values=archive[prefix+'_'+instrument.lower()+'_'+str(days)]
            if instrument=='VIX':
                # The archived fourth return value is raw simulated VIX, not
                # the controlled forward/option matrix. Reconstruct that matrix.
                assert values.ndim==1 and len(values)==512
                source_role='search' if role=='search' else 'fresh'
                donors=pd.read_csv(ROOT/'matteo/model_comparison/option_pricing/runs/m2_joint_curve_fit_002/vix_contracts.csv')
                donors=donors[donors.expiry_days==days].reset_index(drop=True)
                ix=nearest_same_side(g,donors)
                parent=Path(record['parent_run'])
                with np.load(parent/f'vix_state2_{source_role}_t0_d{days}.npz') as saved:
                    beta=saved['beta'][:,np.r_[0,1+ix]]
                    sign=np.where(g.option_type.to_numpy()=='call',1.,-1.)
                    payoff=np.c_[values,np.maximum(sign*(values[:,None]-g.strike.to_numpy()),0.)]
                    controlled=saved['weights'][:,None]*payoff-saved['controls']@beta
                np.testing.assert_allclose(controlled[:,0].mean(),f.loc[f.expiry_days==days,'model_forward'].iloc[0],rtol=1e-12)
                values=controlled[:,1:].reshape(8,64,len(g)).mean(axis=1)
            np.testing.assert_allclose(values.mean(axis=0),g.model_price_U,rtol=1e-12,atol=1e-9)
            np.testing.assert_allclose(values.std(axis=0,ddof=1)/np.sqrt(len(values)),g.mc_se_U,rtol=1e-12,atol=1e-10)
            expiry.append(dict(role=role,candidate=name,instrument=instrument,expiry_days=int(days),quotes=len(g),
                mse=float(np.mean(g.interval_residual**2)),inside_bidask=int(g.interval_residual.eq(0).sum()),
                negative_estimates=int(g.model_price_U.lt(0).sum()),nominal_se_passes=int(g.nominal_se_screen.sum()),
                median_se_over_scale=float(np.median(g.mc_se_U/g.price_scale_U))))
        summary.append(dict(role=role,candidate=name,total=total,components=components,
            quotes=len(q),inside_bidask=int((residual==0).sum()),negative_estimates=int(p.model_price_U.lt(0).sum()),
            nominal_se_passes=int(p.nominal_se_screen.sum())))
        byrole.setdefault(role,{})[name]=(blocks,fb)
    paired={}
    for role,d in byrole.items():
        a,fa=d['baseline'];b,fb=d['selected']
        r=paired_loss(q,targets,a,b,fa,fb,levels['baseline'],levels['selected'])
        r['selected_minus_baseline']=r.pop('selected_minus_anchor');paired[role]=r
    seeds=[];new=[]
    for role in ['search','check']:
        for days in [119,666,1030]:
            r=json.loads((run/f'new_spx_{role}_d{days}/completed.json').read_text());new.append(r);seeds+=r['seeds']
    assert len(seeds)==24 and len(set(seeds))==24
    parent=json.loads((Path(record['parent_run'])/'run.json').read_text())
    assert not set(seeds)&{r['seed'] for r in parent['seed_records']}
    pd.DataFrame(expiry).to_csv(run/'expiry_scores.csv',index=False)
    result=dict(status='Expanded membership, all four score records, archived sample means/SEs, new seeds and input/source hashes replayed.',
        population=dict(spx=len(spx),vix=len(vix),spx_expiry_days=list(SPX_DAYS),vix_expiry_days=list(VIX_DAYS),
            market_out_of_sample=False,original_splits_retained_as_metadata=True),
        scores=summary,paired=paired,new_maturity_generation=new,
        uncertainty='Conditional sampling SE with frozen controls. Eight VIX scramble blocks; SPX iid path contributions. Stock timestep, normalizer interpolation/GH and VIX time/conditional quadrature bias remain separate. Reused prior numerical ensembles do not supply a new independent calibration validation.')
    (run/'score_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='new_maturity_generation'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    audit(p.parse_args().run)
