"""Saved-draw/affine-search replay for the point1 acceleration pilot."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT
from .full_readout_fit import dump
from .calibration_policy import validate, LOWER, UPPER
from .audit_parity_parameter_fit import replay
from .parity_parameter_core import accept_proposal
from .parity_checkpoint import load_result
from .joint_price_local import JointPriceLocal


def audit(run):
    r=json.loads((run/'run.json').read_text());cache=Path(r['cache'])
    original=pd.concat([pd.read_csv(cache/'spx_contracts.csv'),pd.read_csv(cache/'vix_contracts.csv')],ignore_index=True)
    q=pd.read_csv(run/'prices_canonical_baseline.csv');f=pd.read_csv(run/'forwards_baseline.csv')
    baseline=np.array(r['baseline']);selected=np.array(r['selected']);validate(baseline);validate(selected)
    np.testing.assert_array_equal(r['bounds_lower'],LOWER);np.testing.assert_array_equal(r['bounds_upper'],UPPER)
    anchor=replay(run,'baseline',baseline[7:10],original);p=baseline.copy();accepted=0
    for row in r['outer_proposals']:
        x=validate(row['parameters']);trial=replay(run,f"proposal{row['cycle']}",x[7:10],original)
        decision=accept_proposal(anchor,trial,row['decision']['predicted_reduction'])
        assert decision['accepted']==row['decision']['accepted']
        if decision['accepted']:anchor=trial;p=x.copy();accepted+=1
    np.testing.assert_array_equal(p,selected)
    final=replay(run,'selected',selected[7:10],original)
    local=json.loads((run/'local_searches.json').read_text());anchor=load_result(run,'baseline');p=baseline.copy();evaluations=0
    for step in local:
        np.testing.assert_array_equal(step['anchor_parameters'],p)
        model=JointPriceLocal(p,anchor,q,original,f)
        assert step['full_pricer_calls']==0
        for h in step['history']:
            v=model.evaluate(h['parameters']);np.testing.assert_allclose(v['criteria'],h['criteria'],rtol=1e-12,atol=1e-10);evaluations+=1
        row=next((x for x in r['outer_proposals'] if x['cycle']==step['cycle']),None)
        if row is not None and row['decision']['accepted']:
            anchor=load_result(run,f"proposal{step['cycle']}");p=np.array(row['parameters'])
    for path,digest in r['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,('Changed source',path)
    assert r['timing']['full_pricer_calls']==1+len(r['outer_proposals'])
    config=json.loads((run/'selected_model.json').read_text())
    assert config['variance_floor']==selected[10] and config['forward_variance']['levels']==selected[7:10].tolist()
    assert config['initial_state']==r['initial_state']
    for name,x in zip(r['parameter_names'][:7],selected[:7]):assert config['readout'][name]==x
    result=dict(status='point1 price/SE/parity/true-acceptance/local-history/source/policy replay passed',sources_verified=len(r['source_sha256']),
                local_evaluations_replayed=evaluations,accepted_steps=accepted,all_strikes_paths=4096,full_pricer_calls=r['timing']['full_pricer_calls'],
                minimum_volatility=float(np.sqrt(selected[10])),new_independent_checking=False,final=final)
    dump(run/'audit.json',result);print(json.dumps({k:v for k,v in result.items() if k!='final'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();audit(a.run.resolve())
