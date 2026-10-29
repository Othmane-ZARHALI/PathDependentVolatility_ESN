"""Integrity/statistics replay for the VIX identity diagnostics, not a certificate."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def check_run(root):
    record=json.loads((root/'run.json').read_text())
    for relative,digest in record['implementation_sha256'].items():
        assert hashlib.sha256((root/'source'/relative).read_bytes()).hexdigest()==digest,relative
    assert record['scope']['vix_fit_expiry_days']==[27,55,118,153]
    assert not record['scope']['vix_deferred_expiry_days']
    summary=pd.read_csv(root/'summary.csv');rep=pd.read_csv(root/'replications.csv')
    if 'multiplier_mean' in summary:
        for row in summary.itertuples():
            thresholds=[100.,1000.,10000.,100000.];j=thresholds.index(row.threshold)
            group=rep[rep.days==row.days]
            for metric,column,se_column in ((1,'multiplier_mean','mean_se'),(2+2*j,'probability','probability_se'),(3+2*j,'multiplier_contribution','contribution_se')):
                x=group[group.metric==metric].estimate.to_numpy()
                np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x))],[getattr(row,column),getattr(row,se_column)],rtol=2e-12,atol=1e-12)
            np.testing.assert_allclose(row.fraction_of_M128,row.multiplier_contribution/row.M128,rtol=2e-12)
    else:
        keys=[k for k in ['days','time_order','normalizer_order','rule','power'] if k in summary]
        for row in summary.itertuples():
            group=rep
            for key in keys:group=group[group[key]==getattr(row,key)]
            x=group.estimate.to_numpy()
            np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x))],[row.estimate,row.replication_se],rtol=2e-12,atol=1e-12)
            if hasattr(row,'paired_change_se'):
                d=group.change.to_numpy();np.testing.assert_allclose([d.mean(),d.std(ddof=1)/np.sqrt(len(d))],[row.change_from_base,row.paired_change_se],rtol=2e-12,atol=1e-12)
        if (root/'conditional_changes.csv').exists():
            changes=pd.read_csv(root/'conditional_changes.csv')
            for row in changes.itertuples():
                q=rep[(rep.days==row.days)&(rep.time_order==row.time_order)&(rep.normalizer_order==row.normalizer_order)&(rep.power==row.power)]
                a=q[q.rule==row.rule].sort_values('replicate').estimate.to_numpy();b=q[q.rule=='continuous'].sort_values('replicate').estimate.to_numpy();d=a-b
                np.testing.assert_allclose([d.mean(),d.std(ddof=1)/np.sqrt(len(d))],[row.GH_minus_continuous,row.paired_se],rtol=2e-12,atol=1e-12)
        if (root/'time_changes.csv').exists():
            for row in pd.read_csv(root/'time_changes.csv').itertuples():
                q=rep[(rep.days==row.days)&(rep.normalizer_order==row.normalizer_order)&(rep.power==row.power)]
                d=q[q.time_order==16].sort_values('replicate').estimate.to_numpy()-q[q.time_order==8].sort_values('replicate').estimate.to_numpy()
                np.testing.assert_allclose([d.mean(),d.std(ddof=1)/np.sqrt(len(d))],[row.time16_minus8,row.paired_se],rtol=2e-12,atol=1e-12)
    result=dict(status='source snapshots and summary statistics verified',source_hashes=len(record['implementation_sha256']),summary_rows=len(summary),all_four_maturities=True,
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limitation='Artifact correspondence only; does not certify numerical accuracy or population inference.')
    (root/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result,record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('roots',type=Path,nargs='+');args=parser.parse_args();seeds=set()
    for root in args.roots:
        result,record=check_run(root)
        current=[s['seed'] for s in record.get('seed_records',[])]
        assert len(current)==len(set(current)) and not seeds.intersection(current)
        seeds.update(current);print(json.dumps(dict(run=root.name,**result)))
