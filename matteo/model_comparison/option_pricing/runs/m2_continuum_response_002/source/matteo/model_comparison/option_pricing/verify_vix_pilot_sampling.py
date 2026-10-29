"""Artifact correspondence checks, not pricing-error or statistical certificates."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from .m2_spx import ROOT
from esn_eq.market_targets import MarketTargets
from esn_eq.importance_sampling import DefensiveGaussianMixture


def run(folder):
    record=json.loads((folder/'run.json').read_text())
    proposals=json.loads((folder/'proposals.json').read_text())
    assert record['proposal_sha256']==hashlib.sha256((folder/'proposals.json').read_bytes()).hexdigest()
    for relative,digest in record['implementation_sha256'].items():
        assert hashlib.sha256((folder/'source'/relative).read_bytes()).hexdigest()==digest,relative
    old=json.loads((Path(__file__).with_name('runs')/'m2_vix_importance_sampling_001/run.json').read_text())
    old_seeds={x['seed'] for x in old['seed_records']}
    new_seeds=[x['seed'] for x in record['seed_records']]
    assert len(set(new_seeds))==len(new_seeds)
    assert not old_seeds.intersection(new_seeds)
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    assert targets.artifact_sha256==record['target_sha256']
    raw=pd.read_csv(folder/'outer_values.csv'); summary=pd.read_csv(folder/'summary.csv')
    changes=pd.read_csv(folder/'prefix_changes.csv')
    powers=record['arguments']['powers']; R=record['arguments']['replications']; N=2**max(powers)
    assert set(raw.days)=={27,55,118,153} and len(raw)==4*2*R*N
    tested=0
    for days in (27,55,118,153):
        quotes=targets.select('VIX',expiry_days=[days]).reset_index(drop=True)
        sign=np.where(quotes.option_type.to_numpy()=='call',1.,-1.)
        prop=DefensiveGaussianMixture(proposals[str(days)]['shifts'],proposals[str(days)]['probabilities'])
        for method in ('sobol','pilot_mixture'):
            matrices={power:[] for power in powers}
            for rep in range(R):
                x=raw[(raw.days==days)&(raw.method==method)&(raw.replicate==rep)].sort_values('point')
                v=x.vix.to_numpy(); w=x.weight.to_numpy()
                assert len(x)==N
                seed=next(s['seed'] for s in record['seed_records'] if s['days']==days and s['method']==method and s['replicate']==rep)
                expected=prop.draw(max(powers),seed)[1] if method=='pilot_mixture' else np.ones(N)
                np.testing.assert_allclose(w,expected,rtol=1e-12,atol=1e-14)
                h=np.column_stack((v,v**2,np.maximum(sign*(v[:,None]-quotes.strike.to_numpy()),0.)))
                for power in powers:
                    matrices[power].append((w[:2**power,None]*h[:2**power]).mean(axis=0))
            for power,matrix in matrices.items():
                matrix=np.array(matrix); mean=matrix.mean(axis=0); se=matrix.std(axis=0,ddof=1)/np.sqrt(R)
                rows=summary[(summary.days==days)&(summary.method==method)&(summary.power==power)]
                np.testing.assert_allclose(rows.estimate,mean,rtol=2e-12,atol=1e-12)
                np.testing.assert_allclose(rows.replication_se,se,rtol=2e-12,atol=1e-12)
                tested+=len(rows)
            diff=np.array(matrices[max(powers)])-np.array(matrices[min(powers)])
            rows=changes[(changes.days==days)&(changes.method==method)]
            np.testing.assert_allclose(rows.change,diff.mean(axis=0),rtol=2e-12,atol=1e-12)
            np.testing.assert_allclose(rows.paired_replication_se,diff.std(axis=0,ddof=1)/np.sqrt(R),rtol=2e-12,atol=1e-12)
    result=dict(status='artifact correspondence checks passed',source_snapshot_hashes=len(record['implementation_sha256']),
        checking_seeds=len(new_seeds),checking_seeds_disjoint_from_pilot=True,raw_observations=len(raw),
        summary_rows_recomputed=tested,proposal_weights_regenerated=True,paired_prefix_changes_recomputed=True,
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limitation='Checks saved-artifact correspondence only; does not certify price accuracy, convergence or proposal superiority.')
    (folder/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    run(parser.parse_args().folder)
