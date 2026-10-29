"""Recompute fixed-allocation sampling artifacts; no price-error certificate."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT
from .vix_stratified_sampling import prefix_mask
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.stratified_sampling import StratifiedGaussianMixture
from esn_eq.market_targets import MarketTargets


def run(folder):
    record=json.loads((folder/'run.json').read_text());proposal_bytes=(folder/'proposals.json').read_bytes()
    assert hashlib.sha256(proposal_bytes).hexdigest()==record['proposal_sha256']
    proposals=json.loads(proposal_bytes)
    for relative,digest in record['implementation_sha256'].items():
        assert hashlib.sha256((folder/'source'/relative).read_bytes()).hexdigest()==digest,relative
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    assert targets.artifact_sha256==record['target_sha256']
    history=Path(__file__).with_name('runs');old_seeds=set()
    for name in ('m2_vix_pilot_sampling_001','m2_vix_importance_sampling_001'):
        prior=json.loads((history/name/'run.json').read_text());old_seeds.update(s['seed'] for s in prior['seed_records'])
    actual_seeds=[]
    for s in record['seed_records']:
        if s['method']=='stratified':
            actual_seeds.extend(int(np.random.SeedSequence([s['seed'],k]).generate_state(1)[0]) for k in range(7))
        else:
            actual_seeds.append(s['seed'])
    assert len(actual_seeds)==len(set(actual_seeds))
    assert not set(actual_seeds).intersection(old_seeds)
    raw=pd.read_csv(folder/'outer_values.csv');summary=pd.read_csv(folder/'summary.csv');changes=pd.read_csv(folder/'prefix_changes.csv')
    R=record['arguments']['replications'];powers=record['arguments']['powers'];N=2**max(powers)
    assert set(raw.days)=={27,55,118,153} and len(raw)==4*2*R*N
    tested=0;weight_stats=[]
    for days in (27,55,118,153):
        q=targets.select('VIX',expiry_days=[days]).reset_index(drop=True)
        sign=np.where(q.option_type.to_numpy()=='call',1.,-1.)
        mixture=DefensiveGaussianMixture(proposals[str(days)]['shifts'],proposals[str(days)]['probabilities'])
        sampler=StratifiedGaussianMixture(mixture,record['base_counts'],record['base_power'])
        for method in ('random_mixture','stratified'):
            matrices={p:[] for p in powers}
            for replicate in range(R):
                frame=raw[(raw.days==days)&(raw.method==method)&(raw.replicate==replicate)].sort_values('point')
                assert len(frame)==N
                seed=next(s['seed'] for s in record['seed_records'] if s['days']==days and s['method']==method and s['replicate']==replicate)
                if method=='stratified':
                    draw=sampler.draw(max(powers),seed)
                    likelihood,allocation,component,points=draw.likelihood_weights,draw.allocation_factors,draw.components,draw.component_points
                else:
                    z,likelihood,component=mixture.draw(max(powers),seed)
                    allocation=np.ones(N);points=np.full(N,-1)
                np.testing.assert_array_equal(component,frame.component)
                np.testing.assert_array_equal(points,frame.component_point)
                for expected,column in ((likelihood,'likelihood_weight'),(allocation,'allocation_factor'),(likelihood*allocation,'weight')):
                    np.testing.assert_allclose(expected,frame[column],rtol=1e-12,atol=1e-14)
                assert np.max(likelihood)<=record['likelihood_weight_bound']+1e-13
                v=frame.vix.to_numpy();w=frame.weight.to_numpy()
                h=np.column_stack((v,v**2,np.maximum(sign*(v[:,None]-q.strike.to_numpy()),0.)))
                for power in powers:
                    mask=prefix_mask(frame,power)
                    assert int(mask.sum())==2**power
                    if method=='stratified':
                        np.testing.assert_array_equal(np.bincount(frame.component.to_numpy()[mask]),sampler.counts(power))
                        assert abs(frame.allocation_factor.to_numpy()[mask].mean()-1)<1e-13
                    matrices[power].append((w[mask,None]*h[mask]).mean(axis=0))
                    weight_stats.append(dict(days=days,method=method,replicate=replicate,power=power,
                        effective_weight_mean=float(w[mask].mean()),effective_weight_max=float(w[mask].max())))
            for power,matrix in matrices.items():
                x=np.array(matrix);rows=summary[(summary.days==days)&(summary.method==method)&(summary.power==power)]
                np.testing.assert_allclose(rows.estimate,x.mean(axis=0),rtol=2e-12,atol=1e-12)
                np.testing.assert_allclose(rows.replication_se,x.std(axis=0,ddof=1)/np.sqrt(R),rtol=2e-12,atol=1e-12);tested+=len(rows)
            diff=np.array(matrices[max(powers)])-np.array(matrices[min(powers)])
            rows=changes[(changes.days==days)&(changes.method==method)]
            np.testing.assert_allclose(rows.change,diff.mean(axis=0),rtol=2e-12,atol=1e-12)
            np.testing.assert_allclose(rows.paired_replication_se,diff.std(axis=0,ddof=1)/np.sqrt(R),rtol=2e-12,atol=1e-12)
    report=dict(status='artifact correspondence checks passed',source_hashes=len(record['implementation_sha256']),
        actual_scramble_seeds=len(actual_seeds),checking_seeds_disjoint_from_prior=True,raw_observations=len(raw),
        summary_rows_recomputed=tested,weights_regenerated=True,component_allocations_and_prefixes_checked=True,
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limitation='Checks correspondence, not numerical convergence or statistical superiority.')
    (folder/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    pd.DataFrame(weight_stats).to_csv(folder/'weight_diagnostics.csv',index=False)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path)
    run(parser.parse_args().folder)
