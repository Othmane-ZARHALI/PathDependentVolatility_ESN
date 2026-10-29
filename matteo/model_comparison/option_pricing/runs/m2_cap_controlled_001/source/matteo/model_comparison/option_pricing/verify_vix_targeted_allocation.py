"""Recompute fixed-allocation sampling artifacts; no price-error certificate."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT,frozen_stage3
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.pricing import DatedM2Pricer
from .vix_targeted_allocation import prefix_mask
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
    for name in ('m2_vix_pilot_sampling_001','m2_vix_importance_sampling_001','m2_vix_stratified_sampling_001','m2_vix_score_sampling_001','m2_vix_targeted_check_001'):
        prior=json.loads((history/name/'run.json').read_text())
        for s in prior['seed_records']:
            old_seeds.add(s['seed'])
            if s['method']=='stratified':
                old_seeds.update(int(np.random.SeedSequence([s['seed'],k]).generate_state(1)[0]) for k in range(7))
    actual_seeds=[]
    for s in record['seed_records']:
        if s['method']=='fixed_targeted':
            actual_seeds.extend(int(np.random.SeedSequence([s['seed'],k]).generate_state(1)[0]) for k in range(17))
        else:
            actual_seeds.append(s['seed'])
    assert len(actual_seeds)==len(set(actual_seeds))
    assert not set(actual_seeds).intersection(old_seeds)
    raw=pd.read_csv(folder/'outer_values.csv');summary=pd.read_csv(folder/'summary.csv');changes=pd.read_csv(folder/'prefix_changes.csv')
    R=record['arguments']['replications'];powers=record['arguments']['powers'];N=2**max(powers)
    assert set(raw.days)=={27,55,118,153} and len(raw)==4*2*R*N
    tested=0;weight_stats=[];archive=np.load(folder/'states.npz')
    arch,params,premium,_=frozen_stage3();dated=DatedM2Pricer(arch,premium);initial=np.array(record['initial_state'])
    for key,value in [('architecture',arch),('parameters',params),('premium',premium)]:
        assert json.loads(json.dumps(asdict(value)))==record[key]
    assert record['scope']['vix_fit_expiry_days']==[27,55,118,153] and not record['scope']['vix_deferred_expiry_days']
    rep_table=pd.read_csv(folder/'replications.csv')
    for days in (27,55,118,153):
        q=targets.select('VIX',expiry_days=[days]).reset_index(drop=True)
        sign=np.where(q.option_type.to_numpy()=='call',1.,-1.)
        outer=GaussianOuterSampler(dated,float(q.T_years.iloc[0]),initial)
        mixture=DefensiveGaussianMixture(proposals[str(days)]['shifts'],proposals[str(days)]['probabilities'])
        sampler=StratifiedGaussianMixture(mixture,record['base_counts'],record['base_power'])
        for method_id,method in enumerate(('random_targeted','fixed_targeted')):
            matrices={p:[] for p in powers};contributions={p:[] for p in powers};hits={p:[] for p in powers}
            for replicate in range(R):
                frame=raw[(raw.days==days)&(raw.method==method)&(raw.replicate==replicate)].sort_values('point')
                assert len(frame)==N
                seed=next(s['seed'] for s in record['seed_records'] if s['days']==days and s['method']==method and s['replicate']==replicate)
                assert seed==int(np.random.SeedSequence([record['arguments']['seed'],days,method_id,replicate]).generate_state(1)[0])
                if method=='fixed_targeted':
                    draw=sampler.draw(max(powers),seed)
                    z=draw.normals
                    likelihood,allocation,component,points=draw.likelihood_weights,draw.allocation_factors,draw.components,draw.component_points
                else:
                    z,likelihood,component=mixture.draw(max(powers),seed)
                    allocation=np.ones(N);points=np.full(N,-1)
                np.testing.assert_allclose(archive[f'{method}_{days}'][replicate*N:(replicate+1)*N],outer.mean+z@outer.root.T,rtol=1e-12,atol=1e-14)
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
                    if method=='fixed_targeted':
                        np.testing.assert_array_equal(np.bincount(frame.component.to_numpy()[mask]),sampler.counts(power))
                        assert abs(frame.allocation_factor.to_numpy()[mask].mean()-1)<1e-13
                    assert np.max(w[mask])<=record['effective_weight_bound']+1e-12
                    c=w[mask,None]*h[mask];matrices[power].append(c.mean(axis=0));contributions[power].append(c);hits[power].append((h[mask,2:]>0).sum(axis=0))
                    rep_rows=rep_table[(rep_table.days==days)&(rep_table.method==method)&(rep_table.power==power)&(rep_table.replicate==replicate)]
                    np.testing.assert_allclose(rep_rows.estimate,c.mean(axis=0),rtol=2e-12,atol=1e-12)
                    weight_stats.append(dict(days=days,method=method,replicate=replicate,power=power,
                        effective_weight_mean=float(w[mask].mean()),effective_weight_max=float(w[mask].max())))
            for power,matrix in matrices.items():
                x=np.array(matrix);rows=summary[(summary.days==days)&(summary.method==method)&(summary.power==power)]
                np.testing.assert_allclose(rows.estimate,x.mean(axis=0),rtol=2e-12,atol=1e-12)
                np.testing.assert_allclose(rows.replication_se,x.std(axis=0,ddof=1)/np.sqrt(R),rtol=2e-12,atol=1e-12);tested+=len(rows)
                pooled=np.concatenate(contributions[power]);total=pooled.sum(axis=0);positive=total>0
                np.testing.assert_allclose(rows.max_pooled_share.to_numpy()[positive],pooled.max(axis=0)[positive]/total[positive],rtol=2e-12,atol=1e-12)
                np.testing.assert_allclose(rows.payoff_contribution_ess.to_numpy()[positive],total[positive]**2/(pooled*pooled).sum(axis=0)[positive],rtol=2e-12,atol=1e-12)
                np.testing.assert_array_equal(rows.positive_payoff_count.to_numpy()[2:],np.sum(hits[power],axis=0))
                np.testing.assert_array_equal(rows.zero_hit_replications.to_numpy()[2:],np.sum(np.array(hits[power])==0,axis=0))
            diff=np.array(matrices[max(powers)])-np.array(matrices[min(powers)])
            rows=changes[(changes.days==days)&(changes.method==method)]
            np.testing.assert_allclose(rows.change,diff.mean(axis=0),rtol=2e-12,atol=1e-12)
            np.testing.assert_allclose(rows.paired_replication_se,diff.std(axis=0,ddof=1)/np.sqrt(R),rtol=2e-12,atol=1e-12)
    report=dict(status='artifact correspondence checks passed',source_hashes=len(record['implementation_sha256']),
        actual_scramble_seeds=len(actual_seeds),checking_seeds_disjoint_from_prior=True,raw_observations=len(raw),
        summary_rows_recomputed=tested,weights_regenerated=True,full_states_regenerated=True,component_allocations_and_prefixes_checked=True,
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limitation='Checks correspondence, not numerical convergence or statistical superiority.')
    (folder/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    pd.DataFrame(weight_stats).to_csv(folder/'weight_diagnostics.csv',index=False)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path)
    run(parser.parse_args().folder)
