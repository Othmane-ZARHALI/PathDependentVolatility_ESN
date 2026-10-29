"""Replay score-sampling artifacts; no statistical or pricing-accuracy certificate."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT,frozen_stage3
from esn_eq.market_targets import MarketTargets
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.score_sampling import LowRankInflatedGaussian
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.pricing import DatedM2Pricer


def run(folder):
    record=json.loads((folder/'run.json').read_text())
    proposals=json.loads((folder/'proposals.json').read_text())
    old=json.loads((folder/'pilot_proposals.json').read_text())
    assert record['proposal_sha256']==hashlib.sha256((folder/'proposals.json').read_bytes()).hexdigest()
    for relative,digest in record['implementation_sha256'].items():
        assert hashlib.sha256((folder/'source'/relative).read_bytes()).hexdigest()==digest,relative
    prior=json.loads((folder.parent/'m2_vix_stratified_sampling_001/run.json').read_text())
    assert prior['proposal_sha256']==hashlib.sha256((folder/'pilot_proposals.json').read_bytes()).hexdigest()
    new_seeds=[x['seed'] for x in record['seed_records']]
    assert len(set(new_seeds))==len(new_seeds)
    prior_seeds=set()
    for name in ('m2_vix_importance_sampling_001','m2_vix_pilot_sampling_001','m2_vix_stratified_sampling_001'):
        path=folder.parent/name/'run.json'
        if path.exists():
            for x in json.loads(path.read_text()).get('seed_records',[]):
                prior_seeds.add(x['seed'])
                if x['method']=='stratified':
                    prior_seeds.update(int(np.random.SeedSequence([x['seed'],k]).generate_state(1)[0]) for k in range(7))
    assert not prior_seeds.intersection(new_seeds)
    architecture,parameters,premium,_=frozen_stage3()
    assert json.loads(json.dumps(asdict(architecture)))==record['architecture']
    assert json.loads(json.dumps(asdict(parameters)))==record['parameters']
    assert json.loads(json.dumps(asdict(premium)))==record['premium']
    dated=DatedM2Pricer(architecture,premium);initial=np.array(record['initial_state'])
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    assert targets.artifact_sha256==record['target_sha256']
    raw=pd.read_csv(folder/'outer_values.csv');summary=pd.read_csv(folder/'summary.csv')
    replications=pd.read_csv(folder/'replications.csv');changes=pd.read_csv(folder/'prefix_changes.csv')
    archive=np.load(folder/'states.npz')
    powers=record['arguments']['powers'];R=record['arguments']['replications'];N=2**max(powers)
    assert set(raw.days)=={27,55,118,153} and len(raw)==4*2*R*N
    assert record['scope']['vix_fit_expiry_days']==[27,55,118,153] and not record['scope']['vix_deferred_expiry_days']
    tested=0;state_error=0.;weight_error=0.
    for days in (27,55,118,153):
        quotes=targets.select('VIX',expiry_days=[days]).reset_index(drop=True)
        sign=np.where(quotes.option_type.to_numpy()=='call',1.,-1.)
        sampler=GaussianOuterSampler(dated,float(quotes.T_years.iloc[0]),initial)
        base=DefensiveGaussianMixture(old[str(days)]['shifts'],old[str(days)]['probabilities'])
        packet=proposals[str(days)];score=LowRankInflatedGaussian(packet['basis'],packet['scale'])
        np.testing.assert_allclose(score.rotation,packet['rotation'],rtol=0,atol=2e-14)
        for method_id,method in enumerate(('learned_mean','score_inflated')):
            matrices={p:[] for p in powers};contributions={p:[] for p in powers};hits={p:[] for p in powers}
            saved=archive[f'{method}_{days}'];assert saved.shape==(R*N,len(initial))
            for rep in range(R):
                x=raw[(raw.days==days)&(raw.method==method)&(raw.replicate==rep)].sort_values('point')
                v=x.vix.to_numpy();w=x.weight.to_numpy();assert len(x)==N
                seed=next(s['seed'] for s in record['seed_records'] if s['days']==days and s['method']==method and s['replicate']==rep)
                assert seed==int(np.random.SeedSequence([record['arguments']['seed'],days,method_id,rep]).generate_state(1)[0])
                if method=='learned_mean':z,expected,component=base.draw(max(powers),seed)
                else:
                    d=score.draw(max(powers),seed);z,expected,component=d.normals,d.weights,d.broad.astype(int)
                weight_error=max(weight_error,float(np.max(abs(w-expected))))
                np.testing.assert_allclose(w,expected,rtol=1e-12,atol=1e-14)
                np.testing.assert_array_equal(x.component,component)
                full=sampler.mean+z@sampler.root.T
                state_error=max(state_error,float(np.max(abs(saved[rep*N:(rep+1)*N]-full))))
                np.testing.assert_allclose(saved[rep*N:(rep+1)*N],full,rtol=1e-12,atol=1e-14)
                assert np.all(w>0) and np.all(w<=2+1e-12)
                h=np.column_stack((v,v*v,np.maximum(sign*(v[:,None]-quotes.strike.to_numpy()),0.)))
                for p in powers:
                    c=w[:2**p,None]*h[:2**p];matrices[p].append(c.mean(axis=0));contributions[p].append(c);hits[p].append((h[:2**p,2:]>0).sum(axis=0))
                    rows=replications[(replications.days==days)&(replications.method==method)&(replications.power==p)&(replications.replicate==rep)]
                    np.testing.assert_allclose(rows.estimate,c.mean(axis=0),rtol=2e-12,atol=1e-12)
            for p,matrix in matrices.items():
                matrix=np.array(matrix);mean=matrix.mean(axis=0);se=matrix.std(axis=0,ddof=1)/np.sqrt(R)
                rows=summary[(summary.days==days)&(summary.method==method)&(summary.power==p)]
                np.testing.assert_allclose(rows.estimate,mean,rtol=2e-12,atol=1e-12)
                np.testing.assert_allclose(rows.replication_se,se,rtol=2e-12,atol=1e-12)
                pooled=np.concatenate(contributions[p]);total=pooled.sum(axis=0);positive=total>0
                np.testing.assert_allclose(rows.max_pooled_share.to_numpy()[positive],pooled.max(axis=0)[positive]/total[positive],rtol=2e-12,atol=1e-12)
                np.testing.assert_allclose(rows.payoff_contribution_ess.to_numpy()[positive],total[positive]**2/(pooled*pooled).sum(axis=0)[positive],rtol=2e-12,atol=1e-12)
                np.testing.assert_array_equal(rows.positive_payoff_count.to_numpy()[2:],np.sum(hits[p],axis=0))
                np.testing.assert_array_equal(rows.zero_hit_replications.to_numpy()[2:],np.sum(np.array(hits[p])==0,axis=0))
                tested+=len(rows)
            diff=np.array(matrices[max(powers)])-np.array(matrices[min(powers)])
            rows=changes[(changes.days==days)&(changes.method==method)]
            np.testing.assert_allclose(rows.change,diff.mean(axis=0),rtol=2e-12,atol=1e-12)
            np.testing.assert_allclose(rows.paired_replication_se,diff.std(axis=0,ddof=1)/np.sqrt(R),rtol=2e-12,atol=1e-12)
    result=dict(status='artifact correspondence checks passed',source_snapshot_hashes=len(record['implementation_sha256']),checking_seeds=len(new_seeds),checking_seeds_disjoint_from_previous_outer_studies=True,raw_observations=len(raw),summary_rows_recomputed=tested,proposal_weights_regenerated=True,full_states_regenerated=True,max_state_absolute_error=state_error,max_weight_absolute_error=weight_error,paired_prefix_changes_recomputed=True,verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limitation='Saved-artifact correspondence only; does not reprice the conditional functional or certify accuracy, convergence, or proposal superiority.')
    (folder/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path);run(parser.parse_args().folder)
