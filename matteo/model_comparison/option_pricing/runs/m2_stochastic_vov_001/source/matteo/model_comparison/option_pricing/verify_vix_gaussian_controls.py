"""Replay frozen Gaussian-control fits and paired checking statistics."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT,frozen_stage3
from .vix_gaussian_controls import METHODS,fit_variants,payoffs
from esn_eq.gaussian_controls import gaussian_control_rows
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.pricing import DatedM2Pricer
from esn_eq.market_targets import MarketTargets


def run(folder):
    r=json.loads((folder/'run.json').read_text());packets=json.loads((folder/'coefficients.json').read_text());props=json.loads((folder/'proposals.json').read_text())
    assert hashlib.sha256((folder/'coefficients.json').read_bytes()).hexdigest()==r['coefficient_sha256']
    assert hashlib.sha256((folder/'proposals.json').read_bytes()).hexdigest()==r['proposal_sha256']
    for name,digest in r['implementation_sha256'].items():assert hashlib.sha256((folder/'source'/name).read_bytes()).hexdigest()==digest,name
    assert r['scope']['vix_fit_expiry_days']==[27,55,118,153] and not r['scope']['vix_deferred_expiry_days']
    arch,params,premium,_=frozen_stage3()
    for name,obj in [('architecture',arch),('parameters',params),('premium',premium)]:assert json.loads(json.dumps(asdict(obj)))==r[name]
    dated=DatedM2Pricer(arch,premium);target=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');assert target.artifact_sha256==r['target_sha256']
    seeds=[s['seed'] for s in r['seed_records']];assert len(seeds)==len(set(seeds))
    old_seeds=set(r['pilot_seed_records'])
    for name in ('m2_vix_importance_sampling_001','m2_vix_pilot_sampling_001','m2_vix_stratified_sampling_001','m2_vix_score_sampling_001','m2_vix_targeted_check_001','m2_vix_targeted_allocation_001'):
        prior=json.loads((folder.parent/name/'run.json').read_text())
        for s in prior['seed_records']:
            old_seeds.add(s['seed'])
            if s.get('method') in ('stratified','fixed_targeted'):
                for k in range(7 if s['method']=='stratified' else 17):old_seeds.add(int(np.random.SeedSequence([s['seed'],k]).generate_state(1)[0]))
    assert not old_seeds.intersection(seeds)
    pilot=np.load(folder/'pilot_arrays.npz');states=np.load(folder/'states.npz');controls=np.load(folder/'control_rows.npz')
    raw=pd.read_csv(folder/'outer_values.csv');summary=pd.read_csv(folder/'summary.csv');rep_table=pd.read_csv(folder/'replications.csv');changes=pd.read_csv(folder/'prefix_changes.csv')
    powers=r['arguments']['powers'];N=2**max(powers);R=r['arguments']['replications'];tested=0;pilot_record=json.loads((Path(r['pilot_source'])/'run.json').read_text())
    for days in (27,55,118,153):
        q=target.select('VIX',expiry_days=[days]).reset_index(drop=True);U=np.array(packets[str(days)]['basis']);proposal=DefensiveGaussianMixture(props[str(days)]['shifts'],props[str(days)]['probabilities'])
        Y=[];C=[]
        for rep in range(packets[str(days)]['pilot_replications']):
            rec=next(s for s in pilot_record['seed_records'] if s['days']==days and s['method']=='random_targeted' and s['replicate']==rep)
            z,w,_=proposal.draw(max(pilot_record['arguments']['powers']),rec['seed'])
            np.testing.assert_array_equal(z,pilot[f'normals_{days}'][rep]);np.testing.assert_array_equal(w,pilot[f'weights_{days}'][rep])
            Y.append((w[:,None]*payoffs(pilot[f'vix_{days}'][rep],q)).mean(axis=0));C.append(gaussian_control_rows(z,w,U).mean(axis=0))
        np.testing.assert_allclose(pilot[f'outcomes_{days}'],Y,rtol=1e-13,atol=1e-13);np.testing.assert_allclose(pilot[f'controls_{days}'],C,rtol=1e-13,atol=1e-13)
        betas=fit_variants(np.array(Y),np.array(C),r['arguments']['ridge'])
        for method,b in betas.items():np.testing.assert_allclose(packets[str(days)]['coefficients'][method],b,rtol=1e-12,atol=1e-12)
        outer=GaussianOuterSampler(dated,float(q.T_years.iloc[0]),np.array(r['initial_state']));matrices={m:{p:[] for p in powers} for m in METHODS};hits={p:[] for p in powers}
        for rep in range(R):
            seed=next(s['seed'] for s in r['seed_records'] if s['days']==days and s['replicate']==rep)
            assert seed==int(np.random.SeedSequence([r['arguments']['seed'],days,rep]).generate_state(1)[0])
            z,w,component=proposal.draw(max(powers),seed);frame=raw[(raw.days==days)&(raw.replicate==rep)].sort_values('point');assert len(frame)==N
            np.testing.assert_allclose(frame.weight,w,rtol=1e-12,atol=1e-14);np.testing.assert_array_equal(frame.component,component)
            np.testing.assert_allclose(states[str(days)][rep*N:(rep+1)*N],outer.mean+z@outer.root.T,rtol=1e-12,atol=1e-14)
            # Independently spell out the exact zero-mean row formulas.
            y=z@U;C0=np.column_stack((w-1,w[:,None]*y,w[:,None]*(y*y-1)))
            np.testing.assert_allclose(controls[str(days)][rep*N:(rep+1)*N],C0,rtol=1e-13,atol=1e-13)
            h=payoffs(frame.vix.to_numpy(),q)
            for p in powers:
                n=2**p;Y0=(w[:n,None]*h[:n]).mean(axis=0);Cmean=C0[:n].mean(axis=0);hits[p].append((h[:n,2:]>0).sum(axis=0))
                for method in METHODS:
                    adjusted=Y0-Cmean@betas[method];matrices[method][p].append(adjusted)
                    rows=rep_table[(rep_table.days==days)&(rep_table.method==method)&(rep_table.power==p)&(rep_table.replicate==rep)]
                    np.testing.assert_allclose(rows.estimate,adjusted,rtol=2e-12,atol=1e-12)
        for method in METHODS:
            for p in powers:
                x=np.array(matrices[method][p]);plain=np.array(matrices['raw'][p]);d=x-plain
                rows=summary[(summary.days==days)&(summary.method==method)&(summary.power==p)]
                se=x.std(axis=0,ddof=1)/np.sqrt(R)
                for expected,column in ((x.mean(axis=0),'estimate'),(se,'replication_se'),(d.mean(axis=0),'change_from_raw'),(d.std(axis=0,ddof=1)/np.sqrt(R),'paired_change_se'),(se/(plain.std(axis=0,ddof=1)/np.sqrt(R)),'se_ratio_to_raw')):
                    np.testing.assert_allclose(rows[column],expected,rtol=2e-12,atol=1e-12)
                np.testing.assert_array_equal(rows.positive_payoff_count.to_numpy()[2:],np.sum(hits[p],axis=0));np.testing.assert_array_equal(rows.zero_hit_replications.to_numpy()[2:],np.sum(np.array(hits[p])==0,axis=0));np.testing.assert_array_equal(rows.negative_estimate,x.mean(axis=0)<0)
                tested+=len(rows)
            d=np.array(matrices[method][max(powers)])-np.array(matrices[method][min(powers)])
            rows=changes[(changes.days==days)&(changes.method==method)]
            np.testing.assert_allclose(rows.change,d.mean(axis=0),rtol=2e-12,atol=1e-12);np.testing.assert_allclose(rows.paired_replication_se,d.std(axis=0,ddof=1)/np.sqrt(R),rtol=2e-12,atol=1e-12)
    result=dict(status='pilot fits, checking states/controls and paired statistics replay passed',source_snapshot_hashes=len(r['implementation_sha256']),checking_seeds=len(seeds),checking_seeds_disjoint_from_prior=True,raw_observations=len(raw),summary_rows_recomputed=tested,pilot_coefficients_recomputed=True,full_states_and_zero_mean_controls_regenerated=True,negative_estimates_retained=True,verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limitation='Artifact/implementation correspondence, not conditional quadrature or statistical price accuracy certification.')
    (folder/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path);run(parser.parse_args().folder)
