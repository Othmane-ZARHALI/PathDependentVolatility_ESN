"""Historical draws promoted to design: independent-pilot mean-shift screening."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from .vix_pilot_sampling import pilot_targets
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.score_sampling import LowRankInflatedGaussian
from esn_eq.pilot_proposal import fit_payoff_mixture


def blend(base,learned):
    """Equal mixture, merging identical unshifted components; exact bound two."""
    return DefensiveGaussianMixture(np.r_[base.shifts,learned.shifts[1:]],
        np.r_[.5*(base.probabilities[0]+learned.probabilities[0]),.5*base.probabilities[1:],.5*learned.probabilities[1:]])


def fit_candidate(z,w,v,quotes,basis,base,kind,cap):
    target=pilot_targets(v,quotes)
    if kind.endswith('balanced'):target=np.column_stack((target,v,v*v))
    if kind.startswith('score4'):
        learned,diagnostic=fit_payoff_mixture(z@basis,w,target,coordinates=4,max_shift=cap)
        learned=DefensiveGaussianMixture(learned.shifts@basis.T,learned.probabilities)
    else:learned,diagnostic=fit_payoff_mixture(z,w,target,coordinates=8,max_shift=cap)
    return blend(base,learned),diagnostic


def second_moment_contributions(z,source_weights,h,candidate):
    """E_source[(f/source)(f/candidate) h^2] = integral f^2 h^2/candidate.

    This is an iid importance-estimator second moment, not its variance and
    not a scrambled-Sobol variance prediction. No sample re-normalization.
    """
    return source_weights[:,None]*candidate.weights(z)[:,None]*h*h


def metric_arrays(v,quotes):
    calls=quotes.option_type.to_numpy()=='call';strikes=quotes.strike.to_numpy()
    sign=np.where(calls,1.,-1.);h=np.maximum(sign*(v[:,None]-strikes),0.)
    far=calls&(strikes>=np.quantile(strikes[calls],.75))
    masks={'puts':np.r_[False,False,~calls],'calls':np.r_[False,False,calls],
        'far_calls':np.r_[False,False,far],'forward':np.r_[True,False,np.zeros(len(calls),bool)],
        'second_moment':np.r_[False,True,np.zeros(len(calls),bool)]}
    return np.column_stack((v,v*v,h)),masks


def run(output):
    folder=Path(__file__).parent;train=folder/'runs/m2_vix_score_sampling_001';hold=folder/'runs/m2_vix_stratified_sampling_001'
    scope_path=folder/'active_scope.json'
    paths=[Path(__file__),folder/'m2_spx.py',folder/'vix_pilot_sampling.py',scope_path,
        *[train/name for name in ('run.json','proposals.json','pilot_proposals.json','outer_values.csv')],
        *[hold/name for name in ('run.json','proposals.json','outer_values.csv')],*sorted((M2/'src/esn_eq').glob('*.py'))]
    sources=capture_sources(paths);scope=json.loads(scope_path.read_text())
    if scope['vix_fit_expiry_days']!=[27,55,118,153] or scope['vix_deferred_expiry_days']:raise ValueError('Require full scope')
    training=json.loads((train/'run.json').read_text());checking=json.loads((hold/'run.json').read_text())
    from esn_eq.market_targets import MarketTargets
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');arch,params,premium,_=frozen_stage3()
    expected=json.loads(json.dumps(dict(architecture=asdict(arch),parameters=asdict(params),premium=asdict(premium),target_sha256=targets.artifact_sha256)))
    for r in (training,checking):
        if any(r[k]!=x for k,x in expected.items()):raise ValueError('Frozen model/targets changed')
    for key in ('initial_state','curve_levels','conditional_order','normalizer_order','time_order'):
        if training[key]!=checking[key]:raise ValueError('Design inputs differ')
    if (training['conditional_order'],training['normalizer_order'],training['time_order'])!=(32,48,8):raise ValueError('Unexpected rule')
    prior=json.loads((train/'pilot_proposals.json').read_text());score=json.loads((train/'proposals.json').read_text())
    for path,digest in ((train/'pilot_proposals.json',checking['proposal_sha256']),(hold/'proposals.json',checking['proposal_sha256']),(train/'proposals.json',training['proposal_sha256'])):
        if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:raise ValueError('Proposal hash mismatch')
    train_seeds={x['seed'] for x in training['seed_records']};hold_seeds={x['seed'] for x in checking['seed_records'] if x['method']=='random_mixture'}
    if train_seeds&hold_seeds:raise ValueError('Pilot and screening randomness overlaps')
    tr=pd.read_csv(train/'outer_values.csv');ho=pd.read_csv(hold/'outer_values.csv')
    configs=[(kind,cap) for kind in ('pca8_options','pca8_balanced','score4_balanced') for cap in (3.,4.)]
    rows=[];packets={};arrays={};fit_records={}
    for days in scope['vix_fit_expiry_days']:
        quotes=targets.select('VIX',expiry_days=[days]).reset_index(drop=True)
        base=DefensiveGaussianMixture(prior[str(days)]['shifts'],prior[str(days)]['probabilities'])
        inflated=LowRankInflatedGaussian(score[str(days)]['basis'],score[str(days)]['scale']);basis=np.array(score[str(days)]['basis'])
        batches=[]
        for rep in range(16):
            zrows=[];wrows=[];vrows=[]
            for method in ('learned_mean','score_inflated'):
                rec=next(x for x in training['seed_records'] if x['days']==days and x['method']==method and x['replicate']==rep)
                frame=tr[(tr.days==days)&(tr.method==method)&(tr.replicate==rep)].sort_values('point')
                if method=='learned_mean':z,w,_=base.draw(6,rec['seed'])
                else:
                    d=inflated.draw(6,rec['seed']);z,w=d.normals,d.weights
                np.testing.assert_allclose(frame.weight,w,rtol=1e-12,atol=1e-14)
                zrows.append(z);wrows.append(w);vrows.append(frame.vix.to_numpy())
            batches.append((np.concatenate(zrows),np.concatenate(wrows),np.concatenate(vrows)))
        z=np.concatenate([x[0] for x in batches]);w=np.concatenate([x[1] for x in batches]);v=np.concatenate([x[2] for x in batches])
        arrays[f'train_z_{days}']=z;arrays[f'train_w_{days}']=w;arrays[f'train_v_{days}']=v
        held=[]
        for rep in range(checking['arguments']['replications']):
            rec=next(x for x in checking['seed_records'] if x['days']==days and x['method']=='random_mixture' and x['replicate']==rep)
            frame=ho[(ho.days==days)&(ho.method=='random_mixture')&(ho.replicate==rep)].sort_values('point')
            zz,ww,_=base.draw(max(checking['arguments']['powers']),rec['seed'])
            np.testing.assert_allclose(frame.weight,ww,rtol=1e-12,atol=1e-14)
            held.append((zz,ww,frame.vix.to_numpy()))
        arrays[f'hold_z_{days}']=np.array([x[0] for x in held]);arrays[f'hold_w_{days}']=np.array([x[1] for x in held]);arrays[f'hold_v_{days}']=np.array([x[2] for x in held])
        packets[str(days)]={'baseline':dict(shifts=base.shifts.tolist(),probabilities=base.probabilities.tolist())};fit_records[str(days)]={}
        for kind,cap in configs:
            name=f'{kind}_cap{int(cap)}';packets[str(days)][name]={};fit_records[str(days)][name]={}
            for pilot in range(5):
                # Four disjoint 512-point pilots plus one final pooled 2048-point fit.
                ix=slice(pilot*512,(pilot+1)*512) if pilot<4 else slice(None)
                candidate,diagnostic=fit_candidate(z[ix],w[ix],v[ix],quotes,basis,base,kind,cap)
                packets[str(days)][name][str(pilot)]=dict(shifts=candidate.shifts.tolist(),probabilities=candidate.probabilities.tolist())
                fit_records[str(days)][name][str(pilot)]=diagnostic
                for rep,(zz,ww,vv) in enumerate(held):
                    h,masks=metric_arrays(vv,quotes)
                    trial=second_moment_contributions(zz,ww,h,candidate).mean(axis=0)
                    reference=second_moment_contributions(zz,ww,h,base).mean(axis=0)
                    for family,mask in masks.items():
                        # Ratios formed after pooling complete held scrambles below.
                        for j in np.flatnonzero(mask):rows.append(dict(days=days,candidate=name,pilot=pilot,replicate=rep,family=family,metric_index=int(j),candidate_second_moment=float(trial[j]),baseline_second_moment=float(reference[j])))
    table=pd.DataFrame(rows)
    pooled=table.groupby(['days','candidate','pilot','family','metric_index'],as_index=False)[['candidate_second_moment','baseline_second_moment']].mean()
    pooled['second_moment_ratio']=pooled.candidate_second_moment/pooled.baseline_second_moment
    families=pooled.groupby(['days','candidate','pilot','family'],as_index=False).second_moment_ratio.mean()
    ranks=families.groupby(['candidate','pilot']).second_moment_ratio.agg(['mean','max']).reset_index()
    check_sources(sources);output.mkdir(parents=True,exist_ok=False);save_sources(output,sources,ROOT)
    table.to_csv(output/'screening_replications.csv',index=False);pooled.to_csv(output/'screening_metrics.csv',index=False);families.to_csv(output/'screening_families.csv',index=False);ranks.to_csv(output/'screening_ranks.csv',index=False)
    np.savez_compressed(output/'design_arrays.npz',**arrays)
    (output/'proposals.json').write_text(json.dumps(packets,indent=2,allow_nan=False)+'\n');(output/'pilot_diagnostics.json').write_text(json.dumps(fit_records,indent=2,allow_nan=False)+'\n')
    record=dict(status='historical design-only pilot screen; no fresh proposal validation',scope=scope,
        architecture=training['architecture'],parameters=training['parameters'],premium=training['premium'],initial_state=training['initial_state'],curve_levels=training['curve_levels'],target_sha256=training['target_sha256'],conditional_order=32,normalizer_order=48,time_order=8,
        design_source=str(train),screening_source=str(hold),training_seed_count=len(train_seeds),screening_seed_count=len(hold_seeds),disjoint_seed_sets=True,
        pilot_points_per_maturity=512,independent_design_pilots=4,pooled_pilot_points_per_maturity=2048,screening_points_per_maturity=len(held)*len(held[0][1]),
        configurations=[dict(kind=k,cap=c) for k,c in configs],implementation_sha256=source_hashes(sources,ROOT),proposal_sha256=hashlib.sha256((output/'proposals.json').read_bytes()).hexdigest(),
        pilot_policy='Previous checking observations explicitly promoted to proposal-design data. Screening is independent of new centroid fits, but has been seen historically and is used for design selection; it is not fresh validation.',
        limitations=['IID estimator second moment proxy is not variance and does not predict scrambled-Sobol replication SE.',
            'Four pilot fits are disjoint, but their shared screening set makes performance assessments dependent.',
            'Finite screening draws can miss tails and induce selection optimism; a pooled fit is not an independently validated fifth pilot.',
            'No new conditional pricing, normalizer/time refinement, market fit, SPX change, or Lean certificate.'])
    (output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');print(ranks.to_string(index=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():parser.error('Output exists')
    run(args.output.resolve())
