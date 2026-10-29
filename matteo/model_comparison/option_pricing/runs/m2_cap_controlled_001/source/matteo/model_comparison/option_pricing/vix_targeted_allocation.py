"""Frozen targeted VIX density: fresh random versus fixed allocation."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd
import scipy

from .m2_spx import ROOT,M2,frozen_stage3
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.stratified_sampling import StratifiedGaussianMixture
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.outer_sampling import GaussianOuterSampler,replication_statistics
from esn_eq.pricing import DatedM2Pricer
from esn_eq.vix_pricing import ConditionalVIXPricer


def prefix_mask(frame,power):
    if frame.method.iloc[0]=='fixed_targeted':
        counts=np.array([32]+[2]*16)*2**(power-6)
        return frame.component_point.to_numpy()<counts[frame.component.to_numpy()]
    return frame.point.to_numpy()<2**power


def run(args):
    start=time.perf_counter();folder=Path(__file__).parent
    previous=folder/'runs/m2_vix_targeted_check_001'
    scope_path=folder/'active_scope.json';proposal_path=previous/'proposals.json'
    sources=capture_sources([Path(__file__),folder/'m2_spx.py',scope_path,previous/'run.json',proposal_path,
        *sorted((M2/'src/esn_eq').glob('*.py'))])
    prior=json.loads((previous/'run.json').read_text());scope=json.loads(scope_path.read_text())
    if scope['vix_fit_expiry_days']!=[27,55,118,153] or scope['vix_deferred_expiry_days']:
        raise ValueError('Require full four-maturity scope')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    architecture,parameters,premium,frozen_hashes=frozen_stage3()
    dated=DatedM2Pricer(architecture,premium);initial,_=dated.initial_state([0.,0.,0.])
    expected=json.loads(json.dumps(dict(architecture=asdict(architecture),parameters=asdict(parameters),
        premium=asdict(premium),initial_state=initial.tolist(),target_sha256=targets.artifact_sha256)))
    if any(prior[k]!=v for k,v in expected.items()) or (prior['conditional_order'],prior['normalizer_order'],prior['time_order'])!=(32,48,8):
        raise ValueError('Frozen model, target or quadrature mismatch')
    proposal_bytes=proposal_path.read_bytes();declaration=json.loads(proposal_bytes)
    if hashlib.sha256(proposal_bytes).hexdigest()!=prior['proposal_sha256']:
        raise ValueError('Frozen proposal changed')
    curve=NuisanceCurveFamily.from_targets(targets).build(prior['curve_levels'])
    engine=ConditionalVIXPricer(dated,parameters,curve,initial,normalizer_order=48)
    args.output.mkdir(parents=True,exist_ok=False);save_sources(args.output,sources,ROOT)
    (args.output/'proposals.json').write_bytes(proposal_bytes)
    (args.output/'status.json').write_text(json.dumps(dict(status='running',design='frozen same proposal; fixed allocations versus random component selection'))+'\n')
    rows,rep_rows,raw,seeds,changes,refinements=[],[],[],[],[],[];state_archive={}
    for days in scope['vix_fit_expiry_days']:
        quotes=targets.select('VIX',expiry_days=[days]).reset_index(drop=True)
        expiry=float(quotes.T_years.iloc[0]);sampler=GaussianOuterSampler(dated,expiry,initial)
        mixture=DefensiveGaussianMixture(declaration[str(days)]['shifts'],declaration[str(days)]['probabilities'])
        stratifier=StratifiedGaussianMixture(mixture,[32]+[2]*16)
        window=engine.gaussian_window(expiry,time_order=8)
        sign=np.where(quotes.option_type.to_numpy()=='call',1.,-1.);strikes=quotes.strike.to_numpy()
        descriptors=[dict(metric='forward',strike=None,option_type=None,scale=.5),dict(metric='second_moment',strike=None,option_type=None,scale=None)]
        descriptors += [dict(metric='option',strike=float(q.strike),option_type=q.option_type,scale=float(q.price_scale_U)) for q in quotes.itertuples()]
        for method_id,method in enumerate(('random_targeted','fixed_targeted')):
            estimates={p:[] for p in args.powers};hits={p:[] for p in args.powers};saved_states=[]
            squares={p:[] for p in args.powers};maxima={p:[] for p in args.powers}
            for replicate in range(args.replications):
                tick=time.perf_counter()
                seed=int(np.random.SeedSequence([args.seed,days,method_id,replicate]).generate_state(1)[0])
                if method=='fixed_targeted':
                    draw=stratifier.draw(max(args.powers),seed)
                    z,likelihood,allocation,component,component_point=draw.normals,draw.likelihood_weights,draw.allocation_factors,draw.components,draw.component_points
                else:
                    z,likelihood,component=mixture.draw(max(args.powers),seed)
                    allocation=np.ones(len(z));component_point=np.full(len(z),-1,int)
                states=sampler.mean+z@sampler.root.T;saved_states.append(states)
                v=engine.evaluate(states,window,order=32).vix
                weight=likelihood*allocation
                frame=pd.DataFrame(dict(days=days,method=method,replicate=replicate,point=np.arange(len(v)),
                    component=component,component_point=component_point,likelihood_weight=likelihood,allocation_factor=allocation,weight=weight,vix=v))
                h=np.maximum(sign*(v[:,None]-strikes),0.)
                contribution=weight[:,None]*np.column_stack((v,v**2,h))
                for power in args.powers:
                    mask=prefix_mask(frame,power)
                    if int(mask.sum())!=2**power:
                        raise ValueError('Incomplete component prefix')
                    c=contribution[mask]
                    estimates[power].append(c.mean(axis=0));hits[power].append((h[mask]>0).sum(axis=0))
                    squares[power].append((c*c).sum(axis=0));maxima[power].append(c.max(axis=0))
                    for j,desc in enumerate(descriptors):
                        rep_rows.append(dict(days=days,method=method,power=power,replicate=replicate,**desc,estimate=float(c[:,j].mean())))
                raw.append(frame);seeds.append(dict(days=days,method=method,replicate=replicate,seed=seed))
                if method=='fixed_targeted' and replicate==0:
                    indices=np.argsort(v)[[0,len(v)//2,-1]]
                    high=engine.evaluate(states[indices],window,order=48).vix
                    refinements.append(dict(days=days,base_vix=v[indices].tolist(),conditional48_vix=high.tolist(),max_change=float(np.max(np.abs(high-v[indices])))))
                print(json.dumps(dict(days=days,method=method,replicate=replicate+1,seconds=time.perf_counter()-tick)),flush=True)
            state_archive[f'{method}_{days}']=np.concatenate(saved_states)
            for power in args.powers:
                matrix=np.array(estimates[power]);mean,se=replication_statistics(matrix)
                total=matrix.sum(axis=0)*2**power;sq=np.sum(squares[power],axis=0);maximum=np.max(maxima[power],axis=0)
                for j,desc in enumerate(descriptors):
                    counts=np.array(hits[power])[:,j-2] if j>=2 else None
                    rows.append(dict(days=days,method=method,power=power,model_moment_target=engine.moment_target(window) if j==1 else None,**desc,estimate=float(mean[j]),replication_se=float(se[j]),
                        se_over_scale=float(se[j]/desc['scale']) if desc['scale'] else None,
                        positive_payoff_count=int(counts.sum()) if counts is not None else None,
                        zero_hit_replications=int(np.sum(counts==0)) if counts is not None else None,
                        max_pooled_share=float(maximum[j]/total[j]) if total[j]>0 else None,
                        payoff_contribution_ess=float(total[j]**2/sq[j]) if sq[j]>0 else None))
            diff=np.array(estimates[max(args.powers)])-np.array(estimates[min(args.powers)])
            mean,se=replication_statistics(diff)
            for j,desc in enumerate(descriptors):
                changes.append(dict(days=days,method=method,**desc,change=float(mean[j]),paired_replication_se=float(se[j])))
            pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False)
            pd.DataFrame(rep_rows).to_csv(args.output/'replications.csv',index=False)
            pd.concat(raw,ignore_index=True).to_csv(args.output/'outer_values.csv',index=False)
            np.savez_compressed(args.output/'states.npz',**state_archive)
    check_sources(sources)
    if (args.output/'proposals.json').read_bytes()!=proposal_bytes:
        raise RuntimeError('Proposal changed during evaluation')
    report=dict(status='complete same-proposal allocation diagnostic; accuracy gates unresolved',
        arguments=vars(args)|dict(output=str(args.output)),scope=scope,architecture=asdict(architecture),parameters=asdict(parameters),
        premium=asdict(premium),initial_state=initial.tolist(),curve_levels=prior['curve_levels'],target_sha256=targets.artifact_sha256,
        frozen_model_sha256=frozen_hashes,implementation_sha256=source_hashes(sources,ROOT),proposal_sha256=prior['proposal_sha256'],
        seed_records=seeds,base_power=6,base_counts=[32]+[2]*16,conditional_order=32,normalizer_order=48,time_order=8,
        estimator='random: mean[(f/g)H]; fixed: sum_k pi_k mean_k[(f/g)H]',
        effective_weight='(N*pi_k/n_k)*(f/g), no sample normalization',likelihood_weight_bound=2.,effective_weight_bound=8/3,
        uncertainty_unit='independent complete replication, including all independent component scrambles',
        component_seed_rule='SeedSequence([replication_seed,component_index]).generate_state(1)[0]',
        prefix_rule='First n_k(coarse) points within every component; global first N points is incorrect for stratified draws.',
        conditional_refinement=refinements,input_compatibility=targets.audit['variance_audit'],
        limitations=['Only one pilot proposal; no refit or pilot-to-pilot assessment.','Twelve-replication SEs are exploratory, not price certificates.',
            'All current component counts are powers of two, with only 2/4 tail-component points at coarse/fine counts.',
            'Separate method streams are independent; cross-method SEs are not paired.',
            'Conditional/time/normalizer biases remain; selected-state GH checks do not bound complete-price errors.',
            'Zero hits or observed SEs and contribution ESS do not establish tail convergence.'],
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),wall_seconds=time.perf_counter()-start)
    pd.DataFrame(changes).to_csv(args.output/'prefix_changes.csv',index=False)
    (args.output/'run.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    (args.output/'status.json').write_text('{"status":"complete"}\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--powers',type=int,nargs=2,default=[6,7])
    parser.add_argument('--replications',type=int,default=12)
    parser.add_argument('--seed',type=int,default=2026102203)
    args=parser.parse_args()
    if args.output.exists() or args.replications<4 or not 6<=args.powers[0]<args.powers[1]:
        parser.error('Require new output, >=4 replications and increasing powers >=6')
    run(args)
