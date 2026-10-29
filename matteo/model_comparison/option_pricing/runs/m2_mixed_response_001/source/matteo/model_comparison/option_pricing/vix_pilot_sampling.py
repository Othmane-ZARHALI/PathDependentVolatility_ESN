"""Pilot-learned VIX proposal, frozen before fresh all-maturity evaluation."""
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
from esn_eq.pilot_proposal import fit_payoff_mixture
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.outer_sampling import GaussianOuterSampler,replication_statistics
from esn_eq.pricing import DatedM2Pricer
from esn_eq.vix_pricing import ConditionalVIXPricer


def pilot_targets(values,quotes):
    calls=quotes.option_type.to_numpy()=='call'
    strikes=quotes.strike.to_numpy()
    far=calls & (strikes>=np.quantile(strikes[calls],.75))
    masks=[~calls,calls,far]
    sign=np.where(calls,1.,-1.)
    payoffs=np.maximum(sign*(values[:,None]-strikes),0.)/quotes.price_scale_U.to_numpy()
    if any(not mask.any() for mask in masks):
        raise ValueError('Missing pilot payoff family')
    return np.column_stack([payoffs[:,mask].mean(axis=1) for mask in masks])


def run(args):
    start=time.perf_counter()
    folder=Path(__file__).parent
    scope_path=folder/'active_scope.json'
    pilot_folder=folder/'runs/m2_vix_importance_sampling_001'
    old_path=folder/'runs/m2_vix_curve_fit_001/run.json'
    paths=[Path(__file__),folder/'m2_spx.py',scope_path,old_path,pilot_folder/'run.json',
        pilot_folder/'outer_values.csv',*sorted((M2/'src/esn_eq').glob('*.py'))]
    sources=capture_sources(paths)
    scope=json.loads(scope_path.read_text())
    if scope['vix_fit_expiry_days'] != [27,55,118,153] or scope['vix_deferred_expiry_days']:
        raise ValueError('Require restored all-four VIX scope')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    architecture,parameters,premium,frozen_hashes=frozen_stage3()
    dated=DatedM2Pricer(architecture,premium)
    initial,_=dated.initial_state([0.,0.,0.])
    old=json.loads(old_path.read_text()); pilot=json.loads((pilot_folder/'run.json').read_text())
    expected=json.loads(json.dumps(dict(architecture=asdict(architecture),parameters=asdict(parameters),
        premium=asdict(premium),initial_state=initial.tolist(),target_sha256=targets.artifact_sha256)))
    for record in (old,pilot):
        if any(record[k]!=v for k,v in expected.items()):
            raise ValueError('Pilot model/state/targets changed')
    if pilot['curve_levels']!=old['selected_levels'] or (pilot['conditional_order'],pilot['normalizer_order'],pilot['time_order'])!=(32,48,8):
        raise ValueError('Pilot curve or quadrature changed')
    curve=NuisanceCurveFamily.from_targets(targets).build(old['selected_levels'])
    engine=ConditionalVIXPricer(dated,parameters,curve,initial,normalizer_order=48)
    raw_pilot=pd.read_csv(pilot_folder/'outer_values.csv')
    proposals,proposal_records,pilot_arrays={}, {}, {}
    # Freeze every proposal before generating any fresh evaluation draw.
    for days in scope['vix_fit_expiry_days']:
        q=targets.select('VIX',expiry_days=[days]).reset_index(drop=True)
        previous=DefensiveGaussianMixture.symmetric_leading(len(initial))
        zs,ws,vs=[],[],[]
        for seed_record in pilot['seed_records']:
            if seed_record['days']!=days or seed_record['method']!='mixture':
                continue
            frame=raw_pilot[(raw_pilot.days==days)&(raw_pilot.method=='mixture')&
                (raw_pilot.replicate==seed_record['replicate'])].sort_values('point')
            z,w,_=previous.draw(max(pilot['arguments']['powers']),seed_record['seed'])
            np.testing.assert_allclose(w,frame.weight.to_numpy(),rtol=1e-12,atol=1e-14)
            zs.append(z); ws.append(w); vs.append(frame.vix.to_numpy())
        z,w,v=np.concatenate(zs),np.concatenate(ws),np.concatenate(vs)
        h=pilot_targets(v,q)
        proposal,diagnostic=fit_payoff_mixture(z,w,h,coordinates=8,max_shift=3.)
        proposals[days]=proposal
        proposal_records[str(days)]=dict(shifts=proposal.shifts.tolist(),probabilities=proposal.probabilities.tolist(),
            target_families=['mean scaled put payoffs','mean scaled call payoffs','mean scaled upper-strike-quartile call payoffs'],
            pilot_points=len(v),diagnostics=diagnostic)
        pilot_arrays[f'z_{days}']=z; pilot_arrays[f'w_{days}']=w; pilot_arrays[f'h_{days}']=h
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    (args.output/'proposals.json').write_text(json.dumps(proposal_records,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(args.output/'pilot_arrays.npz',**pilot_arrays)
    proposal_hash=hashlib.sha256((args.output/'proposals.json').read_bytes()).hexdigest()
    (args.output/'status.json').write_text(json.dumps(dict(status='running; proposals frozen before evaluation',proposal_sha256=proposal_hash)))
    rows,rep_rows,raw,weights,seeds,changes,refinements=[],[],[],[],[],[],[]
    for days in scope['vix_fit_expiry_days']:
        quotes=targets.select('VIX',expiry_days=[days]).reset_index(drop=True)
        expiry=float(quotes.T_years.iloc[0])
        sampler=GaussianOuterSampler(dated,expiry,initial)
        window=engine.gaussian_window(expiry,time_order=8)
        strikes=quotes.strike.to_numpy(); sign=np.where(quotes.option_type.to_numpy()=='call',1.,-1.)
        descriptors=[dict(metric='forward',strike=None,option_type=None,scale=.5),
            dict(metric='second_moment',strike=None,option_type=None,scale=None)]
        descriptors += [dict(metric='option',strike=float(q.strike),option_type=q.option_type,scale=float(q.price_scale_U)) for q in quotes.itertuples()]
        for method_id,method in enumerate(('sobol','pilot_mixture')):
            estimates={p:[] for p in args.powers}; count={p:[] for p in args.powers}
            shares={p:[] for p in args.powers}; squares={p:[] for p in args.powers}
            maxima={p:[] for p in args.powers}
            for replicate in range(args.replications):
                tick=time.perf_counter()
                seed=int(np.random.SeedSequence([args.seed,days,method_id,replicate]).generate_state(1)[0])
                if method=='pilot_mixture':
                    z,w,component=proposals[days].draw(max(args.powers),seed)
                    states=sampler.mean+z@sampler.root.T
                else:
                    states=sampler.draw(max(args.powers),seed)
                    w=np.ones(len(states)); component=np.zeros(len(states),int)
                v=engine.evaluate(states,window,order=32).vix
                payoff=np.maximum(sign*(v[:,None]-strikes),0.)
                contribution=w[:,None]*np.column_stack((v,v**2,payoff))
                for power in args.powers:
                    n=2**power; c=contribution[:n]; sums=c.sum(axis=0)
                    estimates[power].append(c.mean(axis=0)); count[power].append((payoff[:n]>0).sum(axis=0))
                    squares[power].append((c*c).sum(axis=0)); maxima[power].append(c.max(axis=0))
                    shares[power].append(np.divide(c.max(axis=0),sums,out=np.full(len(sums),np.nan),where=sums>0))
                    weights.append(dict(days=days,method=method,replicate=replicate,power=power,
                        mean=float(w[:n].mean()),minimum=float(w[:n].min()),maximum=float(w[:n].max())))
                for j,desc in enumerate(descriptors):
                    for power in args.powers:
                        rep_rows.append(dict(days=days,method=method,replicate=replicate,power=power,**desc,
                            estimate=float(estimates[power][-1][j])))
                raw.append(pd.DataFrame(dict(days=days,method=method,replicate=replicate,point=np.arange(len(v)),vix=v,weight=w,component=component)))
                seeds.append(dict(days=days,method=method,replicate=replicate,seed=seed))
                # Check conditional quadrature at predetermined first-replication extreme values.
                # Diagnostic selection only; these results never modify the frozen proposals.
                if method=='pilot_mixture' and replicate==0:
                    indices=np.argsort(v)[[0,len(v)//2,-1]]
                    high=engine.evaluate(states[indices],window,order=48).vix
                    refinements.append(dict(days=days,base_vix=v[indices].tolist(),conditional48_vix=high.tolist(),
                        max_change=float(np.max(np.abs(high-v[indices])))))
                print(json.dumps(dict(days=days,method=method,replicate=replicate+1,seconds=time.perf_counter()-tick)),flush=True)
            for power in args.powers:
                matrix=np.array(estimates[power]); mean,se=replication_statistics(matrix)
                total=matrix.sum(axis=0)*(2**power)
                sq=np.sum(squares[power],axis=0); maximum=np.max(maxima[power],axis=0)
                for j,desc in enumerate(descriptors):
                    hits=np.array(count[power])[:,j-2] if j>=2 else None
                    share=np.array(shares[power])[:,j]
                    rows.append(dict(days=days,method=method,power=power,**desc,estimate=float(mean[j]),replication_se=float(se[j]),
                        se_over_scale=float(se[j]/desc['scale']) if desc['scale'] else None,
                        positive_payoff_count=int(hits.sum()) if hits is not None else None,
                        zero_hit_replications=int(np.sum(hits==0)) if hits is not None else None,
                        max_within_replication_share=float(np.nanmax(share)) if np.isfinite(share).any() else None,
                        max_pooled_share=float(maximum[j]/total[j]) if total[j]>0 else None,
                        payoff_contribution_ess=float(total[j]**2/sq[j]) if sq[j]>0 else None))
            change=np.array(estimates[max(args.powers)])-np.array(estimates[min(args.powers)])
            cm,cs=replication_statistics(change)
            for j,desc in enumerate(descriptors):
                changes.append(dict(days=days,method=method,**desc,change=float(cm[j]),paired_replication_se=float(cs[j])))
            # Persist completed blocks without treating them as a complete study.
            pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False)
            pd.DataFrame(rep_rows).to_csv(args.output/'replications.csv',index=False)
            pd.concat(raw,ignore_index=True).to_csv(args.output/'outer_values.csv',index=False)
    check_sources(sources)
    if hashlib.sha256((args.output/'proposals.json').read_bytes()).hexdigest()!=proposal_hash:
        raise RuntimeError('Proposal changed during evaluation')
    record=dict(status='complete fixed-candidate pilot-proposal check; numerical gates unresolved',
        arguments=vars(args)|dict(output=str(args.output)),scope=scope,architecture=asdict(architecture),
        parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),curve_levels=old['selected_levels'],
        target_sha256=targets.artifact_sha256,frozen_model_sha256=frozen_hashes,implementation_sha256=source_hashes(sources,ROOT),
        proposal_sha256=proposal_hash,seed_records=seeds,conditional_order=32,normalizer_order=48,time_order=8,
        conditional_refinement=refinements,input_compatibility=targets.audit['variance_audit'],
        uncertainty_unit='independent complete scrambles, conditional on one frozen pilot proposal',
        pilot_policy='Only old mixture draws promoted to pilot; old Sobol results remain historical, not fresh validation. All new evaluation seeds are independent.',
        limitations=['No refit or held-out market evaluation.','Only one pilot was fitted; pilot-to-pilot variation unmeasured.',
            'Eight-scramble SEs are exploratory; zero payoffs/SEs do not certify rare-event prices.',
            'Outer SE excludes GH, time and normalizer biases; three-state conditional refinements are diagnostics only.',
            'Leading-eight-coordinate heuristic may miss tail directions; no universal variance reduction claim.',
            'No self-normalized evaluation weights or forward rescaling; contribution ESS is not an accuracy certificate.'],
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),wall_seconds=time.perf_counter()-start)
    pd.DataFrame(weights).to_csv(args.output/'weight_diagnostics.csv',index=False)
    pd.DataFrame(changes).to_csv(args.output/'prefix_changes.csv',index=False)
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    (args.output/'status.json').write_text(json.dumps(dict(status='complete',proposal_sha256=proposal_hash))+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--powers',type=int,nargs=2,default=[6,7])
    parser.add_argument('--replications',type=int,default=8)
    parser.add_argument('--seed',type=int,default=20261027)
    args=parser.parse_args()
    if args.output.exists() or args.replications<4 or not 3<=args.powers[0]<args.powers[1]:
        parser.error('Require new output, >=4 replications and increasing powers >=3')
    run(args)
