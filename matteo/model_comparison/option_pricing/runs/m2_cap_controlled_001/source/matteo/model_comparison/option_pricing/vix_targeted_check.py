"""Fresh exploratory VIX check of independently designed targeted mean shifts."""
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
from esn_eq.outer_sampling import GaussianOuterSampler,replication_statistics
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import DatedM2Pricer
from esn_eq.vix_pricing import ConditionalVIXPricer


def run(args):
    start=time.perf_counter();folder=Path(__file__).parent;prior_folder=folder/'runs/m2_vix_stratified_sampling_001'
    reference_folder=folder/'runs/m2_vix_collapsed_inflated_001';scope_path=folder/'active_scope.json'
    design_folder=folder/'runs/m2_vix_targeted_pilots_001'
    sources=capture_sources([design_folder/'run.json',design_folder/'proposals.json',Path(__file__),folder/'m2_spx.py',scope_path,prior_folder/'run.json',prior_folder/'proposals.json',reference_folder/'run.json',reference_folder/'summary.csv',*sorted((M2/'src/esn_eq').glob('*.py'))])
    scope=json.loads(scope_path.read_text());prior=json.loads((prior_folder/'run.json').read_text());reference=json.loads((reference_folder/'run.json').read_text())
    if scope['vix_fit_expiry_days']!=[27,55,118,153] or scope['vix_deferred_expiry_days']:raise ValueError('Require full restored scope')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');architecture,parameters,premium,hashes=frozen_stage3()
    dated=DatedM2Pricer(architecture,premium);initial,_=dated.initial_state([0.,0.,0.])
    expected=json.loads(json.dumps(dict(architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),target_sha256=targets.artifact_sha256)))
    for record in (prior,reference):
        if any(record[k]!=v for k,v in expected.items()):raise ValueError('Frozen model or target changed')
    if prior['curve_levels']!=reference['curve_levels']:raise ValueError('Moment-reference curve differs')
    curve=NuisanceCurveFamily.from_targets(targets).build(prior['curve_levels']);engine=ConditionalVIXPricer(dated,parameters,curve,initial,normalizer_order=48)
    old_bytes=(prior_folder/'proposals.json').read_bytes()
    if hashlib.sha256(old_bytes).hexdigest()!=prior['proposal_sha256']:raise ValueError('Old proposal changed')
    old=json.loads(old_bytes);packets={};samplers={};windows={};proposals={}
    design=json.loads((design_folder/'run.json').read_text())
    design_bytes=(design_folder/'proposals.json').read_bytes()
    if hashlib.sha256(design_bytes).hexdigest()!=design['proposal_sha256']:raise ValueError('Design proposal changed')
    for k in (*expected.keys(),'curve_levels','conditional_order','normalizer_order','time_order'):
        if design[k]!=prior[k]:raise ValueError('Design and checking configuration differ')
    fitted=json.loads(design_bytes)
    # Single candidate and pooled pilot fixed before any fresh pricing draw.
    for days in scope['vix_fit_expiry_days']:
        expiry=float(targets.select('VIX',expiry_days=[days]).T_years.iloc[0])
        samplers[days]=GaussianOuterSampler(dated,expiry,initial)
        windows[days]=engine.gaussian_window(expiry,time_order=8)
        packet=fitted[str(days)]['score4_balanced_cap4']['4'];packets[str(days)]=packet
        proposals[days]=DefensiveGaussianMixture(packet['shifts'],packet['probabilities'])
    args.output.mkdir(parents=True,exist_ok=False);save_sources(args.output,sources,ROOT)
    (args.output/'proposals.json').write_text(json.dumps(packets,indent=2,allow_nan=False)+'\n');(args.output/'pilot_proposals.json').write_bytes(old_bytes)
    proposal_sha=hashlib.sha256((args.output/'proposals.json').read_bytes()).hexdigest();(args.output/'status.json').write_text('{"status":"running; proposals frozen"}\n')
    rows,rep_rows,raw,seeds,changes,refinements=[],[],[],[],[],[];state_archive={}
    for days in scope['vix_fit_expiry_days']:
        quotes=targets.select('VIX',expiry_days=[days]).reset_index(drop=True);window=windows[days];sampler=samplers[days]
        base=DefensiveGaussianMixture(old[str(days)]['shifts'],old[str(days)]['probabilities'])
        sign=np.where(quotes.option_type.to_numpy()=='call',1.,-1.);strikes=quotes.strike.to_numpy()
        target=engine.moment_target(window)
        descriptors=[dict(metric='forward',strike=None,option_type=None,scale=.5),dict(metric='second_moment',strike=None,option_type=None,scale=None)]
        descriptors += [dict(metric='option',strike=float(q.strike),option_type=q.option_type,scale=float(q.price_scale_U)) for q in quotes.itertuples()]
        for method_id,method in enumerate(('learned_mean','targeted_mean')):
            estimates={p:[] for p in args.powers};hits={p:[] for p in args.powers};maxima={p:[] for p in args.powers};squares={p:[] for p in args.powers};states_saved=[]
            for replicate in range(args.replications):
                tick=time.perf_counter();seed=int(np.random.SeedSequence([args.seed,days,method_id,replicate]).generate_state(1)[0])
                if method=='learned_mean':z,w,component=base.draw(max(args.powers),seed)
                else:
                    z,w,component=proposals[days].draw(max(args.powers),seed)
                states=sampler.mean+z@sampler.root.T;states_saved.append(states)
                v=engine.evaluate(states,window,order=32).vix;h=np.maximum(sign*(v[:,None]-strikes),0.)
                contribution=w[:,None]*np.column_stack((v,v*v,h))
                for power in args.powers:
                    n=2**power;c=contribution[:n];estimates[power].append(c.mean(axis=0));hits[power].append((h[:n]>0).sum(axis=0));maxima[power].append(c.max(axis=0));squares[power].append((c*c).sum(axis=0))
                    for j,desc in enumerate(descriptors):rep_rows.append(dict(days=days,method=method,power=power,replicate=replicate,**desc,estimate=float(c[:,j].mean())))
                raw.append(pd.DataFrame(dict(days=days,method=method,replicate=replicate,point=np.arange(len(v)),vix=v,weight=w,component=component)))
                seeds.append(dict(days=days,method=method,replicate=replicate,seed=seed))
                if replicate==0:
                    ix=np.argsort(v)[[0,len(v)//2,-1]];high=engine.evaluate(states[ix],window,order=48).vix
                    refinements.append(dict(days=days,method=method,base_vix=v[ix].tolist(),conditional48_vix=high.tolist(),max_change=float(np.max(np.abs(high-v[ix])))))
                print(json.dumps(dict(days=days,method=method,replicate=replicate+1,replications=args.replications,seconds=time.perf_counter()-tick)),flush=True)
            state_archive[f'{method}_{days}']=np.concatenate(states_saved)
            for power,matrix in estimates.items():
                x=np.array(matrix);mean,se=replication_statistics(x);total=x.sum(axis=0)*2**power;maximum=np.max(maxima[power],axis=0);sq=np.sum(squares[power],axis=0)
                for j,desc in enumerate(descriptors):
                    count=np.array(hits[power])[:,j-2] if j>=2 else None
                    rows.append(dict(days=days,method=method,power=power,**desc,estimate=float(mean[j]),replication_se=float(se[j]),se_over_scale=float(se[j]/desc['scale']) if desc['scale'] else None,
                        positive_payoff_count=int(count.sum()) if count is not None else None,zero_hit_replications=int(np.sum(count==0)) if count is not None else None,
                        max_pooled_share=float(maximum[j]/total[j]) if total[j]>0 else None,payoff_contribution_ess=float(total[j]**2/sq[j]) if sq[j]>0 else None,
                        model_moment_target=target if j==1 else None))
            diff=np.array(estimates[max(args.powers)])-np.array(estimates[min(args.powers)]);dm,ds=replication_statistics(diff)
            for j,desc in enumerate(descriptors):changes.append(dict(days=days,method=method,**desc,change=float(dm[j]),paired_replication_se=float(ds[j])))
            pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False);pd.DataFrame(rep_rows).to_csv(args.output/'replications.csv',index=False);pd.concat(raw,ignore_index=True).to_csv(args.output/'outer_values.csv',index=False)
            np.savez_compressed(args.output/'states.npz',**state_archive)
    check_sources(sources)
    if hashlib.sha256((args.output/'proposals.json').read_bytes()).hexdigest()!=proposal_sha:raise RuntimeError('Proposal changed during evaluation')
    report=dict(status='fixed-model fresh outer-proposal comparison, no refit or moment rescaling',arguments=vars(args)|dict(output=str(args.output)),scope=scope,
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),curve_levels=prior['curve_levels'],target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,
        implementation_sha256=source_hashes(sources,ROOT),proposal_sha256=proposal_sha,seed_records=seeds,conditional_order=32,normalizer_order=48,time_order=8,conditional_refinement=refinements,
        moment_reference=str(reference_folder),design_source=str(design_folder),selected_candidate='score4_balanced_cap4',selected_pilot=4,design_proposal_sha256=design['proposal_sha256'],input_compatibility=targets.audit['variance_audit'],weight_upper_bound=2.,
        uncertainties='independent complete scrambles; nested count differences paired within each method, cross-method streams independent',
        limitations=['Exploratory R8x64 comparison, not a high-budget convergence or superiority certificate.',
            'Independent checking streams; count prefixes paired within each method only.',
            'Pilot IID second moments do not predict randomized Sobol uncertainty.',
            'One pooled pilot proposal checked; pilot-to-pilot variation was screened historically, not repriced on new draws.',
            'Moment agreement does not certify option/forward prices; no identity is imposed on estimates.',
            'Outer SE excludes conditional/time/normalizer bias; selected-state GH checks are limited diagnostics.',
            'Only one frozen model; no held-out market test, refit, or model ranking.'],
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),wall_seconds=time.perf_counter()-start)
    pd.DataFrame(changes).to_csv(args.output/'prefix_changes.csv',index=False);(args.output/'run.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');(args.output/'status.json').write_text('{"status":"complete"}\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--powers',type=int,nargs=2,default=[5,6]);parser.add_argument('--replications',type=int,default=8)
    parser.add_argument('--seed',type=int,default=2026102202);args=parser.parse_args()
    if args.output.exists() or args.replications<4 or not 4<=args.powers[0]<args.powers[1]:parser.error('Require new output, >=4 replications and increasing powers>=4')
    run(args)
