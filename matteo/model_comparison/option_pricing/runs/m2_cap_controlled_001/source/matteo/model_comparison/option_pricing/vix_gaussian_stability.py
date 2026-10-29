"""Larger independent Gaussian-control pilots and fresh 128/256 VIX checking."""
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
from esn_eq.gaussian_controls import gaussian_control_rows,fit_replication_controls
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.outer_sampling import GaussianOuterSampler,replication_statistics
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import DatedM2Pricer
from esn_eq.vix_pricing import ConditionalVIXPricer


METHODS=('raw','mass','gaussian','gaussian_batch_a','gaussian_batch_b','legacy12')


def payoffs(v,quotes):
    sign=np.where(quotes.option_type.to_numpy()=='call',1.,-1.)
    return np.column_stack((v,v*v,np.maximum(sign*(v[:,None]-quotes.strike.to_numpy()),0.)))


def descriptors(quotes):
    return [dict(metric='forward',strike=None,option_type=None,scale=.5),dict(metric='second_moment',strike=None,option_type=None,scale=None)]+[dict(metric='option',strike=float(q.strike),option_type=q.option_type,scale=float(q.price_scale_U)) for q in quotes.itertuples()]


def balanced_pilot_split(counts):
    if any(not isinstance(n,int) or n<2 or n%2 for n in counts):raise ValueError('Require positive even batch sizes')
    first=[];second=[];offset=0
    for n in counts:
        first.extend(range(offset,offset+n//2));second.extend(range(offset+n//2,offset+n));offset+=n
    return np.array(first),np.array(second)


def fit_variants(Y,C,legacy,ridge,counts=(12,16)):
    if len(Y)!=sum(counts) or len(C)!=sum(counts):raise ValueError('Pilot batch sizes differ')
    a,b=balanced_pilot_split(counts)
    mass=np.zeros((C.shape[1],Y.shape[1]));mass[:1]=fit_replication_controls(Y,C[:,:1],ridge)
    return dict(raw=np.zeros_like(mass),mass=mass,gaussian=fit_replication_controls(Y,C,ridge),gaussian_batch_a=fit_replication_controls(Y[a],C[a],ridge),gaussian_batch_b=fit_replication_controls(Y[b],C[b],ridge),legacy12=np.array(legacy,float,copy=True))


def run(args):
    start=time.perf_counter();folder=Path(__file__).parent;pilot=folder/'runs/m2_vix_targeted_allocation_001';recent=folder/'runs/m2_vix_gaussian_controls_001';basis_path=folder/'runs/m2_vix_score_sampling_001/proposals.json';scope_path=folder/'active_scope.json'
    paths=[*[recent/n for n in ('run.json','proposals.json','coefficients.json','outer_values.csv','states.npz')],Path(__file__),folder/'m2_spx.py',scope_path,basis_path,*[pilot/n for n in ('run.json','proposals.json','outer_values.csv','states.npz')],*sorted((M2/'src/esn_eq').glob('*.py'))]
    sources=capture_sources(paths);scope=json.loads(scope_path.read_text());prior=json.loads((pilot/'run.json').read_text());basis_packets=json.loads(basis_path.read_text())
    if scope['vix_fit_expiry_days']!=[27,55,118,153] or scope['vix_deferred_expiry_days']:raise ValueError('Require all four maturities')
    target=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');arch,params,premium,hashes=frozen_stage3();dated=DatedM2Pricer(arch,premium);initial,_=dated.initial_state([0.,0.,0.])
    expected=json.loads(json.dumps(dict(architecture=asdict(arch),parameters=asdict(params),premium=asdict(premium),initial_state=initial.tolist(),target_sha256=target.artifact_sha256)))
    if any(prior[k]!=v for k,v in expected.items()) or (prior['conditional_order'],prior['normalizer_order'],prior['time_order'])!=(32,48,8):raise ValueError('Pilot configuration mismatch')
    proposal_bytes=(pilot/'proposals.json').read_bytes()
    if hashlib.sha256(proposal_bytes).hexdigest()!=prior['proposal_sha256']:raise ValueError('Pilot proposal changed')
    proposals=json.loads(proposal_bytes)
    recent_record=json.loads((recent/'run.json').read_text());legacy_packets=json.loads((recent/'coefficients.json').read_text())
    for key in (*expected.keys(),'curve_levels','conditional_order','normalizer_order','time_order','proposal_sha256'):
        if recent_record[key]!=prior[key]:raise ValueError('Pilot sources differ')
    if (recent/'proposals.json').read_bytes()!=proposal_bytes:raise ValueError('Pilot proposals differ')
    if hashlib.sha256((recent/'coefficients.json').read_bytes()).hexdigest()!=recent_record['coefficient_sha256']:raise ValueError('Legacy coefficients changed')
    raw_sources=[pd.read_csv(src/'outer_values.csv') for src in (pilot,recent)]
    state_sources=[np.load(src/'states.npz') for src in (pilot,recent)]
    source_records=[prior,recent_record]
    counts=[prior['arguments']['replications'],recent_record['arguments']['replications']]
    if counts!=[12,16] or any(max(rec['arguments']['powers'])!=7 for rec in source_records):raise ValueError('Require 12+16 compatible 128-point pilot replications')
    curve=NuisanceCurveFamily.from_targets(target).build(prior['curve_levels']);engine=ConditionalVIXPricer(dated,params,curve,initial,normalizer_order=48)
    packets={};pilot_arrays={};samplers={};windows={};coefficients={};pilot_seeds=[]
    # Every direction, ridge, coefficient and alternative is frozen before checking.
    for days in scope['vix_fit_expiry_days']:
        q=target.select('VIX',expiry_days=[days]).reset_index(drop=True);expiry=float(q.T_years.iloc[0]);sampler=GaussianOuterSampler(dated,expiry,initial)
        proposal=DefensiveGaussianMixture(proposals[str(days)]['shifts'],proposals[str(days)]['probabilities']);U=np.array(basis_packets[str(days)]['basis']);Y=[];C=[];zs=[];ws=[];vs=[]
        for source_index,(rec0,raw0,states0) in enumerate(zip(source_records,raw_sources,state_sources)):
            for rep in range(rec0['arguments']['replications']):
                rec=next(s for s in rec0['seed_records'] if s['days']==days and s['replicate']==rep and (source_index==1 or s['method']=='random_targeted'))
                z,w,_=proposal.draw(7,rec['seed']);pilot_seeds.append(rec['seed'])
                frame=raw0[(raw0.days==days)&(raw0.replicate==rep)]
                if source_index==0:frame=frame[frame.method=='random_targeted']
                frame=frame.sort_values('point');np.testing.assert_allclose(w,frame.weight,rtol=1e-12,atol=1e-14)
                n=len(z);saved=states0[f'random_targeted_{days}'] if source_index==0 else states0[str(days)]
                np.testing.assert_allclose(saved[rep*n:(rep+1)*n],sampler.mean+z@sampler.root.T,rtol=1e-12,atol=1e-14)
                v=frame.vix.to_numpy();Y.append((w[:,None]*payoffs(v,q)).mean(axis=0));C.append(gaussian_control_rows(z,w,U).mean(axis=0));zs.append(z);ws.append(w);vs.append(v)
        Y,C=np.array(Y),np.array(C);np.testing.assert_array_equal(U,legacy_packets[str(days)]['basis']);legacy=np.array(legacy_packets[str(days)]['coefficients']['gaussian']);np.testing.assert_allclose(legacy,fit_replication_controls(Y[:12],C[:12],args.ridge),rtol=1e-12,atol=1e-12);beta=fit_variants(Y,C,legacy,args.ridge,counts)
        packets[str(days)]=dict(basis=U.tolist(),coefficients={k:v.tolist() for k,v in beta.items()},pilot_replications=len(Y),pilot_batch_sizes=counts,pilot_split_indices=[a.tolist() for a in balanced_pilot_split(counts)],pilot_points_per_replication=n,ridge=args.ridge)
        coefficients[days]=beta;samplers[days]=sampler;windows[days]=engine.gaussian_window(expiry,time_order=8)
        for key,value in dict(outcomes=Y,controls=C,normals=np.array(zs),weights=np.array(ws),vix=np.array(vs)).items():pilot_arrays[f'{key}_{days}']=value
    args.output.mkdir(parents=True,exist_ok=False);save_sources(args.output,sources,ROOT);(args.output/'proposals.json').write_bytes(proposal_bytes)
    (args.output/'coefficients.json').write_text(json.dumps(packets,indent=2,allow_nan=False)+'\n');coef_sha=hashlib.sha256((args.output/'coefficients.json').read_bytes()).hexdigest();np.savez_compressed(args.output/'pilot_arrays.npz',**pilot_arrays)
    (args.output/'status.json').write_text('{"status":"running; controls frozen"}\n')
    rows=[];rep_rows=[];changes=[];raw=[];seeds=[];state_archive={};control_archive={};refinements=[]
    for days in scope['vix_fit_expiry_days']:
        q=target.select('VIX',expiry_days=[days]).reset_index(drop=True);desc=descriptors(q);sampler=samplers[days];window=windows[days]
        proposal=DefensiveGaussianMixture(proposals[str(days)]['shifts'],proposals[str(days)]['probabilities']);U=np.array(packets[str(days)]['basis']);beta=coefficients[days]
        matrices={m:{p:[] for p in args.powers} for m in METHODS};hits={p:[] for p in args.powers};saved_states=[];saved_controls=[]
        for rep in range(args.replications):
            tick=time.perf_counter();seed=int(np.random.SeedSequence([args.seed,days,rep]).generate_state(1)[0])
            if seed in pilot_seeds:raise ValueError('Checking seed reuses pilot')
            z,w,component=proposal.draw(max(args.powers),seed);states=sampler.mean+z@sampler.root.T;v=engine.evaluate(states,window,order=32).vix
            h=payoffs(v,q);weighted=w[:,None]*h;C=gaussian_control_rows(z,w,U);saved_states.append(states);saved_controls.append(C)
            for p in args.powers:
                n=2**p;hits[p].append((h[:n,2:]>0).sum(axis=0));y=weighted[:n].mean(axis=0);c=C[:n].mean(axis=0)
                for method in METHODS:
                    adjusted=y-c@beta[method];matrices[method][p].append(adjusted)
                    for j,d in enumerate(desc):rep_rows.append(dict(days=days,method=method,power=p,replicate=rep,**d,estimate=float(adjusted[j])))
            raw.append(pd.DataFrame(dict(days=days,replicate=rep,point=np.arange(len(v)),vix=v,weight=w,component=component)))
            seeds.append(dict(days=days,replicate=rep,seed=seed))
            if rep==0:
                ix=np.argsort(v)[[0,len(v)//2,-1]];high=engine.evaluate(states[ix],window,order=48).vix
                refinements.append(dict(days=days,base_vix=v[ix].tolist(),conditional48_vix=high.tolist(),max_change=float(np.max(abs(high-v[ix])))))
            print(json.dumps(dict(days=days,replicate=rep+1,replications=args.replications,seconds=time.perf_counter()-tick)),flush=True)
        state_archive[str(days)]=np.concatenate(saved_states);control_archive[str(days)]=np.concatenate(saved_controls)
        for method in METHODS:
            for p in args.powers:
                x=np.array(matrices[method][p]);mean,se=replication_statistics(x);raw_matrix=np.array(matrices['raw'][p]);dm,ds=replication_statistics(x-raw_matrix)
                for j,d in enumerate(desc):
                    count=np.array(hits[p])[:,j-2] if j>=2 else None
                    rows.append(dict(days=days,method=method,power=p,**d,estimate=float(mean[j]),replication_se=float(se[j]),se_over_scale=float(se[j]/d['scale']) if d['scale'] else None,change_from_raw=float(dm[j]),paired_change_se=float(ds[j]),se_ratio_to_raw=float(se[j]/raw_matrix[:,j].std(ddof=1)*np.sqrt(args.replications)) if raw_matrix[:,j].std(ddof=1)>0 else None,positive_payoff_count=int(count.sum()) if count is not None else None,zero_hit_replications=int(np.sum(count==0)) if count is not None else None,negative_estimate=bool(mean[j]<0),model_moment_target=engine.moment_target(window) if j==1 else None))
            delta=np.array(matrices[method][max(args.powers)])-np.array(matrices[method][min(args.powers)]);dm,ds=replication_statistics(delta)
            for j,d in enumerate(desc):changes.append(dict(days=days,method=method,**d,change=float(dm[j]),paired_replication_se=float(ds[j])))
        pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False);pd.DataFrame(rep_rows).to_csv(args.output/'replications.csv',index=False);pd.concat(raw,ignore_index=True).to_csv(args.output/'outer_values.csv',index=False);np.savez_compressed(args.output/'states.npz',**state_archive);np.savez_compressed(args.output/'control_rows.npz',**control_archive)
    check_sources(sources)
    if hashlib.sha256((args.output/'coefficients.json').read_bytes()).hexdigest()!=coef_sha:raise RuntimeError('Controls changed during checking')
    record=dict(status='fresh paired VIX Gaussian-control diagnostic; no refit, rescaling or clipping',arguments=vars(args)|dict(output=str(args.output)),scope=scope,**expected,curve_levels=prior['curve_levels'],frozen_model_sha256=hashes,proposal_sha256=prior['proposal_sha256'],coefficient_sha256=coef_sha,implementation_sha256=source_hashes(sources,ROOT),pilot_source=str(pilot),pilot_seed_records=pilot_seeds,seed_records=seeds,conditional_order=32,normalizer_order=48,time_order=8,conditional_refinement=refinements,methods=list(METHODS),pilot_sources=[str(pilot),str(recent)],legacy_coefficient_sha256=recent_record['coefficient_sha256'],primary_method='gaussian',control_definitions=['w-1','w*(U transpose z), four coordinates','w*((U transpose z)^2-1), four coordinates'],pilot_policy='Historical R12 and R16 raw 128-point checks explicitly promoted to design data. Pooled R28 primary, disjoint balanced R14 batches diagnostic; original R12 fit retained unchanged. Fixed ridge and dictionary, no checking-driven selection.',uncertainty_unit='independent complete checking scrambles, conditional on frozen pilot coefficients; variants and prefixes paired on identical paths',limitations=['Gaussian identities do not force VIX squared to its exact-model mean or remove finite conditional/time/normalizer bias.','Pooled and balanced pilot fits share training observations; original R12 legacy fit uses part of this pilot. Shared checking paths make comparisons paired and dependent.','No clipping of negative controlled estimates; small SE alone does not certify rare-event accuracy.','Same coefficients used at 128/256 counts; pilots are 128-point replications. Effectiveness need not transfer to the 256-point net. R12 checking has fewer uncertainty degrees of freedom than prior R16 study.','No model parameter fit, maturity exclusion, market held-out claim or Lean certificate.'],versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),wall_seconds=time.perf_counter()-start)
    pd.DataFrame(changes).to_csv(args.output/'prefix_changes.csv',index=False);(args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');(args.output/'status.json').write_text('{"status":"complete"}\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--powers',type=int,nargs=2,default=[7,8]);parser.add_argument('--replications',type=int,default=12);parser.add_argument('--ridge',type=float,default=1.);parser.add_argument('--seed',type=int,default=2026102205);args=parser.parse_args()
    if args.output.exists() or args.replications<4 or not 4<=args.powers[0]<args.powers[1] or args.ridge<=0:parser.error('Require new output, >=4 replications, increasing counts and positive ridge')
    run(args)
