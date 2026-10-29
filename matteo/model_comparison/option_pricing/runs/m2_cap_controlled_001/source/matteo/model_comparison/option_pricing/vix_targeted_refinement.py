"""Exact saved-state attribution of the targeted 118-day VIX moment warning."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from .vix_matched_moment_refinement import expected_nodes,values_from_nodes
from esn_eq.pricing import DatedM2Pricer
from esn_eq.vix_pricing import ConditionalVIXPricer
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.outer_sampling import replication_statistics


def run(output):
    start=time.perf_counter();folder=Path(__file__).parent;old=folder/'runs/m2_vix_targeted_check_001'
    paths=[Path(__file__),folder/'m2_spx.py',folder/'vix_matched_moment_refinement.py',folder/'active_scope.json',*[old/n for n in ('run.json','states.npz','outer_values.csv')],*sorted((M2/'src/esn_eq').glob('*.py'))]
    sources=capture_sources(paths);prior=json.loads((old/'run.json').read_text());scope=json.loads((folder/'active_scope.json').read_text())
    if scope['vix_fit_expiry_days']!=[27,55,118,153] or scope['vix_deferred_expiry_days']:raise ValueError('Require retained full scope')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');arch,params,premium,_=frozen_stage3();dated=DatedM2Pricer(arch,premium);initial,_=dated.initial_state([0.,0.,0.])
    expected=json.loads(json.dumps(dict(architecture=asdict(arch),parameters=asdict(params),premium=asdict(premium),initial_state=initial.tolist(),target_sha256=targets.artifact_sha256)))
    if any(prior[k]!=v for k,v in expected.items()):raise ValueError('Frozen configuration changed')
    curve=NuisanceCurveFamily.from_targets(targets).build(prior['curve_levels']);engine=ConditionalVIXPricer(dated,params,curve,initial,normalizer_order=48)
    q=targets.select('VIX',expiry_days=[118]).reset_index(drop=True);expiry=float(q.T_years.iloc[0]);R=prior['arguments']['replications'];N=2**max(prior['arguments']['powers'])
    frame=pd.read_csv(old/'outer_values.csv');frame=frame[(frame.days==118)&(frame.method=='targeted_mean')].sort_values(['replicate','point'])
    states=np.load(old/'states.npz')['targeted_mean_118'];w=frame.weight.to_numpy();base=frame.vix.to_numpy()
    if len(states)!=R*N or len(frame)!=R*N:raise ValueError('Incomplete saved block')
    output.mkdir(parents=True,exist_ok=False);save_sources(output,sources,ROOT);(output/'status.json').write_text('{"status":"running"}\n')
    window=engine.gaussian_window(expiry,time_order=8);norm64=dated.normalizer(window.times,initial,params,64)
    gh48=expected_nodes(engine,states,window,48,'118d_all512_GH48')
    rules={'GH32_M48_time8':base,'GH48_M48_time8':values_from_nodes(engine,window,gh48,window.normalizer),'GH48_M64_time8':values_from_nodes(engine,window,gh48,norm64)}
    h=lambda v:np.column_stack((v,v*v,np.maximum(np.where(q.option_type.to_numpy()=='call',1.,-1.)*(v[:,None]-q.strike.to_numpy()),0.)))
    descriptors=[dict(metric='forward',strike=None,option_type=None,scale=.5),dict(metric='second_moment',strike=None,option_type=None,scale=None)]
    descriptors += [dict(metric='option',strike=float(x.strike),option_type=x.option_type,scale=float(x.price_scale_U)) for x in q.itertuples()]
    rows=[];reps=[];target=engine.moment_target(window)
    for rule,v in rules.items():
        values=(w[:,None]*h(v)).reshape(R,N,-1).mean(axis=1);delta=(w[:,None]*(h(v)-h(base))).reshape(R,N,-1).mean(axis=1)
        mean,se=replication_statistics(values);dm,ds=replication_statistics(delta)
        for j,desc in enumerate(descriptors):
            rows.append(dict(rule=rule,**desc,estimate=float(mean[j]),replication_se=float(se[j]),change_from_base=float(dm[j]),paired_change_se=float(ds[j]),moment_target=target if j==1 else None,residual=float(mean[j]-target) if j==1 else None))
            for rep in range(R):reps.append(dict(rule=rule,replicate=rep,**desc,estimate=float(values[rep,j]),change=float(delta[rep,j])))
    checks=[];archive=dict(states=states,weight=w,GH48=gh48,times=window.times,time_weights=window.weights,xi=window.xi,normalizer48=window.normalizer,normalizer64=norm64,**rules)
    for label,points,time_order,order in [('time16',8,16,48),('conditional64',4,8,64)]:
        ix=np.arange(R*N).reshape(R,N)[:,:points].ravel();window2=engine.gaussian_window(expiry,time_order=time_order);m64=dated.normalizer(window2.times,initial,params,64)
        e=expected_nodes(engine,states[ix],window2,order,f'118d_{label}')
        v=values_from_nodes(engine,window2,e,m64);ref=rules['GH48_M64_time8'][ix]
        delta=(w[ix]*(v*v-ref*ref)).reshape(R,points).mean(axis=1);dm,ds=replication_statistics(delta)
        checks.append(dict(check=label,points_per_replication=points,weighted_moment_change=float(dm),paired_change_se=float(ds),max_pointwise_vix_change=float(np.max(abs(v-ref)))))
        archive[label+'_vix']=v;archive[label+'_indices']=ix
    check_sources(sources);pd.DataFrame(rows).to_csv(output/'summary.csv',index=False);pd.DataFrame(reps).to_csv(output/'replications.csv',index=False);np.savez_compressed(output/'nodes_and_rules.npz',**archive)
    record=dict(status='saved 118-day targeted state refinement, no refit or rescaling',scope=scope,**expected,curve_levels=prior['curve_levels'],outer_source=str(old),points_per_replication=N,replications=R,investigated_expiry_days=118,checks=checks,delta=window.delta,variance_floor=engine.variance_floor,implementation_sha256=source_hashes(sources,ROOT),limitations=['All 512 offending 118d states retained with exact original weights; no new outer observations.', 'Time16 uses first eight states/scramble; conditional64 first four only. These cannot certify global bias.', 'Other three maturities remain active but are not newly refined in this attribution study.'],wall_seconds=time.perf_counter()-start)
    (output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');(output/'status.json').write_text('{"status":"complete"}\n');print(pd.DataFrame(rows).query("metric == 'second_moment'").to_string(index=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():parser.error('Output exists')
    run(args.output.resolve())
