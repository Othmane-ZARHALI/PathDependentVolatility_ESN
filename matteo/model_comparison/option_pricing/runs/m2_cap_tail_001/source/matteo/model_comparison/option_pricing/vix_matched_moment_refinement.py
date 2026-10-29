"""Matched-state attribution of VIX moment warning, without changing samples."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.outer_sampling import GaussianOuterSampler,replication_statistics
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import DatedM2Pricer
from esn_eq.vix_pricing import ConditionalVIXPricer


def expected_nodes(engine,states,window,order,label):
    out=np.empty((len(states),len(window.times)))
    for node,t in enumerate(window.times):
        mean,cov=engine.dated.gaussian_projection_law(float(t-window.expiry),states)
        out[:,node]=engine.dated.gaussian_expectation(mean,cov,engine.parameters,order)
        if (node+1)%20==0:print(json.dumps(dict(label=label,node=node+1,nodes=len(window.times))),flush=True)
    return out


def values_from_nodes(engine,window,expected,normalizer):
    variance=(engine.variance_floor+(window.xi-engine.variance_floor)*expected/normalizer)@window.weights
    return 100*np.sqrt(variance/window.delta)


def run(args):
    start=time.perf_counter();folder=Path(__file__).parent;old=folder/'runs/m2_vix_stratified_sampling_001'
    source_paths=[Path(__file__),folder/'m2_spx.py',folder/'active_scope.json',old/'run.json',old/'proposals.json',old/'outer_values.csv',*sorted((M2/'src/esn_eq').glob('*.py'))]
    sources=capture_sources(source_paths);prior=json.loads((old/'run.json').read_text());scope=json.loads((folder/'active_scope.json').read_text())
    if scope['vix_fit_expiry_days']!=[27,55,118,153] or scope['vix_deferred_expiry_days']:raise ValueError('Require full scope')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');architecture,parameters,premium,hashes=frozen_stage3();dated=DatedM2Pricer(architecture,premium);initial,_=dated.initial_state([0.,0.,0.])
    expected=json.loads(json.dumps(dict(architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),target_sha256=targets.artifact_sha256)))
    if any(prior[k]!=v for k,v in expected.items()):raise ValueError('Frozen inputs changed')
    curve=NuisanceCurveFamily.from_targets(targets).build(prior['curve_levels']);engine=ConditionalVIXPricer(dated,parameters,curve,initial,normalizer_order=48)
    props=json.loads((old/'proposals.json').read_text());raw=pd.read_csv(old/'outer_values.csv');R=prior['arguments']['replications']
    rows,rep_rows,checks=[],[],[];stored={}
    args.output.mkdir(parents=True,exist_ok=False);save_sources(args.output,sources,ROOT);(args.output/'status.json').write_text('{"status":"running"}\n')
    for days in scope['vix_fit_expiry_days']:
        N=128 if days==27 else 8
        expiry=float(targets.select('VIX',expiry_days=[days]).T_years.iloc[0]);sampler=GaussianOuterSampler(dated,expiry,initial)
        proposal=DefensiveGaussianMixture(props[str(days)]['shifts'],props[str(days)]['probabilities'])
        xs,ws,vs=[],[],[]
        for replicate in range(R):
            seed=next(s['seed'] for s in prior['seed_records'] if s['days']==days and s['method']=='random_mixture' and s['replicate']==replicate)
            z,w,_=proposal.draw(max(prior['arguments']['powers']),seed)
            frame=raw[(raw.days==days)&(raw.method=='random_mixture')&(raw.replicate==replicate)].sort_values('point')
            np.testing.assert_allclose(w,frame.weight.to_numpy(),rtol=1e-12,atol=1e-14)
            xs.append(sampler.mean+z[:N]@sampler.root.T);ws.append(w[:N]);vs.append(frame.vix.to_numpy()[:N])
        states=np.concatenate(xs);weight=np.concatenate(ws);base=np.concatenate(vs)
        window=engine.gaussian_window(expiry,time_order=8);norm64=dated.normalizer(window.times,initial,parameters,64)
        gh48=expected_nodes(engine,states,window,48,f'{days}d_GH48')
        rules={'GH32_M48_time8':base,'GH48_M48_time8':values_from_nodes(engine,window,gh48,window.normalizer),'GH48_M64_time8':values_from_nodes(engine,window,gh48,norm64)}
        exact=10000*curve.integral(expiry,expiry+window.delta)/window.delta
        for rule,v in rules.items():
            per=(weight*v*v).reshape(R,N).mean(axis=1)
            change=(weight*(v*v-base*base)).reshape(R,N).mean(axis=1)
            mean,se=replication_statistics(per);dm,ds=replication_statistics(change)
            rows.append(dict(days=days,rule=rule,points_per_replication=N,estimate=float(mean),target=exact,residual=float(mean-exact),replication_se=float(se),change_from_base=float(dm),paired_change_se=float(ds)))
            for replicate in range(R):rep_rows.append(dict(days=days,rule=rule,replicate=replicate,estimate=float(per[replicate]),change=float(change[replicate])))
        stored[f'days{days}_GH48']=gh48;stored[f'days{days}_states']=states;stored[f'days{days}_weight']=weight
        if days==27:
            for test,points,time_order,cond_order in [('time16',16,16,48),('conditional64',4,8,64)]:
                indices=np.arange(R*N).reshape(R,N)[:,:points].ravel()
                selected_states=states[indices];window2=engine.gaussian_window(expiry,time_order=time_order);m64=dated.normalizer(window2.times,initial,parameters,64)
                e=expected_nodes(engine,selected_states,window2,cond_order,f'27d_{test}')
                v=values_from_nodes(engine,window2,e,m64);reference=rules['GH48_M64_time8'][indices]
                delta=(weight[indices]*(v*v-reference*reference)).reshape(R,points).mean(axis=1);dm,ds=replication_statistics(delta)
                checks.append(dict(days=days,check=test,points_per_replication=points,weighted_moment_change=float(dm),paired_se=float(ds),max_pointwise_vix_change=float(np.max(np.abs(v-reference)))))
                stored[f'27_{test}_expected']=e
        pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False)
        print(json.dumps(dict(completed_days=days,elapsed=time.perf_counter()-start)),flush=True)
    check_sources(sources)
    report=dict(status='matched original-state quadrature attribution, no rescaling or refit',scope=scope,architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),curve_levels=prior['curve_levels'],
        target_sha256=targets.artifact_sha256,implementation_sha256=source_hashes(sources,ROOT),checks=checks,outer_source=str(old),
        limits=['27d uses all 8x128 original random-mixture states; others use first 8 states per replication, a limited diagnostic subset.',
            '27d time16 uses first 16 states per replication; conditional64 uses first 4, not a complete surface or population check.',
            'Paired changes isolate numerical rules on fixed states but cannot validate unseen outer tails.',
            'All prior prices, scope decisions and warning records remain unchanged.'],wall_seconds=time.perf_counter()-start)
    (args.output/'run.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');pd.DataFrame(rep_rows).to_csv(args.output/'replications.csv',index=False);np.savez_compressed(args.output/'expected_nodes.npz',**stored)
    (args.output/'status.json').write_text('{"status":"complete"}\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():parser.error('Output exists')
    run(args)
