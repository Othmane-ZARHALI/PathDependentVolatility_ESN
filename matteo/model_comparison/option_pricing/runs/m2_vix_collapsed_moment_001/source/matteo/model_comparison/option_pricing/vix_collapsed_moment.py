"""Separate Gaussian-law, normalizer and time-rule checks for E[VIX^2]."""
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
from scipy.special import ndtri
from scipy.stats import qmc
from .m2_spx import ROOT,M2,frozen_stage3
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.outer_sampling import GaussianOuterSampler,replication_statistics
from esn_eq.pricing import DatedM2Pricer
from esn_eq.states import symmetric_square_root
from esn_eq.vix_pricing import ConditionalVIXPricer


def tower_errors(dated,initial,expiry,lag):
    sampler=GaussianOuterSampler(dated,expiry,initial)
    decay=np.exp(-dated.layout.rates*lag)
    shift=dated.q_drift*(-np.expm1(-dated.layout.rates*lag))/dated.layout.rates
    mean,cov=dated.gaussian_projection_law(expiry+lag,initial)
    _,innovation=dated.gaussian_projection_law(lag,np.zeros(len(initial)))
    composed_mean=dated.projection@(decay*sampler.mean+shift)
    composed_cov=dated.projection@(decay[:,None]*sampler.covariance*decay[None,:])@dated.projection.T+innovation
    return float(np.max(np.abs(composed_mean-mean))),float(np.max(np.abs(composed_cov-cov))),float(np.max(np.abs(sampler.root@sampler.root.T-sampler.covariance)))


def collapsed_integrands(engine,window,normals,denominators):
    """Only for mean VIX squared, not VIX draws, futures or option payoffs.

    A common 3D normal is reused across future dates. This need not have the
    correct path law: linearity of expectation requires only correct marginals.
    Averaging future variance directly bypasses conditional VIX integration.
    """
    z=np.asarray(normals,float)
    den=np.asarray(denominators,float)
    if z.ndim!=2 or z.shape[1]!=3 or not np.isfinite(z).all() or den.ndim!=2 or den.shape[1]!=len(window.times) or np.any(den<=0) or not np.isfinite(den).all():
        raise ValueError('Invalid Gaussian normals or node normalizers')
    if window.owner is not engine or window.configuration!=engine._configuration():
        raise ValueError('Stale or foreign window')
    out=np.zeros((len(z),len(den)))
    for node,(t,w,xi) in enumerate(zip(window.times,window.weights,window.xi)):
        mean,cov=engine.dated.gaussian_projection_law(float(t),engine.initial_state)
        g=engine.dated.multiplier(mean+z@symmetric_square_root(cov).T,engine.parameters)
        out += (10000*w/window.delta)*(engine.variance_floor+(xi-engine.variance_floor)*g[:,None]/den[:,node])
    return out


def run(args):
    start=time.perf_counter();folder=Path(__file__).parent
    prior_path=folder/'runs/m2_vix_stratified_sampling_001/run.json';scope_path=folder/'active_scope.json'
    sources=capture_sources([Path(__file__),folder/'m2_spx.py',scope_path,prior_path,*sorted((M2/'src/esn_eq').glob('*.py'))])
    prior=json.loads(prior_path.read_text());scope=json.loads(scope_path.read_text())
    if scope['vix_fit_expiry_days']!=[27,55,118,153] or scope['vix_deferred_expiry_days']:
        raise ValueError('Require restored full scope')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    architecture,parameters,premium,hashes=frozen_stage3();dated=DatedM2Pricer(architecture,premium)
    initial,_=dated.initial_state([0.,0.,0.])
    expected=json.loads(json.dumps(dict(architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),target_sha256=targets.artifact_sha256)))
    if any(prior[k]!=v for k,v in expected.items()):raise ValueError('Frozen model or target changed')
    curve=NuisanceCurveFamily.from_targets(targets).build(prior['curve_levels'])
    engine=ConditionalVIXPricer(dated,parameters,curve,initial,normalizer_order=48)
    args.output.mkdir(parents=True,exist_ok=False);save_sources(args.output,sources,ROOT)
    (args.output/'status.json').write_text('{"status":"running"}\n')
    summaries,rep_rows,node_rows,law_rows,ref_rows,seed_rows,time_changes=[],[],[],[],[],[],[]
    stored={}
    for days in scope['vix_fit_expiry_days']:
        expiry=float(targets.select('VIX',expiry_days=[days]).T_years.iloc[0]);all_estimates={}
        exact=10000*curve.integral(expiry,expiry+30/365)/(30/365)
        for time_order in (8,16):
            window=engine.gaussian_window(expiry,time_order=time_order)
            normal64=dated.normalizer(window.times,initial,parameters,64)
            denominators=np.array([window.normalizer,normal64])
            stored[f'days{days}_time{time_order}_times']=window.times
            stored[f'days{days}_time{time_order}_weights']=window.weights
            stored[f'days{days}_time{time_order}_normalizers']=denominators
            for node,t in enumerate(window.times):
                me,ce,re=tower_errors(dated,initial,expiry,float(t-expiry))
                law_rows.append(dict(days=days,time_order=time_order,node=node,time=float(t),mean_error=me,covariance_error=ce,root_error=re))
                node_rows.append(dict(days=days,time_order=time_order,node=node,time=float(t),normalizer48=float(denominators[0,node]),normalizer64=float(normal64[node])))
            if time_order==8:
                for node in np.unique(np.linspace(0,len(window.times)-1,5).astype(int)):
                    for order in (48,64,96,128):
                        value=denominators[0 if order==48 else 1,node] if order in (48,64) else dated.normalizer([window.times[node]],initial,parameters,order)[0]
                        ref_rows.append(dict(days=days,node=int(node),time=float(window.times[node]),order=order,value=float(value)))
            estimates={power:[] for power in args.powers}
            for replicate in range(args.replications):
                seed=int(np.random.SeedSequence([args.seed,days,replicate]).generate_state(1)[0])
                u=qmc.Sobol(3,scramble=True,bits=52,rng=seed).random_base2(max(args.powers))
                if np.any((u<=0)|(u>=1)):raise FloatingPointError('Sobol endpoint')
                integral=collapsed_integrands(engine,window,ndtri(u),denominators)
                for power in args.powers:
                    mean=integral[:2**power].mean(axis=0);estimates[power].append(mean)
                    for j,order in enumerate((48,64)):
                        rep_rows.append(dict(days=days,time_order=time_order,normalizer_order=order,power=power,replicate=replicate,estimate=float(mean[j]),target=exact))
                if time_order==8:seed_rows.append(dict(days=days,replicate=replicate,seed=seed))
            all_estimates[time_order]=estimates
            for power,matrix in estimates.items():
                mean,se=replication_statistics(matrix)
                for j,order in enumerate((48,64)):
                    summaries.append(dict(days=days,time_order=time_order,normalizer_order=order,power=power,estimate=float(mean[j]),target=exact,residual=float(mean[j]-exact),replication_se=float(se[j])))
            print(json.dumps(dict(days=days,time_order=time_order,nodes=len(window.times),elapsed=time.perf_counter()-start)),flush=True)
            pd.DataFrame(summaries).to_csv(args.output/'summary.csv',index=False)
        for power in args.powers:
            diff=np.array(all_estimates[16][power])-np.array(all_estimates[8][power]);mean,se=replication_statistics(diff)
            for j,order in enumerate((48,64)):time_changes.append(dict(days=days,power=power,normalizer_order=order,time16_minus8=float(mean[j]),paired_se=float(se[j])))
    check_sources(sources)
    record=dict(status='collapsed unconditional Gaussian moment diagnostic; no price correction or fit',arguments=vars(args)|dict(output=str(args.output)),scope=scope,
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),curve_levels=prior['curve_levels'],
        target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,implementation_sha256=source_hashes(sources,ROOT),seed_records=seed_rows,
        estimand='10000/delta integral E[V_t] dt, using direct 3D future Gaussian marginals; compare with conditional VIX squared identity',
        randomness='same 3D scrambled normals across future times and time/normalizer rules; complete independent scrambles define SE',
        limitations=['Not simulated VIX values; cannot take square roots to price VIX futures/options.','Gaussian algebra errors are floating-point checks, not formal proofs.',
            'Direct-marginal uncertainty excludes remaining time and normalizer approximation biases.','Sparse GH96/128 node references do not bound all-node errors.',
            'Diagnoses normalization/tower identity, not joint-fit feasibility or physical compatibility.'],versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),wall_seconds=time.perf_counter()-start)
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    for name,rows in [('replications',rep_rows),('normalizer_nodes',node_rows),('gaussian_law_checks',law_rows),('normalizer_reference',ref_rows),('time_changes',time_changes)]:pd.DataFrame(rows).to_csv(args.output/(name+'.csv'),index=False)
    np.savez_compressed(args.output/'nodes.npz',**stored)
    (args.output/'status.json').write_text('{"status":"complete"}\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--powers',type=int,nargs='+',default=[12,14,16]);parser.add_argument('--replications',type=int,default=8);parser.add_argument('--seed',type=int,default=20261030)
    args=parser.parse_args()
    if args.output.exists() or args.replications<4 or min(args.powers)<4:parser.error('Require new output, >=4 replications and powers>=4')
    run(args)
