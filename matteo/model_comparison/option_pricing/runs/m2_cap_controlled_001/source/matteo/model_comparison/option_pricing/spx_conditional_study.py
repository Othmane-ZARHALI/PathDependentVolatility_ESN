"""Target-Q partial conditional Gaussian SPX prototype and frozen-panel check."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .fine_grid import CachedNormalizerPricer
from .spx_hedging_study import refinement
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.conditional_spx import simulate_conditional_spx,independent_equity_fraction
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import dated_grid


FOLDER=Path(__file__).resolve().parent
PARENT=FOLDER/'runs/m2_spx_central_panel_001'


def run(args):
    args.output.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    paths=[Path(__file__).resolve(),FOLDER/'active_scope.json',PARENT/'run.json',PARENT/'normalizers.npz',PARENT/'normalizer_prefix_check.json',PARENT/'contracts.csv',
        FOLDER/'runs/m2_vix_curve_fit_001/run.json',*[FOLDER/(name+'.py') for name in ['m2_spx','fine_grid','spx_hedging_study']],*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(paths);save_sources(args.output,captured,ROOT)
    prior=json.loads((PARENT/'run.json').read_text());targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    architecture,parameters,premium,hashes=frozen_stage3();p=CachedNormalizerPricer(architecture,premium);initial,_=p.initial_state([0.,0.,0.])
    fit=json.loads((FOLDER/'runs/m2_vix_curve_fit_001/run.json').read_text())
    if prior['target_sha256']!=targets.artifact_sha256 or fit['target_sha256']!=targets.artifact_sha256 or prior['initial_state']!=initial.tolist():raise ValueError('Cache target/state mismatch')
    for key,value in dict(architecture=architecture,parameters=parameters,premium=premium).items():
        if prior[key]!=json.loads(json.dumps(asdict(value))):raise ValueError('Frozen model mismatch')
    prefix=str((M2/'src/esn_eq').relative_to(ROOT))+'/'
    for name,digest in prior['implementation_sha256'].items():
        if name.startswith(prefix) and hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Prior normalizer source changed')
    provenance=json.loads((PARENT/'normalizer_prefix_check.json').read_text())
    if hashlib.sha256((PARENT/'normalizers.npz').read_bytes()).hexdigest()!=provenance['new_sha256'] or prior['normalizer_order']!=48:raise ValueError('Cache hash/order mismatch')
    curve=NuisanceCurveFamily.from_targets(targets).build(fit['selected_levels'])
    if list(curve.levels)!=prior['curve_levels']:raise ValueError('Frozen curve mismatch')
    panel=pd.read_csv(PARENT/'contracts.csv');panel.to_csv(args.output/'contracts.csv',index=False)
    scope=json.loads((FOLDER/'active_scope.json').read_text())
    if scope['spx_research_priority']['criterion']['cutoff']!=3:raise ValueError('Scope changed')
    events=np.r_[panel.T_years.unique(),7/365,21/365,[t for t in curve.ends if t<panel.T_years.max()]]
    times=dated_grid(events,1/(252*64))
    for _ in range(3):times=refinement(times)
    with np.load(PARENT/'normalizers.npz') as old:
        np.testing.assert_array_equal(times[:-1],old['times']);values=old['values'].copy()
    values.setflags(write=False);p.normalizers[(times[:-1].tobytes(),initial.tobytes(),parameters,48)]=values
    fraction=independent_equity_fraction(p.layout)
    design=dict(arguments=vars(args)|dict(output=str(args.output)),architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),
        initial_state=initial.tolist(),curve_levels=list(curve.levels),target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,
        implementation_sha256=source_hashes(captured,ROOT),steps_per_trading_day=512,normalizer_order=48,independent_equity_variance_fraction=fraction,
        estimator='Uncontrolled direct-Q raw payoff vs positive partial conditional payoff on coupled retained OU-driving Gaussian shocks.',
        conditioning='Full retained driving shocks determine the OU state path. Independent primitive equity direction integrated analytically. Not a three-projection Markov or full-state pseudoinverse approximation.',
        limitations=['Same intended Gaussian transition law; numerical covariance reconstruction versus legacy sampler is audited, not bitwise same states or a formal payoff-error bound.',
            'No importance sampling: old full-path likelihoods cannot be reused after marginalizing equity noise.',
            'Comparison is against raw target-Q Monte Carlo, not the controlled/importance-sampled reference pricer.',
            'Only 0.64% independent equity variance in this frozen model; positivity does not guarantee a major variance gain or rare-tail accuracy.',
            '24 representatives of707 priority quotes; all original quote/VIX scopes remain; no fit, timestep or normalizer refinement.'])
    (args.output/'design.json').write_text(json.dumps(design,indent=2,allow_nan=False)+'\n')
    rows=[];paired=[];archive={};costs=[];seed_records=[];audits=None
    for trial in range(args.trials):
        pieces={};seeds=[]
        for batch,start in enumerate(range(0,args.paths,args.batch_size)):
            n=min(args.batch_size,args.paths-start);seed=int(np.random.SeedSequence([args.seed,trial,batch]).generate_state(1)[0]);seeds.append(seed)
            (args.output/'status.json').write_text(json.dumps(dict(status='direct-Q raw/conditional checking',trial=trial,paths=start))+'\n')
            tick=time.perf_counter();result=simulate_conditional_spx(p,times,initial,parameters,curve,panel.T_years.to_numpy(),panel.K_over_F.to_numpy(),
                panel.option_type.to_numpy(),n,seed,order=48)
            costs.append(dict(trial=trial,batch=batch,paths=n,joint_seconds=time.perf_counter()-tick))
            if audits is None:audits=result.transition_audit
            data=dict(raw=result.raw_payoffs,conditional=result.conditional_payoffs,raw_stock=result.raw_stock,conditional_stock=result.conditional_stock,
                log_mean=result.log_means,log_variance=result.log_variances)
            for key,value in data.items():pieces.setdefault(key,[]).append(value)
            print(f'trial={trial} paths={start+n}/{args.paths}; independent equity variance={fraction:.6%}',flush=True)
        data={key:np.concatenate(value) for key,value in pieces.items()};seed_records.append(dict(trial=trial,seeds=seeds))
        for key,value in data.items():archive[f't{trial}_{key}']=value
        for j,q in enumerate(panel.itertuples()):
            for method in ['raw','conditional']:
                x=data[method][:,j]*q.forward;center=x-x.mean();squares=center*center
                rows.append(dict(trial=trial,contract=j,method=method,expiry_days=q.expiry_days,strike=q.strike,option_type=q.option_type,band=q.band,z=q.z_mid,
                    estimate=float(x.mean()),se=float(x.std(ddof=1)/np.sqrt(len(x))),variance=float(x.var(ddof=1)),scale=q.price_scale_U,
                    raw_positive_hits=int((data['raw'][:,j]>0).sum()),conditional_positive_paths=int((data['conditional'][:,j]>0).sum()),
                    largest_residual_share=float(squares.max()/squares.sum()) if squares.sum() else 0.))
            difference=(data['conditional'][:,j]-data['raw'][:,j])*q.forward
            paired.append(dict(trial=trial,contract=j,expiry_days=q.expiry_days,band=q.band,difference=float(difference.mean()),se=float(difference.std(ddof=1)/np.sqrt(len(difference)))))
    check_sources(captured)
    for name,data in [('estimates',rows),('paired_differences',paired),('costs',costs)]:pd.DataFrame(data).to_csv(args.output/(name+'.csv'),index=False)
    np.savez_compressed(args.output/'draws.npz',**archive)
    (args.output/'transition_audit.json').write_text(json.dumps(list(audits),indent=2,allow_nan=False)+'\n')
    record=design|dict(status='completed partial conditional SPX investigation; no pricing convergence clearance',seed_records=seed_records,grid_intervals=len(times)-1,
        wall_seconds=time.perf_counter()-begin,uncertainty_unit='iid direct-Q outer paths; raw and conditional paired, strikes/maturities correlated')
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');(args.output/'status.json').write_text('{"status":"completed"}\n')
    print(f'Completed partial conditional study in {record["wall_seconds"]:.1f}s',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--paths',type=int,default=256)
    parser.add_argument('--batch-size',type=int,default=128);parser.add_argument('--trials',type=int,default=2);parser.add_argument('--seed',type=int,default=20261123)
    args=parser.parse_args()
    if args.output.exists() or min(args.paths,args.batch_size)<8 or args.trials!=2 or args.paths%args.batch_size==1:parser.error('Fresh output, valid budget and two trials required')
    run(args)
