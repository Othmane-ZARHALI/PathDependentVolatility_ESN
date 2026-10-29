"""Frozen market-selected cutoff-3 panel: independent direct512 sampling checks."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time
from types import SimpleNamespace
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .fine_grid import CachedNormalizerPricer
from .spx_hedging_study import refinement
from .spx_innovation_study import collect,datasets,train,adjusted,METHODS
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import dated_grid
from esn_eq.volatility_path_importance import volatility_path_basis


FOLDER=Path(__file__).resolve().parent
MEMBERSHIP=FOLDER/'runs/spx_scaled_moneyness_003/training_membership.csv'


def choose_panel(quotes):
    if set(quotes.split)!={'train'} or set(quotes.instrument)!={'SPX'}:
        raise ValueError('Panel source must be SPX training only')
    rows=[]
    for _,group in quotes.groupby('expiry',sort=True):
        for kind in ('put','call'):
            for band,center in enumerate((.5,1.5,2.5)):
                eligible=group[(group.option_type==kind)&(group.z_mid>=band)&(group.z_mid<band+1 if band<2 else group.z_mid<=3)].copy()
                if eligible.empty:raise ValueError('Missing maturity/side/band representative')
                eligible['selection_distance']=abs(eligible.z_mid-center)
                row=eligible.sort_values(['selection_distance','strike','option_type']).iloc[0].copy()
                row['band']=band;row['target_z']=center;rows.append(row)
    panel=pd.DataFrame(rows).sort_values(['T_years','strike']).reset_index(drop=True)
    if panel.duplicated(['expiry','strike','option_type']).any():raise ValueError('Panel duplicates')
    if not np.all(panel.z_mid<=3):raise ValueError('Panel outside priority band')
    return panel


def cache_normalizer(p,times,initial,parameters,output,chunk=2048):
    """Compute the unchanged pointwise GH48 rule in bounded progress chunks."""
    pieces=[]
    for start in range(0,len(times),chunk):
        pieces.append(p.normalizer(times[start:start+chunk],initial,parameters,48))
        print(f'normalizer {min(start+chunk,len(times))}/{len(times)}',flush=True)
    values=np.concatenate(pieces);values.setflags(write=False)
    p.normalizers[(times.tobytes(),initial.tobytes(),parameters,48)]=values
    np.savez_compressed(output/'normalizers.npz',times=times,values=values)
    return values


def record_trial(data,coefficients,quotes,trial):
    y,b,v=datasets(data)['price512'];x=adjusted(y,b,v,coefficients,quotes.forward.to_numpy());rows=[];paired=[]
    for j,row in enumerate(quotes.itertuples()):
        for k,method in enumerate(METHODS):
            values=x[:,j,k];centered=values-values.mean();squares=centered**2
            rows.append(dict(trial=trial,contract=j,method=method,expiry_days=row.expiry_days,strike=row.strike,
                option_type=row.option_type,band=row.band,z=row.z_mid,target_z=row.target_z,
                estimate=float(values.mean()),se=float(values.std(ddof=1)/np.sqrt(len(values))),variance=float(values.var(ddof=1)),
                scale=row.price_scale_U,nonzero_raw=int((y[:,j]>0).sum()),
                largest_share_centered_squares=float(squares.max()/squares.sum()) if squares.sum()>0 else 0.))
        delta=x[:,j,3]-x[:,j,2]
        paired.append(dict(trial=trial,contract=j,difference=float(delta.mean()),se=float(delta.std(ddof=1)/np.sqrt(len(delta)))))
    return x,rows,paired


def run(args):
    args.output.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    helpers=['m2_spx','fine_grid','spx_hedging_study','spx_innovation_study','spx_path_importance','spx_volatility_sampling','spx_path_score_controls']
    fit_path=FOLDER/'runs/m2_vix_curve_fit_001/run.json'
    captured=capture_sources([Path(__file__).resolve(),fit_path,MEMBERSHIP,FOLDER/'active_scope.json',
        *[FOLDER/(name+'.py') for name in helpers],*sorted((M2/'src/esn_eq').glob('*.py'))])
    save_sources(args.output,captured,ROOT)
    panel=choose_panel(pd.read_csv(MEMBERSHIP));panel.to_csv(args.output/'contracts.csv',index=False)
    scope=json.loads((FOLDER/'active_scope.json').read_text())
    if scope['spx_research_priority']['criterion']['cutoff']!=3:raise ValueError('Priority changed')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    # Confirm membership rows still correspond exactly to unchanged training inputs.
    current=targets.select('SPX').set_index(['expiry','strike','option_type'])
    for row in panel.itertuples():
        actual=current.loc[(row.expiry,row.strike,row.option_type)]
        np.testing.assert_allclose(actual[['T_years','forward','price_scale_U']].to_numpy(float),[row.T_years,row.forward,row.price_scale_U],rtol=0,atol=1e-12)
    architecture,parameters,premium,hashes=frozen_stage3();p=CachedNormalizerPricer(architecture,premium)
    initial,_=p.initial_state([0.,0.,0.]);fit=json.loads(fit_path.read_text())
    if fit['target_sha256']!=targets.artifact_sha256 or fit['initial_state']!=initial.tolist():raise ValueError('Frozen targets/state changed')
    for key,value in dict(architecture=architecture,parameters=parameters,premium=premium).items():
        if fit[key]!=json.loads(json.dumps(asdict(value))):raise ValueError('Frozen model changed')
    curve=NuisanceCurveFamily.from_targets(targets).build(fit['selected_levels'])
    events=np.r_[panel.T_years.unique(),7/365,21/365,[t for t in curve.ends if t<panel.T_years.max()]]
    master=[dated_grid(events,1/(252*64))]
    for _ in range(3):master.append(refinement(master[-1]))
    design=dict(arguments=vars(args)|dict(output=str(args.output)),contracts=len(panel),maturities=sorted(panel.expiry_days.unique().tolist()),
        selection='Nearest z=.5/1.5/2.5 within each unit band, per put/call side and expiry; training market membership only.',
        primary='ridge8 linear vs ridge5 baseline, penalty 1; legacy OLS5 and both pilot halves retained',
        uncertainty_unit='iid checking paths conditional on independently frozen pilot; paired methods; independent seeds per maturity/trial',
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),
        curve_levels=list(curve.levels),target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,
        implementation_sha256=source_hashes(captured,ROOT),normalizer_order=48,grid_max_steps_per_trading_day=512,
        limitations=['24 representatives, not all 707 priority quotes; no model/curve/state refit or held-out evaluation.',
            'Per-maturity horizon-specific vol7 proposals: short-maturity selection evidence does not establish long-maturity efficiency.',
            'Master event grid includes all panel expiries and curve knots; each maturity uses its exact prefix.',
            'Single512 grid only; no stock timestep or GH refinement clearance.',
            'Small pilot/checks are diagnostics; extrapolated budgets are conditional, not convergence guarantees.',
            'Negative controlled estimates, zero payoff hits, leverage and control-mean warnings are retained without clipping.'])
    (args.output/'design.json').write_text(json.dumps(design,indent=2,allow_nan=False)+'\n')
    (args.output/'status.json').write_text('{"status":"normalizer setup"}\n')
    normalizers=cache_normalizer(p,master[3][:-1],initial,parameters,args.output)
    setup=time.perf_counter()-begin;all_costs=[];all_estimates=[];all_paired=[];records=[]
    for days,quotes in panel.groupby('expiry_days',sort=True):
        tick=time.perf_counter();quotes=quotes.reset_index(drop=True);horizon=float(quotes.T_years.iloc[0]);out=args.output/f'expiry_{int(days)}';out.mkdir()
        quotes.to_csv(out/'contracts.csv',index=False)
        grids=[g[:np.searchsorted(g,horizon)+1] for g in master]
        n=len(grids[3])-1;values=normalizers[:n];values.setflags(write=False)
        p.normalizers[(grids[3][:-1].tobytes(),initial.tobytes(),parameters,48)]=values
        basis=volatility_path_basis(p,grids[3],7/365)
        local=SimpleNamespace(pilot=args.pilot,paths=args.paths,batch_size=args.batch_size,seed=args.seed+int(days),direct_only=True,linear_only=True)
        (args.output/'status.json').write_text(json.dumps(dict(status='pilot',expiry_days=int(days)))+'\n')
        training,seeds,costs=collect(p,parameters,curve,initial,quotes,grids,basis,local,0,0)
        coefficients=train(*datasets(training)['price512'])
        frozen=dict(coefficients={k:v.tolist() for k,v in coefficients.items()},checking_not_used=True,
            args=vars(local),pilot_seeds=seeds)
        (out/'frozen.json').write_text(json.dumps(frozen,indent=2,allow_nan=False)+'\n');np.savez_compressed(out/'pilot_draws.npz',**training)
        seed_records=[dict(role=0,trial=0,seeds=seeds)];archive={}
        for trial in range(args.trials):
            (args.output/'status.json').write_text(json.dumps(dict(status='checking; pilot frozen',expiry_days=int(days),trial=trial))+'\n')
            data,seeds,cc=collect(p,parameters,curve,initial,quotes,grids,basis,local,1,trial);costs+=cc
            seed_records.append(dict(role=1,trial=trial,seeds=seeds));x,rows,paired=record_trial(data,coefficients,quotes,trial)
            for key,value in data.items():archive[f't{trial}_'+key]=value
            archive[f't{trial}_price512_adjusted']=x
            all_estimates+=rows;all_paired +=[r|dict(expiry_days=int(days)) for r in paired]
        np.savez_compressed(out/'checking_draws.npz',**archive)
        all_costs +=[row|dict(expiry_days=int(days)) for row in costs]
        record=dict(expiry_days=int(days),seed_records=seed_records,grid_intervals=n,wall_seconds=time.perf_counter()-tick,
            normalizer_seconds_excluded=True,proposal_horizon=horizon)
        (out/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');records.append(record)
        print(f'Completed expiry {days}: {record["wall_seconds"]:.1f}s',flush=True)
    check_sources(captured)
    for name,data in [('estimates',all_estimates),('paired_differences',all_paired),('costs',all_costs)]:pd.DataFrame(data).to_csv(args.output/(name+'.csv'),index=False)
    record=design|dict(status='completed bounded sampling diagnostic; accuracy gates remain open',setup_seconds=setup,
        wall_seconds=time.perf_counter()-begin,expiry_records=records)
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');(args.output/'status.json').write_text('{"status":"completed"}\n')
    print(f'Completed central panel in {record["wall_seconds"]:.1f}s',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--pilot',type=int,default=128);parser.add_argument('--paths',type=int,default=256)
    parser.add_argument('--batch-size',type=int,default=128);parser.add_argument('--trials',type=int,default=2);parser.add_argument('--seed',type=int,default=20261101)
    args=parser.parse_args()
    if args.output.exists() or min(args.pilot,args.paths,args.batch_size)<8 or args.pilot%2 or args.trials<2:parser.error('Require fresh output, even pilot, valid budgets and at least two trials')
    if any(n%args.batch_size==1 for n in [args.pilot,args.paths]):parser.error('One-path final batch')
    run(args)
