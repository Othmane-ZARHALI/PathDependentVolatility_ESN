"""Fresh fixed +3 equity-path tilt vs original sampler, five unchanged controls."""
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
from .spx_hedging_study import refinement,controls
from .spx_path_importance import weighted_controls
from .spx_path_score_controls import ridge_fit
from .spx_sampling_design import call_proposal
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import dated_grid
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.variance_reduction import deterministic_stock_controls
from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths


FOLDER=Path(__file__).resolve().parent
PARENT=FOLDER/'runs/m2_spx_central_panel_001'
DESIGN=FOLDER/'runs/spx_sampling_design_001'
PROPOSALS=('original','call3')
METHODS=('raw','ridge','ridge_A','ridge_B','old128')


def fit_controls(data,old):
    y,b=data['y'],data['b'];n=len(y);fits={}
    for name,index in [('ridge',np.arange(n)),('ridge_A',np.arange(n//2)),('ridge_B',np.arange(n//2,n))]:
        fits[name]=np.stack([ridge_fit(y[index,j],b[index,j]) for j in range(y.shape[1])])
    fits['old128']=np.asarray(old);return fits


def adjusted(data,fits,forwards):
    values=[data['y']]
    for name in METHODS[1:]:values.append(data['y']-np.einsum('njk,jk->nj',data['b'],np.asarray(fits[name])))
    return np.stack(values,axis=-1)*forwards[None,:,None]


def collect(p,parameters,curve,initial,quotes,times,basis,args,day,role,trial):
    count=args.pilot if role==0 else args.paths;pieces={name:{} for name in PROPOSALS};costs=[];seeds=[]
    maturities=quotes.T_years.to_numpy();strikes=quotes.K_over_F.to_numpy();kinds=quotes.option_type.to_numpy()
    for batch,start in enumerate(range(0,count,args.batch_size)):
        n=min(args.batch_size,count-start);seed=int(np.random.SeedSequence([args.seed,int(day),role,trial,batch]).generate_state(1)[0]);seeds.append(seed)
        order=PROPOSALS if batch%2==0 else PROPOSALS[::-1]
        for label in order:
            tick=time.perf_counter();mixture=call_proposal(2. if label=='original' else 3.)
            sample=prepare_volatility_paths(p,times,initial,n,seed,basis,mixture)
            cache,w=sample.cache,sample.weights
            r=deterministic_stock_controls(cache,curve,np.unique(maturities));reference=np.stack([r[t] for t in maturities],axis=1)
            h=spx_hedging_controls(p,cache,parameters,curve,maturities,strikes,kinds,order=48)
            stock=np.stack([h.terminal_ratios[t] for t in maturities],axis=1)
            pairs=controls(dict(stock=stock,hedges=h.hedges),reference,quotes,curve)
            weighted=[weighted_controls(y,b,w) for y,b in pairs]
            data=dict(y=np.stack([a for a,b in weighted],axis=1),b=np.stack([b for a,b in weighted],axis=1),
                w=w,z=sample.score,components=sample.components,monitors=np.c_[w,w*stock[:,0],w*reference[:,0]])
            seconds=time.perf_counter()-tick
            for key,value in data.items():pieces[label].setdefault(key,[]).append(value)
            costs.append(dict(expiry_days=int(day),role=role,trial=trial,batch=batch,proposal=label,paths=n,seconds=seconds,first=bool(label==order[0])))
            del sample,cache,h,stock,reference
        print(f'expiry={day} role={role} trial={trial} paired paths={start+n}/{count}',flush=True)
    return {label:{key:np.concatenate(v) for key,v in data.items()} for label,data in pieces.items()},costs,seeds


def summarize(values,data,quotes,day,trial,label):
    rows=[]
    for j,row in enumerate(quotes.itertuples()):
        for k,method in enumerate(METHODS):
            x=values[:,j,k];squares=(x-x.mean())**2
            rows.append(dict(expiry_days=int(day),trial=trial,proposal=label,contract=j,method=method,band=row.band,
                strike=row.strike,option_type=row.option_type,z=row.z_mid,estimate=float(x.mean()),se=float(x.std(ddof=1)/np.sqrt(len(x))),
                variance=float(x.var(ddof=1)),scale=row.price_scale_U,nonzero_raw=int((data['y'][:,j]>0).sum()),
                largest_residual_share=float(squares.max()/squares.sum()) if squares.sum() else 0.))
    return rows


def run(args):
    args.output.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    files=[Path(__file__).resolve(),FOLDER/'spx_sampling_design.py',FOLDER/'active_scope.json',PARENT/'run.json',PARENT/'normalizers.npz',PARENT/'contracts.csv',
        DESIGN/'design.json',DESIGN/'pilot_components.csv',DESIGN/'pilot_reweighting.csv',FOLDER/'runs/m2_vix_curve_fit_001/run.json']
    files += [FOLDER/(name+'.py') for name in ['m2_spx','fine_grid','spx_hedging_study','spx_path_importance','spx_path_score_controls','spx_volatility_sampling']]
    files += sorted((M2/'src/esn_eq').glob('*.py'))
    files += [folder/name for folder in sorted(PARENT.glob('expiry_*')) for name in ['pilot_draws.npz','frozen.json','contracts.csv']]
    captured=capture_sources(files);save_sources(args.output,captured,ROOT)
    prior=json.loads((PARENT/'run.json').read_text());selection=json.loads((DESIGN/'design.json').read_text())
    if selection['selected_positive_equity_shift']!=3:raise ValueError('Design changed')
    if json.loads((FOLDER/'active_scope.json').read_text())['spx_research_priority']['criterion']['cutoff']!=3:raise ValueError('Scope changed')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');architecture,parameters,premium,hashes=frozen_stage3()
    p=CachedNormalizerPricer(architecture,premium);initial,_=p.initial_state([0.,0.,0.]);fit=json.loads((FOLDER/'runs/m2_vix_curve_fit_001/run.json').read_text())
    if prior['target_sha256']!=targets.artifact_sha256 or fit['target_sha256']!=targets.artifact_sha256:raise ValueError('Target cache mismatch')
    if prior['initial_state']!=initial.tolist() or prior['normalizer_order']!=48:raise ValueError('Cache state/order mismatch')
    for key,value in dict(architecture=architecture,parameters=parameters,premium=premium).items():
        expected=json.loads(json.dumps(asdict(value)))
        if prior[key]!=expected or fit[key]!=expected:raise ValueError('Cached model mismatch')
    for path in sorted((M2/'src/esn_eq').glob('*.py')):
        if hashlib.sha256(path.read_bytes()).hexdigest()!=prior['implementation_sha256'][str(path.relative_to(ROOT))]:raise ValueError('Cache package source mismatch')
    curve=NuisanceCurveFamily.from_targets(targets).build(fit['selected_levels'])
    if list(curve.levels)!=prior['curve_levels']:raise ValueError('Curve mismatch')
    panel=pd.read_csv(PARENT/'contracts.csv');panel.to_csv(args.output/'contracts.csv',index=False)
    events=np.r_[panel.T_years.unique(),7/365,21/365,[t for t in curve.ends if t<panel.T_years.max()]]
    grid=dated_grid(events,1/(252*64))
    for _ in range(3):grid=refinement(grid)
    with np.load(PARENT/'normalizers.npz') as old:
        np.testing.assert_array_equal(grid[:-1],old['times']);normalizers=old['values'].copy()
    if not np.isfinite(normalizers).all() or np.any(normalizers<=0):raise ValueError('Invalid cached normalizers')
    normalizers.setflags(write=False)
    setup=time.perf_counter()-begin
    design=dict(arguments=vars(args)|dict(output=str(args.output)),parent=str(PARENT.relative_to(ROOT)),design_source=str(DESIGN.relative_to(ROOT)),
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),curve_levels=list(curve.levels),
        target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,implementation_sha256=source_hashes(captured,ROOT),
        primary_comparison='Original +2 vs call +3 equity shift; own 256-path ridge5 pilots at equal budget; no new controls.',
        fixed_coefficient_comparison='Identical old128 ridge5 beta under both proposals isolates sampling law; original new256 vs old128 isolates learning.',
        normalizer_order=48,steps_per_trading_day=512,cache_validation='Exact full master time grid, model/state/order/target/curve and all current package sources match previous normalizers; cached bytes frozen.',
        limitations=['24 representatives; full cutoff3/all VIX scope preserved; no model/state/curve fit.',
            'Old pilots used for design; no old checking archive read or used to tune shift. Fresh checking remains unseen until fits frozen.',
            'Hypothesis targets call coverage; aggregate pilot proxies do not establish a robust gain.',
            'Same seeds pair proposals via identical component labels and underlying normals; paths/payoffs differ and cannot be treated as identical.',
            'Five zero-mean controls unchanged; mass compensation retained; no self-normalization or clipping.',
            'GH48 cached approximation and single512 stock discretization remain open bias gates.'])
    (args.output/'design.json').write_text(json.dumps(design,indent=2,allow_nan=False)+'\n')
    all_rows=[];all_costs=[];records=[]
    for day,quotes in panel.groupby('expiry_days',sort=True):
        tick=time.perf_counter();quotes=quotes.reset_index(drop=True);out=args.output/f'expiry_{int(day)}';out.mkdir();quotes.to_csv(out/'contracts.csv',index=False)
        times=grid[:np.searchsorted(grid,float(quotes.T_years.iloc[0]))+1];values=normalizers[:len(times)-1];values.setflags(write=False)
        p.normalizers[(times[:-1].tobytes(),initial.tobytes(),parameters,48)]=values
        basis=volatility_path_basis(p,times,7/365)
        (args.output/'status.json').write_text(json.dumps(dict(status='pilot',expiry_days=int(day)))+'\n')
        pilot,costs,seeds=collect(p,parameters,curve,initial,quotes,times,basis,args,day,0,0);all_costs+=costs
        old=json.loads((PARENT/f'expiry_{int(day)}/frozen.json').read_text())['coefficients']['base']
        fits={label:fit_controls(data,old) for label,data in pilot.items()}
        frozen=dict(coefficients={label:{name:value.tolist() for name,value in f.items()} for label,f in fits.items()},pilot_seeds=seeds,checking_not_used=True)
        (out/'frozen.json').write_text(json.dumps(frozen,indent=2,allow_nan=False)+'\n')
        np.savez_compressed(out/'pilot_draws.npz',**{label+'_'+key:value for label,data in pilot.items() for key,value in data.items()})
        seeds_record=[dict(role=0,trial=0,seeds=seeds)];archive={}
        for trial in range(args.trials):
            (args.output/'status.json').write_text(json.dumps(dict(status='checking; pilots frozen',expiry_days=int(day),trial=trial))+'\n')
            data,costs,seeds=collect(p,parameters,curve,initial,quotes,times,basis,args,day,1,trial);all_costs+=costs;seeds_record.append(dict(role=1,trial=trial,seeds=seeds))
            for label in PROPOSALS:
                x=adjusted(data[label],fits[label],quotes.forward.to_numpy());all_rows+=summarize(x,data[label],quotes,day,trial,label)
                for key,value in data[label].items():archive[f't{trial}_{label}_{key}']=value
                archive[f't{trial}_{label}_adjusted']=x
        np.savez_compressed(out/'checking_draws.npz',**archive)
        record=dict(expiry_days=int(day),seed_records=seeds_record,grid_intervals=len(times)-1,wall_seconds=time.perf_counter()-tick)
        (out/'run.json').write_text(json.dumps(record,indent=2)+'\n');records.append(record)
        print(f'Completed expiry {day}: {record["wall_seconds"]:.1f}s',flush=True)
    check_sources(captured)
    pd.DataFrame(all_rows).to_csv(args.output/'estimates.csv',index=False);pd.DataFrame(all_costs).to_csv(args.output/'costs.csv',index=False)
    record=design|dict(status='completed fixed call-tilt sampling diagnostic',setup_seconds=setup,wall_seconds=time.perf_counter()-begin,expiry_records=records)
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');(args.output/'status.json').write_text('{"status":"completed"}\n')
    print(f'Completed call-tilt study in {record["wall_seconds"]:.1f}s',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--pilot',type=int,default=256);parser.add_argument('--paths',type=int,default=256);parser.add_argument('--batch-size',type=int,default=128)
    parser.add_argument('--trials',type=int,default=2);parser.add_argument('--seed',type=int,default=20261122);args=parser.parse_args()
    if args.output.exists() or min(args.pilot,args.paths,args.batch_size)<8 or args.pilot%2 or args.trials!=2:parser.error('Require fresh output, even pilot, valid budgets, two trials')
    if any(n%args.batch_size==1 for n in [args.pilot,args.paths]):parser.error('No one-path final batches')
    run(args)
