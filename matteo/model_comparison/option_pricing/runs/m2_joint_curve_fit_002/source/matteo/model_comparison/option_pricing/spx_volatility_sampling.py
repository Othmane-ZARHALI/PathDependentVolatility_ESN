"""Design-only screening of OU-excursion path mixtures, then fresh SPX checks."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT, M2, frozen_stage3
from .fine_grid import CachedNormalizerPricer
from .spx_hedging_study import choose_quotes, refinement, controls, beta_fit
from .spx_path_importance import weighted_controls
from .provenance import capture_sources, source_hashes, save_sources, check_sources
from esn_eq.volatility_path_importance import volatility_path_basis, prepare_volatility_paths
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import dated_grid, coarsen_paths
from esn_eq.variance_reduction import deterministic_stock_controls


NAMES = ('equity', 'vol7', 'vol21')


def proposal(method):
    shifts = np.zeros((3 if method == 0 else 9, 4))
    shifts[1,0], shifts[2,0] = 2., -2.
    if method == 0:
        return DefensiveGaussianMixture(shifts,[.5,.25,.25])
    for axis in range(1,4):
        shifts[2*axis+1,axis], shifts[2*axis+2,axis] = 2.,-2.
    return DefensiveGaussianMixture(shifts,[.5,.125,.125]+[1/24]*6)


def fit_variants(y, c, w):
    fits = []
    for indices in (np.arange(len(y)), np.arange(len(y)//2), np.arange(len(y)//2,len(y))):
        fits.append(np.stack([beta_fit(*weighted_controls(y[indices,j],c[indices,j],w[indices])) for j in range(y.shape[1])]))
    return np.stack(fits)


def adjusted(y, c, w, fits, forwards):
    result = [w[:,None]*y]
    for fit in fits:
        result.append(np.stack([weighted_controls(y[:,j],c[:,j],w)[0]-weighted_controls(y[:,j],c[:,j],w)[1]@fit[j] for j in range(y.shape[1])],axis=1))
    return np.stack(result,axis=-1)*forwards[None,:,None]


def collect(p, parameters, curve, initial, quotes, grids, basis, mixture, args, method, role, trial):
    budget = (args.pilot,args.selection,args.paths)[role]
    used = [2] if role < 2 else [0,1,2]
    pieces = {level:[] for level in used}
    scores, weights, components, monitors, seeds, costs = [],[],[],[],[],[]
    maturities = quotes.T_years.to_numpy()
    for batch,start in enumerate(range(0,budget,args.batch_size)):
        n = min(args.batch_size,budget-start)
        seed = int(np.random.SeedSequence([args.seed,method,role,trial,batch]).generate_state(1)[0]); seeds.append(seed)
        tick = time.perf_counter()
        proposal_paths = prepare_volatility_paths(p,grids[2],initial,n,seed,basis,mixture)
        cache,w = proposal_paths.cache, proposal_paths.weights
        state_seconds = time.perf_counter()-tick
        ref = deterministic_stock_controls(cache,curve,np.unique(maturities))
        reference = np.stack([ref[t] for t in maturities],axis=1)
        scores.append(proposal_paths.score); weights.append(w); components.append(proposal_paths.components)
        for level in used:
            tick = time.perf_counter()
            c = cache if level == 2 else coarsen_paths(cache,grids[level])
            h = spx_hedging_controls(p,c,parameters,curve,maturities,quotes.K_over_F.to_numpy(),quotes.option_type.to_numpy(),order=48)
            stock = np.stack([h.terminal_ratios[t] for t in maturities],axis=1)
            pairs = controls(dict(stock=stock,hedges=h.hedges),reference,quotes,curve)
            pieces[level].append((np.stack([a for a,b in pairs],axis=1),np.stack([b for a,b in pairs],axis=1)))
            costs.append(dict(method=NAMES[method],role=role,trial=trial,batch=batch,paths=n,steps=64*2**level,
                              state_seconds=state_seconds if level==2 else 0., stock_control_seconds=time.perf_counter()-tick))
            if level == 2:
                monitors.append(np.c_[w,w*stock[:,0],w*stock[:,-1],w*reference[:,0],w*reference[:,-1]])
        print(f'{NAMES[method]} role={role} trial={trial} paths={start+n}/{budget}',flush=True)
    return ({level:tuple(np.concatenate([part[j] for part in parts]) for j in range(2)) for level,parts in pieces.items()},
            np.concatenate(weights),np.concatenate(scores),np.concatenate(components),np.concatenate(monitors),seeds,costs)


def run(args):
    args.output.mkdir(parents=True,exist_ok=False)
    tick = time.perf_counter()
    folder = Path(__file__).resolve().parent
    fit_path = folder/'runs/m2_vix_curve_fit_001/run.json'
    captured = capture_sources([Path(__file__).resolve(),fit_path,folder/'active_scope.json',
        *[folder/(name+'.py') for name in ('m2_spx','fine_grid','spx_hedging_study','spx_path_importance')],
        *sorted((M2/'src/esn_eq').glob('*.py'))])
    save_sources(args.output,captured,ROOT)
    targets = MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    quotes = choose_quotes(targets)
    architecture,parameters,premium,hashes = frozen_stage3()
    p = CachedNormalizerPricer(architecture,premium); initial,_=p.initial_state([0.,0.,0.])
    fit = json.loads(fit_path.read_text())
    if fit['target_sha256'] != targets.artifact_sha256 or fit['initial_state'] != initial.tolist():
        raise ValueError('Frozen inputs changed')
    for key,value in dict(architecture=architecture,parameters=parameters,premium=premium).items():
        if fit[key] != json.loads(json.dumps(asdict(value))):
            raise ValueError('Frozen model changed')
    curve = NuisanceCurveFamily.from_targets(targets).build(fit['selected_levels'])
    events=np.r_[quotes.T_years.unique(),7/365,21/365,[t for t in curve.ends if t < quotes.T_years.max()]]
    g=dated_grid(events,1/(252*64)); grids=[g,refinement(g),refinement(refinement(g))]
    bases={7:volatility_path_basis(p,grids[2],7/365),21:volatility_path_basis(p,grids[2],21/365)}
    # Warm deterministic quadrature separately: timing comparisons exclude first-call setup.
    for grid in grids:
        p.normalizer(grid[:-1],initial,parameters,48)
    setup_seconds=time.perf_counter()-tick
    (args.output/'status.json').write_text(json.dumps(dict(status='design',arguments=vars(args)|dict(output=str(args.output)))))
    archive, fitted, seed_records, costs, design = {},{},[],[],[]
    forwards=quotes.forward.to_numpy()
    central=[2,3,8,9]  # Two nearest-ATM representatives at each expiry, fixed before data.
    for method,name in enumerate(NAMES):
        basis=bases[21 if method==2 else 7]; mix=proposal(method)
        for role in (0,1):
            data,w,z,comp,monitor,seeds,timing=collect(p,parameters,curve,initial,quotes,grids,basis,mix,args,method,role,0)
            y,c=data[2]
            prefix=f'{name}_r{role}_t0'
            for key,value in dict(y=y,c=c,w=w,z=z,components=comp,monitors=monitor).items():
                archive[prefix+'_'+key]=value
            seed_records.append(dict(method=name,role=role,trial=0,seeds=seeds));costs+=timing
            if role == 0:
                fitted[name]=fit_variants(y,c,w)
            else:
                x=adjusted(y,c,w,fitted[name],forwards)[:,:,1]
                var=x.var(axis=0,ddof=1)
                # Four liquid representatives avoid selecting on unobserved tail payoffs.
                score=float(np.log(np.maximum(var[central]/quotes.price_scale_U.to_numpy()[central]**2,1e-24)).mean())
                design.append(dict(method=name,score=score,central_contracts=central,variance=var.tolist(),positive_hits=(y>0).sum(axis=0).tolist()))
    selected=min(design,key=lambda row:row['score'])['method']
    selection=dict(selected=selected,objective='mean log controlled variance / quote_scale² on fixed near-ATM contracts 2,3,8,9',
                   candidates=design,selection_paths=args.selection,checking_not_used=True)
    (args.output/'selection.json').write_text(json.dumps(selection,indent=2,allow_nan=False)+'\n')
    # Freeze this record before any role-2 draw. All candidates remain in fresh comparisons.
    np.savez_compressed(args.output/'design_draws.npz',**archive);archive={}
    estimates,paired,diagnostics=[],[],[]
    for trial in range(args.trials):
        for method,name in enumerate(NAMES):
            basis=bases[21 if method==2 else 7];mix=proposal(method)
            data,w,z,comp,monitor,seeds,timing=collect(p,parameters,curve,initial,quotes,grids,basis,mix,args,method,2,trial)
            prefix=f'{name}_r2_t{trial}'
            for key,value in dict(w=w,z=z,components=comp,monitors=monitor).items():
                archive[prefix+'_'+key]=value
            seed_records.append(dict(method=name,role=2,trial=trial,seeds=seeds)); costs+=timing
            for j,label in enumerate(('mass','stock29','stock57','reference29','reference57')):
                diagnostics.append(dict(trial=trial,method=name,statistic=label,mean=float(monitor[:,j].mean()),se=float(monitor[:,j].std(ddof=1)/np.sqrt(args.paths))))
            adjusted_levels={}
            for level,(y,c) in data.items():
                x=adjusted(y,c,w,fitted[name],forwards); adjusted_levels[level]=x
                archive[prefix+f'_g{level}_y']=y;archive[prefix+f'_g{level}_c']=c
                archive[prefix+f'_g{level}_adjusted']=x
                for j,row in enumerate(quotes.itertuples()):
                    for k,estimator in enumerate(('raw','pooled','halfA','halfB')):
                        estimates.append(dict(trial=trial,method=name,contract=j,expiry_days=row.expiry_days,strike=row.strike,
                            option_type=row.option_type,steps=64*2**level,estimator=estimator,price=float(x[:,j,k].mean()),
                            se=float(x[:,j,k].std(ddof=1)/np.sqrt(args.paths)),scale=row.price_scale_U,positive_payoffs=int((y[:,j]>0).sum())))
            for low,high in ((0,1),(1,2)):
                diff=adjusted_levels[high]-adjusted_levels[low]
                for j,row in enumerate(quotes.itertuples()):
                    for k,estimator in enumerate(('raw','pooled','halfA','halfB')):
                        paired.append(dict(trial=trial,method=name,contract=j,estimator=estimator,coarse_steps=64*2**low,fine_steps=64*2**high,
                            change=float(diff[:,j,k].mean()),se=float(diff[:,j,k].std(ddof=1)/np.sqrt(args.paths)),
                            variance=float(diff[:,j,k].var(ddof=1)),scale=row.price_scale_U))
    check_sources(captured)
    record=dict(status='completed exploratory design/check diagnostic; no convergence clearance',arguments=vars(args)|dict(output=str(args.output)),
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),curve_levels=list(curve.levels),
        target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,implementation_sha256=source_hashes(captured,ROOT),
        proposals={name:dict(shifts=proposal(i).shifts.tolist(),probabilities=proposal(i).probabilities.tolist(),anchor_days=21 if i==2 else 7) for i,name in enumerate(NAMES)},
        fitted_controls={name:value.tolist() for name,value in fitted.items()},selected=selected,seed_records=seed_records,
        grid_intervals=[len(grid)-1 for grid in grids],normalizer_order=48,setup_seconds=setup_seconds,wall_seconds=time.perf_counter()-tick,
        uncertainty_unit='iid checking paths conditional on frozen proposal selection and pilot fits; shared paths across estimator variants and grids',
        limitations=['Twelve predeclared training contracts only; all original scopes retained.',
            'No refitting of model/curve/state; GH48 fixed, earlier GH48/64 check remains sample-conditional.',
            'Selection objective covers four near-ATM representatives, not tails or the full surface.',
            'All candidate checking results retained; only design-selected proposal is primary.',
            'Independent candidate draws; empirical variance ratios are descriptive.',
            'Pilot halves share checking paths; checking trials share frozen pilot fits.',
            'Zero hits/zero empirical SE never certify tail accuracy; no negative-price clipping.',
            'Cost rows exclude warmed normalizers; state cost common to all coarsenings, not standalone coarse simulation cost.',
            'No MLMC optimal allocation or integration/normalizer population bound.'])
    for filename,rows in (('estimates',estimates),('paired_changes',paired),('diagnostics',diagnostics),('costs',costs)):
        pd.DataFrame(rows).to_csv(args.output/(filename+'.csv'),index=False)
    quotes.to_csv(args.output/'contracts.csv',index=False)
    np.savez_compressed(args.output/'checking_draws.npz',**archive)
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    (args.output/'status.json').write_text('{"status":"completed"}\n')
    print(f'Completed, selected={selected}, wall={record["wall_seconds"]:.1f}s',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--pilot',type=int,default=512)
    parser.add_argument('--selection',type=int,default=512)
    parser.add_argument('--paths',type=int,default=1024)
    parser.add_argument('--batch-size',type=int,default=256)
    parser.add_argument('--trials',type=int,default=2)
    parser.add_argument('--seed',type=int,default=20261026)
    args=parser.parse_args()
    if args.output.exists() or min(args.pilot,args.selection,args.paths,args.batch_size)<8 or args.trials<1:
        parser.error('Require fresh output and valid budgets')
    if args.pilot%2 or any(n%args.batch_size==1 for n in (args.pilot,args.selection,args.paths)):
        parser.error('Require even pilot count; avoid single-path final batches')
    run(args)
