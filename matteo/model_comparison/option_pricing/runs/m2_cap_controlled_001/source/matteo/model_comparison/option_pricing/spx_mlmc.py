"""Frozen SPX MLMC pilot, standalone level costs, and independent direct checks."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .fine_grid import CachedNormalizerPricer
from .spx_hedging_study import choose_quotes,refinement,controls,beta_fit
from .spx_path_importance import weighted_controls
from .spx_volatility_sampling import proposal
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.mlmc import allocate_levels,independent_level_statistics
from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.pricing import dated_grid,coarsen_paths
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.variance_reduction import deterministic_stock_controls


def fit_controls(y,c):
    return np.stack([beta_fit(y[:,j],c[:,j]) for j in range(y.shape[1])])


def draw_estimators(y,c,fits,forwards):
    return np.stack([y]+[y-np.einsum('njk,jk->nj',c,beta) for beta in fits],axis=-1)*forwards[None,:,None]


def collect(p,parameters,curve,initial,quotes,grids,bases,level,paths,args,role,trial,direct=False):
    parts={key:[] for key in ('y','c','fine_y','fine_c','w','z','components','monitors','hits')}
    seeds,costs=[],[];maturities=quotes.T_years.to_numpy()
    for batch,start in enumerate(range(0,paths,args.batch_size)):
        n=min(args.batch_size,paths-start)
        seed=int(np.random.SeedSequence([args.seed,role,trial,level,batch,int(direct)]).generate_state(1)[0]);seeds.append(seed)
        begin=time.perf_counter()
        sampled=prepare_volatility_paths(p,grids[level],initial,n,seed,bases[level],proposal(1))
        cache,w=sampled.cache,sampled.weights
        reference_paths=deterministic_stock_controls(cache,curve,np.unique(maturities))
        reference=np.stack([reference_paths[t] for t in maturities],axis=1)
        def price(cached):
            h=spx_hedging_controls(p,cached,parameters,curve,maturities,quotes.K_over_F.to_numpy(),quotes.option_type.to_numpy(),order=48)
            stock=np.stack([h.terminal_ratios[t] for t in maturities],axis=1)
            pairs=controls(dict(stock=stock,hedges=h.hedges),reference,quotes,curve)
            return np.stack([a for a,b in pairs],axis=1),np.stack([b for a,b in pairs],axis=1),stock
        fine_y,fine_c,stock=price(cache)
        fine_seconds=time.perf_counter()-begin
        if level and not direct:
            coarse_y,coarse_c,_=price(coarsen_paths(cache,grids[level-1]))
            y,c=fine_y-coarse_y,fine_c-coarse_c
            hits=np.stack([fine_y>0,coarse_y>0,y!=0],axis=-1)
        else:
            y,c=fine_y,fine_c
            hits=np.stack([fine_y>0,fine_y>0,y!=0],axis=-1)
        weighted=[weighted_controls(y[:,j],c[:,j],w) for j in range(len(quotes))]
        wy,wc=np.stack([a for a,b in weighted],axis=1),np.stack([b for a,b in weighted],axis=1)
        fine_weighted=[weighted_controls(fine_y[:,j],fine_c[:,j],w) for j in range(len(quotes))]
        fy,fc=np.stack([a for a,b in fine_weighted],axis=1),np.stack([b for a,b in fine_weighted],axis=1)
        seconds=time.perf_counter()-begin
        values=dict(y=wy,c=wc,fine_y=fy,fine_c=fc,w=w,z=sampled.score,components=sampled.components,
            monitors=np.c_[w,w*stock[:,0],w*stock[:,-1],w*reference[:,0],w*reference[:,-1]],hits=hits)
        for key,value in values.items(): parts[key].append(value)
        costs.append(dict(role=role,trial=trial,level=level,direct=direct,batch=batch,paths=n,
                          seconds=seconds,fine_seconds=fine_seconds))
        print(f'role={role} trial={trial} level={level} direct={direct} paths={start+n}/{paths}',flush=True)
    return {key:np.concatenate(value) for key,value in parts.items()},seeds,costs


def run(args):
    args.output.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();folder=Path(__file__).resolve().parent
    fit_path=folder/'runs/m2_vix_curve_fit_001/run.json'
    captured=capture_sources([Path(__file__).resolve(),fit_path,folder/'active_scope.json',
        *[folder/(s+'.py') for s in ('m2_spx','fine_grid','spx_hedging_study','spx_path_importance','spx_volatility_sampling')],
        *sorted((M2/'src/esn_eq').glob('*.py'))])
    save_sources(args.output,captured,ROOT)
    (args.output/'status.json').write_text('{"status":"setup"}\n')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');quotes=choose_quotes(targets)
    architecture,parameters,premium,hashes=frozen_stage3();p=CachedNormalizerPricer(architecture,premium)
    initial,_=p.initial_state([0.,0.,0.]);fit=json.loads(fit_path.read_text())
    if fit['target_sha256']!=targets.artifact_sha256 or fit['initial_state']!=initial.tolist():raise ValueError('Frozen inputs changed')
    for k,v in dict(architecture=architecture,parameters=parameters,premium=premium).items():
        if fit[k]!=json.loads(json.dumps(asdict(v))):raise ValueError('Frozen model changed')
    curve=NuisanceCurveFamily.from_targets(targets).build(fit['selected_levels'])
    events=np.r_[quotes.T_years.unique(),7/365,21/365,[t for t in curve.ends if t<quotes.T_years.max()]]
    grids=[dated_grid(events,1/(252*64))]
    for _ in range(3):grids.append(refinement(grids[-1]))
    bases=[volatility_path_basis(p,g,7/365) for g in grids]
    for g in grids:p.normalizer(g[:-1],initial,parameters,48)
    setup=time.perf_counter()-begin
    archive,coefficients,seeds,costs,allocation_variances,allocation_costs={},{},[],[],[],[]
    forwards=quotes.forward.to_numpy();scales=quotes.price_scale_U.to_numpy();central=[2,3,8,9]
    direct_coefficients=None;direct_cost=None
    for level in range(4):
        training,ss,cc=collect(p,parameters,curve,initial,quotes,grids,bases,level,args.pilot,args,0,0)
        seeds.append(dict(role=0,level=level,trial=0,direct=False,seeds=ss));costs+=cc
        splits=[np.arange(args.pilot),np.arange(args.pilot//2),np.arange(args.pilot//2,args.pilot)]
        coefficients[level]=np.stack([fit_controls(training['y'][i],training['c'][i]) for i in splits])
        if level==3:
            direct_coefficients=np.stack([fit_controls(training['fine_y'][i],training['fine_c'][i]) for i in splits])
        for key,value in training.items():archive[f'train_l{level}_{key}']=value
        selection,ss,cc=collect(p,parameters,curve,initial,quotes,grids,bases,level,args.allocation,args,1,0)
        seeds.append(dict(role=1,level=level,trial=0,direct=False,seeds=ss));costs+=cc
        for key,value in selection.items():archive[f'alloc_l{level}_{key}']=value
        x=draw_estimators(selection['y'],selection['c'],coefficients[level],forwards)[:,:,1]
        allocation_variances.append(float((x[:,central]/scales[central]).var(axis=0,ddof=1).mean()))
        # Use standalone generation/coupled pricing times from all design batches.
        timing=[r for r in costs if r['level']==level]
        allocation_costs.append(sum(r['seconds'] for r in timing)/sum(r['paths'] for r in timing))
        if level==3:direct_cost=sum(r['fine_seconds'] for r in timing)/sum(r['paths'] for r in timing)
    budget=args.direct_paths*direct_cost
    counts=allocate_levels(allocation_variances,allocation_costs,budget,args.minimum,64)
    allocation=dict(objective='mean controlled sampling variance / quote_scale² on fixed contracts 2,3,8,9',
        variances=allocation_variances,cost_per_path=allocation_costs,direct_cost_per_path=direct_cost,
        direct_paths=args.direct_paths,budget_seconds=budget,counts=counts.tolist(),
        predicted_cost_seconds=float(counts@allocation_costs),minimum=args.minimum,quantum=64,
        checking_not_used=True,coefficients={str(k):v.tolist() for k,v in coefficients.items()},direct_coefficients=direct_coefficients.tolist())
    (args.output/'allocation.json').write_text(json.dumps(allocation,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(args.output/'design_draws.npz',**archive);archive={}
    (args.output/'status.json').write_text('{"status":"checking; allocation frozen"}\n')
    estimates,level_rows,diagnostics=[],[],[]
    for trial in range(args.trials):
        levels=[]
        for level,n in enumerate(counts):
            data,ss,cc=collect(p,parameters,curve,initial,quotes,grids,bases,level,int(n),args,2,trial)
            seeds.append(dict(role=2,level=level,trial=trial,direct=False,seeds=ss));costs+=cc
            prefix=f'mlmc_t{trial}_l{level}'
            for key,value in data.items():archive[prefix+'_'+key]=value
            x=draw_estimators(data['y'],data['c'],coefficients[level],forwards);levels.append(x);archive[prefix+'_adjusted']=x
            for j,row in enumerate(quotes.itertuples()):
                for k,estimator in enumerate(('raw','pooled','halfA','halfB')):
                    level_rows.append(dict(trial=trial,level=level,steps=64*2**level,contract=j,estimator=estimator,paths=int(n),
                        mean=float(x[:,j,k].mean()),se=float(x[:,j,k].std(ddof=1)/np.sqrt(n)),variance=float(x[:,j,k].var(ddof=1)),
                        scale=row.price_scale_U,fine_hits=int(data['hits'][:,j,0].sum()),coarse_hits=int(data['hits'][:,j,1].sum()),nonzero_corrections=int(data['hits'][:,j,2].sum())))
            for j,label in enumerate(('mass','stock29','stock57','reference29','reference57')):
                m=data['monitors'][:,j];diagnostics.append(dict(trial=trial,method='mlmc',level=level,statistic=label,mean=float(m.mean()),se=float(m.std(ddof=1)/np.sqrt(len(m)))))
        mlmc_mean,mlmc_se=independent_level_statistics(levels)
        direct,ss,cc=collect(p,parameters,curve,initial,quotes,grids,bases,3,args.direct_paths,args,3,trial,direct=True)
        seeds.append(dict(role=3,level=3,trial=trial,direct=True,seeds=ss));costs+=cc
        prefix=f'direct_t{trial}'
        for key,value in direct.items():archive[prefix+'_'+key]=value
        dx=draw_estimators(direct['y'],direct['c'],direct_coefficients,forwards);archive[prefix+'_adjusted']=dx
        for j,label in enumerate(('mass','stock29','stock57','reference29','reference57')):
            m=direct['monitors'][:,j];diagnostics.append(dict(trial=trial,method='direct',level=3,statistic=label,mean=float(m.mean()),se=float(m.std(ddof=1)/np.sqrt(len(m)))))
        for j,row in enumerate(quotes.itertuples()):
            for k,estimator in enumerate(('raw','pooled','halfA','halfB')):
                for method,mean,se in (('mlmc',mlmc_mean[j,k],mlmc_se[j,k]),('direct',dx[:,j,k].mean(),dx[:,j,k].std(ddof=1)/np.sqrt(len(dx)))):
                    estimates.append(dict(trial=trial,method=method,contract=j,expiry_days=row.expiry_days,strike=row.strike,
                        option_type=row.option_type,estimator=estimator,price=float(mean),se=float(se),scale=row.price_scale_U,
                        sampling_screen=bool(se<=row.price_scale_U/8),note='Empirical screen only; inspect level/tail hits'))
    check_sources(captured)
    record=dict(status='completed MLMC prototype; targets discretized 512-step price, no continuous-model clearance',
        arguments=vars(args)|dict(output=str(args.output)),architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),
        initial_state=initial.tolist(),curve_levels=list(curve.levels),target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,
        implementation_sha256=source_hashes(captured,ROOT),grid_intervals=[len(g)-1 for g in grids],max_lambda_dt=[float(p.layout.rates.max()*np.diff(g).max()) for g in grids],
        proposal='vol7 fixed from prior design; basis rebuilt at each independently sampled fine level',normalizer_order=48,seed_records=seeds,
        uncertainty_unit='iid within level, independent levels/direct and trials, conditional on frozen training/allocation',
        setup_seconds=setup,wall_seconds=time.perf_counter()-begin,
        limitations=['Twelve training representatives; all original SPX/VIX scopes retained.',
            'MLMC telescopes target expectations, not individual independently drawn paths.',
            'Each level may use a different proposal; exact full-path weights recover the same target coarse marginal.',
            'Allocation minimizes four central contract sampling variances, not all tails or finest-grid bias.',
            'Budget is matched by pilot-predicted cost; realized times must be compared and reported.',
            'Setup/pilot cost excluded from checking budget and reported separately; amortization is not assumed.',
            '512 steps and GH48 are approximations; no extrapolation, bias-order claim or population accuracy certificate.',
            'Zero observed payoff/variance is not a tail-price accuracy guarantee; signed corrections and negative estimates retained.'])
    for name,rows in (('estimates',estimates),('levels',level_rows),('diagnostics',diagnostics),('costs',costs)):
        pd.DataFrame(rows).to_csv(args.output/(name+'.csv'),index=False)
    quotes.to_csv(args.output/'contracts.csv',index=False)
    np.savez_compressed(args.output/'checking_draws.npz',**archive)
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    (args.output/'status.json').write_text('{"status":"completed"}\n')
    print(f'Completed MLMC counts={counts.tolist()} in {record["wall_seconds"]:.1f}s',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--pilot',type=int,default=512)
    parser.add_argument('--allocation',type=int,default=256)
    parser.add_argument('--direct-paths',type=int,default=1024)
    parser.add_argument('--batch-size',type=int,default=256)
    parser.add_argument('--minimum',type=int,default=64)
    parser.add_argument('--trials',type=int,default=2)
    parser.add_argument('--seed',type=int,default=20261027)
    args=parser.parse_args()
    if args.output.exists() or min(args.pilot,args.allocation,args.direct_paths,args.batch_size)<8 or args.trials<1 or args.pilot%2 or args.minimum<64 or args.minimum%64:
        parser.error('Require new output, even pilot and positive budgets/minimum in 64-path blocks')
    if any(n%args.batch_size==1 for n in (args.pilot,args.allocation,args.direct_paths)):
        parser.error('Avoid one-path final batches')
    run(args)
