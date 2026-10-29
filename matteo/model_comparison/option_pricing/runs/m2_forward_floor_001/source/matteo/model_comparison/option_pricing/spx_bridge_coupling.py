"""Bounded reflected-Gaussian correction experiment with actual ordinary timings."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .fine_grid import CachedNormalizerPricer
from .spx_hedging_study import choose_quotes,refinement,controls
from .spx_path_importance import weighted_controls
from .spx_volatility_sampling import proposal
from .spx_mlmc import fit_controls,draw_estimators
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.bridge_coupling import EndpointBridge
from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import dated_grid,coarsen_paths
from esn_eq.variance_reduction import deterministic_stock_controls


def collect(p,parameters,curve,initial,quotes,grids,basis,kernel,args,level,role,trial):
    count=args.pilot if role==0 else args.paths
    pieces={key:[] for key in ('plain_y','plain_c','bridge_y','bridge_c','w','z','anti_z','components','fine_payoffs','anti_payoffs','coarse_payoffs','monitors')}
    seeds,costs=[],[];maturities=quotes.T_years.to_numpy()
    def price(cache,reference):
        h=spx_hedging_controls(p,cache,parameters,curve,maturities,quotes.K_over_F.to_numpy(),quotes.option_type.to_numpy(),order=48)
        stock=np.stack([h.terminal_ratios[t] for t in maturities],axis=1)
        pairs=controls(dict(stock=stock,hedges=h.hedges),reference,quotes,curve)
        return np.stack([a for a,b in pairs],axis=1),np.stack([b for a,b in pairs],axis=1),stock
    def ref(cache):
        r=deterministic_stock_controls(cache,curve,np.unique(maturities))
        return np.stack([r[t] for t in maturities],axis=1)
    def weight(y,c,w):
        pairs=[weighted_controls(y[:,j],c[:,j],w) for j in range(len(quotes))]
        return np.stack([a for a,b in pairs],axis=1),np.stack([b for a,b in pairs],axis=1)
    for batch,start in enumerate(range(0,count,args.batch_size)):
        n=min(args.batch_size,count-start);seed=int(np.random.SeedSequence([args.seed,level,role,trial,batch]).generate_state(1)[0]);seeds.append(seed)
        def ordinary():
            begin=time.perf_counter();sample=prepare_volatility_paths(p,grids[level],initial,n,seed,basis,proposal(1))
            reference=ref(sample.cache);f,fc,_=price(sample.cache,reference)
            c,cc,_=price(coarsen_paths(sample.cache,grids[level-1]),reference)
            y,control=weight(f-c,fc-cc,sample.weights)
            return y,control,time.perf_counter()-begin
        def paired():
            begin=time.perf_counter();sample=kernel.draw(initial,n,seed,proposal(1));original=sample.original;w=original.weights
            reference=ref(original.cache);f,fc,stock=price(original.cache,reference)
            a,ac,anti_stock=price(sample.reflected,reference)
            c,cc,_=price(coarsen_paths(original.cache,grids[level-1]),reference)
            py,pc=weight(f-c,fc-cc,w);by,bc=weight(.5*(f+a)-c,.5*(fc+ac)-cc,w)
            eq_error=float(np.max(np.abs(original.cache.equity_increments[::2]+original.cache.equity_increments[1::2]-sample.reflected.equity_increments[::2]-sample.reflected.equity_increments[1::2])))
            if eq_error>5e-10:raise FloatingPointError('Coarse equity increments not preserved')
            values=dict(plain_y=py,plain_c=pc,bridge_y=by,bridge_c=bc,w=w,z=original.score,anti_z=sample.reflected_scores,
                components=original.components,fine_payoffs=f,anti_payoffs=a,coarse_payoffs=c,
                monitors=np.c_[w,w*stock[:,0],w*stock[:,-1],w*anti_stock[:,0],w*anti_stock[:,-1],w*reference[:,0],w*reference[:,-1]])
            return values,time.perf_counter()-begin,sample.endpoint_error,eq_error
        # Alternate execution order to limit a systematic timing-order effect.
        if batch%2==0:
            ordinary_y,ordinary_c,ordinary_time=ordinary();data,paired_time,endpoint_error,equity_error=paired()
        else:
            data,paired_time,endpoint_error,equity_error=paired();ordinary_y,ordinary_c,ordinary_time=ordinary()
        np.testing.assert_allclose(ordinary_y,data['plain_y'],rtol=1e-8,atol=1e-10)
        np.testing.assert_allclose(ordinary_c,data['plain_c'],rtol=1e-8,atol=1e-10)
        for key,value in data.items():pieces[key].append(value)
        costs.append(dict(role=role,trial=trial,level=level,batch=batch,pairs=n,ordinary_seconds=ordinary_time,bridge_seconds=paired_time,
            endpoint_error=endpoint_error,equity_error=equity_error,execution_order='ordinary_first' if batch%2==0 else 'bridge_first'))
        print(f'level={level} role={role} trial={trial} pairs={start+n}/{count}',flush=True)
    return {key:np.concatenate(value) for key,value in pieces.items()},seeds,costs


def fit_variants(data,method):
    y,c=data[method+'_y'],data[method+'_c'];n=len(y)
    return np.stack([fit_controls(y[i],c[i]) for i in (np.arange(n),np.arange(n//2),np.arange(n//2,n))])


def run(args):
    args.output.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();folder=Path(__file__).resolve().parent
    fit_path=folder/'runs/m2_vix_curve_fit_001/run.json'
    captured=capture_sources([Path(__file__).resolve(),fit_path,folder/'active_scope.json',
        *[folder/(name+'.py') for name in ('m2_spx','fine_grid','spx_hedging_study','spx_path_importance','spx_volatility_sampling','spx_mlmc')],
        *sorted((M2/'src/esn_eq').glob('*.py'))]);save_sources(args.output,captured,ROOT)
    (args.output/'status.json').write_text('{"status":"setup"}\n')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');quotes=choose_quotes(targets)
    architecture,parameters,premium,hashes=frozen_stage3();p=CachedNormalizerPricer(architecture,premium);initial,_=p.initial_state([0.,0.,0.])
    fit=json.loads(fit_path.read_text())
    if fit['target_sha256']!=targets.artifact_sha256 or fit['initial_state']!=initial.tolist():raise ValueError('Frozen inputs changed')
    for key,value in dict(architecture=architecture,parameters=parameters,premium=premium).items():
        if fit[key]!=json.loads(json.dumps(asdict(value))):raise ValueError('Frozen model changed')
    curve=NuisanceCurveFamily.from_targets(targets).build(fit['selected_levels'])
    events=np.r_[quotes.T_years.unique(),7/365,21/365,[t for t in curve.ends if t<quotes.T_years.max()]]
    grids=[dated_grid(events,1/(252*64))]
    for _ in range(3):grids.append(refinement(grids[-1]))
    bases={level:volatility_path_basis(p,grids[level],7/365) for level in (2,3)}
    kernels={level:EndpointBridge(p,grids[level],bases[level]) for level in (2,3)}
    for level in (1,2,3):p.normalizer(grids[level][:-1],initial,parameters,48)
    setup=time.perf_counter()-begin;archive,coefficients,seeds,costs={},{},[],[]
    for level in (2,3):
        data,ss,cc=collect(p,parameters,curve,initial,quotes,grids,bases[level],kernels[level],args,level,0,0)
        coefficients[level]={method:fit_variants(data,method) for method in ('plain','bridge')}
        for key,value in data.items():archive[f'pilot_l{level}_'+key]=value
        seeds.append(dict(level=level,role=0,trial=0,seeds=ss));costs+=cc
    frozen=dict(coefficients={str(level):{method:value.tolist() for method,value in fits.items()} for level,fits in coefficients.items()},
        checking_not_used=True,uncertainty_unit='independent original/reflected PAIRS, not individual fine paths',
        comparisons='no method selection; both correction couplings and all pilot fits retained')
    (args.output/'frozen.json').write_text(json.dumps(frozen,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(args.output/'pilot_draws.npz',**archive);archive={}
    (args.output/'status.json').write_text('{"status":"checking; coefficients frozen"}\n')
    rows,differences,diagnostics=[],[],[];forwards=quotes.forward.to_numpy()
    for trial in range(args.trials):
        for level in (2,3):
            data,ss,cc=collect(p,parameters,curve,initial,quotes,grids,bases[level],kernels[level],args,level,1,trial)
            prefix=f't{trial}_l{level}';seeds.append(dict(level=level,role=1,trial=trial,seeds=ss));costs+=cc
            for key,value in data.items():archive[prefix+'_'+key]=value
            adjusted={}
            for method in ('plain','bridge'):
                x=draw_estimators(data[method+'_y'],data[method+'_c'],coefficients[level][method],forwards)
                archive[prefix+'_'+method+'_adjusted']=x;adjusted[method]=x
                for j,row in enumerate(quotes.itertuples()):
                    for k,estimator in enumerate(('raw','pooled','halfA','halfB')):
                        rows.append(dict(trial=trial,level=level,coarse_steps=64*2**(level-1),fine_steps=64*2**level,
                            method=method,contract=j,expiry_days=row.expiry_days,strike=row.strike,option_type=row.option_type,estimator=estimator,pairs=args.paths,
                            correction=float(x[:,j,k].mean()),se=float(x[:,j,k].std(ddof=1)/np.sqrt(args.paths)),variance=float(x[:,j,k].var(ddof=1)),
                            scale=row.price_scale_U,fine_hits=int((data['fine_payoffs'][:,j]>0).sum()),
                            reflected_hits=int((data['anti_payoffs'][:,j]>0).sum()),coarse_hits=int((data['coarse_payoffs'][:,j]>0).sum())))
            diff=adjusted['bridge']-adjusted['plain']
            for j,row in enumerate(quotes.itertuples()):
                for k,estimator in enumerate(('raw','pooled','halfA','halfB')):
                    differences.append(dict(trial=trial,level=level,contract=j,estimator=estimator,scale=row.price_scale_U,
                        difference=float(diff[:,j,k].mean()),se=float(diff[:,j,k].std(ddof=1)/np.sqrt(args.paths))))
            for j,label in enumerate(('mass','stock29','stock57','reflected_stock29','reflected_stock57','reference29','reference57')):
                x=data['monitors'][:,j];diagnostics.append(dict(trial=trial,level=level,statistic=label,mean=float(x.mean()),se=float(x.std(ddof=1)/np.sqrt(len(x)))))
    check_sources(captured)
    record=dict(status='completed bridge-coupling diagnostic; no continuous-model or market-fit clearance',arguments=vars(args)|dict(output=str(args.output)),
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),curve_levels=list(curve.levels),
        frozen_model_sha256=hashes,target_sha256=targets.artifact_sha256,implementation_sha256=source_hashes(captured,ROOT),seed_records=seeds,
        grid_intervals=[len(g)-1 for g in grids],normalizer_order=48,kernel_diagnostics={str(level):kernel.diagnostics for level,kernel in kernels.items()},
        setup_seconds=setup,wall_seconds=time.perf_counter()-begin,
        uncertainty_unit='iid reflected pairs; methods/fits share original paths; trials independent conditional on frozen pilots',
        likelihood='same original f/q weight for original and reflected fine payoffs; target Gaussian reflection invariant',
        limitations=['12 training representatives; scopes and frozen model/curve/state unchanged.',
            'Corrections only: no new complete MLMC estimate or fine-price accuracy certificate.',
            'One shared-path contribution per pair, not two independent observations.',
            'Full ordinary comparator is rerun for actual timing and path correspondence; no checking-based selection.',
            'SVD completion may retain extra null directions if numerical rank deficient; variance gain not guaranteed.',
            'Costs are measured at equal pair counts; variance-cost gains are descriptive, not exact matched-runtime trials.',
            '512/GH48 remain approximate; zero observed payoffs and negative corrections are preserved.'])
    for name,values in (('estimates',rows),('paired_differences',differences),('diagnostics',diagnostics),('costs',costs)):
        pd.DataFrame(values).to_csv(args.output/(name+'.csv'),index=False)
    quotes.to_csv(args.output/'contracts.csv',index=False);np.savez_compressed(args.output/'checking_draws.npz',**archive)
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');(args.output/'status.json').write_text('{"status":"completed"}\n')
    print(f'Completed bridge coupling in {record["wall_seconds"]:.1f}s',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--pilot',type=int,default=256);parser.add_argument('--paths',type=int,default=512)
    parser.add_argument('--batch-size',type=int,default=256);parser.add_argument('--trials',type=int,default=2)
    parser.add_argument('--seed',type=int,default=20261028);args=parser.parse_args()
    if args.output.exists() or min(args.pilot,args.paths,args.batch_size)<8 or args.pilot%2 or args.trials<1:
        parser.error('Require fresh output, even pilot and positive budgets')
    if any(n%args.batch_size==1 for n in (args.pilot,args.paths)):parser.error('Avoid one-pair final batches')
    run(args)
