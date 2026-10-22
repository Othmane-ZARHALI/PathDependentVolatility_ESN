"""Frozen SPX price/correction tests of exact OU-innovation martingale controls."""
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
from .spx_path_score_controls import ridge_fit
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.innovation_controls import InnovationRecorder,spx_innovation_controls
from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import dated_grid,coarsen_paths
from esn_eq.variance_reduction import deterministic_stock_controls


METHODS=('raw','legacy','base','augmented','base_A','base_B','aug_A','aug_B')


def datasets(data):
    result={}
    for level in (1,2,3):
        result['price'+str(64*2**level)]=(data[f'g{level}_y'],data[f'g{level}_base'],data[f'g{level}_vol'])
    for low,high in ((1,2),(2,3)):
        # Mass control is w-1, not a difference of identical mass columns.
        b=data[f'g{high}_base']-data[f'g{low}_base'];b=b.copy();b[:,:,-1]=data['w'][:,None]-1.
        result[f'corr{64*2**low}_{64*2**high}']=(data[f'g{high}_y']-data[f'g{low}_y'],b,data[f'g{high}_vol']-data[f'g{low}_vol'])
    return result


def train(y,base,vol):
    n=len(y);full=np.arange(n);a=np.arange(n//2);b=np.arange(n//2,n);aug=np.concatenate([base,vol],axis=2)
    fits={}
    for name,c,index,regularized in (('legacy',base,full,False),('base',base,full,True),('augmented',aug,full,True),
            ('base_A',base,a,True),('base_B',base,b,True),('aug_A',aug,a,True),('aug_B',aug,b,True)):
        fitter=ridge_fit if regularized else beta_fit
        fits[name]=np.stack([fitter(y[index,j],c[index,j]) for j in range(y.shape[1])])
    return fits


def adjusted(y,base,vol,fits,forwards):
    aug=np.concatenate([base,vol],axis=2);values=[y]
    for name in METHODS[1:]:
        c=aug if name in ('augmented','aug_A','aug_B') else base
        values.append(y-np.einsum('njk,jk->nj',c,np.asarray(fits[name])))
    return np.stack(values,axis=-1)*forwards[None,:,None]


def collect(p,parameters,curve,initial,quotes,grids,basis,args,role,trial):
    count=args.pilot if role==0 else args.paths;pieces={};seeds=[];costs=[];maturities=quotes.T_years.to_numpy()
    for batch,start in enumerate(range(0,count,args.batch_size)):
        n=min(args.batch_size,count-start);seed=int(np.random.SeedSequence([args.seed,role,trial,batch]).generate_state(1)[0]);seeds.append(seed)
        def sample(augmented):
            begin=time.perf_counter()
            recorder=InnovationRecorder(p,grids[3],{level:grids[level] for level in (1,2,3)},n) if augmented else None
            sampled=prepare_volatility_paths(p,grids[3],initial,n,seed,basis,proposal(1),state_observer=recorder)
            cache,w=sampled.cache,sampled.weights
            r=deterministic_stock_controls(cache,curve,np.unique(maturities));reference=np.stack([r[t] for t in maturities],axis=1)
            result={};grid_costs={}
            for level in (1,2,3):
                tick=time.perf_counter();c=cache if level==3 else coarsen_paths(cache,grids[level])
                if augmented:
                    innovations=recorder.result(level)
                    h=spx_innovation_controls(p,c,innovations,parameters,curve,maturities,quotes.K_over_F.to_numpy(),quotes.option_type.to_numpy(),order=48)
                    hedges=h.stock_hedges
                else:
                    h=spx_hedging_controls(p,c,parameters,curve,maturities,quotes.K_over_F.to_numpy(),quotes.option_type.to_numpy(),order=48);hedges=h.hedges
                stock=np.stack([h.terminal_ratios[t] for t in maturities],axis=1)
                pairs=controls(dict(stock=stock,hedges=hedges),reference,quotes,curve)
                weighted=[weighted_controls(y,control,w) for y,control in pairs]
                result[f'g{level}_y']=np.stack([a for a,b in weighted],axis=1);result[f'g{level}_base']=np.stack([b for a,b in weighted],axis=1)
                if augmented:
                    result[f'g{level}_vol']=w[:,None,None]*h.volatility_hedges
                    eps=innovations.values;dw=c.equity_increments
                    raw=np.concatenate([eps.sum(axis=0),(eps*eps-np.diagonal(innovations.covariance,axis1=1,axis2=2)[:,None,:]).sum(axis=0),
                        (eps*dw[:,:,None]-innovations.equity_covariance[:,None,:]).sum(axis=0)],axis=1)
                    result[f'g{level}_innovation_moments']=w[:,None]*raw
                if level==3:result['monitors']=np.c_[w,w*stock[:,0],w*stock[:,-1],w*reference[:,0],w*reference[:,-1]]
                grid_costs[level]=time.perf_counter()-tick
            result.update(w=w,z=sampled.score,components=sampled.components)
            return result,time.perf_counter()-begin,grid_costs
        if batch%2==0:
            ordinary,old_time,old_grids=sample(False);data,new_time,new_grids=sample(True)
        else:
            data,new_time,new_grids=sample(True);ordinary,old_time,old_grids=sample(False)
        for level in (1,2,3):
            np.testing.assert_allclose(data[f'g{level}_y'],ordinary[f'g{level}_y'],rtol=1e-10,atol=1e-12)
            np.testing.assert_allclose(data[f'g{level}_base'],ordinary[f'g{level}_base'],rtol=1e-10,atol=1e-12)
        for key,value in data.items():pieces.setdefault(key,[]).append(value)
        costs.append(dict(role=role,trial=trial,batch=batch,paths=n,ordinary_seconds=old_time,augmented_seconds=new_time,
            execution_order='ordinary_first' if batch%2==0 else 'augmented_first',
            **{f'ordinary_grid{level}_seconds':old_grids[level] for level in (1,2,3)},**{f'augmented_grid{level}_seconds':new_grids[level] for level in (1,2,3)}))
        print(f'role={role} trial={trial} paths={start+n}/{count}',flush=True)
    return {key:np.concatenate(value) for key,value in pieces.items()},seeds,costs


def run(args):
    args.output.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();folder=Path(__file__).resolve().parent;fit_path=folder/'runs/m2_vix_curve_fit_001/run.json'
    captured=capture_sources([Path(__file__).resolve(),fit_path,folder/'active_scope.json',
        *[folder/(s+'.py') for s in ('m2_spx','fine_grid','spx_hedging_study','spx_path_importance','spx_volatility_sampling','spx_path_score_controls')],
        *sorted((M2/'src/esn_eq').glob('*.py'))]);save_sources(args.output,captured,ROOT)
    (args.output/'status.json').write_text('{"status":"setup"}\n')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');quotes=choose_quotes(targets);architecture,parameters,premium,hashes=frozen_stage3()
    p=CachedNormalizerPricer(architecture,premium);initial,_=p.initial_state([0.,0.,0.]);fit=json.loads(fit_path.read_text())
    if fit['target_sha256']!=targets.artifact_sha256 or fit['initial_state']!=initial.tolist():raise ValueError('Frozen inputs changed')
    for key,value in dict(architecture=architecture,parameters=parameters,premium=premium).items():
        if fit[key]!=json.loads(json.dumps(asdict(value))):raise ValueError('Frozen model changed')
    curve=NuisanceCurveFamily.from_targets(targets).build(fit['selected_levels'])
    events=np.r_[quotes.T_years.unique(),7/365,21/365,[t for t in curve.ends if t<quotes.T_years.max()]];grids=[dated_grid(events,1/(252*64))]
    for _ in range(3):grids.append(refinement(grids[-1]))
    basis=volatility_path_basis(p,grids[3],7/365)
    normalizers={}
    for level in (1,2,3):
        normalizers[f'times{level}']=grids[level][:-1];normalizers[f'values{level}']=p.normalizer(grids[level][:-1],initial,parameters,48)
    np.savez_compressed(args.output/'normalizers.npz',**normalizers);setup=time.perf_counter()-begin
    training,ss,costs=collect(p,parameters,curve,initial,quotes,grids,basis,args,0,0)
    coefficients={name:train(*triple) for name,triple in datasets(training).items()}
    frozen=dict(coefficients={name:{method:value.tolist() for method,value in fits.items()} for name,fits in coefficients.items()},checking_not_used=True,
        primary_comparison='augmented ridge14 vs baseline ridge5; ridge=1 for both',legacy_comparison='original OLS5 retained separately',
        control_dictionary='linear eps, eps²-Rii, eps*dW-Cov(eps,dW), each times predictable log-variance Black vega/readout derivative')
    (args.output/'frozen.json').write_text(json.dumps(frozen,indent=2,allow_nan=False)+'\n');np.savez_compressed(args.output/'pilot_draws.npz',**training)
    seeds=[dict(role=0,trial=0,seeds=ss)];archive={};rows=[];paired=[];diagnostics=[];forwards=quotes.forward.to_numpy()
    (args.output/'status.json').write_text('{"status":"checking; coefficients frozen"}\n')
    for trial in range(args.trials):
        data,ss,cc=collect(p,parameters,curve,initial,quotes,grids,basis,args,1,trial);costs+=cc;seeds.append(dict(role=1,trial=trial,seeds=ss))
        for key,value in data.items():archive[f't{trial}_'+key]=value
        for name,triple in datasets(data).items():
            x=adjusted(*triple,coefficients[name],forwards);archive[f't{trial}_{name}_adjusted']=x
            for j,row in enumerate(quotes.itertuples()):
                for k,method in enumerate(METHODS):
                    rows.append(dict(trial=trial,target=name,method=method,contract=j,expiry_days=row.expiry_days,strike=row.strike,
                        option_type=row.option_type,estimate=float(x[:,j,k].mean()),se=float(x[:,j,k].std(ddof=1)/np.sqrt(len(x))),
                        variance=float(x[:,j,k].var(ddof=1)),scale=row.price_scale_U,nonzero_raw=int((triple[0][:,j]!=0).sum())))
                difference=x[:,j,3]-x[:,j,2];paired.append(dict(trial=trial,target=name,contract=j,
                    difference=float(difference.mean()),se=float(difference.std(ddof=1)/np.sqrt(len(difference))),scale=row.price_scale_U))
        for j,label in enumerate(('mass','stock29','stock57','reference29','reference57')):
            values=data['monitors'][:,j];diagnostics.append(dict(trial=trial,target='forward',statistic=label,mean=float(values.mean()),se=float(values.std(ddof=1)/np.sqrt(len(values))),expectation=1.))
        for level in (1,2,3):
            for j in range(9):
                values=data[f'g{level}_innovation_moments'][:,j]
                diagnostics.append(dict(trial=trial,target=f'innovations{64*2**level}',statistic=str(j),mean=float(values.mean()),se=float(values.std(ddof=1)/np.sqrt(len(values))),expectation=0.))
    check_sources(captured)
    record=dict(status='completed innovation-control diagnostic; no continuous-model/calibration clearance',arguments=vars(args)|dict(output=str(args.output)),
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),curve_levels=list(curve.levels),
        target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,implementation_sha256=source_hashes(captured,ROOT),seed_records=seeds,
        grid_intervals=[len(g)-1 for g in grids],normalizer_order=48,setup_seconds=setup,wall_seconds=time.perf_counter()-begin,
        uncertainty_unit='iid checking paths conditional on frozen pilots; same paths across controls, grids and targets',
        limitations=['12 training representatives; all original scopes retained; no parameter/state/curve fitting.',
            'Control means are under target Q; multiply whole path gains by original likelihood, never recenter proposal drift.',
            'Primary ridge comparison holds penalty fixed; legacy OLS separates regularization effects.',
            'Actual cost covers one joint 128/256/512-grid diagnostic, including recording and innovation diagnostics; not standalone price-only cost.',
            'Direct-price/correction samples share finest paths; no new MLMC allocation or complete estimator is constructed.',
            'GH48 and 512-step discretization remain approximate; variance reduction does not remove bias.',
            'Zero observed payoffs and signed/negative adjusted estimates retained; no clipping.'])
    for name,values in (('estimates',rows),('paired_differences',paired),('diagnostics',diagnostics),('costs',costs)):
        pd.DataFrame(values).to_csv(args.output/(name+'.csv'),index=False)
    quotes.to_csv(args.output/'contracts.csv',index=False);np.savez_compressed(args.output/'checking_draws.npz',**archive)
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');(args.output/'status.json').write_text('{"status":"completed"}\n')
    print(f'Completed OU-innovation controls in {record["wall_seconds"]:.1f}s',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--pilot',type=int,default=512);parser.add_argument('--paths',type=int,default=512)
    parser.add_argument('--batch-size',type=int,default=128);parser.add_argument('--trials',type=int,default=2);parser.add_argument('--seed',type=int,default=20261029)
    args=parser.parse_args()
    if args.output.exists() or min(args.pilot,args.paths,args.batch_size)<8 or args.pilot%2 or args.trials<1:
        parser.error('Require fresh output, even pilot and valid budgets')
    if any(n%args.batch_size==1 for n in (args.pilot,args.paths)):parser.error('Avoid one-path final batches')
    run(args)
