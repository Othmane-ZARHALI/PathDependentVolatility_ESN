"""Seven SPX maturities and four VIX maturities, all observations in sample.

Only the three variance levels are refitted. State/readout/Q stay at the prior
working reference. Existing compatible paths/exposures are reused; previously
unpriced SPX expiries receive new primitives and fixed known-zero controls.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .m2_spx import ROOT,M2,frozen_stage3
from .all_maturity_core import populations,nearest_same_side,AllMaturityData
from .joint_curve_fit import SPXData,VIXData,dump
from .joint_curve_core import CurvePathBatch,donor_indices
from .joint_state_core import AdaptiveDatedNormalizer,StateDatedPricer,projected_initial_shift
from .spx_conditional_control_study import setup
from .spx_volatility_sampling import proposal
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from .audit_joint_curve_fit import score
from esn_eq.pricing import DatedM2Pricer,ForwardVarianceCurve,dated_grid
from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.variance_reduction import deterministic_stock_controls
from esn_eq.conditional_spx import lognormal_payoffs

HERE=Path(__file__).resolve().parent
PARENT=HERE/'runs/m2_joint_state_fit_003'
OLD=HERE/'runs/m2_joint_curve_fit_002'
NEW_DAYS=(119,666,1030)
OLD_DAYS=(29,57,211,302)


def representatives(q):
    indices=[]
    for side in ['put','call']:
        subset=q[q.option_type==side]
        for z in [.5,1.5,2.5]:
            indices.append((subset.scaled_moneyness-z).abs().idxmin())
    return q.loc[list(dict.fromkeys(indices))].reset_index(drop=True)


def worker(output,days,role):
    design=json.loads((output/'design.json').read_text())
    destination=output/f'new_spx_{role}_d{days}';destination.mkdir(exist_ok=False)
    captured=capture_sources([ROOT/p for p in design['implementation_sha256']])
    q=pd.read_csv(output/'spx_contracts.csv');q=q[q.expiry_days==days].reset_index(drop=True)
    architecture,parameters,premium,_=frozen_stage3()
    p=StateDatedPricer(architecture,premium);zero=np.zeros(len(p.layout.rates))
    initial=np.asarray(design['initial_state']);T=float(q.T_years.iloc[0])
    ends=tuple(design['ends_years']);anchor=ForwardVarianceCurve(ends,tuple(design['control_anchor_levels']))
    events=np.r_[T,7/365,21/365,[t for t in ends if t<T]]
    grid=dated_grid(events,1/(252*64))
    for _ in range(3):grid=np.r_[0.,np.concatenate([np.linspace(a,b,3)[1:] for a,b in zip(grid[:-1],grid[1:])])]
    g0=AdaptiveDatedNormalizer(p,parameters,zero,T)
    g1=AdaptiveDatedNormalizer(p,parameters,initial,T)
    p.state_normalizer=g0
    basis=volatility_path_basis(p,grid,7/365)
    reps=representatives(q);donor=donor_indices(q,reps)
    frozen=np.load(output/f'new_beta_d{days}.npz')['beta']
    np.save(destination/'grid.npy',grid)
    controls=[];weights=[];stocks=[];seeds=[]
    begin=time.perf_counter()
    for batch in range(design['paths']//design['batch_size']):
        seed=int(np.random.SeedSequence([design['seed'],days,0 if role=='search' else 1,batch]).generate_state(1)[0])
        r=prepare_volatility_paths(p,grid,zero,design['batch_size'],seed,basis,proposal(1));cache=r.cache
        h=spx_hedging_controls(p,cache,parameters,anchor,reps.T_years.to_numpy(),reps.K_over_F.to_numpy(),reps.option_type.to_numpy())
        stock=h.terminal_ratios[T]
        ref=deterministic_stock_controls(cache,anchor,[T])[T]
        total=anchor.integral(0.,T)
        known=lognormal_payoffs(-.5*total,total,q.K_over_F.to_numpy(),q.option_type.to_numpy())
        sign=np.where(q.option_type.to_numpy()=='call',1.,-1.)
        c=np.concatenate([np.broadcast_to((stock-1)[:,None,None],(len(stock),len(q),1)),
            (np.maximum(sign*(ref[:,None]-q.K_over_F.to_numpy()),0.)-known)[:,:,None],h.hedges[:,donor]],axis=2)
        c=np.concatenate([r.weights[:,None,None]*c,np.broadcast_to((r.weights-1)[:,None,None],(len(stock),len(q),1))],axis=2)
        controls.append(c);weights.append(r.weights);stocks.append(stock);seeds.append(seed)
        stem=Path(design['scratch'])/f'{role}_d{days}_b{batch}'
        np.save(str(stem)+'_equity.npy',cache.equity_increments)
        ratio=np.lib.format.open_memmap(str(stem)+'_ratio.npy',mode='w+',dtype='float64',shape=cache.equity_increments.shape)
        for k in range(0,len(grid)-1,2048):
            times=grid[k:min(k+2048,len(grid)-1)]
            projections=cache.projections[k:k+len(times)]+projected_initial_shift(p,times,initial)[:,None,:]
            ratio[k:k+len(times)]=p.multiplier(projections,parameters)/g1(times)[:,None]
        ratio.flush();del ratio,cache,r,h
        dump(destination/'status.json',dict(batch=batch+1,paths=(batch+1)*design['batch_size'],seconds=time.perf_counter()-begin))
        print(f'{role} SPX {days} days: {(batch+1)*design["batch_size"]} paths',flush=True)
    np.savez_compressed(destination/'controls.npz',controls=np.concatenate(controls),weights=np.concatenate(weights),
        anchor_stock=np.concatenate(stocks),beta=frozen)
    check_sources(captured)
    dump(destination/'completed.json',dict(days=days,role=role,seeds=seeds,grid_steps=len(grid)-1,
        normalizer_zero=g0.audit,normalizer_selected=g1.audit,seconds=time.perf_counter()-begin,
        controls='New maturity zero-state anchor stock/reference/representative adapted hedges; old independent frozen beta transferred with forward scaling, not refitted on searching/checking paths.'))


def run(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    begin=time.perf_counter()
    targets,p,parameters,zero,anchor,panel,master,_,hashes=setup()
    parent=json.loads((PARENT/'run.json').read_text());old=json.loads((OLD/'run.json').read_text())
    initial=parent['selected']['state']['initial'];baseline=np.array(parent['selected']['levels']);ends=list(anchor.ends)
    spx,vix,membership,proxies=populations(targets.quotes)
    forwards=targets.vix_forwards.sort_values('expiry').reset_index(drop=True)
    for name,table in [('spx_contracts',spx),('vix_contracts',vix),('spx_membership',membership),('atm_proxies',proxies),('forward_targets',forwards)]:table.to_csv(out/(name+'.csv'),index=False)
    sources=[Path(__file__),HERE/'all_maturity_core.py',HERE/'joint_state_core.py',HERE/'joint_curve_core.py',
        HERE/'joint_curve_fit.py',HERE/'audit_joint_curve_fit.py',HERE/'spx_conditional_control_study.py',
        HERE/'spx_volatility_sampling.py',HERE/'active_scope.json',HERE/'ALL_MATURITY_FIT_PLAN.md',
        HERE/'scope_decisions/2026-10-23-before-all-maturity-fit.json',PARENT/'run.json',OLD/'run.json',
        OLD/'spx_contracts.csv',OLD/'vix_contracts.csv',*sorted((M2/'src/esn_eq').glob('*.py'))]
    for role,oldrole in [('search','search'),('check','fresh')]:
        for days in OLD_DAYS:sources.append(PARENT/f'controls_{oldrole}_t0_d{days}.npz')
        for days in [27,55,118,153]:sources.append(PARENT/f'vix_state2_{oldrole}_t0_d{days}.npz')
    captured=capture_sources(sources);save_sources(out,captured,ROOT)
    design=dict(scope='ALL retained market observations in sample; no out-of-sample market evaluation.',
        spx_expiry_days=[29,57,119,211,302,666,1030],excluded_spx_days=[1,1401,1765,2129],cutoff=3,
        spx_quotes=len(spx),vix_quotes=len(vix),vix_expiry_days=[27,55,118,153],
        architecture=asdict(p.architecture),parameters=asdict(parameters),premium=asdict(p.premium),initial_state=initial,
        state_penalty=.1*parent['selected']['state']['prior_squared_norm']/3,
        baseline_levels=baseline.tolist(),ends_years=ends,control_anchor_levels=list(anchor.levels),
        implementation_sha256=source_hashes(captured,ROOT),target_sha256=targets.artifact_sha256,
        seed=args.seed,paths=512,batch_size=128,scratch=tempfile.mkdtemp(prefix='m2-all-maturity-',dir='/private/tmp'),
        old_scratch=parent['scratch'],parent_run=str(PARENT),old_master_steps=len(master)-1,
        budget=dict(curve_starts=3,max_nfev_per_start=30,optimization_seconds=1800,new_maturity_workers=2,
            stock_steps_per_trading_day=512,normalizer_GH=48,normalizer_interpolation_screen=2e-7,
            vix_conditional_GH=32,vix_time_order=8,vix_nets=8,vix_points=64,numerical_checks=1),
        check='Same in-sample market quotes. Reused prior numerical search/fresh-t0 streams at four old SPX and all VIX dates; new independent search/check streams at added SPX dates. Not new independent complete validation or market out-of-sample results.',
        controls='Existing fixed controls/beta at original dates, nearest same-side market-key transfer for added strikes. New maturities use original five-control family and frozen nearest-302-day beta scaled by market forward ratio; valid known-zero controls, extrapolated efficiency unqualified.',
        limitations=['Exploratory path budgets, not price convergence.','Initial state/readout/Q fixed; only three curve levels fitted.','All seven SPX expiry weights balanced within family; VIX/forward family weights and price scales unchanged.','Do not compare absolute loss with older four-maturity objectives.'])
    dump(out/'design.json',design)
    donor_tables=[]
    oldspx=pd.read_csv(OLD/'spx_contracts.csv');oldvix=pd.read_csv(OLD/'vix_contracts.csv')
    for days in NEW_DAYS:
        q=spx[spx.expiry_days==days].reset_index(drop=True)
        donors=oldspx[oldspx.expiry_days==302].reset_index(drop=True);ix=nearest_same_side(q,donors)
        beta=np.load(PARENT/'controls_search_t0_d302.npz')['beta'][ix]
        beta=beta*(q.forward.to_numpy()/donors.forward.to_numpy()[ix])[:,None]
        np.savez_compressed(out/f'new_beta_d{days}.npz',beta=beta)
        t=q[['expiry_days','strike','option_type','original_split']].copy()
        t['donor_expiry_days']=302;t['donor_strike']=donors.strike.to_numpy()[ix];donor_tables.append(t)
    pd.concat(donor_tables).to_csv(out/'new_maturity_beta_donors.csv',index=False)
    scratch=Path(design['scratch']);oldscratch=Path(parent['scratch'])
    def status(stage,**kw):dump(out/'status.json',dict(stage=stage,seconds=time.perf_counter()-begin,**kw));print(stage,json.dumps(kw),flush=True)
    def load(role):
        ds={};dv={};oldrole='search' if role=='search' else 'fresh'
        for days,g in spx.groupby('expiry_days',sort=True):
            q=g.reset_index(drop=True);batches=[]
            if days in OLD_DAYS:
                donors=oldspx[oldspx.expiry_days==days].reset_index(drop=True);ix=nearest_same_side(q,donors)
                with np.load(PARENT/f'controls_{oldrole}_t0_d{days}.npz') as a:
                    c,w,beta,stock=a['controls'][:,ix].copy(),a['weights'].copy(),a['beta'][ix].copy(),a['anchor_stock'].copy()
                grid=master[:np.searchsorted(master,q.T_years.iloc[0])+1]
                for batch in range(4):
                    stem=oldscratch/f'{oldrole}_t0_d{days}_b{batch}'
                    batches.append(CurvePathBatch(np.load(str(stem)+'_state2_ratio.npy',mmap_mode='r'),np.load(str(stem)+'_equity.npy',mmap_mode='r'),
                        np.diff(grid),np.minimum(np.searchsorted(ends,grid[:-1],side='right'),2)))
            else:
                folder=out/f'new_spx_{role}_d{days}'
                grid=np.load(folder/'grid.npy')
                with np.load(folder/'controls.npz') as a:c,w,beta,stock=[a[k].copy() for k in ['controls','weights','beta','anchor_stock']]
                for batch in range(4):
                    stem=scratch/f'{role}_d{days}_b{batch}'
                    batches.append(CurvePathBatch(np.load(str(stem)+'_ratio.npy',mmap_mode='r'),np.load(str(stem)+'_equity.npy',mmap_mode='r'),
                        np.diff(grid),np.minimum(np.searchsorted(ends,grid[:-1],side='right'),2)))
            d=SPXData(q,batches,w,c,stock);d.beta=beta;ds[int(days)]=d
        for days,g in vix.groupby('expiry_days',sort=True):
            q=g.reset_index(drop=True);donors=oldvix[oldvix.expiry_days==days].reset_index(drop=True);ix=nearest_same_side(q,donors)
            with np.load(PARENT/f'vix_state2_{oldrole}_t0_d{days}.npz') as a:
                dv[int(days)]=VIXData(q,a['exposures'].copy(),float(a['floor_integral']),float(a['delta']),
                    a['weights'].copy(),a['controls'].copy(),a['beta'][:,np.r_[0,1+ix]].copy(),8)
        return AllMaturityData(ds,dv,spx,vix,forwards)
    def task(days,role):
        log=out/f'worker_{role}_d{days}.log'
        with log.open('w') as stream:
            subprocess.run([sys.executable,'-m','matteo.model_comparison.option_pricing.all_maturity_fit',
                '--output',str(out),'--worker-days',str(days),'--worker-role',role],stdout=stream,stderr=subprocess.STDOUT,check=True)
    archive={};records=[];evaluation=[]
    def record(data,name,xi,role):
        r=data.predict(xi,False);q,f=data.loss.tables(r['prices'],r['errors'],r['forwards'],r['forward_errors'])
        q.to_csv(out/f'prices_{role}_{name}.csv',index=False);f.to_csv(out/f'forwards_{role}_{name}.csv',index=False)
        independent,components,_=score(data.all_quotes,forwards,r['prices'],r['forwards'],xi)
        np.testing.assert_allclose(independent,r['total'],rtol=1e-12)
        components['state_penalty']=design['state_penalty']
        records.append(dict(role=role,candidate=name,total=independent+design['state_penalty'],components=components,
            inside_bidask=int(q.interval_residual.eq(0).sum()),negatives=int(q.negative_estimate.sum()),nominal_se_passes=int(q.nominal_se_screen.sum())))
        for k,a in r['archive'].items():archive[f'{role}_{name}_{k}']=a
        archive[f'{role}_{name}_blocks']=r['blocks'];archive[f'{role}_{name}_forward_blocks']=r['forward_blocks']
        status('Comparison',**records[-1]);return r
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures={(role,days):pool.submit(task,days,role) for role in ['search','check'] for days in NEW_DAYS}
        for days in NEW_DAYS:futures[('search',days)].result()
        search=load('search');baseprediction=record(search,'baseline',baseline,'search')
        # All previously priced common training contracts must replay exactly.
        prior=pd.read_csv(PARENT/'prices_search_t0_selected.csv');new=pd.read_csv(out/'prices_search_baseline.csv')
        merged=prior.merge(new,on=['instrument','expiry','strike','option_type'],suffixes=('_old','_new'))
        assert len(merged)==882
        replay=float(np.max(np.abs(merged.model_price_U_old-merged.model_price_U_new)))
        assert replay<1e-8
        dump(out/'prior_price_replay.json',dict(common_contracts=len(merged),maximum_absolute_change=replay))
        best=[baseprediction['total'],baseline.copy()];ledger=[];optstart=time.perf_counter()
        def evaluate(x):
            if time.perf_counter()-optstart>1800:raise TimeoutError('Declared optimization budget exhausted')
            r=search.predict(x);evaluation.append(dict(levels=x.tolist(),total=r['total'],negatives=r['negative_count']))
            if r['total']<best[0]:best[:]=[r['total'],x.copy()]
            return r
        for start in [baseline,np.array(old['anchor_levels']),np.full(3,.06)]:
            cache={}
            def current(x):
                key=tuple(x)
                if key not in cache:cache.clear();cache[key]=evaluate(x)
                return cache[key]
            try:
                fit=least_squares(lambda x:current(x)['residuals'],start,jac=lambda x:current(x)['jacobian'],
                    bounds=(np.full(3,.00010001),np.full(3,.25)),max_nfev=30,ftol=1e-7,xtol=1e-7,gtol=1e-7)
                ledger.append(dict(start=start.tolist(),levels=fit.x.tolist(),nfev=fit.nfev,success=fit.success,message=fit.message))
            except TimeoutError as e:ledger.append(dict(start=start.tolist(),status='budget exhausted',message=str(e)));break
        selected=best[1];dump(out/'selected.json',dict(levels=selected.tolist(),selection='search numerical ensemble, all market quotes in sample',
            checking_not_used=True,optimizer_seconds=time.perf_counter()-optstart))
        dump(out/'evaluations.json',evaluation);dump(out/'optimizer_ledger.json',ledger)
        record(search,'selected',selected,'search');del search
        for days in NEW_DAYS:futures[('check',days)].result()
        check=load('check');record(check,'baseline',baseline,'check');record(check,'selected',selected,'check')
    np.savez_compressed(out/'comparison_paths.npz',**archive)
    check_sources(captured)
    dump(out/'run.json',design|dict(status='completed in-sample all-maturity curve fit and separate numerical ensemble check',
        selected_levels=selected.tolist(),scores=records,optimizer_ledger=ledger,seconds=time.perf_counter()-begin))
    status('Completed')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seed',type=int,default=2026102307);parser.add_argument('--worker-days',type=int)
    parser.add_argument('--worker-role',choices=['search','check'])
    args=parser.parse_args()
    try:
        if args.worker_days:worker(args.output,args.worker_days,args.worker_role)
        else:run(args)
    except Exception as e:
        if args.worker_days:
            dump(args.output/f'worker_failure_{args.worker_role}_d{args.worker_days}.json',dict(type=type(e).__name__,message=str(e)))
        elif args.output.exists():dump(args.output/'failure.json',dict(type=type(e).__name__,message=str(e)))
        raise
