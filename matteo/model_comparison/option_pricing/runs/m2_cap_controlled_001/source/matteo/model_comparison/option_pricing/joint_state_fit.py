"""Seven restricted initial states, joint curve fitting, three NEW fresh checks.

Reuses old SEARCHING primitive streams and frozen anchor controls; no previous
checking stream becomes fresh evidence. Every changed state gets new G/M0 and
conditional VIX exposures. The seven-point state grid is NOT continuous state
optimization or a global optimum claim. Readout and risk prices remain fixed.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import tempfile
import time
import shutil

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .m2_spx import ROOT, M2
from .spx_conditional_control_study import setup
from .joint_curve_fit import SPXData, VIXData, FitData, dump
from .joint_curve_core import CurvePathBatch, donor_indices
from .joint_state_core import AdaptiveDatedNormalizer, StateDatedPricer, projected_initial_shift
from .spx_volatility_sampling import proposal
from .provenance import capture_sources, source_hashes, save_sources, check_sources
from esn_eq.volatility_path_importance import volatility_path_basis, prepare_volatility_paths
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.variance_reduction import deterministic_stock_controls
from esn_eq.conditional_spx import lognormal_payoffs
from esn_eq.vix_pricing import ConditionalVIXPricer
from esn_eq.vix_curve_profile import prepare_curve_exposures
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.gaussian_controls import gaussian_control_rows

HERE = Path(__file__).resolve().parent
OLD = HERE/'runs/m2_joint_curve_fit_002'
VIX_PARENT = HERE/'runs/m2_vix_gaussian_stability_001'


def run(args):
    resume=Path(args.resume_search).resolve() if args.resume_search else None
    args.output.mkdir(parents=True,exist_ok=False)
    begin=time.perf_counter()
    targets,p,parameters,zero,anchor_curve,panel,master,exact_zero_normalizer,hashes=setup()
    old=json.loads((OLD/'run.json').read_text())
    oldscratch=Path(old['scratch'])
    baseline=np.array(old['selected_levels']); anchor=np.array(old['anchor_levels'])
    lo=np.array(old['family']['lower']); hi=np.array(old['family']['upper'])
    spx=pd.read_csv(OLD/'spx_contracts.csv');vix=pd.read_csv(OLD/'vix_contracts.csv')
    membership=pd.read_csv(OLD/'spx_membership.csv');primary=membership.priority_selected.to_numpy(bool)
    forwards=targets.vix_forwards.sort_values('expiry').reset_index(drop=True)
    assert primary.sum()==707 and len(vix)==175 and list(forwards.expiry_days)==[27,55,118,153]
    np.testing.assert_array_equal(spx[['expiry','strike','option_type']],membership[['expiry','strike','option_type']])
    cov=p.projection@p.stationary@p.projection.T;L=np.linalg.cholesky(cov)
    zgrid=[np.zeros(3)]+[sign*np.eye(3)[j] for j in range(3) for sign in [-1,1]]
    states=[]
    for j,z in enumerate(zgrid):
        u=L@z;initial,norm=p.initial_state(u)
        np.testing.assert_allclose(norm,z@z,atol=1e-12)
        states.append(dict(index=j,whitened=z.tolist(),projections=u.tolist(),initial=initial.tolist(),prior_squared_norm=norm))
    scratch=Path(tempfile.mkdtemp(prefix='m2-joint-state-',dir='/private/tmp'))
    sources=[Path(__file__),HERE/'joint_state_core.py',HERE/'joint_curve_fit.py',HERE/'joint_curve_core.py',
             HERE/'spx_conditional_control_study.py',HERE/'spx_volatility_sampling.py',HERE/'active_scope.json',
             OLD/'run.json',OLD/'selected.json',OLD/'spx_membership.csv',OLD/'frozen_spx_controls.json',
             *[OLD/f'{kind}_r1_t0_d{days}.npz' for kind,days_list in [('spx',[29,57,211,302]),('vix',[27,55,118,153])] for days in days_list],
             VIX_PARENT/'proposals.json',VIX_PARENT/'coefficients.json',*sorted((M2/'src/esn_eq').glob('*.py'))]
    if resume:
        previous_design=json.load(open(resume/'design.json'))
        for key,value in [('architecture',asdict(p.architecture)),('parameters',asdict(parameters)),('premium',asdict(p.premium)),('target_sha256',targets.artifact_sha256)]:
            assert previous_design[key]==json.loads(json.dumps(value)),key
        sources += [resume/'design.json',resume/'selected.json',resume/'state_profiles.json',resume/'optimizer_ledger.json',resume/'failure.json']
    captured=capture_sources(sources);save_sources(args.output,captured,ROOT)
    design=dict(states=states,baseline_levels=baseline.tolist(),anchor_levels=anchor.tolist(),
        architecture=asdict(p.architecture),parameters=asdict(parameters),premium=asdict(p.premium),
        target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,implementation_sha256=source_hashes(captured,ROOT),
        budgets=dict(state_nodes=7,curve_starts=2,max_nfev_per_start=args.max_nfev,search_wall_seconds=3600,
                     spx_paths=512,vix_nets=8,vix_points_per_net=64,fresh_trials=3,stock_steps_per_trading_day=512,
                     conditional_GH=32,normalizer_GH=48,time_order=8,normalizer_interpolation_relative_screen=2e-7),
        state_penalty='0.1*minimum_stationary_prior_squared_norm/3',
        primitive_reuse='Old searching draws only. OU paths translate by exp(-rates*t)*X0; recompute G, dated M0, stock and VIX exposures. Fixed zero-state controls/likelihoods remain fixed functions of the SAME primitive Gaussian law.',
        normalizer='NEW log-GH48 CubicSpline on log1p(time/tau); adaptive midpoint and independent-node screens. Separate approximation error, not GH convergence.',
        objective='707 cutoff-3 SPX +175 VIX training options +4 VIX forwards; same frozen scales/weights/curve penalty; signed estimates retained.',
        scope='SPX fit 29/57/211/302 only. 119/666/1030 not evaluated; 1/1401/1765/2129 deferred. All four VIX maturities.',
        seed=args.seed,scratch=str(scratch),old_scratch=str(oldscratch),
        limitations=['Seven-point state grid is restricted search, not continuous/global calibration.',
                     'Exploratory low-path budgets; unresolved sampling, stock timestep and GH/time bias.',
                     'New normalizer interpolation screening is not a rigorous uniform bound.',
                     'Fresh draws check the same training quotes, not held-out market performance.'])
    if resume:
        design.update(resumed_search_from=str(resume),search_repeated=False,
                      checking_seed=args.seed,checking_reason='Fresh checking only after retained pre-score last-chunk broadcast failure. All search selections frozen and unchanged.')
    dump(args.output/'design.json',design)
    membership.to_csv(args.output/'spx_membership.csv',index=False);vix.to_csv(args.output/'vix_contracts.csv',index=False)
    forwards.to_csv(args.output/'forward_targets.csv',index=False)
    seeds=[];normalizer_audit=[];costs=[];basis={};grids={}
    for days in [29,57,211,302]:
        q=spx[spx.expiry_days==days];grid=master[:np.searchsorted(master,q.T_years.iloc[0])+1]
        grids[days]=grid;basis[days]=volatility_path_basis(p,grid,7/365)
        values=exact_zero_normalizer[:len(grid)-1];values.setflags(write=False)
        p.normalizers[(grid[:-1].tobytes(),zero.tobytes(),parameters,48)]=values

    def status(stage,**details):
        dump(args.output/'status.json',dict(stage=stage,elapsed_seconds=time.perf_counter()-begin,**details))
        print(stage,json.dumps(details),flush=True)

    def primitives(role,trial):
        result={}
        for days,grid in grids.items():
            q=spx[spx.expiry_days==days].reset_index(drop=True)
            reps=panel[panel.expiry_days==days].reset_index(drop=True);donor=donor_indices(q,reps)
            if role=='search':
                saved=np.load(OLD/f'spx_r1_t0_d{days}.npz')
                controls,weights,beta,stocks=saved['controls'],saved['weights'],saved['beta'],saved['anchor_stock']
            else:
                saved=np.load(OLD/f'spx_r1_t0_d{days}.npz');beta=saved['beta']
                controls=[];weights=[];stocks=[]
            batches=[]
            for batch in range(4):
                if role=='search':
                    seed=next(r['seed'] for r in old['seed_records'] if r['instrument']=='SPX' and r['role']==1 and r['trial']==0 and r['days']==days and r['batch']==batch)
                else:
                    seed=int(np.random.SeedSequence([args.seed,100,trial,days,batch]).generate_state(1)[0])
                tick=time.perf_counter();r=prepare_volatility_paths(p,grid,zero,128,seed,basis[days],proposal(1))
                cache=r.cache;stem=scratch/f'{role}_t{trial}_d{days}_b{batch}'
                np.save(str(stem)+'_projections.npy',cache.projections)
                np.save(str(stem)+'_equity.npy',cache.equity_increments)
                if role=='search':
                    np.testing.assert_array_equal(r.weights,weights[batch*128:(batch+1)*128])
                    oldstem=oldscratch/f'r1_t0_d{days}_b{batch}'
                    np.testing.assert_array_equal(cache.equity_increments,np.load(str(oldstem)+'_equity.npy',mmap_mode='r'))
                    ix=np.unique(np.linspace(0,len(grid)-2,37).astype(int))
                    np.testing.assert_allclose(p.multiplier(cache.projections[ix],parameters)/exact_zero_normalizer[ix,None],
                        np.load(str(oldstem)+'_ratio.npy',mmap_mode='r')[ix],rtol=2e-12,atol=2e-12)
                    ratiofile=Path(str(oldstem)+'_ratio.npy')
                else:
                    h=spx_hedging_controls(p,cache,parameters,anchor_curve,reps.T_years.to_numpy(),reps.K_over_F.to_numpy(),reps.option_type.to_numpy())
                    stock=h.terminal_ratios[float(q.T_years.iloc[0])]
                    ref=deterministic_stock_controls(cache,anchor_curve,[float(q.T_years.iloc[0])])[float(q.T_years.iloc[0])]
                    total=anchor_curve.integral(0.,float(q.T_years.iloc[0]))
                    known=lognormal_payoffs(-.5*total,total,q.K_over_F.to_numpy(),q.option_type.to_numpy())
                    sign=np.where(q.option_type.to_numpy()=='call',1.,-1.)
                    c=np.concatenate([np.broadcast_to((stock-1)[:,None,None],(128,len(q),1)),
                        (np.maximum(sign*(ref[:,None]-q.K_over_F.to_numpy()),0.)-known)[:,:,None],h.hedges[:,donor]],axis=2)
                    c=np.concatenate([r.weights[:,None,None]*c,np.broadcast_to((r.weights-1)[:,None,None],(128,len(q),1))],axis=2)
                    controls.append(c);weights.append(r.weights);stocks.append(stock)
                    ratiofile=Path(str(stem)+'_zero_ratio.npy')
                    out=np.lib.format.open_memmap(ratiofile,mode='w+',dtype='float64',shape=cache.equity_increments.shape)
                    local_normalizer=exact_zero_normalizer[:len(grid)-1]
                    for k in range(0,len(out),2048):out[k:k+2048]=p.multiplier(cache.projections[k:k+2048],parameters)/local_normalizer[k:k+2048,None]
                    out.flush();del out
                batches.append(dict(stem=stem,zero_ratio=ratiofile))
                seeds.append(dict(instrument='SPX',role=role,trial=trial,days=days,batch=batch,seed=seed,reused_old_search=role=='search'))
                costs.append(dict(stage='primitive_paths',role=role,trial=trial,days=days,seconds=time.perf_counter()-tick))
                status('Primitive SPX paths',role=role,trial=trial,days=days,batch=batch)
                del cache,r
            if role!='search':controls=np.concatenate(controls);weights=np.concatenate(weights);stocks=np.concatenate(stocks)
            np.savez_compressed(args.output/f'controls_{role}_t{trial}_d{days}.npz',controls=controls,weights=weights,beta=beta,anchor_stock=stocks)
            result[days]=dict(quotes=q,batches=batches,controls=controls,weights=weights,beta=beta,stocks=stocks)
        return result

    normalizers={}
    def candidate_normalizer(state):
        j=state['index']
        if j not in normalizers:
            tick=time.perf_counter();g=AdaptiveDatedNormalizer(p,parameters,np.array(state['initial']),float(master[-1]))
            normalizers[j]=g;normalizer_audit.append(dict(state=j,**g.audit,seconds=time.perf_counter()-tick))
            if j==0:
                errors=np.abs(g(master[:-1])/exact_zero_normalizer-1)
                normalizer_audit[-1]['zero_full_stock_grid_max_relative_error']=float(errors.max())
                assert errors.max()<2e-7
            np.savez_compressed(args.output/f'normalizer_state{j}.npz',knots=g.knots,log_values=g.logvalues,initial=g.initial)
            dump(args.output/'normalizer_audit.json',normalizer_audit)
        return normalizers[j]

    def spx_for(state,base,role,trial):
        result={};j=state['index'];initial=np.array(state['initial'])
        g=candidate_normalizer(state)
        for days,entry in base.items():
            grid=grids[days];ratio=[]
            delta=projected_initial_shift(p,grid[:-1],initial)
            m=g(grid[:-1]) if j else exact_zero_normalizer[:len(grid)-1]
            for batch,item in enumerate(entry['batches']):
                stem=item['stem'];equity=np.load(str(stem)+'_equity.npy',mmap_mode='r')
                if j==0:ratiofile=item['zero_ratio']
                else:
                    ratiofile=Path(str(stem)+f'_state{j}_ratio.npy')
                    projections=np.load(str(stem)+'_projections.npy',mmap_mode='r')
                    out=np.lib.format.open_memmap(ratiofile,mode='w+',dtype='float64',shape=equity.shape)
                    for k in range(0,len(out),2048):
                        out[k:k+2048]=p.multiplier(projections[k:k+2048]+delta[k:k+2048,None,:],parameters)/m[k:k+2048,None]
                    out.flush();del out,projections
                ratio.append(CurvePathBatch(np.load(ratiofile,mmap_mode='r'),equity,np.diff(grid),
                    np.minimum(np.searchsorted(anchor_curve.ends,grid[:-1],side='right'),2)))
            d=SPXData(entry['quotes'],ratio,entry['weights'],entry['controls'],entry['stocks']);d.beta=entry['beta'];result[days]=d
        return result

    proposals=json.loads((VIX_PARENT/'proposals.json').read_text());coeff=json.loads((VIX_PARENT/'coefficients.json').read_text())
    def vix_for(state,role,trial):
        result={};j=state['index'];initial=np.array(state['initial'])
        if j:
            dated=StateDatedPricer(p.architecture,p.premium);dated.state_normalizer=candidate_normalizer(state)
        else:dated=p
        engine=ConditionalVIXPricer(dated,parameters,anchor_curve,initial,normalizer_order=48)
        for days in [27,55,118,153]:
            q=vix[vix.expiry_days==days].reset_index(drop=True);T=float(q.T_years.iloc[0])
            if resume and role=='search':
                path=resume/f'vix_state{j}_search_t0_d{days}.npz'
                with np.load(path) as a:
                    result[days]=VIXData(q,a['exposures'].copy(),float(a['floor_integral']),float(a['delta']),
                        a['weights'].copy(),a['controls'].copy(),a['beta'].copy(),8)
                shutil.copyfile(path,args.output/path.name)
                continue
            if role=='search':
                saved=np.load(OLD/f'vix_r1_t0_d{days}.npz');x=saved['states'].copy();w=saved['weights'];c=saved['controls'];beta=saved['beta']
            else:
                mix=DefensiveGaussianMixture(proposals[str(days)]['shifts'],proposals[str(days)]['probabilities'])
                sampler=GaussianOuterSampler(p,T,zero);U=np.asarray(coeff[str(days)]['basis'])
                xx=[];ww=[];cc=[]
                for rep in range(8):
                    seed=int(np.random.SeedSequence([args.seed,200,trial,days,rep]).generate_state(1)[0])
                    z,weight,_=mix.draw(6,seed);xx.append(sampler.mean+z@sampler.root.T);ww.append(weight);cc.append(gaussian_control_rows(z,weight,U))
                    if j==0 and not any(r['seed']==seed for r in seeds):
                        seeds.append(dict(instrument='VIX',role=role,trial=trial,days=days,rep=rep,seed=seed,reused_old_search=False))
                x=np.concatenate(xx);w=np.concatenate(ww);c=np.concatenate(cc)
                beta=np.load(OLD/f'vix_r1_t0_d{days}.npz')['beta']
            if j:x=x+np.exp(-p.layout.rates*T)*initial
            if role=='search' and j==0:
                exposure=saved['exposures'];floor=float(saved['floor_integral']);delta=float(saved['delta'])
            else:
                tick=time.perf_counter();status('Candidate VIX exposures',state=j,role=role,trial=trial,days=days)
                window=engine.gaussian_window(T,time_order=8)
                cache=prepare_curve_exposures(engine,x,window,order=32)
                exposure,floor,delta=cache.coefficients,cache.floor_integral,cache.delta
                costs.append(dict(stage='VIX_exposures',state=j,role=role,trial=trial,days=days,seconds=time.perf_counter()-tick))
            d=VIXData(q,exposure,floor,delta,w,c,beta,8);result[days]=d
            np.savez_compressed(args.output/f'vix_state{j}_{role}_t{trial}_d{days}.npz',states=x,exposures=exposure,floor_integral=floor,delta=delta,weights=w,controls=c,beta=beta)
        return result

    def data_for(state,base,role,trial):
        status('Build candidate',state=state['index'],role=role,trial=trial)
        return FitData(spx_for(state,base,role,trial),vix_for(state,role,trial),spx,vix,forwards,primary)

    search_begin=time.perf_counter()
    if resume:
        base={};prior_scratch=Path(previous_design['scratch'])
        for days in [29,57,211,302]:
            q=spx[spx.expiry_days==days].reset_index(drop=True)
            with np.load(resume/f'controls_search_t0_d{days}.npz') as a:
                controls,weights,beta,stocks=[a[k].copy() for k in ['controls','weights','beta','anchor_stock']]
            batches=[]
            for batch in range(4):
                source_stem=prior_scratch/f'search_t0_d{days}_b{batch}'
                stem=scratch/f'search_t0_d{days}_b{batch}'
                for suffix in ['_projections.npy','_equity.npy']:
                    Path(str(stem)+suffix).symlink_to(Path(str(source_stem)+suffix))
                batches.append(dict(stem=stem,zero_ratio=Path(str(oldscratch/f'r1_t0_d{days}_b{batch}')+'_ratio.npy')))
            base[days]=dict(quotes=q,batches=batches,controls=controls,weights=weights,beta=beta,stocks=stocks)
            shutil.copyfile(resume/f'controls_search_t0_d{days}.npz',args.output/f'controls_search_t0_d{days}.npz')
        for item in old['seed_records']:
            if item['role']==1:seeds.append(item|dict(role='search',reused_old_search=True))
        for filename in ['state_profiles.json','optimizer_ledger.json','evaluations.jsonl']:
            shutil.copyfile(resume/filename,args.output/filename)
    else:base=primitives('search',0)
    ledger=[];evaluations=[];candidates=[];search_records={}
    for state in ([] if resume else states):
        if time.perf_counter()-search_begin>3600:break
        data=data_for(state,base,'search',0);memo={};penalty=.1*state['prior_squared_norm']/3
        def evaluate(x):
            key=np.asarray(x).tobytes()
            if key not in memo:
                if time.perf_counter()-search_begin>3600:raise TimeoutError('Declared search wall budget exhausted')
                r=data.predict(x);memo.clear();memo[key]=r
                row=dict(evaluation=len(evaluations),state=state['index'],levels=np.asarray(x).tolist(),total=r['total']+penalty,
                         components=r['components']|dict(state_penalty=penalty),negative_count=r['negative_count'])
                evaluations.append(row)
                with (args.output/'evaluations.jsonl').open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n')
            return memo[key]
        best=(float('inf'),None)
        for number,start in enumerate([baseline,anchor]):
            try:
                fit=least_squares(lambda x:evaluate(x)['residuals'],start,jac=lambda x:evaluate(x)['jacobian'],bounds=(lo,hi),
                    max_nfev=args.max_nfev,ftol=1e-6,xtol=1e-6,gtol=1e-5)
                r=evaluate(fit.x);score=r['total']+penalty
                ledger.append(dict(state=state['index'],start=number,levels=fit.x.tolist(),total=score,success=bool(fit.success),nfev=fit.nfev,message=fit.message))
                if score<best[0]:best=(score,fit.x.copy())
            except TimeoutError as e:
                ledger.append(dict(state=state['index'],start=number,status='interrupted',message=str(e)));break
        if best[1] is not None:
            candidates.append((best[0],state,best[1]));search_records[state['index']]=dict(state=state,levels=best[1].tolist(),total=best[0])
            status('Completed state profile',state=state['index'],total=best[0],levels=best[1].tolist())
        dump(args.output/'optimizer_ledger.json',ledger);dump(args.output/'state_profiles.json',list(search_records.values()))
        del data
    if resume:
        selected_record=json.load(open(resume/'selected.json'))
        selected_state=selected_record['state'];selected=np.array(selected_record['levels'])
        ledger=json.load(open(resume/'optimizer_ledger.json'))
        selected_record=selected_record|dict(resumed_search_from=str(resume),selection_repeated=False)
    else:
        if not candidates:raise RuntimeError('No completed state profile')
        _,selected_state,selected=min(candidates,key=lambda x:x[0])
        selected_record=dict(state=selected_state,levels=selected.tolist(),selection='search only',fresh_checking_not_used=True,
                             completed_state_profiles=len(candidates),search_seconds=time.perf_counter()-search_begin)
    dump(args.output/'selected.json',selected_record)
    records=[];blocks={};paths={}
    def record(name,state,levels,data,role,trial):
        r=data.predict(levels,derivative=False);penalty=.1*state['prior_squared_norm']/3
        q,f=data.loss.tables(r['prices'][data.primary],r['errors'][data.primary],r['forwards'],r['forward_errors'])
        q.to_csv(args.output/f'prices_{role}_t{trial}_{name}.csv',index=False);f.to_csv(args.output/f'forwards_{role}_t{trial}_{name}.csv',index=False)
        full=data.all_quotes.copy();full['model_price_U']=r['prices'];full['mc_se_U']=r['errors'];full['priority_objective_member']=data.primary
        full.to_csv(args.output/f'full_training_{role}_t{trial}_{name}.csv',index=False)
        blocks[f'{role}_t{trial}_{name}_prices']=r['blocks'];blocks[f'{role}_t{trial}_{name}_forwards']=r['forward_blocks']
        for key,array in r['archive'].items():paths[f'{role}_t{trial}_{name}_{key}']=array
        records.append(dict(candidate=name,role=role,trial=trial,state_index=state['index'],levels=levels.tolist(),
                            total=r['total']+penalty,components=r['components']|dict(state_penalty=penalty),negative_count=r['negative_count']))
        status('Comparison score',candidate=name,role=role,trial=trial,total=r['total']+penalty)
    for name,state,levels in [('baseline',states[0],baseline),('selected',selected_state,selected)]:
        record(name,state,levels,data_for(state,base,'search',0),'search',0)
    # Selection is frozen BEFORE any new checking stream is generated.
    for trial in range(3):
        basefresh=primitives('fresh',trial)
        for name,state,levels in [('baseline',states[0],baseline),('selected',selected_state,selected)]:
            record(name,state,levels,data_for(state,basefresh,'fresh',trial),'fresh',trial)
        del basefresh
    check_sources(captured)
    np.savez_compressed(args.output/'comparison_blocks.npz',**blocks)
    np.savez_compressed(args.output/'comparison_paths.npz',**paths)
    pd.DataFrame(costs).to_csv(args.output/'costs.csv',index=False)
    dump(args.output/'scores.json',records)
    freshseeds=[r['seed'] for r in seeds if not r['reused_old_search']]
    assert len(freshseeds)==len(set(freshseeds))
    assert not set(freshseeds)&{r['seed'] for r in old['seed_records']}
    result=design|dict(status='completed restricted seven-state joint calibration and three fresh numerical checks',
                       selected=selected_record,optimizer_ledger=ledger,seed_records=seeds,scores=records,
                       normalizer_audit=normalizer_audit,wall_seconds=time.perf_counter()-begin)
    dump(args.output/'run.json',result);status('completed',selected_state=selected_state['index'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--max-nfev',type=int,default=25)
    parser.add_argument('--seed',type=int,default=2026102401)
    parser.add_argument('--resume-search',type=Path)
    args=parser.parse_args()
    try:run(args)
    except Exception as e:
        if args.output.exists():dump(args.output/'failure.json',dict(type=type(e).__name__,message=str(e)))
        raise
