"""Forward-constrained VIX fits: fixed floor versus independently fitted floor."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
from scipy.optimize import least_squares,minimize
from .forward_floor_core import ForwardFloorFit,pack,unpack
from .mixed_response_core import MixedFit
from .mixed_response_study import ZSHIFTS,ZPROB,seed_for
from .vix_only_curve_fit import PRIOR,ARCHIVE,ACTIVE,DAYS,load_draws,paired
from .full_readout_fit import HERE,dump
from .m2_spx import ROOT,M2,frozen_stage3
from .postcap_vov_core import PostKernel,exposures
from .cap_controlled_core import shape_metrics
from .controlled_readout_fast_fit import ControlledEngine
from .plot_joint_curve_fit import implied_vol
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture

PARENT=HERE/'runs/m2_mixed_response_001'
CHECK_SEED=20261026083


def metadata(x,p):
    levels,floor,theta=unpack(x,p)
    return dict(b=float(theta[0]),weight=float(theta[1]),half_life_days=float(np.exp(theta[2])),floor_variance=float(floor),floor_volatility=float(np.sqrt(floor)),active_levels=levels[ACTIVE].tolist())


def record(e,x,folder,name):
    r=e.predict(x,True);folder.mkdir(parents=True,exist_ok=True)
    native=np.concatenate([d['native'] for d in r['raw'].values()],axis=1);call=np.concatenate([d['call'] for d in r['raw'].values()],axis=1);fw=np.column_stack([d['forward'] for d in r['raw'].values()])
    bn=native.reshape(8,-1,native.shape[1]).mean(1);bc=call.reshape(8,-1,call.shape[1]).mean(1);bf=fw.reshape(8,-1,4).mean(1)
    for typ,q,val,block in [('otm',e.market,r['prices'],bn),('canonical',e.calls,r['calls'],bc)]:
        q=q.copy();q['model_price_U']=val;q['mc_se_U']=block.std(0,ddof=1)/np.sqrt(8);q['interval_residual']=np.where(val<q.price_bid_U,val-q.price_bid_U,np.where(val>q.price_ask_U,val-q.price_ask_U,0.))/q.price_scale_U
        q.to_csv(folder/f'prices_{typ}_{name}.csv',index=False)
    f=e.forwards.copy();f['model_forward']=r['forwards'];f['mc_se']=bf.std(0,ddof=1)/np.sqrt(8);f['forward_error']=r['forwards']-f.forward_vix_points.to_numpy();f.to_csv(folder/f'forwards_{name}.csv',index=False)
    np.savez_compressed(folder/f'blocks_{name}.npz',native=bn,call=bc,forward=bf)
    np.savez_compressed(folder/f'paths_{name}.npz',**{f'd{day}_{k}':v for day,d in r['raw'].items() for k,v in d.items()})
    q=e.market.copy();q['model_price_U']=r['prices'];q['model_forward']=q.expiry_days.map(dict(zip(f.expiry_days,f.model_forward)))
    q['market_native_IV']=[implied_vol(a.price_mid_U,a.forward,a.strike,a.T_years,a.option_type) for a in q.itertuples()]
    q['model_native_IV']=[implied_vol(a.model_price_U,a.model_forward,a.strike,a.T_years,a.option_type) for a in q.itertuples()]
    q.to_csv(folder/f'native_iv_{name}.csv',index=False)
    puts=q[q.option_type.eq('put')].sort_values(['expiry_days','strike']);low=puts.groupby('expiry_days').head(1).copy();low['price_error']=low.model_price_U-low.price_mid_U;low.to_csv(folder/f'lowest_puts_{name}.csv',index=False)
    per=[]
    for day,g in puts.groupby('expiry_days'):
        per.append(dict(day=int(day),put_price_MAE=float(np.mean(abs(g.model_price_U-g.price_mid_U))),put_IV_MAE=float(np.nanmean(abs(g.model_native_IV-g.market_native_IV))),lowest_strike=float(g.strike.min())))
    summary=dict(name=name,parameters=x.tolist(),model=metadata(x,e.p),loss=r['loss'],option_loss=r['option_loss'],inside_bidask=int(np.sum((r['prices']>=e.market.price_bid_U)&(r['prices']<=e.market.price_ask_U))),negative_native_prices=int(np.sum(r['prices']<0)),maximum_forward_error=float(abs(f.forward_error).max()),forward_errors=f.forward_error.to_list(),forward_SEs=f.mc_se.to_list(),shape=shape_metrics(e.calls,r['calls']),native_IV_unavailable=int(np.sum(~np.isfinite(q.model_native_IV))),put_diagnostics=per,moments=r['moments'])
    dump(folder/f'score_{name}.json',summary);return summary


def fresh_inputs(out,dated,initial):
    cache=out/'fresh_cache';prod=cache/'production';prod.mkdir(parents=True);seeds=[];proposals=json.loads((ARCHIVE/'m2_vix_targeted_check_001/proposals.json').read_text())
    for day in DAYS:
        dec=proposals[str(day)];shifts=np.array([np.r_[old,z] for old in dec['shifts'] for z in ZSHIFTS]);prob=np.array([p*q for p in dec['probabilities'] for q in ZPROB]);mix=DefensiveGaussianMixture(shifts,prob);sampler=GaussianOuterSampler(dated,day/365,initial);zz=[];ww=[];cc=[]
        for net in range(8):
            seed=seed_for(CHECK_SEED,day,net);z,w,c=mix.draw(12,seed);zz.append(z);ww.append(w);cc.append(c);seeds.append(dict(day=day,net=net,seed=seed))
        z=np.concatenate(zz);np.savez_compressed(prod/f'vix_d{day}.npz',states=sampler.mean+z[:,:-1]@sampler.root.T,normals=z[:,:-1],z_normals=z[:,-1],weights=np.concatenate(ww),components=np.concatenate(cc),nets=np.repeat(np.arange(8),4096))
    dump(cache/'design.json',dict(seed_base=CHECK_SEED,seeds=seeds,points_per_expiry=32768,old_proposals=proposals,Z_shifts=ZSHIFTS.tolist(),Z_probabilities=ZPROB.tolist(),created_after_selection=True))
    return cache


def run(out):
    out.mkdir(parents=True,exist_ok=False);start=time.perf_counter();d=json.loads((PRIOR/'run.json').read_text());p=np.array(d['selected']);ends=np.array(d['ends_years']);parent=json.loads((PARENT/'selection.json').read_text())
    prior=parent['candidates']['additive_price']['parameters'];levels=p[7:-1].copy();levels[ACTIVE]=prior[:7];x0=pack(levels,prior[7],prior[8],np.exp(prior[9]),p[-1])
    market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv');forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv');controls=json.loads((PRIOR/'vix_controls.json').read_text());nodes=dict(np.load(PARENT/'fresh_conditional_nodes.npz'));draws=load_draws(PARENT/'fresh_cache')
    sources=[*sorted(HERE.glob('*.py')),HERE/'FORWARD_FLOOR_PLAN.md',HERE/'postcap_vov_kernel.cpp',HERE/'stochastic_vov_kernel.cpp',HERE/'cap_tail_kernel.cpp',PRIOR/'run.json',PRIOR/'selected_model.json',PRIOR/'vix_controls.json',PARENT/'selection.json',PARENT/'experimental_model_additive_price.json',ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv',ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv',ARCHIVE/'m2_vix_targeted_check_001/proposals.json',M2/'experiments/stage3_001/manifest.json',M2/'experiments/stage3_001/top_configurations.csv',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(sources);save_sources(out,captured,ROOT);inputs=[PARENT/'fresh_conditional_nodes.npz',*(PARENT/'fresh_cache/production'/f'vix_d{day}.npz' for day in DAYS)]
    design=dict(status='running',source_sha256=source_hashes(captured,ROOT),input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},parent=str(PARENT),training_note='Previous mixture fresh checking draws explicitly reused as training;16384points/expiry',curve_piece_count=16,active_curve_indices_zero_based=ACTIVE.tolist(),SPX_priced=False,baseline_replaced=False,floor_volatility_bounds=[.05,.12],forward_constraints='Four equality constraints on training estimator; selection feasibility <=1e-5VIXpoint',fresh_check_seed_base=CHECK_SEED)
    dump(out/'design.json',design);e=ForwardFloorFit(nodes,draws,market,forwards,controls,p,ends);replay=e.predict(x0);saved=pd.read_csv(PARENT/'fresh_check/prices_otm_additive_price.csv');error=float(np.max(abs(replay['prices']-saved.model_price_U.to_numpy())));assert error<1e-9,error;dump(out/'baseline_replay.json',dict(maximum_native_price_error=error,loss=replay['loss']))
    check=x0.copy();check[:7]=np.clip(check[:7],.05,.8);check[7:]=[.6,.65,.45,.35];rr=e.predict(check);checks=[]
    for j in range(11):
        dx=np.zeros(11);dx[j]=1e-6;rp=e.predict(check+dx);rm=e.predict(check-dx);fdp=(rp['prices']-rm['prices'])/2e-6;fdf=(rp['forwards']-rm['forwards'])/2e-6
        pe=float(np.linalg.norm(fdp-rr['native_jacobian'][:,j])/max(1,np.linalg.norm(fdp)));fe=float(np.linalg.norm(fdf-rr['forward_jacobian'][:,j])/max(1,np.linalg.norm(fdf)));assert pe<.003 and fe<1e-6,(j,pe,fe);checks.append(dict(parameter=j,native_price_relative_error=pe,forward_relative_error=fe))
    dump(out/'jacobian_checks.json',checks);print('REPLAY AND 11-DIRECTION JACOBIAN PASS',flush=True)
    target=forwards.forward_vix_points.to_numpy();ledger=[];evaluations=[];candidates={'previous_mixture':x0.copy()}
    def fit(label,initial,free_floor):
        initial=np.clip(initial,0.,1.);memo={};free=np.arange(11 if free_floor else 10)
        def forward_prefit(y):
            x=initial.copy();x[:7]=y;return e.predict(x)
        cf=least_squares(lambda y:forward_prefit(y)['forwards']-target,initial[:7],jac=lambda y:forward_prefit(y)['forward_jacobian'][:,:7],bounds=(np.zeros(7),np.ones(7)),max_nfev=50,ftol=1e-10,xtol=1e-10,gtol=1e-9)
        initial[:7]=cf.x;fixed=initial.copy()
        def current(xx):
            key=tuple(xx)
            if key not in memo:
                x=fixed.copy();x[free]=xx;r=e.predict(x);memo.clear();memo[key]=r;evaluations.append(dict(label=label,parameters=x.tolist(),option_loss=r['option_loss'],maximum_forward_error=float(np.max(abs(r['forwards']-target)))))
            return memo[key]
        def callback(xx):
            if len(evaluations)%20==0:
                r=current(xx);print('ITER',label,'option loss',r['option_loss'],'forward err',max(abs(r['forwards']-target)),flush=True)
        opt=minimize(lambda xx:current(xx)['option_loss']/10,initial[free],jac=lambda xx:(2*current(xx)['residuals'][:-4]@current(xx)['jacobian'][:-4,free])/10,method='SLSQP',bounds=[(0.,1.) for _ in free],constraints=[dict(type='eq',fun=lambda xx:current(xx)['forwards']-target,jac=lambda xx:current(xx)['forward_jacobian'][:,free])],callback=callback,options=dict(maxiter=120,ftol=1e-9))
        x=fixed.copy();x[free]=opt.x;r=e.predict(x);gap=float(np.max(abs(r['forwards']-target)));rec=dict(label=label,free_floor=free_floor,parameters=x.tolist(),model=metadata(x,p),option_loss=r['option_loss'],maximum_forward_error=gap,eligible=gap<=1e-5,success=bool(opt.success),message=opt.message,iterations=opt.nit,nfev=opt.nfev,forward_prefit_nfev=cf.nfev)
        ledger.append(rec);dump(out/'optimizer_ledger.json',ledger);dump(out/'evaluations.json',evaluations);print('FIT',label,'option loss',r['option_loss'],'floor vol',metadata(x,p)['floor_volatility'],'F gap',gap,'success',opt.success,flush=True);return x,rec
    fixed_starts=[x0.copy()]
    for b,w,half in [(4.,.25,365),(5.,.9,5),(4.,.8,15),(3.,.5,100),(5.,.95,30)]:
        x=x0.copy();x[7:10]=[b/5,w,np.log(half)/np.log(365.)];fixed_starts.append(x)
    for j,x in enumerate(fixed_starts):fit(f'fixed_start{j}',x,False)
    eligible=[r for r in ledger if not r['free_floor'] and r['eligible']]
    if not eligible:raise RuntimeError('No forward-feasible fixed-floor candidate; retain all diagnostics')
    fixed=np.array(min(eligible,key=lambda r:r['option_loss'])['parameters']);candidates['tight_forwards']=fixed
    free_starts=[]
    for sigma in [.05,np.sqrt(p[-1]),.065,.08,.10,.115]:
        x=fixed.copy();x[10]=(sigma-.05)/.07;free_starts.append(x)
    x=x0.copy();x[7:]=[.8,.25,1.,(.08-.05)/.07];free_starts.append(x)
    for j,x in enumerate(free_starts):fit(f'free_start{j}',x,True)
    eligible=[r for r in ledger if r['free_floor'] and r['eligible']];best=min(eligible,key=lambda r:r['option_loss']) if eligible else None
    candidates['tight_forwards_free_floor']=np.array(best['parameters']) if best and best['option_loss']<=e.predict(fixed)['option_loss'] else fixed.copy()
    dump(out/'selection.json',dict(parameters={name:x.tolist() for name,x in candidates.items()},frozen_before_fresh_check=True,baseline_replaced=False,market_out_of_sample=False))
    train=[record(e,x,out/'training',name) for name,x in candidates.items()];dump(out/'training/scores.json',train)
    table=pd.DataFrame(dict(start_day=np.r_[0,ends[:-1]]*365,end_day=ends*365,active_in_VIX_fit=np.isin(np.arange(16),ACTIVE)))
    for name,x in candidates.items():
        levels,floor,_=unpack(x,p);pp=p.copy();pp[7:-1]=levels;pp[-1]=floor;table[name]=levels;assert np.all(levels>=floor)
        dump(out/f'experimental_model_{name}.json',dict(parameters=pp.tolist(),ends_years=ends.tolist(),readout=metadata(x,p),architecture=d['architecture'],premium=d['premium'],initial_state=d['initial_state'],adopted=False))
    table.to_csv(out/'curve_levels.csv',index=False)
    a,_,premium,_=frozen_stage3();assert json.loads(json.dumps(asdict(a)))==d['architecture'];assert json.loads(json.dumps(asdict(premium)))==d['premium'];dated=DatedM2Pricer(a,premium);initial=np.array(d['initial_state']);cache=fresh_inputs(out,dated,initial)
    lib=out/'postcap.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'postcap_vov_kernel.cpp'),'-o',str(lib)],check=True);k=PostKernel(lib,ends)
    ctx=object.__new__(ControlledEngine);ctx.dated=dated;ctx.initial=initial;ctx.k=k;ctx.ends=ends;ctx.tau=1/max(dated.layout.rates);ctx.horizon=183/365;ctx.pool=ThreadPoolExecutor(max_workers=4);ctx.normal_laws={};ctx.last_readout=None;ctx.normal_cached=None
    try:
        norm,na=ctx.normalizer(p);dump(out/'normalizer.json',na);_,fresh,_,nn=exposures(ctx,k,p,norm,cache,specs=[('baseline',0.,1.,1.,1.)]);fn={f'd{day}_{key}':v for day,r in nn.items() for key,v in r.items()};np.savez_compressed(out/'fresh_conditional_nodes.npz',**fn)
        fe=ForwardFloorFit(fn,fresh,market,forwards,controls,p,ends);scores=[record(fe,x,out/'fresh_check',name) for name,x in candidates.items()];dump(out/'fresh_check/scores.json',scores)
        pairs=[paired(fe,out/'fresh_check',name,ref) for name,ref in [('tight_forwards','previous_mixture'),('tight_forwards_free_floor','previous_mixture'),('tight_forwards_free_floor','tight_forwards')]];dump(out/'fresh_check/paired_losses.json',pairs)
        probeprod=out/'quadrature_cache/production';probeprod.mkdir(parents=True);ids_by_day={}
        for day,aa in fresh.items():
            ids=list(np.r_[0,1,np.argmax(abs(aa['states']),axis=0),np.argmin(aa['z_normals']),np.argmax(aa['z_normals'])])
            for name in candidates:
                with np.load(out/'fresh_check'/f'paths_{name}.npz') as paths:
                    V=paths[f'd{day}_V'];w=paths[f'd{day}_weights'];ids.extend(np.argsort(w*np.maximum(V-market[market.expiry_days.eq(day)].strike.max(),0))[-3:])
            ids=np.unique(ids);ids_by_day[day]=ids;np.savez_compressed(probeprod/f'vix_d{day}.npz',**{key:val[ids] for key,val in aa.items()})
        _,qa,_,qn=exposures(ctx,k,p,norm,probeprod.parent,specs=[('baseline',0.,1.,1.,1.)],time_order=16,order=48);coarse={};fine={}
        for day,ids in ids_by_day.items():
            for key in ['times','weights']:coarse[f'd{day}_{key}']=fn[f'd{day}_{key}']
            coarse[f'd{day}_ratios']=fn[f'd{day}_ratios'][:,ids]
            for key,val in qn[str(day)].items():fine[f'd{day}_{key}']=val
        ce=ForwardFloorFit(coarse,qa,market,forwards,controls,p,ends);qe=ForwardFloorFit(fine,qa,market,forwards,controls,p,ends);qc=[]
        for name,x in candidates.items():
            cr=ce.predict(x,True);qr=qe.predict(x,True)
            for day in DAYS:qc.append(dict(candidate=name,day=day,maximum_VIX_difference=float(np.max(abs(cr['raw'][day]['V']-qr['raw'][day]['V']))),probe_count=len(ids_by_day[day]),coarse='GH32/time8',fine='GH48/time16'))
        dump(out/'quadrature_checks.json',qc)
    finally:ctx.close()
    check_sources(captured);dump(out/'run.json',design|dict(status='completed',training_scores=train,fresh_scores=scores,paired_losses=pairs,wall_seconds=time.perf_counter()-start,limitations=['Equality constraints apply to the training estimator; fresh forward residuals and SE reported.','Finite bounded multistart search, not a global optimum.','Same market data for numerical checking; not market OOS.','Sampled quadrature refinements do not establish uniform or population error bounds.','No new component/bank/readout, SPX assessment, baseline replacement or Lean theorem.']))
    print('COMPLETED',time.perf_counter()-start,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    try:run(args.output.resolve())
    except Exception as error:
        if args.output.exists():dump(args.output/'failure.json',dict(type=type(error).__name__,message=str(error)))
        raise
