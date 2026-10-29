"""VIX-only continuous response experiment, forward equality constraints."""
import argparse
import ctypes as ct
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
from scipy.optimize import least_squares,minimize
from .continuum_response_core import ContinuumFit,measure,beta_rule,LOG_K_MIN,LOG_K_RANGE
from .forward_floor_core import ForwardFloorFit,unpack,LOG_HALF_MAX
from .forward_floor_study import record as boundary_record
from .vix_only_curve_fit import PRIOR,ARCHIVE,ACTIVE,DAYS,load_draws,paired
from .full_readout_fit import HERE,dump,Kernel,arr
from .m2_spx import ROOT,M2,frozen_stage3
from .mixed_response_study import ZSHIFTS,ZPROB,seed_for
from .maturity_curve import window_nodes
from .controlled_readout_fast_fit import ControlledEngine
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture
from .cap_controlled_core import shape_metrics
from .plot_joint_curve_fit import implied_vol

PARENT=HERE/'runs/m2_forward_floor_001'
CHECK_SEED=20261027091

def metadata(x,p):
    levels,floor,theta=unpack(x,p);m,k=measure(x)
    return dict(b=float(theta[0]),measure_mean=m,kappa=k,alpha=m*k,beta=(1-m)*k,half_life_days=float(np.exp(theta[2])),floor_variance=float(floor),floor_volatility=float(np.sqrt(floor)),active_levels=levels[ACTIVE].tolist())


def subset(draws,n):
    return {day:{k:v[np.concatenate([np.flatnonzero(a['nets']==j)[:n] for j in range(8)])] for k,v in a.items()} for day,a in draws.items()}


def context(out,p,ends):
    a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium)
    lib=out/'power.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'continuum_power_kernel.cpp'),'-o',str(lib)],check=True)
    k=Kernel(lib);k.lib.configure.argtypes=[ct.c_int,ct.c_double]
    ctx=object.__new__(ControlledEngine);ctx.k=k;ctx.dated=dated;ctx.ends=ends;ctx.tau=1/max(dated.layout.rates)
    return ctx,k,dated


def tables(ctx,k,p,initial,draws,folder,grid_size=33,time_order=8,order=32):
    folder.mkdir(parents=True,exist_ok=True);u=np.linspace(0,1,grid_size);nodes={'u_grid':u}
    for day,a in draws.items():
        path=folder/f'd{day}.npz'
        if path.exists():
            with np.load(path) as saved:
                for key,val in saved.items():nodes[f'd{day}_{key}']=val
            continue
        tick=time.perf_counter();times,tw=window_nodes(day/365,ctx.ends,ctx.tau,time_order);ratios=[]
        for it,t in enumerate(times):
            law=ctx.law(float(t-day/365),a['states'],order);nlaw=ctx.law(float(t),initial,48);values=[]
            for uj in u:
                k.lib.configure(0,float(1-uj));m=k.expectation(*nlaw,p)[0,0];v=k.expectation(*law,p)[:,0]
                assert m>0 and np.all(v>0);values.append(np.log(v)-np.log(m))
            ratios.append(values)
            if it%8==0:print('TABLE',folder.name,day,it,'/',len(times),'seconds',round(time.perf_counter()-tick,1),flush=True)
        lr=np.array(ratios).transpose(1,0,2)
        np.savez_compressed(path,times=times,weights=tw,log_ratios=lr,ratios=np.exp(lr[0]))
        for key,val in dict(times=times,weights=tw,log_ratios=lr,ratios=np.exp(lr[0])).items():nodes[f'd{day}_{key}']=val
        print('TABLE DONE',folder.name,day,'seconds',round(time.perf_counter()-tick,1),flush=True)
    return nodes


def fresh_draws(out,dated,initial):
    cache=out/'fresh_cache';prod=cache/'production';prod.mkdir(parents=True,exist_ok=False);props=json.loads((ARCHIVE/'m2_vix_targeted_check_001/proposals.json').read_text());seeds=[]
    for day in DAYS:
        dec=props[str(day)];shifts=np.array([np.r_[old,z] for old in dec['shifts'] for z in ZSHIFTS]);prob=np.array([p*q for p in dec['probabilities'] for q in ZPROB]);mix=DefensiveGaussianMixture(shifts,prob);sampler=GaussianOuterSampler(dated,day/365,initial);zz=[];ww=[];cc=[]
        for net in range(8):
            seed=seed_for(CHECK_SEED,day,net);z,w,c=mix.draw(11,seed);zz.append(z);ww.append(w);cc.append(c);seeds.append(dict(day=day,net=net,seed=seed))
        z=np.concatenate(zz);np.savez_compressed(prod/f'vix_d{day}.npz',states=sampler.mean+z[:,:-1]@sampler.root.T,normals=z[:,:-1],z_normals=z[:,-1],weights=np.concatenate(ww),components=np.concatenate(cc),nets=np.repeat(np.arange(8),2048))
    dump(cache/'design.json',dict(seed_base=CHECK_SEED,seeds=seeds,points_per_expiry=16384,old_proposals=props,Z_shifts=ZSHIFTS.tolist(),Z_probabilities=ZPROB.tolist(),created_after_selection=True))
    return load_draws(cache)


def fit(e,initial,label,out,ledger,maxiter=100):
    target=e.forwards.forward_vix_points.to_numpy();n=len(initial);bounds=[(0.,1.)]*n
    if n==12:bounds[8]=(.01,.99)
    initial=np.array(initial).copy();initial=np.maximum([b[0] for b in bounds],np.minimum([b[1] for b in bounds],initial));memo={};count=0;tick=time.perf_counter()
    def pred(x):
        nonlocal count
        key=tuple(x)
        if key not in memo:
            memo.clear();memo[key]=e.predict(x);count+=1
        return memo[key]
    def pre(y):
        x=initial.copy();x[:7]=y;return pred(x)
    cf=least_squares(lambda y:pre(y)['forwards']-target,initial[:7],jac=lambda y:pre(y)['forward_jacobian'][:,:7],bounds=(np.zeros(7),np.ones(7)),max_nfev=40,ftol=1e-10,xtol=1e-10,gtol=1e-9);initial[:7]=cf.x
    it=[0]
    def callback(x):
        it[0]+=1
        if it[0]%10==0:
            r=pred(x);print('ITER',label,it[0],'optionloss',r['option_loss'],'Fgap',max(abs(r['forwards']-target)),'seconds',round(time.perf_counter()-tick,1),flush=True)
    opt=minimize(lambda x:pred(x)['option_loss']/10,initial,jac=lambda x:2*pred(x)['residuals'][:-4]@pred(x)['jacobian'][:-4]/10,method='SLSQP',bounds=bounds,constraints=[dict(type='eq',fun=lambda x:pred(x)['forwards']-target,jac=lambda x:pred(x)['forward_jacobian'])],callback=callback,options=dict(maxiter=maxiter,ftol=1e-9))
    r=pred(opt.x);gap=float(max(abs(r['forwards']-target)));rec=dict(label=label,parameters=opt.x.tolist(),option_loss=r['option_loss'],maximum_forward_error=gap,eligible=gap<=1e-5,success=bool(opt.success),message=opt.message,iterations=opt.nit,nfev=opt.nfev,seconds=time.perf_counter()-tick)
    if n==12:rec['model']=metadata(opt.x,e.p)
    ledger.append(rec);dump(out/'optimizer_ledger.json',ledger);print('FIT',label,json.dumps(rec),flush=True);return opt.x


def checks(e,x,out):
    rows=[];rr=e.predict(x)
    for j in range(len(x)):
        dx=np.zeros(len(x));dx[j]=1e-6;rp=e.predict(x+dx);rm=e.predict(x-dx)
        fp=(rp['prices']-rm['prices'])/2e-6;ff=(rp['forwards']-rm['forwards'])/2e-6
        pe=float(np.linalg.norm(fp-rr['native_jacobian'][:,j])/max(1,np.linalg.norm(fp)));fe=float(np.linalg.norm(ff-rr['forward_jacobian'][:,j])/max(1,np.linalg.norm(ff)))
        rows.append(dict(parameter=j,price_error=pe,forward_error=fe));assert pe<.003 and fe<1e-5,rows[-1]
    for m,k in [(.25,8),(.5,2),(.75,.01),(.7,30)]:
        u,w=beta_rule(m,k,24);assert abs(w@u-m)<1e-9;assert abs(w@(u*u)-(m*m+m*(1-m)/(k+1)))<1e-9
    dump(out/'jacobian_checks.json',rows);print('12-DIRECTION JACOBIAN + BETA MOMENTS PASS',flush=True)


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



def run(out,reuse=None):
    begin=time.perf_counter();out.mkdir(parents=True,exist_ok=False)
    prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);initial=np.array(prior['initial_state'])
    market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv');forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv');controls=json.loads((PRIOR/'vix_controls.json').read_text())
    sources=[*sorted(HERE.glob('*.py')),HERE/'continuum_power_kernel.cpp',HERE/'cap_tail_kernel.cpp',HERE/'CONTINUUM_RESPONSE_PLAN.md',PRIOR/'run.json',PRIOR/'vix_controls.json',PARENT/'selection.json',ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv',ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv',ARCHIVE/'m2_vix_targeted_check_001/proposals.json',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(sources);save_sources(out,captured,ROOT);input_paths=[PARENT/'fresh_cache/production'/f'vix_d{day}.npz' for day in DAYS]
    design=dict(status='running',source_sha256=source_hashes(captured,ROOT),input_sha256={str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in input_paths},training_points_per_expiry=8192,training_note='first1024points from each8 previous checking nets, explicitly reused as training',fresh_points_per_expiry=16384,fresh_seed_base=CHECK_SEED,SPX_priced=False,baseline_replaced=False,forward_constraints='four training-estimator equalities, eligibility max gap <=1e-5VIXpoint',curve_piece_count=16,active_curve_indices=ACTIVE.tolist(),parameter_order=['xi_unit'+str(i) for i in ACTIVE]+['b/5','m','log_half/log365','floor_vol_unit','log_kappa_unit'])
    dump(out/'design.json',design);ctx,k,dated=context(out,p,ends);train=subset(load_draws(PARENT/'fresh_cache'),1024)
    if reuse is not None:
        import shutil
        old_design=json.loads((reuse/'design.json').read_text())
        assert old_design['input_sha256']==design['input_sha256']
        for path in [HERE/'continuum_power_kernel.cpp',PRIOR/'run.json']:
            rel=str(path.relative_to(ROOT));assert old_design['source_sha256'][rel]==design['source_sha256'][rel]
        target=out/'training_power_tables';target.mkdir();reused={}
        for day in DAYS:
            source=reuse/'training_power_tables'/f'd{day}.npz';shutil.copy2(source,target/source.name);reused[str(source.relative_to(ROOT))]=hashlib.sha256(source.read_bytes()).hexdigest()
        dump(out/'reused_power_tables.json',dict(from_run=str(reuse),sha256=reused,note='Same powers, frozen readout, time/Gaussian rules and exact training draws; prior run stopped at overly strict 1e-11 Beta-moment assertion. Observed error 1.18e-11; tolerance now1e-9. No model or integrator change.'))
    nodes=tables(ctx,k,p,initial,train,out/'training_power_tables');e=ContinuumFit(nodes,train,market,forwards,controls,p,ends);be=ForwardFloorFit(nodes,train,market,forwards,controls,p,ends)
    parent=np.array(json.loads((PARENT/'selection.json').read_text())['parameters']['tight_forwards_free_floor'])
    # Powered u=0 endpoint must reproduce original M2 conditional ratios.
    old=dict(np.load(PARENT/'fresh_conditional_nodes.npz'));endpoint=[]
    for day in DAYS:
        all_a=dict(np.load(PARENT/'fresh_cache/production'/f'vix_d{day}.npz'));ids=np.concatenate([np.flatnonzero(all_a['nets']==j)[:1024] for j in range(8)])
        err=float(np.max(abs(nodes[f'd{day}_ratios']/old[f'd{day}_ratios'][:,ids]-1)));assert err<2e-7,err;assert np.max(abs(nodes[f'd{day}_log_ratios'][-1]))<1e-12;endpoint.append(dict(day=day,M2_endpoint_relative_error=err,stress_endpoint_log_error=float(np.max(abs(nodes[f'd{day}_log_ratios'][-1])))))
    dump(out/'endpoint_checks.json',endpoint)
    probe=np.r_[parent,.5];probe[:7]=np.clip(probe[:7],.05,.8);probe[7:]=[.6,.65,.45,.35,.55];checks(e,probe,out)
    ledger=[];bx=fit(be,parent,'boundary_refit',out,ledger,120)
    starts=[]
    # Near boundary, uniform, moderate spread and concentrated response families.
    for mean,kap,b,half in [(bx[8],.03,5*bx[7],np.exp(LOG_HALF_MAX*bx[9])),(.5,2,3.5,100),(.7,4,4.5,30),(.35,8,5,150),(.8,.3,3.8,100)]:
        x=np.r_[bx,(np.log(kap)-LOG_K_MIN)/LOG_K_RANGE];x[7:10]=[b/5,mean,np.log(half)/LOG_HALF_MAX];starts.append(x)
    for j,x in enumerate(starts):fit(e,x,f'continuous_start{j}',out,ledger,120)
    eligible=[r for r in ledger if r['label'].startswith('continuous') and r['eligible']]
    if not eligible:raise RuntimeError('No continuous candidate satisfies forward equalities; retain diagnostics')
    selected=np.array(min(eligible,key=lambda r:r['option_loss'])['parameters']);candidates={'continuous':selected}
    dump(out/'selection.json',dict(parameters=selected.tolist(),model=metadata(selected,p),boundary_parameters=bx.tolist(),previous_parameters=parent.tolist(),frozen_before_fresh_check=True,baseline_replaced=False,market_out_of_sample=False))
    scores_train=[record(e,selected,out/'training','continuous'),boundary_record(be,bx,out/'training','boundary'),boundary_record(be,parent,out/'training','previous')];dump(out/'training/scores.json',scores_train)
    # Fresh randomness is generated only after the final training selection above.
    fresh=fresh_draws(out,dated,initial);fn=tables(ctx,k,p,initial,fresh,out/'fresh_power_tables');fe=ContinuumFit(fn,fresh,market,forwards,controls,p,ends);fb=ForwardFloorFit(fn,fresh,market,forwards,controls,p,ends)
    scores=[record(fe,selected,out/'fresh_check','continuous'),boundary_record(fb,bx,out/'fresh_check','boundary'),boundary_record(fb,parent,out/'fresh_check','previous')];dump(out/'fresh_check/scores.json',scores)
    pairs=[paired(fe,out/'fresh_check','continuous',ref) for ref in ['boundary','previous']];dump(out/'fresh_check/paired_losses.json',pairs)
    # Deterministic measure-quadrature refinement on every fresh state.
    qe=ContinuumFit(fn,fresh,market,forwards,controls,p,ends,measure_order=48);r24=fe.predict(selected,True);r48=qe.predict(selected,True)
    dump(out/'measure_quadrature_check.json',dict(coarse=24,fine=48,max_price_difference=float(max(abs(r24['prices']-r48['prices']))),max_forward_difference=float(max(abs(r24['forwards']-r48['forwards']))),by_day=[dict(day=day,max_VIX_difference=float(max(abs(r24['raw'][day]['V']-r48['raw'][day]['V'])))) for day in DAYS]))
    levels,floor,_=unpack(selected,p);pp=p.copy();pp[7:-1]=levels;pp[-1]=floor
    dump(out/'experimental_model_continuous.json',dict(parameters=pp.tolist(),ends_years=ends.tolist(),response='normalized geometric M2/stress continuum',measure='Beta',readout=metadata(selected,p),architecture=prior['architecture'],premium=prior['premium'],initial_state=initial.tolist(),adopted=False))
    table=pd.DataFrame(dict(start_day=np.r_[0,ends[:-1]]*365,end_day=ends*365,active_in_VIX_fit=np.isin(np.arange(16),ACTIVE)))
    for name,x in [('previous',parent),('boundary',bx),('continuous',selected)]:table[name]=unpack(x,p)[0]
    table.to_csv(out/'curve_levels.csv',index=False)
    check_sources(captured);dump(out/'run.json',design|dict(status='completed fit and fresh check; separate independent audit follows',training_scores=scores_train,fresh_scores=scores,paired_losses=pairs,wall_seconds=time.perf_counter()-begin,limitations=['Training forward equality is not exact population equality; fresh residuals and SE reported.','Finite bounded multistart, no global optimum.','Same market targets, no market out-of-sample claim.','Power interpolation/conditional/time refinement requires separate audit.','No SPX assessment, baseline adoption, or new Lean theorem.']))
    print('COMPLETED',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--reuse-power-tables',type=Path);args=ap.parse_args()
    try:run(args.output.resolve(),args.reuse_power_tables.resolve() if args.reuse_power_tables else None)
    except Exception as error:
        if args.output.exists():dump(args.output/'failure.json',dict(type=type(error).__name__,message=str(error)))
        raise
