"""Bounded joint SPX/VIX calibration of the selected two-speed model."""
import argparse,json,time,hashlib,subprocess,shutil
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize,least_squares
from .full_readout_fit import HERE,dump,arr
from .vix_only_curve_fit import PRIOR,ARCHIVE,DAYS,ACTIVE,load_draws
from .two_speed_study import engine as vix_engine,PREVIOUS as CONTINUUM
from .continuum_response_study import PARENT
from .two_speed_joint_core import StockKernel,SPXFit,JointFit,physical,metadata,vix_x
from .all_maturity_core import populations
from .m2_spx import ROOT
from .plot_joint_curve_fit import implied_vol
from esn_eq.market_targets import MarketTargets

SELECTED=HERE/'runs/m2_two_speed_001'
CACHE=HERE/'runs/m2_two_speed_joint_cache_001'


def load_spx(role):
    levels=np.array(json.loads((PRIOR/'design.json').read_text())['reference_levels']);data={}
    for day in [29,57,119,211,302,666,1030]:
        ds=[]
        for folder in sorted((CACHE/role/f'd{day}').glob('batch*')):
            d={name:np.load(folder/f'{name}.npy',mmap_mode='r') for name in ['old','equity','eps','time','dt','seg','grid','refbrownian','duration','weights']};d['refvariance']=float(d['duration']@levels);d['folder']=folder;ds.append(d)
        data[day]=ds
    return data,levels


def build_vix(out,role,market,forwards,p,ends):
    controls=json.loads((PRIOR/'vix_controls.json').read_text())
    if role=='search':nodes=dict(np.load(PARENT/'fresh_conditional_nodes.npz'));draws=load_draws(PARENT/'fresh_cache');seed=20261027201
    else:nodes=dict(np.load(CONTINUUM/'fresh_base_nodes.npz'));draws=load_draws(CONTINUUM/'fresh_base_cache');seed=20261027237
    return vix_engine(nodes,draws,p,ends,market,forwards,controls,seed,out/f'{role}_vix_fast')


def numpy_stock(d,x):
    from .two_speed_core import LOG_HALF_MAX
    hf=np.exp(LOG_HALF_MAX*x[21]);hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf]);b=5*x[[16,20]];z=np.zeros((d['old'].shape[1],2));logs=np.zeros(len(z));levels,f=physical(x)
    for t,dt in enumerate(d['dt']):
        q=-np.expm1(-2*rates*d['time'][t]);A=np.exp(z*b-.5*b*b*q);Q=(1-x[17])*d['old'][t]+x[17]*((1-x[22])*A[:,0]+x[22]*A[:,1]);v=f+(levels[d['seg'][t]]-f)*Q
        logs+=-.5*v*dt+np.sqrt(v)*d['equity'][t];z=np.exp(-rates*dt)*z+np.sqrt(-np.expm1(-2*rates*dt))*d['eps'][t]
    return np.exp(logs)


def kernel_checks(k,data,x,out):
    original=data[1030][0];ids=np.array([0,7,19,47]);d={name:arr(original[name][:,ids]) for name in ['old','equity','eps']};d.update({name:original[name] for name in ['time','dt','seg']});probe=x.copy();probe[:16]=.22;probe[16:]=[.35,.65,.6,.5,.4,.45,.4]
    actual=k.stock(d,probe);oracle=numpy_stock(d,probe);np.testing.assert_allclose(actual[:,0],oracle,rtol=1e-10,atol=1e-12);checks=[]
    for j in range(23):
        dx=np.zeros(23);dx[j]=1e-6;a=k.stock(d,probe+dx,False)[:,0];b=k.stock(d,probe-dx,False)[:,0];fd=(a-b)/2e-6;error=float(max(abs(fd-actual[:,j+1])/(1+abs(fd))));assert error<2e-5,(j,error);checks.append(dict(parameter=j,error=error))
    dump(out/'stock_kernel_checks.json',dict(numpy_replay_max_error=float(max(abs(actual[:,0]-oracle))),derivatives=checks,probe_days=1030,probe_paths=ids.tolist()));print('NUMPY STOCK REPLAY + ALL23 DERIVATIVES PASS',flush=True)


def record(e,x,out,name,phase):
    folder=out/phase;folder.mkdir(exist_ok=True);r=e.predict(x,True);s=r['spx'];v=r['vix'];records=[]
    for instrument,q,prices,blocks in [('SPX',e.spx.market,s['prices'],s['blocks']),('VIX',e.vix.market,v['prices'],v['blocks']['native'])]:
        q=q.copy();q['model_price_U']=prices;q['mc_se_U']=blocks.std(0,ddof=1)/np.sqrt(8);q['interval_residual']=np.where(prices<q.price_bid_U,prices-q.price_bid_U,np.where(prices>q.price_ask_U,prices-q.price_ask_U,0))/q.price_scale_U
        q['model_forward']=q.forward if instrument=='SPX' else q.expiry_days.map(dict(zip(DAYS,v['forwards'])))
        q['market_native_IV']=[implied_vol(a.price_mid_U,a.forward,a.strike,a.T_years,a.option_type) for a in q.itertuples()];q['model_native_IV']=[implied_vol(a.model_price_U,a.model_forward,a.strike,a.T_years,a.option_type) for a in q.itertuples()]
        q.to_csv(folder/f'prices_{instrument.lower()}_{name}.csv',index=False);np.savez_compressed(folder/f'blocks_{instrument.lower()}_{name}.npz',native=blocks)
        per=[]
        for day,g in q.groupby('expiry_days',sort=True):per.append(dict(day=int(day),inside=int(g.interval_residual.eq(0).sum()),quotes=len(g),loss=float(np.mean(g.interval_residual**2)),IV_MAE=float(np.nanmean(abs(g.model_native_IV-g.market_native_IV))),negative=int((g.model_price_U<0).sum()),median_SE_in_scales=float(np.median(g.mc_se_U/g.price_scale_U)),SE_gate_passes=int((2*g.mc_se_U<=g.price_scale_U/4).sum())))
        records.append(dict(instrument=instrument,inside=int(q.interval_residual.eq(0).sum()),quotes=len(q),negative=int((q.model_price_U<0).sum()),per_maturity=per))
    f=e.vix.forwards.copy();f['model_forward']=v['forwards'];f['forward_error']=f.model_forward-f.forward_vix_points;f['mc_se']=v['blocks']['forward'].std(0,ddof=1)/np.sqrt(8);f.to_csv(folder/f'forwards_{name}.csv',index=False)
    np.savez_compressed(folder/f'stocks_{name}.npz',**{f'd{day}_{key}':val for day,d in s['raw'].items() for key,val in d.items()})
    score=dict(name=name,parameters=x.tolist(),model=metadata(x,e.vix.p),total_option_loss=r['loss'],SPX_loss=r['spx_loss'],VIX_loss=r['vix_loss'],forward_errors=f.forward_error.tolist(),forward_SEs=f.mc_se.tolist(),instruments=records,stock_moments=s['moments'])
    dump(folder/f'score_{name}.json',score);print('RECORD',phase,name,'SPX',r['spx_loss'],'VIX',r['vix_loss'],'Fgap',max(abs(f.forward_error)),flush=True);return score


def fit(e,initial,label,free,out,ledger,maxiter=50):
    target=e.vix.forwards.forward_vix_points.to_numpy();fixed=initial.copy();memo={};ticks=[0];begin=time.perf_counter()
    def predict(u):
        key=tuple(u)
        if key not in memo:
            x=fixed.copy();x[free]=u;memo.clear();memo[key]=e.predict(x)
        return memo[key]
    def callback(u):
        ticks[0]+=1
        if ticks[0]%5==0:
            r=predict(u);print('JOINT ITER',label,ticks[0],'SPX',r['spx_loss'],'VIX',r['vix_loss'],'Fgap',max(abs(r['forwards']-target)),'seconds',round(time.perf_counter()-begin,1),flush=True)
    opt=minimize(lambda u:predict(u)['loss']/10,fixed[free],jac=lambda u:2*predict(u)['residuals']@predict(u)['jacobian'][:,free]/10,method='SLSQP',bounds=[(0,1)]*len(free),constraints=[dict(type='eq',fun=lambda u:predict(u)['forwards']-target,jac=lambda u:predict(u)['forward_jacobian'][:,free])],callback=callback,options=dict(maxiter=maxiter,ftol=1e-8))
    x=fixed.copy();x[free]=opt.x;r=e.predict(x);gap=float(max(abs(r['forwards']-target)));row=dict(label=label,parameters=x.tolist(),SPX_loss=r['spx_loss'],VIX_loss=r['vix_loss'],total=r['loss'],forward_gap=gap,eligible=gap<=1e-5,success=bool(opt.success),message=opt.message,iterations=opt.nit,nfev=opt.nfev,seconds=time.perf_counter()-begin);ledger.append(row);dump(out/'optimizer_ledger.json',ledger);print('JOINT FIT',json.dumps(row),flush=True);return x


def run(out):
    out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);selected=json.loads((SELECTED/'experimental_model_two_speed.json').read_text());y=np.array(json.loads((SELECTED/'selection.json').read_text())['parameters']);floor=selected['parameters'][-1];x=np.r_[(np.array(selected['parameters'][7:-1])-floor)/(.2-floor),y[7:]]
    sourcefiles=[Path(__file__),HERE/'two_speed_joint_core.py',HERE/'two_speed_spx_kernel.cpp',HERE/'prepare_two_speed_joint.py',HERE/'TWO_SPEED_JOINT_PLAN.md',HERE/'active_scope.json',HERE/'two_speed_core.py',HERE/'two_speed_multilevel_core.py',HERE/'two_speed_study.py',HERE/'forward_floor_core.py',HERE/'parity_adaptive.py',PRIOR/'run.json',SELECTED/'selection.json',SELECTED/'experimental_model_two_speed.json',CACHE/'design.json']
    src=out/'source';src.mkdir();hashes={}
    for f in sourcefiles:hashes[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest();shutil.copy2(f,src/f.name)
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');spx,vix,membership,atm=populations(targets.quotes);assert len(spx)==1294 and len(vix)==217;forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv')
    for name,q in [('spx_contracts',spx),('vix_contracts',vix),('spx_membership',membership),('atm_proxies',atm),('forward_targets',forwards)]:q.to_csv(out/f'{name}.csv',index=False)
    lib=out/'stock.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'two_speed_spx_kernel.cpp'),'-o',str(lib)],check=True);k=StockKernel(lib);data,ref=load_spx('search');s=SPXFit(data,spx,k,ref);kernel_checks(k,data,x,out)
    v=build_vix(out,'search',vix,forwards,p,ends);e=JointFit(s,v);replay=v.predict(vix_x(x));old=pd.read_csv(SELECTED/'training/prices_otm_two_speed.csv');assert max(abs(replay['prices']-old.model_price_U))<1e-10
    dump(out/'design.json',dict(source_hashes=hashes,target_hashes=targets.artifact_sha256,baseline=x.tolist(),SPX_primary_quotes=len(spx),SPX_eligible_quotes=len(membership),VIX_quotes=len(vix),stock_steps_per_trading_day=64,stock_paths_per_expiry=512,free_parameters=23,free_curve_pieces=16,original_M2_readout_frozen=True,SPX_paused=False,baseline_replaced=False,checking='Previously stored separate old-state streams plus new stress normals; VIX checking streams reused; same market snapshot',forward_constraints='four VIX training-estimator equalities; max gap<=1e-5'))
    scores=[record(e,x,out,'vix_only','training')];ledger=[]
    curve=fit(e,x,'curve_floor',np.r_[np.arange(16),19],out,ledger,40);scores.append(record(e,curve,out,'curve_floor','training'))
    joint=fit(e,curve,'joint_all23',np.arange(23),out,ledger,60);scores.append(record(e,joint,out,'joint','training'))
    candidates=[dict(label='vix_only',parameters=x.tolist(),eligible=True,total=scores[0]['total_option_loss'])]+ledger
    eligible=[a for a in candidates if a['eligible']];winner=min(eligible,key=lambda a:a['total']);xx=np.array(winner['parameters']);dump(out/'selection.json',dict(label=winner['label'],parameters=xx.tolist(),model=metadata(xx,p),selected_before_checking=True,adopted=False))
    del e,v,s,data
    checkdata,ref=load_spx('check');s=SPXFit(checkdata,spx,k,ref);v=build_vix(out,'check',vix,forwards,p,ends);ce=JointFit(s,v)
    checking=[record(ce,z,out,name,'checking') for name,z in [('vix_only',x),('curve_floor',curve),('joint',joint)]]
    levels,floor=physical(xx);pp=p.copy();pp[7:-1]=levels;pp[-1]=floor;dump(out/'experimental_model_joint.json',dict(parameters=pp.tolist(),ends_years=ends.tolist(),readout=metadata(xx,p),architecture=prior['architecture'],premium=prior['premium'],initial_state=prior['initial_state'],adopted=False))
    for f in sourcefiles:assert hashlib.sha256(f.read_bytes()).hexdigest()==hashes[str(f)]
    dump(out/'run.json',dict(status='bounded joint fit completed; refinement/audit follows',training_scores=scores,checking_scores=checking,selected=winner['label'],wall_seconds=time.perf_counter()-begin,baseline_replaced=False));print('JOINT STUDY COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();run(args.output.resolve())
