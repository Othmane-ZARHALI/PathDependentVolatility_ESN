"""Repair forward estimation using a matching cheap boundary as multilevel control."""
import argparse,json,hashlib,subprocess,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from .continuum_response_study import metadata,fit,checks,subset,PARENT
from .continuum_response_core import ContinuumFit
from .continuum_multilevel_core import MultilevelFit
from .forward_floor_core import ForwardFloorFit,unpack
from .forward_floor_study import record as boundary_record
from .vix_only_curve_fit import PRIOR,ARCHIVE,ACTIVE,DAYS,load_draws,paired
from .full_readout_fit import HERE,dump
from .m2_spx import ROOT,M2,frozen_stage3
from .mixed_response_study import ZSHIFTS,ZPROB,seed_for
from .postcap_vov_core import PostKernel,exposures
from .controlled_readout_fast_fit import ControlledEngine
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from .cap_controlled_core import shape_metrics
from .plot_joint_curve_fit import implied_vol
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture

INITIAL=HERE/'runs/m2_continuum_response_002'
BASE_SEED=20261027117

def power_nodes(folder):
    nodes={'u_grid':np.linspace(0,1,33)}
    for day in DAYS:
        with np.load(folder/f'd{day}.npz') as a:
            for key,v in a.items():nodes[f'd{day}_{key}']=v
    return nodes


def new_base(out,dated,initial):
    cache=out/'fresh_base_cache';prod=cache/'production';prod.mkdir(parents=True);props=json.loads((ARCHIVE/'m2_vix_targeted_check_001/proposals.json').read_text());seeds=[]
    for day in DAYS:
        dec=props[str(day)];mix=DefensiveGaussianMixture(np.array([np.r_[old,z] for old in dec['shifts'] for z in ZSHIFTS]),np.array([p*q for p in dec['probabilities'] for q in ZPROB]));sampler=GaussianOuterSampler(dated,day/365,initial);zz=[];ww=[];cc=[]
        for net in range(8):
            seed=seed_for(BASE_SEED,day,net);z,w,c=mix.draw(13,seed);zz.append(z);ww.append(w);cc.append(c);seeds.append(dict(day=day,net=net,seed=seed))
        z=np.concatenate(zz);np.savez_compressed(prod/f'vix_d{day}.npz',states=sampler.mean+z[:,:-1]@sampler.root.T,normals=z[:,:-1],z_normals=z[:,-1],weights=np.concatenate(ww),components=np.concatenate(cc),nets=np.repeat(np.arange(8),8192))
    dump(cache/'design.json',dict(seed_base=BASE_SEED,seeds=seeds,points_per_expiry=65536,old_proposals=props,Z_shifts=ZSHIFTS.tolist(),Z_probabilities=ZPROB.tolist(),created_after_selection=True));return cache

def record(e,x,folder,name):
    r=e.predict(x,True);folder.mkdir(parents=True,exist_ok=True)
    bn=r['blocks']['native'];bc=r['blocks']['call'];bf=r['blocks']['forward']
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



def run(out):
    out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);initial=np.array(prior['initial_state']);market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv');forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv');controls=json.loads((PRIOR/'vix_controls.json').read_text())
    sources=[*sorted(HERE.glob('*.py')),HERE/'CONTINUUM_MULTILEVEL_PLAN.md',HERE/'continuum_power_kernel.cpp',HERE/'postcap_vov_kernel.cpp',HERE/'stochastic_vov_kernel.cpp',HERE/'cap_tail_kernel.cpp',PRIOR/'run.json',PRIOR/'vix_controls.json',INITIAL/'selection.json',INITIAL/'optimizer_ledger.json',PARENT/'selection.json',*sorted((M2/'src/esn_eq').glob('*.py'))];captured=capture_sources(sources);save_sources(out,captured,ROOT)
    inputs=[PARENT/'fresh_conditional_nodes.npz',*[PARENT/'fresh_cache/production'/f'vix_d{day}.npz' for day in DAYS],*[INITIAL/'training_power_tables'/f'd{day}.npz' for day in DAYS],*[INITIAL/'fresh_power_tables'/f'd{day}.npz' for day in DAYS],*[INITIAL/'fresh_cache/production'/f'vix_d{day}.npz' for day in DAYS]]
    design=dict(status='running',source_sha256=source_hashes(captured,ROOT),input_sha256={str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in inputs},training='boundary32768 + paired continuum-minus-boundary8192; nested8nets',checking='new boundary65536 + reused difference16384; independent base/difference nets; not wholly fresh',model_unchanged=True,SPX_priced=False,baseline_replaced=False,forward_constraints='four equalities on multilevel training estimator')
    dump(out/'design.json',design);large_nodes=dict(np.load(PARENT/'fresh_conditional_nodes.npz'));large=load_draws(PARENT/'fresh_cache');small=subset(large,1024);nodes=power_nodes(INITIAL/'training_power_tables');ce=ContinuumFit(nodes,small,market,forwards,controls,p,ends);bs=ForwardFloorFit(nodes,small,market,forwards,controls,p,ends);bl=ForwardFloorFit(large_nodes,large,market,forwards,controls,p,ends);e=MultilevelFit(ce,bs,bl)
    old=json.loads((INITIAL/'selection.json').read_text());x=np.array(old['parameters']);probe=x.copy();probe[:7]=np.clip(probe[:7],.05,.8);probe[7:]=[.6,.65,.45,.35,.55];checks(e,probe,out)
    ledger=[];parent=np.array(old['previous_parameters']);bx=fit(bl,parent,'boundary_refit',out,ledger,100)
    oldledger=json.loads((INITIAL/'optimizer_ledger.json').read_text());broad=sorted([z for z in oldledger if z['label'].startswith('continuous') and z['eligible'] and z['model']['kappa']>=.1],key=lambda z:z['option_loss']);starts=[x,*[np.array(z['parameters']) for z in broad[:2]]]
    for j,xx in enumerate(starts):fit(e,xx,f'continuous_start{j}',out,ledger,80)
    eligible=[z for z in ledger if z['label'].startswith('continuous') and z['eligible']];assert eligible;winner=min(eligible,key=lambda z:z['option_loss']);selected=np.array(winner['parameters']);broads=[z for z in eligible if z['model']['kappa']>=.1];secondary=min(broads,key=lambda z:z['option_loss']) if broads else None
    dump(out/'selection.json',dict(parameters=selected.tolist(),model=metadata(selected,p),boundary_parameters=bx.tolist(),previous_parameters=parent.tolist(),secondary=secondary,frozen_before_new_base_draws=True,baseline_replaced=False))
    train=[record(e,selected,out/'training','continuous'),boundary_record(bl,bx,out/'training','boundary'),boundary_record(bl,parent,out/'training','previous')]
    if secondary:train.append(record(e,np.array(secondary['parameters']),out/'training','interior'))
    dump(out/'training/scores.json',train)
    # Release large response splines before constructing the checking evaluator.
    del e,ce,bs,bl,nodes,large_nodes,large,small
    a,_,prem,_=frozen_stage3();dated=DatedM2Pricer(a,prem);cache=new_base(out,dated,initial);lib=out/'base.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'postcap_vov_kernel.cpp'),'-o',str(lib)],check=True);k=PostKernel(lib,ends)
    ctx=object.__new__(ControlledEngine);ctx.dated=dated;ctx.initial=initial;ctx.k=k;ctx.ends=ends;ctx.tau=1/max(dated.layout.rates);ctx.horizon=183/365;ctx.pool=ThreadPoolExecutor(max_workers=4);ctx.normal_laws={};ctx.last_readout=None;ctx.normal_cached=None
    try:
        norm,na=ctx.normalizer(p);dump(out/'base_normalizer.json',na);_,big,_,bn=exposures(ctx,k,p,norm,cache,specs=[('baseline',0.,1.,1.,1.)]);base_nodes={f'd{day}_{key}':v for day,d in bn.items() for key,v in d.items()};np.savez_compressed(out/'fresh_base_nodes.npz',**base_nodes)
    finally:ctx.close()
    draws=load_draws(INITIAL/'fresh_cache');nodes=power_nodes(INITIAL/'fresh_power_tables');ce=ContinuumFit(nodes,draws,market,forwards,controls,p,ends);bs=ForwardFloorFit(nodes,draws,market,forwards,controls,p,ends);bl=ForwardFloorFit(base_nodes,big,market,forwards,controls,p,ends);fe=MultilevelFit(ce,bs,bl)
    scores=[record(fe,selected,out/'fresh_check','continuous'),boundary_record(bl,bx,out/'fresh_check','boundary'),boundary_record(bl,parent,out/'fresh_check','previous')]
    if secondary:scores.append(record(fe,np.array(secondary['parameters']),out/'fresh_check','interior'))
    dump(out/'fresh_check/scores.json',scores);pairs=[paired(fe,out/'fresh_check','continuous',ref) for ref in ['boundary','previous']]
    if secondary:pairs.extend([paired(fe,out/'fresh_check','interior',ref) for ref in ['boundary','previous','continuous']])
    dump(out/'fresh_check/paired_losses.json',pairs);table=pd.DataFrame(dict(start_day=np.r_[0,ends[:-1]]*365,end_day=ends*365,active_in_VIX_fit=np.isin(np.arange(16),ACTIVE)))
    for name,xx in [('continuous',selected),('boundary',bx),('previous',parent)]+([('interior',np.array(secondary['parameters']))] if secondary else []):
        levels,floor,_=unpack(xx,p);pp=p.copy();pp[7:-1]=levels;pp[-1]=floor;table[name]=levels;dump(out/f'experimental_model_{name}.json',dict(parameters=pp.tolist(),ends_years=ends.tolist(),readout=metadata(xx,p) if len(xx)==12 else dict(boundary_parameters=xx.tolist(),floor_volatility=float(np.sqrt(floor))),architecture=prior['architecture'],premium=prior['premium'],initial_state=initial.tolist(),adopted=False))
    table.to_csv(out/'curve_levels.csv',index=False);check_sources(captured);dump(out/'run.json',design|dict(status='completed multilevel fit/check; audit follows',training_scores=train,checking_scores=scores,paired_losses=pairs,wall_seconds=time.perf_counter()-begin));print('MULTILEVEL COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    try:run(args.output.resolve())
    except Exception as error:
        if args.output.exists():dump(args.output/'failure.json',dict(type=type(error).__name__,message=str(error)))
        raise
