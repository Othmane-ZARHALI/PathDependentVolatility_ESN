"""Finite VIX-only fits of two-response readouts, with independently free xi pieces."""
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
from scipy.optimize import least_squares
from .mixed_response_core import MixedFit,parameters,response
from .vix_only_curve_fit import PRIOR,POST,ARCHIVE,ACTIVE,DAYS,DELTA,VIXFit,load_draws,paired
from .full_readout_fit import HERE,dump
from .m2_spx import ROOT,M2,frozen_stage3
from .postcap_vov_core import PostKernel,exposures,conditional_multiplier
from .postcap_vov_study import shape_score
from .cap_controlled_core import shape_metrics
from .controlled_readout_fast_fit import ControlledEngine
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture

REFERENCE=HERE/'runs/m2_vix_only_curve_fit_002'
TRAIN_SEED=20261026071
CHECK_SEED=20261026072
ZSHIFTS=np.array([0.,2.,4.,6.]);ZPROB=np.array([.5,.25,.2,.05])
PRODUCT_STARTS=[(0.,1.,.1,7),(0.,2.5,.05,7),(0.,4.5,.01,7),(0.,2.5,.1,90),(0.,4.5,.03,90),(0.,4.5,.1,365),(.25,1.5,.2,7),(.25,3.,.1,30),(.25,4.5,.1,180),(.5,2.5,.25,365),(.1,4.5,.5,7),(0.,4.5,.8,365)]
ADDITIVE_STARTS=[(1.,.2,7),(2.5,.1,7),(4.5,.03,7),(2.5,.3,90),(4.5,.1,90),(4.5,.4,365),(1.5,.8,365),(4.5,.8,30)]


def seed_for(base,*ids):return int(np.random.SeedSequence([base,*ids]).generate_state(1,dtype=np.uint64)[0])


def training_inputs(out):
    old=load_draws(POST/'fresh_cache');saved=dict(np.load(POST/'fresh_check/base_conditional_nodes.npz'));draws={};nodes={};seeds=[]
    mix=DefensiveGaussianMixture(ZSHIFTS[:,None],ZPROB)
    for day in DAYS:
        ids=[];zz=[];ww=[]
        for net in range(8):
            for rep in range(2):
                seed=seed_for(TRAIN_SEED,day,net,rep);z,w,_=mix.draw(9,seed)
                ids.extend(range(net*512,(net+1)*512));zz.extend(z[:,0]);ww.extend(w);seeds.append(dict(day=day,net=net,replica=rep,seed=seed))
        ids=np.array(ids);a={key:value[ids] for key,value in old[day].items()}
        a['z_normals']=np.array(zz);a['weights']=old[day]['weights'][ids]*np.array(ww);a['nets']=np.repeat(np.arange(8),1024)
        draws[day]=a
        for key in ['times','weights']:nodes[f'd{day}_{key}']=saved[f'd{day}_{key}']
        nodes[f'd{day}_ratios']=saved[f'd{day}_ratios'][:,ids]
    dump(out/'training_design.json',dict(seeds=seeds,old_cache=str(POST/'fresh_cache'),old_nodes=str(POST/'fresh_check/base_conditional_nodes.npz'),replicas_per_old_state=2,points_per_expiry=8192,Z_shifts=ZSHIFTS.tolist(),Z_probabilities=ZPROB.tolist(),weights='w_old times exact independent Z target/proposal weight; no self normalization',independent_uncertainty_units=8))
    return nodes,draws


def fresh_inputs(out,dated,initial):
    cache=out/'fresh_cache';prod=cache/'production';prod.mkdir(parents=True);seeds=[]
    proposals=json.loads((ARCHIVE/'m2_vix_targeted_check_001/proposals.json').read_text())
    for day in DAYS:
        dec=proposals[str(day)];shifts=np.array([np.r_[old,z] for old in dec['shifts'] for z in ZSHIFTS]);prob=np.array([p*q for p in dec['probabilities'] for q in ZPROB]);mix=DefensiveGaussianMixture(shifts,prob)
        sampler=GaussianOuterSampler(dated,day/365,initial);zz=[];ww=[];cc=[]
        for net in range(8):
            seed=seed_for(CHECK_SEED,day,net);z,w,c=mix.draw(11,seed);zz.append(z);ww.append(w);cc.append(c);seeds.append(dict(day=day,net=net,seed=seed))
        z=np.concatenate(zz);np.savez_compressed(prod/f'vix_d{day}.npz',states=sampler.mean+z[:,:-1]@sampler.root.T,normals=z[:,:-1],z_normals=z[:,-1],weights=np.concatenate(ww),components=np.concatenate(cc),nets=np.repeat(np.arange(8),2048))
    dump(cache/'design.json',dict(seed_base=CHECK_SEED,seeds=seeds,points_per_expiry=16384,old_proposals=proposals,Z_shifts=ZSHIFTS.tolist(),Z_probabilities=ZPROB.tolist(),proposal='Product of old17-component mixture and four Z components:68 joint37D components; exact weights',created_after_selection=True))
    return cache


def record(e,x,folder,name,market_shape):
    r=e.predict(x,True);folder.mkdir(parents=True,exist_ok=True)
    native=np.concatenate([d['native'] for d in r['raw'].values()],axis=1);call=np.concatenate([d['call'] for d in r['raw'].values()],axis=1);fw=np.column_stack([d['forward'] for d in r['raw'].values()])
    bn=native.reshape(8,-1,native.shape[1]).mean(1);bc=call.reshape(8,-1,call.shape[1]).mean(1);bf=fw.reshape(8,-1,4).mean(1)
    for typ,q,val,block in [('otm',e.market,r['prices'],bn),('canonical',e.calls,r['calls'],bc)]:
        q=q.copy();q['model_price_U']=val;q['mc_se_U']=block.std(0,ddof=1)/np.sqrt(8)
        q['interval_residual']=np.where(val<q.price_bid_U,val-q.price_bid_U,np.where(val>q.price_ask_U,val-q.price_ask_U,0.))/q.price_scale_U
        q.to_csv(folder/f'prices_{typ}_{name}.csv',index=False)
    ft=e.forwards.copy();ft['model_forward']=r['forwards'];ft['mc_se']=bf.std(0,ddof=1)/np.sqrt(8);ft.to_csv(folder/f'forwards_{name}.csv',index=False)
    np.savez_compressed(folder/f'blocks_{name}.npz',native=bn,call=bc,forward=bf)
    np.savez_compressed(folder/f'paths_{name}.npz',**{f'd{day}_{k}':v for day,d in r['raw'].items() for k,v in d.items()})
    shape=shape_metrics(e.calls,r['calls']);summary=dict(name=name,family=e.family,parameters=x.tolist(),response=parameters(e.family,x),loss=r['loss'],inside_bidask=int(np.sum((r['prices']>=e.market.price_bid_U)&(r['prices']<=e.market.price_ask_U))),negative_prices=int(np.sum(r['prices']<0)),shape=shape,shape_score=shape_score(shape,market_shape),moments=r['moments'])
    dump(folder/f'score_{name}.json',summary);return summary


def validation(engines,prior,p,nodes,draws,market,forwards,controls,ends):
    checks=[];old=VIXFit(nodes,draws,market,forwards,controls,p,ends)
    previous=np.array(prior['parameters']['factor_curve']);base=np.r_[p[7:-1][ACTIVE],0.,np.log(7.)]
    for name,x in [('baseline',base),('single_factor',previous)]:
        xx=np.r_[x[:7],x[7],0.,.2,x[8]];a=engines['product'].predict(xx);b=old.predict(x)
        err=float(max(np.max(abs(a['prices']-b['prices'])),np.max(abs(a['forwards']-b['forwards']))));assert err<2e-10
        checks.append(dict(kind='single_factor_nesting',name=name,max_error=err))
    xx=np.r_[base[:7],2.,0.,np.log(30.)];a=engines['additive'].predict(xx);b=old.predict(base);err=float(np.max(abs(a['prices']-b['prices'])));assert err<2e-10;checks.append(dict(kind='additive_zero_weight_nesting',max_error=err))
    for family,e in engines.items():
        x=np.r_[base[:7],.3,.35,.17,np.log(37.)] if family=='product' else np.r_[base[:7],1.7,.3,np.log(37.)]
        a=e.predict(x)
        for j in range(len(x)):
            h=1e-6 if j<7 else 1e-5;dx=np.zeros_like(x);dx[j]=h;fd=(e.predict(x+dx)['residuals']-e.predict(x-dx)['residuals'])/(2*h)
            error=float(np.linalg.norm(fd-a['jacobian'][:,j])/max(1,np.linalg.norm(fd)));assert error<.003,(family,j,error);checks.append(dict(kind='jacobian',family=family,parameter=j,relative_error=error))
    # Independent Gaussian integration over future Z innovations, including b=5.
    z,w=np.polynomial.hermite.hermgauss(96);z=z*np.sqrt(2);w=w/np.sqrt(np.pi)
    for b in [0.,.3,2.,5.]:
        for half in [1.,30.,365.]:
            rate=np.log(2)*365/half;T=55/365;lag=17/365;qt=-np.expm1(-2*rate*T);qu=-np.expm1(-2*rate*(T+lag));decay=np.exp(-rate*lag)
            for eps in [-2.,0.,3.]:
                ZT=np.sqrt(qt)*eps;zu=decay*ZT+np.sqrt(-np.expm1(-2*rate*lag))*z
                dense=float(w@np.exp(b*zu-.5*b*b*qu));analytic=float(conditional_multiplier(ZT,T,lag,b,rate));error=abs(dense/analytic-1);assert error<1e-11
                checks.append(dict(kind='conditional_gaussian_oracle',b=b,half=half,eps=eps,relative_error=error))
            mean=float(w@np.exp(b*np.sqrt(qt)*z-.5*b*b*qt));assert abs(mean-1)<1e-11
            checks.append(dict(kind='mean_one_oracle',b=b,half=half,error=abs(mean-1)))
    return checks


def run(out):
    out.mkdir(parents=True,exist_ok=False);start=time.perf_counter();d=json.loads((PRIOR/'run.json').read_text());p=np.array(d['selected']);ends=np.array(d['ends_years'])
    market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv');forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv');controls=json.loads((PRIOR/'vix_controls.json').read_text());previous=json.loads((REFERENCE/'selection.json').read_text())
    sources=[*sorted(HERE.glob('*.py')),HERE/'MIXED_RESPONSE_PLAN.md',HERE/'postcap_vov_kernel.cpp',HERE/'stochastic_vov_kernel.cpp',HERE/'cap_tail_kernel.cpp',PRIOR/'run.json',PRIOR/'selected_model.json',PRIOR/'vix_controls.json',REFERENCE/'selection.json',ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv',ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv',ARCHIVE/'m2_vix_targeted_check_001/proposals.json',M2/'experiments/stage3_001/manifest.json',M2/'experiments/stage3_001/top_configurations.csv',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(sources);save_sources(out,captured,ROOT)
    inputs=[POST/'fresh_check/base_conditional_nodes.npz',*(POST/'fresh_cache/production'/f'vix_d{day}.npz' for day in DAYS)]
    design=dict(status='running',source_sha256=source_hashes(captured,ROOT),input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},curve_piece_count=16,active_curve_indices_zero_based=ACTIVE.tolist(),SPX_priced=False,baseline_replaced=False,product_starts=PRODUCT_STARTS,additive_starts=ADDITIVE_STARTS,training_seed_base=TRAIN_SEED,checking_seed_base=CHECK_SEED,market_out_of_sample=False)
    dump(out/'design.json',design);nodes,draws=training_inputs(out);engines={family:MixedFit(nodes,draws,market,forwards,controls,p,ends,family=family) for family in ['product','additive']}
    dump(out/'implementation_checks.json',validation(engines,previous,p,nodes,draws,market,forwards,controls,ends));print('IMPLEMENTATION CHECKS PASS',flush=True)
    market_shape=shape_metrics(engines['product'].calls,engines['product'].calls.price_mid_U.to_numpy());dump(out/'market_shape.json',market_shape)
    candidates={'baseline':('product',np.r_[p[7:-1][ACTIVE],0.,0.,0.,np.log(7.)])};ledger=[];evaluations=[];pool=[]
    def fit(family,label,initial,free=None):
        e=engines[family];free=np.arange(len(initial)) if free is None else free;memo={};fixed=initial.copy()
        lo=np.r_[np.full(7,p[-1]),0.,0.,0.,np.log(1.)] if family=='product' else np.r_[np.full(7,p[-1]),0.,0.,np.log(1.)]
        hi=np.r_[np.full(7,.2),3.,1.,1.,np.log(365.)] if family=='product' else np.r_[np.full(7,.2),5.,1.,np.log(365.)]
        scale=np.r_[np.full(7,.06),1.,.3,.2,2.] if family=='product' else np.r_[np.full(7,.06),2.,.3,2.]
        def current(xx):
            key=tuple(xx)
            if key not in memo:
                x=fixed.copy();x[free]=xx;rr=e.predict(x);memo.clear();memo[key]=rr;evaluations.append(dict(label=label,family=family,parameters=x.tolist(),loss=rr['loss']))
            return memo[key]
        opt=least_squares(lambda xx:current(xx)['residuals'],initial[free],jac=lambda xx:current(xx)['jacobian'][:,free],bounds=(lo[free],hi[free]),x_scale=scale[free],max_nfev=120,ftol=1e-7,xtol=1e-7,gtol=1e-6)
        x=fixed.copy();x[free]=opt.x;rr=e.predict(x);shape=shape_metrics(e.calls,rr['calls'])
        rec=dict(label=label,family=family,parameters=x.tolist(),response=parameters(family,x),loss=rr['loss'],shape=shape,shape_score=shape_score(shape,market_shape),success=bool(opt.success),message=opt.message,nfev=opt.nfev,optimality=float(opt.optimality))
        ledger.append(rec);dump(out/'optimizer_ledger.json',ledger);dump(out/'evaluations.json',evaluations);print('FIT',label,'loss',rr['loss'],'response',rec['response'],'nfev',opt.nfev,flush=True);return x,rec
    x,rec=fit('product','curve_only',candidates['baseline'][1],np.arange(7));candidates['curve_only']=('product',x)
    old=np.array(previous['parameters']['factor_curve']);candidates['previous_single_factor']=('product',np.r_[old[:7],old[7],0.,0.,old[8]])
    for family,starts in [('product',PRODUCT_STARTS),('additive',ADDITIVE_STARTS)]:
        for i,start_values in enumerate(starts):
            if family=='product':
                b1,b2,w,half=start_values;theta=[b1,(b2-b1)/(5-b1),w,np.log(half)]
            else:b,w,half=start_values;theta=[b,w,np.log(half)]
            xx=np.r_[candidates['curve_only'][1][:7],theta];xx,rec=fit(family,f'{family}_start{i}',xx);pool.append(rec)
        eligible=[r for r in pool if r['family']==family];winner=min(eligible,key=lambda r:r['loss']);candidates[family+'_price']=(family,np.array(winner['parameters']))
        shape_eligible=[r for r in eligible if r['shape_score'] is not None]
        if shape_eligible:
            winner=min(shape_eligible,key=lambda r:r['shape_score']);candidates[family+'_shape']=(family,np.array(winner['parameters']))
    frozen={name:dict(family=family,parameters=x.tolist(),response=parameters(family,x)) for name,(family,x) in candidates.items()}
    dump(out/'selection.json',dict(candidates=frozen,frozen_before_check=True,market_out_of_sample=False,adopted=False))
    train=[record(engines[family],x,out/'training',name,market_shape) for name,(family,x) in candidates.items()];dump(out/'training/scores.json',train)
    curve=pd.DataFrame(dict(start_day=np.r_[0,ends[:-1]]*365,end_day=ends*365,active_in_VIX_fit=np.isin(np.arange(16),ACTIVE)))
    for name,(family,x) in candidates.items():
        pp=p.copy();pp[7+ACTIVE]=x[:7];curve[name]=pp[7:-1]
        dump(out/f'experimental_model_{name}.json',dict(family=family,readout=parameters(family,x),parameters=pp.tolist(),ends_years=ends.tolist(),architecture=d['architecture'],premium=d['premium'],initial_state=d['initial_state'],adopted=False))
    curve.to_csv(out/'curve_levels.csv',index=False)
    a,_,premium,_=frozen_stage3();assert json.loads(json.dumps(asdict(a)))==d['architecture'];assert json.loads(json.dumps(asdict(premium)))==d['premium'];dated=DatedM2Pricer(a,premium);initial=np.array(d['initial_state'])
    cache=fresh_inputs(out,dated,initial);library=out/'postcap.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'postcap_vov_kernel.cpp'),'-o',str(library)],check=True);k=PostKernel(library,ends)
    ctx=object.__new__(ControlledEngine);ctx.dated=dated;ctx.initial=initial;ctx.k=k;ctx.ends=ends;ctx.tau=1/max(dated.layout.rates);ctx.horizon=183/365;ctx.pool=ThreadPoolExecutor(max_workers=4);ctx.normal_laws={};ctx.last_readout=None;ctx.normal_cached=None
    try:
        norm,na=ctx.normalizer(p);dump(out/'normalizer.json',na);_,fresh,_,nn=exposures(ctx,k,p,norm,cache,specs=[('baseline',0.,1.,1.,1.)])
        fn={f'd{day}_{key}':v for day,r in nn.items() for key,v in r.items()};np.savez_compressed(out/'fresh_conditional_nodes.npz',**fn)
        fe={family:MixedFit(fn,fresh,market,forwards,controls,p,ends,family=family) for family in engines}
        scores=[record(fe[family],x,out/'fresh_check',name,market_shape) for name,(family,x) in candidates.items()];dump(out/'fresh_check/scores.json',scores)
        comparisons=[paired(fe['product'],out/'fresh_check',name,ref) for name in candidates for ref in ['baseline','curve_only','previous_single_factor'] if name!=ref];dump(out/'fresh_check/paired_losses.json',comparisons)
        # Select probe states using exposure coordinates AND large weighted payoffs.
        probeprod=out/'quadrature_cache/production';probeprod.mkdir(parents=True);ids_by_day={}
        for day,aa in fresh.items():
            ids=list(np.r_[0,1,np.argmax(abs(aa['states']),axis=0),np.argmin(aa['z_normals']),np.argmax(aa['z_normals'])])
            for name,(family,x) in candidates.items():
                with np.load(out/'fresh_check'/f'paths_{name}.npz') as paths:
                    v=paths[f'd{day}_V'];w=paths[f'd{day}_weights'];ids.extend(np.argsort(w*np.maximum(v-market[market.expiry_days.eq(day)].strike.max(),0))[-3:]);ids.append(np.argmax(v))
            ids=np.unique(ids);ids_by_day[day]=ids;np.savez_compressed(probeprod/f'vix_d{day}.npz',**{key:val[ids] for key,val in aa.items()})
        _,qa,_,qn=exposures(ctx,k,p,norm,probeprod.parent,specs=[('baseline',0.,1.,1.,1.)],time_order=16,order=48)
        coarse={};fine={}
        for day,ids in ids_by_day.items():
            for key in ['times','weights']:coarse[f'd{day}_{key}']=fn[f'd{day}_{key}']
            coarse[f'd{day}_ratios']=fn[f'd{day}_ratios'][:,ids]
            for key,value in qn[str(day)].items():fine[f'd{day}_{key}']=value
        ce={family:MixedFit(coarse,qa,market,forwards,controls,p,ends,family=family) for family in engines};qe={family:MixedFit(fine,qa,market,forwards,controls,p,ends,family=family) for family in engines};qc=[]
        for name,(family,x) in candidates.items():
            cr=ce[family].predict(x,True);qr=qe[family].predict(x,True)
            for day in DAYS:qc.append(dict(candidate=name,day=day,maximum_VIX_difference=float(np.max(abs(cr['raw'][day]['V']-qr['raw'][day]['V']))),probe_count=len(ids_by_day[day]),coarse='GH32/time8',fine='GH48/time16'))
        dump(out/'quadrature_checks.json',qc)
    finally:ctx.close()
    check_sources(captured);dump(out/'run.json',design|dict(status='completed',training_scores=train,fresh_scores=scores,paired_losses=comparisons,wall_seconds=time.perf_counter()-start,limitations=['Bounded finite multistart search, not a global optimum.','Same market quotes for fresh numerical checking, not market out-of-sample.','One added independent OU direction; old readout/floor/banks fixed.','Sampled quadrature refinement and scramble SEs are not population convergence certificates.','No SPX assessment or baseline replacement.']))
    print('COMPLETED',time.perf_counter()-start,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    try:run(args.output.resolve())
    except Exception as error:
        if args.output.exists():dump(args.output/'failure.json',dict(type=type(error).__name__,message=str(error)))
        raise
