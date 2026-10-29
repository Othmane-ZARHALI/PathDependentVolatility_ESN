"""VIX-only calibration with seven independent active pieces of the 16-level curve."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from .full_readout_fit import HERE, dump
from .m2_spx import ROOT, M2, frozen_stage3
from .parity_adaptive import canonical_market, gaussian_controls
from .cap_controlled_core import shape_metrics
from .controlled_readout_fast_fit import ControlledEngine
from .postcap_vov_core import PostKernel, exposures
from .provenance import capture_sources, save_sources, check_sources, source_hashes
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture

DAYS=(27,55,118,153)
ACTIVE=np.array([1,2,3,4,7,8,10])
PRIOR=HERE/'runs/m2_maturity_curve_fit_001'
POST=HERE/'runs/m2_postcap_vov_001'
ARCHIVE=POST/'source/matteo/model_comparison/option_pricing/runs'
SEED=2026102691
DELTA=30/365


def load_draws(cache):
    return {day:dict(np.load(cache/'production'/f'vix_d{day}.npz')) for day in DAYS}


class VIXFit:
    def __init__(self, nodes, draws, market, forwards, controls, p, ends):
        self.market=market.reset_index(drop=True);self.calls=canonical_market(self.market)
        self.forwards=forwards;self.p=p;self.ends=ends;self.data={}
        for day,g in self.market.groupby('expiry_days',sort=True):
            a=draws[day];times=nodes[f'd{day}_times'];tw=nodes[f'd{day}_weights']
            seg=np.minimum(np.searchsorted(ends,times,side='right'),len(ends)-1)
            c=gaussian_controls(a['normals'],a['weights']);fixed=controls[str(day)]
            b=np.stack([np.interp(np.log(g.strike),np.log(fixed['anchors']),row) for row in fixed['beta']])
            self.data[day]=dict(a=a, times=times, seg=seg,
                weighted_ratios=nodes[f'd{day}_ratios']*tw[:,None],
                correction=c@b, fc=c@np.array(fixed['forward_beta']), q=g)

    def predict(self,x,detail=False):
        levels=self.p[7:-1].copy();levels[ACTIVE]=x[:7];gamma=x[7];half=np.exp(x[8]);rate=np.log(2)*365/half
        prices=[];calls=[];pj=[];fwd=[];fj=[];raw={};moments=[]
        for day,d in self.data.items():
            T=day/365;lag=d['times']-T;q=-np.expm1(-2*rate*T)
            base=np.exp(-rate*lag)*np.sqrt(q);b=gamma*base[:,None]
            eps=d['a']['z_normals'][None,:];L=np.exp(b*eps-.5*b*b)
            W=d['weighted_ratios']*L;coeff=levels[d['seg']]-self.p[-1]
            integral=self.p[-1]*DELTA+coeff@W
            deriv=np.zeros((len(integral),9))
            for j,seg in enumerate(ACTIVE):deriv[:,j]=W[d['seg']==seg].sum(0)
            diff=eps-b
            deriv[:,7]=np.sum(W*coeff[:,None]*diff*base[:,None],axis=0)
            dlogbase=rate*lag-rate*T*np.exp(-2*rate*T)/q
            deriv[:,8]=np.sum(W*coeff[:,None]*diff*b*dlogbase[:,None],axis=0)
            V=100*np.sqrt(integral/DELTA);dv=50/np.sqrt(DELTA*integral)[:,None]*deriv
            a=d['a'];w=a['weights'];g=d['q'];K=g.strike.to_numpy();put=g.option_type.eq('put').to_numpy()
            yf=w*V-d['fc'];yc=w[:,None]*np.maximum(V[:,None]-K,0)-d['correction']
            yn=yc-put[None,:]*(yf[:,None]-K)
            jac=np.einsum('nq,nj->qj',w[:,None]*(V[:,None]>K),dv)/len(V)
            f_jac=np.mean(w[:,None]*dv,axis=0);jac-=put[:,None]*f_jac
            prices.extend(yn.mean(0));calls.extend(yc.mean(0));pj.extend(jac);fwd.append(yf.mean());fj.append(f_jac)
            if detail:
                raw[day]=dict(V=V,weights=w,native=yn,call=yc,forward=yf,nets=a['nets'])
                exact=10000/30*np.sum(np.maximum(0,np.minimum(self.ends*365,day+30)-np.maximum(np.r_[0,self.ends[:-1]]*365,day))*levels)
                block=np.array([np.mean((w*V**2)[a['nets']==j]) for j in range(8)])
                tail=w*np.maximum(V-g.strike.max(),0);total=tail.sum();top=np.sort(tail)[::-1]
                moments.append(dict(day=int(day),exact_second_moment=float(exact),sample_second_moment=float(block.mean()),second_moment_se=float(block.std(ddof=1)/np.sqrt(8)),mass=float(w.mean()),max_VIX=float(V.max()),last_strike=float(g.strike.max()),last_strike_raw_call=float(tail.mean()),last_strike_payoff_ESS=float(total**2/np.sum(tail**2)) if total else 0.,last_strike_top10_fraction=float(top[:10].sum()/total) if total else 0.))
        prices=np.array(prices);pj=np.array(pj);fwd=np.array(fwd);fj=np.array(fj);m=self.market
        dist=np.where(prices<m.price_bid_U,prices-m.price_bid_U,np.where(prices>m.price_ask_U,prices-m.price_ask_U,0.))
        weights=np.concatenate([np.full(len(g),1/(4*len(g))) for _,g in m.groupby('expiry_days',sort=True)])
        r=np.r_[np.sqrt(weights)*dist/m.price_scale_U,(fwd-self.forwards.forward_vix_points.to_numpy())/(.5*2)]
        outside=(prices<m.price_bid_U)|(prices>m.price_ask_U)
        jac=np.concatenate([np.sqrt(weights)[:,None]*outside.to_numpy()[:,None]*pj/m.price_scale_U.to_numpy()[:,None],fj/(.5*2)])
        return dict(residuals=r,jacobian=jac,loss=float(r@r),prices=prices,calls=np.array(calls),forwards=fwd,raw=raw,moments=moments)


def loss_from_prices(e,price,fwd):
    q=e.market;d=np.where(price<q.price_bid_U,price-q.price_bid_U,np.where(price>q.price_ask_U,price-q.price_ask_U,0.))/q.price_scale_U.to_numpy()
    return float(np.mean([np.mean(d[g.index]**2) for _,g in q.groupby('expiry_days')])+np.mean(((fwd-e.forwards.forward_vix_points.to_numpy())/.5)**2))


def record(e,x,folder,name):
    r=e.predict(x,True);folder.mkdir(parents=True,exist_ok=True)
    native=np.concatenate([d['native'] for d in r['raw'].values()],axis=1)
    call=np.concatenate([d['call'] for d in r['raw'].values()],axis=1)
    fw=np.column_stack([d['forward'] for d in r['raw'].values()])
    bn=native.reshape(8,-1,native.shape[1]).mean(1);bc=call.reshape(8,-1,call.shape[1]).mean(1);bf=fw.reshape(8,-1,4).mean(1)
    for typ,q,val,block in [('otm',e.market,r['prices'],bn),('canonical',e.calls,r['calls'],bc)]:
        q=q.copy();q['model_price_U']=val;q['mc_se_U']=block.std(0,ddof=1)/np.sqrt(8)
        q['interval_residual']=np.where(val<q.price_bid_U,val-q.price_bid_U,np.where(val>q.price_ask_U,val-q.price_ask_U,0.))/q.price_scale_U
        q.to_csv(folder/f'prices_{typ}_{name}.csv',index=False)
    ft=e.forwards.copy();ft['model_forward']=r['forwards'];ft['mc_se']=bf.std(0,ddof=1)/np.sqrt(8);ft.to_csv(folder/f'forwards_{name}.csv',index=False)
    np.savez_compressed(folder/f'blocks_{name}.npz',native=bn,call=bc,forward=bf)
    np.savez_compressed(folder/f'paths_{name}.npz',**{f'd{day}_{k}':v for day,d in r['raw'].items() for k,v in d.items()})
    shape=shape_metrics(e.calls,r['calls']);summary=dict(name=name,parameters=x.tolist(),gamma=float(x[7]),half_life_days=float(np.exp(x[8])),loss=r['loss'],inside_bidask=int(np.sum((r['prices']>=e.market.price_bid_U)&(r['prices']<=e.market.price_ask_U))),negative_prices=int(np.sum(r['prices']<0)),shape=shape,moments=r['moments'])
    dump(folder/f'score_{name}.json',summary);return summary


def paired(e,folder,name,reference='baseline'):
    with np.load(folder/f'blocks_{name}.npz') as a, np.load(folder/f'blocks_{reference}.npz') as b:
        full=loss_from_prices(e,a['native'].mean(0),a['forward'].mean(0))-loss_from_prices(e,b['native'].mean(0),b['forward'].mean(0))
        jk=[]
        for j in range(8):
            ids=np.arange(8)!=j
            jk.append(loss_from_prices(e,a['native'][ids].mean(0),a['forward'][ids].mean(0))-loss_from_prices(e,b['native'][ids].mean(0),b['forward'][ids].mean(0)))
    return dict(candidate=name,reference=reference,loss_change=full,paired_jackknife_se=float(np.std(jk,ddof=0)*np.sqrt(7)),leave_one_net_changes=jk)


def run(out):
    out.mkdir(parents=True,exist_ok=False);start=time.perf_counter()
    d=json.loads((PRIOR/'run.json').read_text());p=np.array(d['selected']);ends=np.array(d['ends_years'])
    market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv').sort_values(['expiry_days','strike']).reset_index(drop=True)
    forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv').sort_values('expiry_days').reset_index(drop=True)
    controls=json.loads((PRIOR/'vix_controls.json').read_text())
    # Stored rows are in the original source order. Verify sorted order is unchanged.
    pd.testing.assert_frame_equal(market,pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv'))
    nodes=dict(np.load(POST/'fresh_check/base_conditional_nodes.npz'));draws=load_draws(POST/'fresh_cache')
    sources=[Path(__file__),HERE/'VIX_ONLY_CURVE_PLAN.md',PRIOR/'run.json',PRIOR/'selected_model.json',PRIOR/'vix_controls.json',ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv',ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv',ARCHIVE/'m2_vix_targeted_check_001/proposals.json',HERE/'postcap_vov_core.py',HERE/'postcap_vov_kernel.cpp',HERE/'stochastic_vov_kernel.cpp',HERE/'cap_tail_kernel.cpp',HERE/'maturity_curve.py',HERE/'controlled_readout_fast_fit.py',HERE/'parity_adaptive.py',HERE/'cap_controlled_core.py',HERE/'plot_maturity_curve.py',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(sources);save_sources(out,captured,ROOT)
    inputs=[POST/'fresh_check/base_conditional_nodes.npz',*(POST/'fresh_cache/production'/f'vix_d{day}.npz' for day in DAYS)]
    design=dict(status='running',source_sha256=source_hashes(captured,ROOT),input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},active_curve_indices_zero_based=ACTIVE.tolist(),curve_piece_count=16,floor=float(p[-1]),baseline_replaced=False,SPX_priced=False,fresh_check_seed_base=SEED,training_note='Prior checking draws deliberately reused as training; new draws generated after selection',objective='native OTM bid-ask scaled distances, equal expiry weights, plus forward residuals /0.5; no penalties')
    dump(out/'design.json',design);e=VIXFit(nodes,draws,market,forwards,controls,p,ends)
    x0=np.r_[p[7:-1][ACTIVE],0.,np.log(7.)];replay=e.predict(x0)
    saved=pd.read_csv(POST/'fresh_check/prices_otm_baseline.csv')
    err=float(np.max(abs(saved.model_price_U.to_numpy()-replay['prices'])));assert err<2e-8,err
    dump(out/'baseline_replay.json',dict(maximum_native_price_difference=err,loss=replay['loss']))
    probe=x0.copy();probe[7]=.7;probe[8]=np.log(19.)
    r=e.predict(probe);checks=[]
    for j in range(9):
        h=1e-6 if j<7 else 1e-5;dx=np.zeros(9);dx[j]=h
        rp=e.predict(probe+dx);rm=e.predict(probe-dx)
        # Check smooth price and forward derivatives, away from payoff kinks at this step.
        fd=(rp['residuals']-rm['residuals'])/(2*h)
        error=float(np.linalg.norm(fd-r['jacobian'][:,j])/max(1,np.linalg.norm(fd)))
        checks.append(dict(parameter=j,relative_error=error));assert error<.002,(j,error)
    dump(out/'jacobian_check.json',checks)
    lo=np.r_[np.full(7,p[-1]),0.,np.log(1.)];hi=np.r_[np.full(7,.2),3.,np.log(365.)]
    candidates={'baseline':x0.copy()};ledger=[];evaluations=[]
    def fit(label,initial,free):
        memo={};fixed=initial.copy()
        def current(xx):
            key=tuple(xx)
            if key not in memo:
                x=fixed.copy();x[free]=xx;rr=e.predict(x);memo.clear();memo[key]=rr
                evaluations.append(dict(label=label,parameters=x.tolist(),loss=rr['loss']))
            return memo[key]
        opt=least_squares(lambda xx:current(xx)['residuals'],initial[free],jac=lambda xx:current(xx)['jacobian'][:,free],bounds=(lo[free],hi[free]),x_scale=np.r_[np.full(7,.06),1.,2.][free],max_nfev=120,ftol=1e-7,xtol=1e-7,gtol=1e-6)
        x=fixed.copy();x[free]=opt.x;rr=e.predict(x)
        rec=dict(label=label,parameters=x.tolist(),loss=rr['loss'],success=bool(opt.success),message=opt.message,nfev=opt.nfev,optimality=float(opt.optimality))
        ledger.append(rec);dump(out/'optimizer_ledger.json',ledger);dump(out/'evaluations.json',evaluations)
        print('FIT',label,'loss',rr['loss'],'gamma',x[7],'half',np.exp(x[8]),'nfev',opt.nfev,flush=True);return x
    candidates['curve_only']=fit('curve_only',x0,np.arange(7))
    strong=x0.copy();strong[7]=1.;candidates['fixed_factor_curve']=fit('fixed_factor_curve',strong,np.arange(7))
    fits=[]
    for i,(gamma,half) in enumerate([(.25,7),(1.,7),(1.,30),(1.,90),(2.,7),(2.,90)]):
        seed=candidates['curve_only'].copy();seed[7:]=[gamma,np.log(half)]
        fits.append(fit(f'joint_start{i}',seed,np.arange(9)))
    candidates['factor_curve']=min(fits,key=lambda x:e.predict(x)['loss'])
    dump(out/'selection.json',dict(parameters={name:x.tolist() for name,x in candidates.items()},frozen_before_fresh_check=True,automatic_adoption=False))
    train=[record(e,x,out/'training',name) for name,x in candidates.items()];dump(out/'training/scores.json',train)
    # Preserve all 16 levels in each experimental model; off-window pieces are unchanged.
    left=np.r_[0.,ends[:-1]]*365;table=pd.DataFrame(dict(start_day=left,end_day=ends*365,active_in_VIX_fit=np.isin(np.arange(16),ACTIVE),baseline=p[7:-1]))
    for name,x in candidates.items():
        pp=p.copy();pp[7+ACTIVE]=x[:7];table[name]=pp[7:-1]
        dump(out/f'experimental_model_{name}.json',dict(parameters=pp.tolist(),ends_years=ends.tolist(),gamma=float(x[7]),half_life_days=float(np.exp(x[8])),architecture=d['architecture'],premium=d['premium'],initial_state=d['initial_state'],adopted=False))
    table.to_csv(out/'curve_levels.csv',index=False)
    # New joint37D streams are generated only now, with immutable candidate parameters.
    a,_,premium,_=frozen_stage3();assert json.loads(json.dumps(asdict(a)))==d['architecture'];assert json.loads(json.dumps(asdict(premium)))==d['premium']
    dated=DatedM2Pricer(a,premium);initial=np.array(d['initial_state']);proposals=json.loads((ARCHIVE/'m2_vix_targeted_check_001/proposals.json').read_text())
    cache=out/'fresh_cache';prod=cache/'production';prod.mkdir(parents=True);seeds=[]
    for day in DAYS:
        dec=proposals[str(day)];mix=DefensiveGaussianMixture(np.c_[dec['shifts'],np.zeros(len(dec['shifts']))],dec['probabilities']);sampler=GaussianOuterSampler(dated,day/365,initial)
        zz=[];ww=[];cc=[]
        for net in range(8):
            seed=int(np.random.SeedSequence([SEED,day,net]).generate_state(1,dtype=np.uint64)[0]);z,w,c=mix.draw(10,seed);zz.append(z);ww.append(w);cc.append(c);seeds.append(dict(day=day,net=net,seed=seed))
        z=np.concatenate(zz);np.savez_compressed(prod/f'vix_d{day}.npz',states=sampler.mean+z[:,:-1]@sampler.root.T,normals=z[:,:-1],z_normals=z[:,-1],weights=np.concatenate(ww),components=np.concatenate(cc),nets=np.repeat(np.arange(8),1024))
    dump(cache/'design.json',dict(seeds=seeds,seed_base=SEED,created_after_selection=True,points_per_expiry=8192,proposals=proposals))
    k=PostKernel(POST/'postcap.dylib',ends)
    # VIX-only Gaussian context avoids constructing or reading any SPX paths.
    ctx=object.__new__(ControlledEngine);ctx.dated=dated;ctx.initial=initial;ctx.k=k;ctx.ends=ends;ctx.tau=1/max(dated.layout.rates);ctx.horizon=183/365;ctx.pool=ThreadPoolExecutor(max_workers=4);ctx.normal_laws={};ctx.last_readout=None;ctx.normal_cached=None
    try:
        norm,na=ctx.normalizer(p);dump(out/'normalizer.json',na)
        _,fresh,_,nn=exposures(ctx,k,p,norm,cache,specs=[('baseline',0.,1.,1.,1.)])
        freshnodes={f'd{day}_{key}':v for day,rn in nn.items() for key,v in rn.items()}
        np.savez_compressed(out/'fresh_conditional_nodes.npz',**freshnodes)
        fe=VIXFit(freshnodes,fresh,market,forwards,controls,p,ends)
        scores=[record(fe,x,out/'fresh_check',name) for name,x in candidates.items()];dump(out/'fresh_check/scores.json',scores)
        comparisons=[paired(fe,out/'fresh_check',name,reference) for name,reference in [('curve_only','baseline'),('fixed_factor_curve','baseline'),('factor_curve','baseline'),('factor_curve','curve_only')]]
        dump(out/'fresh_check/paired_losses.json',comparisons)
        probeprod=out/'quadrature_cache/production';probeprod.mkdir(parents=True);ids_by_day={}
        for day,aa in fresh.items():
            ids=np.unique(np.r_[0,1,np.argmax(abs(aa['states']),axis=0),np.argmin(aa['z_normals']),np.argmax(aa['z_normals'])]);ids_by_day[day]=ids
            np.savez_compressed(probeprod/f'vix_d{day}.npz',**{key:val[ids] for key,val in aa.items()})
        _,qa,_,qn=exposures(ctx,k,p,norm,probeprod.parent,specs=[('baseline',0.,1.,1.,1.)],time_order=16,order=48)
        coarse={};fine={}
        for day,ids in ids_by_day.items():
            for key in ['times','weights']:coarse[f'd{day}_{key}']=freshnodes[f'd{day}_{key}']
            coarse[f'd{day}_ratios']=freshnodes[f'd{day}_ratios'][:,ids]
            for key,v in qn[str(day)].items():fine[f'd{day}_{key}']=v
        ce=VIXFit(coarse,qa,market,forwards,controls,p,ends);qe=VIXFit(fine,qa,market,forwards,controls,p,ends);qc=[]
        for name,x in candidates.items():
            cr=ce.predict(x,True);qr=qe.predict(x,True)
            for day in DAYS:qc.append(dict(candidate=name,day=day,maximum_VIX_difference=float(np.max(abs(cr['raw'][day]['V']-qr['raw'][day]['V']))),probe_count=len(ids_by_day[day]),coarse='GH32/time8',fine='GH48/time16'))
        dump(out/'quadrature_checks.json',qc)
    finally:ctx.close()
    check_sources(captured)
    dump(out/'run.json',design|dict(status='completed',training_scores=train,fresh_scores=scores,paired_losses=comparisons,wall_seconds=time.perf_counter()-start,limitations=['Finite multistart bounded fit; no global optimum claim.','Reused training simulation; fresh numerical check uses the same market observations.','Old readout and floor fixed; seven independent VIX-sensitive curve pieces.','No SPX assessment or production model replacement.','Sampled deterministic refinement and eight-net SE do not establish population convergence.']))
    print('COMPLETED',time.perf_counter()-start,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    try:run(args.output.resolve())
    except Exception as error:
        if args.output.exists():dump(args.output/'failure.json',dict(type=type(error).__name__,message=str(error)))
        raise
