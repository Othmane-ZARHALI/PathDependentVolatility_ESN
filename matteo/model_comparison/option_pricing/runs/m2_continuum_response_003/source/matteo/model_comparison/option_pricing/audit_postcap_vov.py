"""Independent saved-price, conditional-factor, seed and source replay."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .full_readout_fit import HERE,arr,dump
from .m2_spx import ROOT,frozen_stage3
from .postcap_vov_core import conditional_multiplier
from .postcap_vov_study import native_loss,shape_score
from .audit_maturity_curve import replay
from .audit_cap_tail import paired
from .stochastic_vov_core import seed_for
from .cap_controlled_core import shape_metrics
from esn_eq.pricing import DatedM2Pricer,_normal_rule
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.states import symmetric_square_root
from esn_eq.importance_sampling import DefensiveGaussianMixture


def audit(run):
    d=json.loads((run/'design.json').read_text());result=json.loads((run/'run.json').read_text());p=np.array(d['baseline']);parent=Path(d['parent']);old=json.loads((parent/'run.json').read_text());cache=Path(d['cache'])
    vix=pd.read_csv(cache/'vix_contracts.csv');spx=pd.read_csv(cache/'spx_contracts.csv');q=pd.concat([spx,vix],ignore_index=True);f=pd.read_csv(cache/'forward_targets.csv')
    selection=json.loads((run/'selection.json').read_text());names=selection['full_candidates'];profiles=json.loads((run/'profiles.json').read_text());market=json.loads((run/'market_shape.json').read_text());records=[];negative=[]
    for folder,rows,quotes in [(run,[r['name'] for r in profiles],vix),(run/'full',names,q),(run/'fresh_check',names,vix)]:
        for name in rows:
            rr=replay(folder,name,p[7:-1],quotes);records.append(dict(folder=str(folder.relative_to(run)),name=name,native_loss=rr['original_otm_total']))
            n=pd.read_csv(folder/f'prices_otm_{name}.csv');bad=n[n.model_price_U.lt(0)]
            negative.append(dict(folder=str(folder.relative_to(run)),name=name,count=len(bad),maximum_magnitude_in_SE=float(np.max(-bad.model_price_U/bad.mc_se_U)) if len(bad) else 0))
            if folder!=run:dump(folder/f'paired_{name}_baseline.json',paired(folder,name,folder,'baseline',quotes,f,p[7:-1],p[7:-1]))
    # Reconstruct each VIX value from archived baseline conditional nodes, then
    # independently integrate conditional OU exponentials on a dense real grid.
    factor_checks=[];x=np.linspace(-10,10,200001);density=np.exp(-x*x/2)/np.sqrt(2*np.pi)
    for folder,rows,drawcache in [(run,[r['name'] for r in profiles],cache),(run/'fresh_check',names,run/'fresh_cache')]:
        nodes=np.load(folder/'base_conditional_nodes.npz')
        for name in rows:
            spec=next(s for s in d['specs'] if s[0]==name);_,gamma,rate,_,_=spec
            with np.load(folder/f'paths_{name}.npz') as paths:
                for day in [27,55,118,153]:
                    if folder==run:
                        from .stochastic_vov_core import outer_z
                        eps,_=outer_z(cache,day)
                    else:
                        with np.load(drawcache/'production'/f'vix_d{day}.npz') as a:eps=a['z_normals']
                    zT=np.sqrt(-np.expm1(-2*rate*day/365))*eps;times=nodes[f'd{day}_times'];weights=nodes[f'd{day}_weights'];ratios=nodes[f'd{day}_ratios'];integ=np.full(len(eps),p[-1]*30/365)
                    for t,w,ratio in zip(times,weights,ratios):
                        lag=t-day/365;decay=np.exp(-rate*lag);qT=-np.expm1(-2*rate*day/365)
                        # Formula written independently of the production helper.
                        factor=np.exp(gamma*decay*zT-.5*gamma**2*decay**2*qT)
                        j=min(np.searchsorted(old['ends_years'],t,side='right'),len(old['ends_years'])-1)
                        integ+=w*(p[7+j]-p[-1])*ratio*factor
                    values=100*np.sqrt(integ/(30/365));np.testing.assert_allclose(values,paths[f'vix_{day}'],rtol=5e-14,atol=1e-12)
                    if day==153 and gamma:
                        T=day/365;lag=.02;z=float(zT[np.argmax(zT)]);mu=np.exp(-rate*lag)*z;sd=np.sqrt(-np.expm1(-2*rate*lag));qU=-np.expm1(-2*rate*(T+lag))
                        direct=np.trapezoid(np.exp(gamma*(mu+sd*x)-.5*gamma**2*qU)*density,x);fast=float(conditional_multiplier(z,T,lag,gamma,rate));error=abs(fast/direct-1);assert error<2e-12
                        factor_checks.append(dict(folder=str(folder.relative_to(run)),name=name,relative_error=error))
    eligible=[r for r in profiles if r['name']!='baseline' and r['shape_score'] is not None]
    assert min(eligible,key=lambda r:(native_loss(r['score']),r['shape_score']))['name']==selection['price_winner']
    assert min(eligible,key=lambda r:(r['shape_score'],native_loss(r['score'])))['name']==selection['shape_winner']
    for r in profiles:
        table=pd.read_csv(run/f"prices_canonical_{r['name']}.csv");sh=shape_metrics(table,table.model_price_U.to_numpy());np.testing.assert_allclose(shape_score(sh,market),r['shape_score'],rtol=1e-11)
    # Paired delete-one-scramble sensitivity of the plotted wing regressions.
    fresh_shapes={};tables={};bs={}
    for name in names:
        table=pd.read_csv(run/'fresh_check'/f'prices_canonical_{name}.csv');tables[name]=table
        with np.load(run/'fresh_check'/f'paths_{name}.npz') as paths:
            bs[name]=np.concatenate([paths[f'vix_blocks_{day}'] for day in [27,55,118,153]],axis=1)
        estimates=[]
        for omit in [None,*range(8)]:
            keep=np.ones(8,dtype=bool) if omit is None else np.arange(8)!=omit
            metric=shape_metrics(table,bs[name][keep].mean(0))
            estimates.append([[r['wing_slope'],r['wing_curvature']] for r in metric['rows']])
        fresh_shapes[name]=np.array(estimates,dtype=float)
    shape_paired=[]
    for name in names:
        changes=fresh_shapes[name]-fresh_shapes['baseline'];se=np.sqrt(7/8*np.sum((changes[1:]-changes[1:].mean(0))**2,axis=0))
        shape_paired.append(dict(name=name,days=[27,55,118,153],order=['wing_slope','wing_curvature'],change=changes[0].tolist(),jackknife_SE=se.tolist(),leave_one_scramble_changes=changes[1:].tolist()))
    dump(run/'fresh_check/paired_wing_shape.json',dict(records=shape_paired,limitation='Eight-scramble descriptive sampling sensitivity; excludes quadrature, pilot and selection uncertainty.'))
    # Exact joint37D draw regeneration and likelihood-weight verification.
    a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial=np.array(old['initial_state']);fd=json.loads((run/'fresh_cache/design.json').read_text());allseeds=[]
    for record in fd['files']:
        day=record['day'];dec=fd['proposals'][str(day)];mix=DefensiveGaussianMixture(np.c_[dec['shifts'],np.zeros(len(dec['shifts']))],dec['probabilities']);sampler=GaussianOuterSampler(dated,day/365,initial)
        path=run/'fresh_cache/production'/f'vix_d{day}.npz';assert hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256']
        with np.load(path) as saved:
            for row in record['seeds']:
                seed=row['seed'];allseeds.append(seed);assert seed==int(np.random.SeedSequence([2026102631,day,row['net']]).generate_state(1,dtype=np.uint64)[0]);z,w,c=mix.draw(9,seed);sl=slice(row['net']*512,(row['net']+1)*512)
                for actual,expected in [(z[:,:-1],saved['normals'][sl]),(z[:,-1],saved['z_normals'][sl]),(w,saved['weights'][sl]),(c,saved['components'][sl]),(sampler.mean+z[:,:-1]@sampler.root.T,saved['states'][sl]),(mix.weights(z),saved['weights'][sl])]:np.testing.assert_array_equal(actual,expected)
    assert len(set(allseeds))==32
    historical=set()
    for oldcache in [cache,Path(old['check_cache'])]:
        for day in [27,55,118,153]:
            historical.update(row['seed'] for row in json.loads((oldcache/'production'/f'vix_d{day}_seeds.json').read_text()))
    for rows in json.loads((run/'training_Z_seeds.json').read_text()).values():historical.update(row['seed'] for row in rows)
    previous=json.loads((HERE/'runs/m2_stochastic_vov_001/fresh_cache/design.json').read_text())
    for rec in previous['files']:historical.update(row['seed'] for row in rec['seeds'])
    assert not historical&set(allseeds)
    selection_time=(run/'selection.json').stat().st_mtime_ns;assert selection_time<(run/'fresh_cache/design.json').stat().st_mtime_ns
    assert all(selection_time<(run/'full'/f'score_{name}.json').stat().st_mtime_ns for name in names)
    for rel,sha in d['source_sha256'].items():
        assert hashlib.sha256((run/'source'/rel).read_bytes()).hexdigest()==sha,rel
        assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==sha,rel
    # Direct NumPy original-state conditional integral, independent of C++.
    node_checks=[]
    from .postcap_vov_core import PostKernel
    k=PostKernel(run/'postcap.dylib',old['ends_years'])
    with np.load(run/'fresh_cache/production/vix_d153.npz') as draw:
        ids=np.unique(np.r_[np.arange(32),int(np.argmax(abs(draw['states']).max(1)))])
        for lag in [.0002,.02,30/365]:
            mean,cov=dated.gaussian_projection_law(lag,draw['states'][ids]);z,w=_normal_rule(32);offset=z@symmetric_square_root(cov).T
            fast=k.expectation(arr(mean),arr(offset),arr(w),p)[:,0]
            states=mean[:,None,:]+offset[None,:,:];eta=p[0]*states[:,:,0]+p[1]*(states[:,:,1]-p[2])**2+p[3]*(states[:,:,2]-p[4])**2
            h=p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6];direct=np.exp(2*h)@w;err=float(np.max(abs(fast/direct-1)));assert err<2e-7
            node_checks.append(dict(lag=lag,ids=ids.tolist(),relative_error=err))
    normal_checks=[]
    for t in [27/365,153/365,1030/365]:
        mean,cov=dated.gaussian_projection_law(t,initial);vv=[]
        for order in [48,96]:
            z,w=_normal_rule(order);offset=z@symmetric_square_root(cov).T
            vv.append(float(k.expectation(arr(mean),arr(offset),arr(w),p)[0,0]))
        normal_checks.append(dict(time=t,GH48=vv[0],GH96=vv[1],relative_change=vv[1]/vv[0]-1))
    # Per-candidate declarative experimental configurations; no baseline replacement.
    configs=run/'experimental_models';configs.mkdir(exist_ok=True)
    for name in names:
        if name=='baseline':continue
        spec=next(s for s in d['specs'] if s[0]==name)
        dump(configs/f'{name}.json',dict(status='Experimental fixed-grid comparison; not adopted',baseline_model=json.loads((parent/'selected_model.json').read_text()),baseline_parameters=p.tolist(),layer=dict(gamma=spec[1],rate=spec[2],half_life_days=spec[3],Z0=0,stationary_variance=1,independent=True,formula='L=exp(gamma*Z-gamma^2*(1-exp(-2*a*t))/2)'),variance='f+(xi-f)*(G_base/M_base)*L',normalizer='Unchanged original M; exact independence identity',outer_cap='Unchanged original cap on G_base only',P_calibrated=False))
    stock=json.loads((run/'full/stock_Z_audit.json').read_text());stockstats=[]
    for name in names:
        spec=next(s for s in d['specs'] if s[0]==name)
        for day in [29,57,119,211,302,666,1030]:
            rows=sorted([r for r in stock if r['candidate']==name and r['day']==day],key=lambda r:r['batch']);assert len(rows)==32
            for r in rows:assert r['seed']==seed_for(2,day,r['batch'])
            z=np.concatenate([r['terminal_z'] for r in rows]);stockstats.append(dict(name=name,day=day,mean=float(z.mean()),variance=float(z.var()),target_variance=float(-np.expm1(-2*spec[2]*day/365))))
    dump(run/'audit.json',dict(status='PASS',price_records=records,conditional_factor_oracles=factor_checks,original_state_NumPy_oracles=node_checks,normalizer_order_checks=normal_checks,checked_historical_seed_count=len(historical),exact_joint_scrambles_replayed=32,source_checks=len(d['source_sha256']),stock_OU=stockstats,negative_estimates=negative,selection_precedes_SPX_and_fresh=True,baseline_replaced=False,limitation='Finite implementation checks, not a population integration/sampling or full calibration certificate.'))
    print('PASS',len(records),'saved price records;',len(factor_checks),'conditional-factor integrals;32 joint scrambles',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);audit(p.parse_args().run.resolve())
