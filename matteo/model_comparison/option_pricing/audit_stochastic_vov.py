"""Saved-draw replay, joint-scramble regeneration and independent 4D oracle."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .full_readout_fit import HERE,arr,dump
from .m2_spx import ROOT,frozen_stage3
from .stochastic_vov_core import LayerKernel,SPECS,outer_z,seed_for,normal_rule,response_table
from .maturity_curve import window_nodes
from .stochastic_vov_study import load_normalizers
from .audit_maturity_curve import replay
from .audit_cap_tail import paired
from .cap_controlled_core import shape_metrics
from .controlled_readout_fast_fit import ControlledEngine
from .test_stochastic_vov import python_innovations
from esn_eq.pricing import DatedM2Pricer,_normal_rule
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.states import symmetric_square_root
from esn_eq.importance_sampling import DefensiveGaussianMixture


def oracle(dated,t,state,p,continuation,sensitivity,mu_z,sd_z,old_order,z_order):
    """Direct numpy 4D quadrature. No cap C++, convolution or Chebyshev table."""
    m,cov=dated.gaussian_projection_law(t,state);z,w=_normal_rule(old_order)
    states=np.asarray(m).reshape(1,3)+z@symmetric_square_root(cov).T
    base=p[0]*states[:,0]+p[1]*(states[:,1]-p[2])**2;spike=p[3]*(states[:,2]-p[4])**2
    gz,gw=normal_rule(z_order);answer=0.
    for zz,ww in zip(mu_z+sd_z*gz,gw):
        eta=base+np.exp(sensitivity*np.tanh(zz))*spike
        h=p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6]
        if continuation:h+=continuation*np.log1p(np.logaddexp(0,p[6]*(eta-p[5]))/p[6])
        answer+=ww*np.dot(w,np.exp(2*h))
    return float(answer)


def audit(run):
    design=json.loads((run/'design.json').read_text());r=json.loads((run/'run.json').read_text());p=np.array(design['baseline']);parent=Path(design['parent']);d=json.loads((parent/'run.json').read_text());cache=Path(design['cache'])
    spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');quotes=pd.concat([spx,vix],ignore_index=True);f=pd.read_csv(cache/'forward_targets.csv')
    profiles=json.loads((run/'profiles.json').read_text());selection=json.loads((run/'selection.json').read_text());names=selection['full_candidates'];records=[];regions=[]
    for row in profiles:
        rr=replay(run,row['name'],p[7:-1],vix);records.append(dict(role='grid',name=row['name'],native_loss=rr['original_otm_total']))
        q=pd.read_csv(run/f"prices_otm_{row['name']}.csv")
        for day,g in q.groupby('expiry_days'):
            x=np.log(g.K_over_F.to_numpy())
            for region,mask in [('near_ATM',abs(x)<=.2),('right_wing',x>np.log(2)),('other',(abs(x)>.2)&(x<=np.log(2)))]:
                v=g.interval_residual.to_numpy()[mask];regions.append(dict(candidate=row['name'],days=int(day),region=region,quotes=int(mask.sum()),native_price_loss_contribution=float(np.sum(v*v)/(4*len(g)))))
        total=sum(x['native_price_loss_contribution'] for x in regions if x['candidate']==row['name'])
        assert abs(total-row['score']['original_otm_components']['VIX'])<1e-8
    pd.DataFrame(regions).to_csv(run/'grid_region_diagnostics.csv',index=False)
    negative_estimates=[]
    for folder,label,qq in [(run/'full','SPX+VIX',quotes),(run/'fresh_check','freshVIX',vix)]:
        for name in names:
            rr=replay(folder,name,p[7:-1],qq);records.append(dict(role=label,name=name,native_loss=rr['original_otm_total']))
            native=pd.read_csv(folder/f'prices_otm_{name}.csv');bad=native[native.model_price_U.lt(0)]
            negative_estimates.append(dict(role=label,candidate=name,count=len(bad),largest_negative_magnitude_in_reported_SE=float(np.max(-bad.model_price_U/bad.mc_se_U)) if len(bad) else 0.))
            dump(folder/f'paired_{name}_baseline.json',paired(folder,name,folder,'cap0_ref',qq,f,p[7:-1],p[7:-1]))
            if not name.endswith('_ref'):
                ref=name.split('_')[0]+'_ref'
                dump(folder/f'paired_{name}_same_cap.json',paired(folder,name,folder,ref,qq,f,p[7:-1],p[7:-1]))
    # Baseline both markets must replay old saved prices after adding the new ABI.
    old=pd.read_csv(parent/'prices_canonical_selected.csv');new=pd.read_csv(run/'full/prices_canonical_cap0_ref.csv')
    maximum=float(np.max(abs(old.model_price_U-new.model_price_U)));assert maximum<2e-6,maximum
    for relative,sha in design['source_sha256'].items():
        assert hashlib.sha256((run/'source'/relative).read_bytes()).hexdigest()==sha,relative
    selected_hash=design['source_sha256'][str((parent/'selected_model.json').relative_to(ROOT))]
    assert hashlib.sha256((parent/'selected_model.json').read_bytes()).hexdigest()==selected_hash
    experimental=run/'experimental_models';experimental.mkdir(exist_ok=True)
    for name in names:
        if name.endswith('_ref'):continue
        row=next(x for x in profiles if x['name']==name);spec=row['spec']
        dump(experimental/f'{name}.json',dict(schema='m2_dedicated_ou_spike_sensitivity_experiment_v1',
             status='Fixed-grid experiment; not an adopted or calibrated production configuration',
             baseline_model=json.loads((parent/'selected_model.json').read_text()),baseline_parameters=p.tolist(),
             independent_layer=dict(initial_Z=0.,Q_mean=0.,rate_per_year=spec[2],half_life_calendar_days=spec[3],stationary_variance=1.,
                 diffusion='sqrt(2*rate)',sensitivity=spec[1],multiplier='exp(sensitivity*tanh(Z))',multiplier_range=[1/spec[4],spec[4]],independent_of_existing_drivers=True),
             outer_transform=dict(logistic_cap_level=float(p[5]),sharpness=float(p[6]),logarithmic_continuation=row['continuation']),
             equations='eta=cC*C+qY*(Y-betaY)^2+qJ*exp(s*tanh(Z))*(F-betaJ)^2; G=exp(2*h(eta)); v=f+(xi-f)*G/E_Q[G]',
             conditioning='Full original 36-dimensional state plus Z_T',normalizer='Candidate-specific Q expectation',
             implementation='stochastic_vov_core.py and stochastic_vov_kernel.cpp; use experimental LayerEngine',
             baseline_replaced=False,P_calibrated=False))
    frozen_time=(run/'selection.json').stat().st_mtime_ns
    assert frozen_time<min((run/'full'/f'score_{name}.json').stat().st_mtime_ns for name in names)
    assert frozen_time<(run/'fresh_cache/design.json').stat().st_mtime_ns
    # Independently replay all 32 newly generated joint mixture scrambles.
    a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial=np.array(d['initial_state']);fd=json.loads((run/'fresh_cache/design.json').read_text());allseeds=[]
    for rec in fd['files']:
        day=rec['days'];dec=fd['proposal'][str(day)];shifts=np.c_[np.array(dec['shifts']),np.zeros(len(dec['shifts']))]
        mix=DefensiveGaussianMixture(shifts,dec['probabilities']);sampler=GaussianOuterSampler(dated,day/365,initial)
        with np.load(run/f'fresh_cache/production/vix_d{day}.npz') as saved:
            for row in rec['seeds']:
                net=row['net'];seed=row['seed'];assert seed==seed_for(3,day,net);allseeds.append(seed)
                z,w,c=mix.draw(9,seed);sl=slice(net*512,(net+1)*512)
                np.testing.assert_array_equal(z[:,:-1],saved['normals'][sl]);np.testing.assert_array_equal(z[:,-1],saved['z_normals'][sl]);np.testing.assert_array_equal(w,saved['weights'][sl]);np.testing.assert_array_equal(c,saved['components'][sl]);np.testing.assert_array_equal(sampler.mean+z[:,:-1]@sampler.root.T,saved['states'][sl])
                np.testing.assert_allclose(mix.weights(z),saved['weights'][sl],rtol=0,atol=0)
        assert hashlib.sha256((run/f'fresh_cache/production/vix_d{day}.npz').read_bytes()).hexdigest()==rec['sha256']
    assert len(set(allseeds))==32
    training_seeds=json.loads((run/'training_Z_seeds.json').read_text());z_audit=[]
    for day in [27,55,118,153]:
        z,rec=outer_z(cache,day);assert rec==training_seeds[str(day)]
        z_audit.append(dict(day=day,mean=float(z.mean()),variance=float(z.var()),seeds=rec));assert not set(x['seed'] for x in rec)&set(allseeds)
    # New seeds are disjoint from every previously saved production/pilot/check seed.
    historical=set()
    for path in HERE.joinpath('runs').rglob('vix_d*_seeds.json'):
        if run in path.parents:continue
        for row in json.loads(path.read_text()):historical.add(int(row['seed']))
    assert not historical&set(allseeds)
    # Recorded seed is parameter-independent and all paired stock layers used the same seed.
    stock=json.loads((run/'full/stock_Z_audit.json').read_text());z_terminal=[]
    for row in stock:assert row['seed']==seed_for(2,row['day'],row['batch'])
    for name in selection['matched_specs']:
        spec=next(s for s in SPECS if s[0]==name)
        for day in [29,57,119,211,302,666,1030]:
            rows=[x for x in stock if x['candidate']==name and x['day']==day]
            # Same candidate under both cap choices repeats exact independent Z paths.
            assert len(rows)==64
            for batch in range(32):
                pair=[x for x in rows if x['batch']==batch];assert len(pair)==2
                np.testing.assert_array_equal(pair[0]['terminal_z'],pair[1]['terminal_z'])
            z=np.concatenate([np.array(next(x for x in rows if x['batch']==b)['terminal_z']) for b in range(32)])
            z_terminal.append(dict(candidate=name,day=day,mean=float(z.mean()),variance=float(z.var()),target_variance=float(-np.expm1(-2*spec[2]*day/365))))
    k=LayerKernel(run/'stochastic_vov.dylib',d['ends_years']);normal_audits=json.loads((run/'normalizer_audits.json').read_text());fourD=[];time_checks=[];normalizer_comparison=[];normalizer_refinements=[]
    tau=1/max(dated.layout.rates)
    for cap,continuation in [('cap0',0.),('cap05',.5)]:
        k.configure(continuation);nn=load_normalizers(run/f'normalizers_{cap}.npz',tau,normal_audits[cap])
        for day in [27,55,118,153]:
            t=day/365;m,cov=dated.gaussian_projection_law(t,initial);z,w=_normal_rule(48)
            original=float(k.expectation(arr(m),arr(z@symmetric_square_root(cov).T),arr(w),p)[0,0])
            for spec in SPECS:
                new=float(np.exp(nn[spec[0]][0](t)[0]));normalizer_comparison.append(dict(cap=cap,days=day,candidate=spec[0],original_M=original,new_M=new,ratio=new/original))
        for t in [27/365,153/365,1030/365]:
            m,cov=dated.gaussian_projection_law(t,initial);z,w=_normal_rule(96)
            c=response_table(k,(arr(m),arr(z@symmetric_square_root(cov).T),arr(w)),p,[],'normalizer-GH96')
            for spec in SPECS:
                s,rate=spec[1:3];fine=k.response(c,0.,np.sqrt(-np.expm1(-2*rate*t)),s,256)[0]
                coarse=float(np.exp(nn[spec[0]][0](t)[0]));normalizer_refinements.append(dict(cap=cap,candidate=spec[0],time=t,old_GH48_to96_and_new192_to256_relative_change=float(coarse/fine-1)))
        spec=next(s for s in SPECS if s[0]==selection['matched_specs'][0]);name,s,rate,_,_=spec
        for t in [27/365,153/365]:
            direct=oracle(dated,t,initial,p,continuation,s,0.,np.sqrt(-np.expm1(-2*rate*t)),48,192)
            fast=float(np.exp(nn[name][0](t)[0]));err=abs(fast/direct-1);assert err<1e-6,err
            fourD.append(dict(cap=cap,candidate=name,kind='unconditional',time=t,relative_error=err))
        # Independent conditional 4D oracle, followed by a time-order refinement
        # at ordinary and extreme new checking states. Full Z_T is conditioned on.
        day=153;draw=np.load(run/f'fresh_cache/production/vix_d{day}.npz')
        paths=np.load(run/'fresh_check'/f'paths_{cap}_{name}.npz');vv=paths[f'vix_{day}'];ids=np.unique([0,int(np.argmax(vv))])
        lag=.02;m,cov=dated.gaussian_projection_law(lag,draw['states'][ids]);z,w=_normal_rule(32);law=arr(m),arr(z@symmetric_square_root(cov).T),arr(w)
        c=response_table(k,law,p,[],'conditional-oracle');zT=np.sqrt(-np.expm1(-2*rate*day/365))*draw['z_normals'][ids];mu=np.exp(-rate*lag)*zT;sd=np.sqrt(-np.expm1(-2*rate*lag))
        fast=k.response(c,mu,sd,s,128)
        for j,index in enumerate(ids):
            direct=oracle(dated,lag,draw['states'][index],p,continuation,s,float(mu[j]),sd,32,128);err=abs(fast[j]/direct-1);assert err<3e-6,err
            fourD.append(dict(cap=cap,candidate=name,kind='conditional',day=day,lag=lag,index=int(index),relative_error=float(err)))
        times,tw=window_nodes(day/365,np.array(d['ends_years']),tau,16);integral=np.full(len(ids),p[-1]*30/365)
        for t,weight in zip(times,tw):
            lag=float(t-day/365);m,cov=dated.gaussian_projection_law(lag,draw['states'][ids]);z,w=_normal_rule(32)
            c=response_table(k,(arr(m),arr(z@symmetric_square_root(cov).T),arr(w)),p,[],'time16')
            value=k.response(c,np.exp(-rate*lag)*zT,np.sqrt(-np.expm1(-2*rate*lag)),s,128)
            j=min(np.searchsorted(d['ends_years'],t,side='right'),len(d['ends_years'])-1)
            integral+=weight*(p[7+j]-p[-1])*value/np.exp(nn[name][0](t)[0])
        fine=100*np.sqrt(integral/(30/365));change=float(np.max(abs(fine-vv[ids])))
        time_checks.append(dict(cap=cap,candidate=name,day=day,indices=ids.tolist(),time_order8_to16_maximum_VIX_change=change))
    pd.DataFrame(normalizer_comparison).to_csv(run/'normalizer_comparison.csv',index=False)
    dump(run/'audit.json',dict(status='PASS saved prices/parity/SE/loss, archived sources, seeds and independent normalizer oracle',score_records=records,
         frozen_sources=len(design['source_sha256']),baseline_full_maximum_price_difference=maximum,exact_joint_scrambles_replayed=32,historical_seed_count=len(historical),
         training_Z=z_audit,stock_Z=z_terminal,negative_estimate_audit=negative_estimates,independent_4D_checks=fourD,time_order_checks=time_checks,normalizer_order_checks=normalizer_refinements,selection_precedes_SPX_and_fresh=True,
         limitation='Finite oracles and sampling diagnostics; no population convergence, global-fit or P certificate.'))
    print('PASS',len(records),'price records;',len(fourD),'independent 4D oracles',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);audit(p.parse_args().run.resolve())
