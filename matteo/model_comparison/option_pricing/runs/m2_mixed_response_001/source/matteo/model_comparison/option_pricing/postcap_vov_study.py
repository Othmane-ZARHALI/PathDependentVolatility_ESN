"""Frozen-baseline comparison of independent normalized post-cap amplitudes."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
from .postcap_vov_core import PostKernel,PostEngine,SPECS,exposures
from .stochastic_vov_study import install
from .full_readout_fit import HERE,dump
from .m2_spx import ROOT,M2,frozen_stage3
from .cap_tail_study import engine,tails
from .cap_controlled_core import shape_metrics
from .maturity_curve import vix_only_loss
from .fit_parity_parameters_v3 import save_result
from .parity_adaptive import canonical_market
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture


def native_loss(score):
    c=score['original_otm_components'];return c['VIX']+c['forwards']


def shape_score(shape,market):
    if not shape['wing_complete']:return None
    errors=[]
    for r,m in zip(shape['rows'],market['rows']):
        for key in ['wing_slope','wing_curvature']:
            if r[key] is None:return None
            errors.append(((r[key]-m[key])/max(10,abs(m[key])))**2)
    return float(np.mean(errors))


def fresh_draws(out,dated,initial,cache):
    dest=out/'fresh_cache';prod=dest/'production';prod.mkdir(parents=True,exist_ok=False)
    proposals=json.loads((HERE/'runs/m2_vix_targeted_check_001/proposals.json').read_text());records=[]
    for day in [27,55,118,153]:
        dec=proposals[str(day)];shifts=np.c_[dec['shifts'],np.zeros(len(dec['shifts']))]
        mix=DefensiveGaussianMixture(shifts,dec['probabilities']);sampler=GaussianOuterSampler(dated,day/365,initial)
        zs=[];ww=[];cc=[];seeds=[]
        for net in range(8):
            seed=int(np.random.SeedSequence([2026102631,day,net]).generate_state(1,dtype=np.uint64)[0]);z,w,c=mix.draw(9,seed)
            zs.append(z);ww.append(w);cc.append(c);seeds.append(dict(net=net,seed=seed))
        z=np.concatenate(zs);states=sampler.mean+z[:,:-1]@sampler.root.T
        path=prod/f'vix_d{day}.npz';np.savez_compressed(path,states=states,normals=z[:,:-1],z_normals=z[:,-1],weights=np.concatenate(ww),components=np.concatenate(cc),nets=np.repeat(np.arange(8),512))
        records.append(dict(day=day,seeds=seeds,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    dump(dest/'design.json',dict(files=records,proposals=proposals,seed_base=2026102631,created_after_selection=True,weights='Exact original joint37D mixture weights; no self-normalization'))
    return dest


def run(out):
    out.mkdir(parents=True,exist_ok=False);start=time.perf_counter();prior=HERE/'runs/m2_maturity_curve_fit_001'
    d=json.loads((prior/'run.json').read_text());p=np.array(d['selected']);cache=Path(d['cache']);ctrl=Path(d['controls'])
    spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv');vc=json.loads((prior/'vix_controls.json').read_text())
    a,_,premium,_=frozen_stage3();assert json.loads(json.dumps(asdict(a)))==d['architecture'];assert json.loads(json.dumps(asdict(premium)))==d['premium']
    dated=DatedM2Pricer(a,premium);initial=np.array(d['initial_state']);np.testing.assert_array_equal(initial,dated.initial_state([1,0,0])[0])
    files=[*sorted(HERE.glob('*.py')),HERE/'postcap_vov_kernel.cpp',HERE/'stochastic_vov_kernel.cpp',HERE/'cap_tail_kernel.cpp',HERE/'POSTCAP_VOV_PLAN.md',prior/'run.json',prior/'selected_model.json',prior/'vix_controls.json',cache/'design.json',cache/'spx_contracts.csv',cache/'vix_contracts.csv',cache/'forward_targets.csv',ctrl/'coefficients.json',HERE/'runs/m2_vix_targeted_check_001/proposals.json',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(files);save_sources(out,captured,ROOT)
    design=dict(status='running',baseline=p.tolist(),parent=str(prior),cache=str(cache),controls=str(ctrl),source_sha256=source_hashes(captured,ROOT),specs=SPECS,paths_per_strike=4096,market_out_of_sample=False,baseline_replaced=False)
    dump(out/'design.json',design);lib=out/'postcap.dylib'
    subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'postcap_vov_kernel.cpp'),'-o',str(lib)],check=True)
    k=PostKernel(lib,d['ends_years']);base=engine(out,cache,k,dated,initial,spx,vix,f,vc,ctrl,d)
    def new_engine(full=False):
        ee=PostEngine(cache,k,dated,initial,spx,vix,f,controls=vc,frozen_controls=ctrl,ends=d['ends_years'],reference_levels=d['reference_levels'])
        if not full:
            ee.spx=spx.iloc[:0];ee.original_quotes=vix.reset_index(drop=True);ee.quotes=canonical_market(ee.original_quotes);ee.loss=vix_only_loss(ee.quotes,f);ee.otm_loss=vix_only_loss(ee.original_quotes,f)
        return ee
    def save_nodes(folder,nodes):
        np.savez_compressed(folder/'base_conditional_nodes.npz',**{f'd{day}_{key}':v for day,r in nodes.items() for key,v in r.items()})
    try:
        norm,na=base.normalizer(p);dump(out/'normalizer.json',na)
        ex,draws,seeds,nodes=exposures(base,k,p,norm,cache);save_nodes(out,nodes);dump(out/'training_Z_seeds.json',seeds)
        np.savez_compressed(out/'exposures.npz',**{f'{name}_d{day}':value for name,dd in ex.items() for day,value in dd.items()})
        profiles=[];tt=[];mm=[];ee=new_engine();market=shape_metrics(ee.quotes,ee.quotes.price_mid_U.to_numpy());dump(out/'market_shape.json',market)
        try:
            for spec in SPECS:
                name=spec[0];install(ee,k,spec,(norm,na),ex[name],draws);r=ee.predict(p,True,False);score=save_result(ee,r,p,name,out);sh=shape_metrics(ee.quotes,r['prices'])
                if name=='baseline':
                    saved=pd.read_csv(prior/'prices_canonical_selected.csv');saved=saved[saved.instrument.eq('VIX')];err=float(np.max(abs(saved.model_price_U.to_numpy()-r['prices'])));assert err<2e-6;dump(out/'baseline_replay.json',dict(maximum_VIX_price_difference=err))
                profiles.append(dict(name=name,spec=spec,score=score,shape=sh,shape_score=shape_score(sh,market)))
                t,m=tails(out,name,vix,p,np.array(d['ends_days']));tt+=t;mm+=m
                dump(out/'profiles.json',profiles);pd.DataFrame(tt).to_csv(out/'tails.csv',index=False);pd.DataFrame(mm).to_csv(out/'moments.csv',index=False)
                print('PROFILE',name,'native',native_loss(score),'shape',shape_score(sh,market),'wing',sh['wing_MAE'],'ATM',sh['ATM_MAE'],flush=True)
        finally:ee.close()
        eligible=[r for r in profiles if r['name']!='baseline' and r['shape_score'] is not None]
        if not eligible:raise RuntimeError('No candidate with complete wing; retain all failures')
        price_best=min(eligible,key=lambda r:(native_loss(r['score']),r['shape_score']))
        shape_best=min(eligible,key=lambda r:(r['shape_score'],native_loss(r['score'])))
        names=list(dict.fromkeys(['baseline',price_best['name'],shape_best['name']]))
        selection=dict(price_winner=price_best['name'],shape_winner=shape_best['name'],full_candidates=names,frozen_before_SPX_and_checking=True,automatic_adoption=False)
        dump(out/'selection.json',selection);print('SELECTION',selection,flush=True)
        full=out/'full';full.mkdir();fullscores=[]
        for name in names:
            spec=next(s for s in SPECS if s[0]==name);ee=new_engine(True);install(ee,k,spec,(norm,na),ex[name],draws)
            try:
                r=ee.predict(p,True,False);score=save_result(ee,r,p,name,full);fullscores.append(score);dump(full/'scores.json',fullscores);dump(full/'stock_Z_audit.json',k.stock_audit)
                if name=='baseline':
                    saved=pd.read_csv(prior/'prices_canonical_selected.csv');err=float(np.max(abs(saved.model_price_U.to_numpy()-r['prices'])));assert err<2e-6;dump(full/'baseline_replay.json',dict(maximum_price_difference=err))
                print('FULL DONE',name,score['original_otm_components'],flush=True)
            finally:ee.close()
        fresh=fresh_draws(out,dated,initial,cache);check=out/'fresh_check';check.mkdir();specs=[s for s in SPECS if s[0] in names]
        fe,fa,_,fn=exposures(base,k,p,norm,fresh,specs=specs);save_nodes(check,fn)
        np.savez_compressed(check/'exposures.npz',**{f'{name}_d{day}':value for name,dd in fe.items() for day,value in dd.items()})
        fs=[];ft=[];fm=[];ee=new_engine()
        try:
            for spec in specs:
                name=spec[0];install(ee,k,spec,(norm,na),fe[name],fa);r=ee.predict(p,True,False);sc=save_result(ee,r,p,name,check);sh=shape_metrics(ee.quotes,r['prices']);fs.append(dict(name=name,score=sc,shape=sh,shape_score=shape_score(sh,market)))
                t,m=tails(check,name,vix,p,np.array(d['ends_days']));ft+=t;fm+=m;print('FRESH',name,'native',native_loss(sc),'shape',shape_score(sh,market),flush=True)
        finally:ee.close()
        dump(check/'profiles.json',fs);pd.DataFrame(ft).to_csv(check/'tails.csv',index=False);pd.DataFrame(fm).to_csv(check/'moments.csv',index=False)
        # Refine representative conditional values and complete selected VIX windows.
        # Use a small explicit cache for deterministic refinements, preserving extreme states.
        probe=out/'quadrature_cache'/'production';probe.mkdir(parents=True)
        probe_ids={}
        for day,aa in fa.items():
            ids=np.unique(np.r_[0,1,np.argmax(abs(aa['states']),axis=0),np.argmax(abs(aa['z_normals']))]);probe_ids[day]=ids
            np.savez_compressed(probe/f'vix_d{day}.npz',**{key:value[ids] for key,value in aa.items() if isinstance(value,np.ndarray) and len(value)==4096})
        fine,_,_,_=exposures(base,k,p,norm,probe.parent,specs=specs,time_order=16,order=48)
        quadrature=[]
        for spec in specs:
            name=spec[0]
            for day in fa:
                ids=probe_ids[day];v0=100*np.sqrt(p[-1]+fe[name][day][ids]@(p[7:-1]-p[-1])/(30/365));v1=100*np.sqrt(p[-1]+fine[name][day]@(p[7:-1]-p[-1])/(30/365))
                quadrature.append(dict(name=name,day=day,ids=ids.tolist(),maximum_VIX_point_difference=float(np.max(abs(v0-v1))),coarse='GH32,time8',fine='GH48,time16'))
        dump(out/'quadrature_checks.json',quadrature)
        check_sources(captured);dump(out/'run.json',design|dict(status='completed',selection=selection,wall_seconds=time.perf_counter()-start,full_scores=fullscores,limitations=['Fixed-curve finite grid; no global joint calibration.','Fresh numerical VIX streams use same market data; no fresh SPX sample.','Sampled quadrature checks and Monte Carlo errors are not population convergence certificates.']))
        print('COMPLETED',time.perf_counter()-start,flush=True)
    finally:base.close()

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    try:run(args.output.resolve())
    except Exception as error:
        if args.output.exists():dump(args.output/'failure.json',dict(type=type(error).__name__,message=str(error)))
        raise
