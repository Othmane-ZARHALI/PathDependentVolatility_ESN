"""Bounded two-cap, two-strength, three-persistence dedicated-layer study."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

from .stochastic_vov_core import (LayerKernel,LayerEngine,SPECS,SEED,seed_for,
    make_normalizers,make_exposures,response_table,normal_rule,arr)
from .full_readout_fit import HERE,dump
from .m2_spx import ROOT,M2,frozen_stage3
from .cap_tail_study import engine,tails
from .cap_controlled_core import shape_metrics
from .maturity_curve import vix_only_loss
from .fit_parity_parameters_v3 import save_result
from .parity_adaptive import canonical_market
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture

CAPS=[('cap0',0.),('cap05',.5)]


def load_normalizers(file,tau,audit):
    with np.load(file) as a:spline=CubicSpline(a['knots'],a['log_values'],axis=0,extrapolate=False)
    def single(j):
        def norm(t):
            val=spline(np.log1p(np.asarray(t)/tau))[...,j];out=np.zeros(np.shape(val)+(8,));out[...,0]=val;return out
        return norm
    return {spec[0]:(single(j),audit) for j,spec in enumerate(SPECS)}


def install(e,k,spec,normalizer,exposures,draws):
    k.spec=spec;e.layer_normalizer,e.layer_normalizer_audit=normalizer
    e.layer_exposures=exposures;e.layer_draws=draws


def fresh_cache(out,dated,initial,cache):
    dest=out/'fresh_cache';production=dest/'production';production.mkdir(parents=True,exist_ok=False)
    proposals=json.loads((HERE/'runs/m2_vix_targeted_check_001/proposals.json').read_text());records=[]
    for day in [27,55,118,153]:
        dec=proposals[str(day)];shifts=np.c_[np.array(dec['shifts']),np.zeros(len(dec['shifts']))]
        mixture=DefensiveGaussianMixture(shifts,dec['probabilities']);sampler=GaussianOuterSampler(dated,day/365,initial)
        zs=[];ww=[];cc=[];seeds=[]
        for net in range(8):
            seed=seed_for(3,day,net);z,w,c=mixture.draw(9,seed);zs.append(z);ww.append(w);cc.append(c);seeds.append(dict(net=net,seed=seed,power=9))
        z=np.concatenate(zs);states=sampler.mean+z[:,:-1]@sampler.root.T
        source=np.load(cache/'production'/f'vix_d{day}.npz');path=production/f'vix_d{day}.npz'
        np.savez_compressed(path,states=states,normals=z[:,:-1],z_normals=z[:,-1],weights=np.concatenate(ww),components=np.concatenate(cc),nets=np.repeat(np.arange(8),512),times=source['times'],time_weights=source['time_weights'])
        dump(production/f'vix_d{day}_seeds.json',seeds);records.append(dict(days=day,seeds=seeds,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    dump(dest/'design.json',dict(seed_base=SEED,seed_formula='SeedSequence([base,3,day,net]) uint64',proposal=proposals,
         augmented_shift='Append one zero coordinate for independent Z; joint 37D mixture and exact original weights',files=records,created_after_selection=True))
    return dest


def quadrature_checks(e,k,p,normalizers,cache):
    """Representative 4D integral and original-state/order refinements."""
    checks=[]
    for t in [.001,27/365,153/365,1030/365]:
        c=response_table(k,e.law(t,e.initial,48),p,[],f'normal-refine{t}')
        for name,s,rate,_,_ in SPECS:
            sd=np.sqrt(-np.expm1(-2*rate*t));coarse=k.response(c,0.,sd,s,192)[0];fine=k.response(c,0.,sd,s,256)[0]
            checks.append(dict(kind='new_Z_normalizer_192_to_256',time=t,candidate=name,relative_error=float(abs(coarse/fine-1))))
    for day in [27,153]:
        a=np.load(cache/'production'/f'vix_d{day}.npz');ids=np.unique(np.r_[0,1,np.argmax(abs(a['states']),axis=0)[:4]])
        if 'z_normals' in a:eps=a['z_normals']
        else:
            from .stochastic_vov_core import outer_z
            eps,_=outer_z(cache,day)
        for lag in [.0002,.02,30/365]:
            c=response_table(k,e.law(lag,a['states'][ids],32),p,[],f'conditional-refine{day}_{lag}')
            fine_c=response_table(k,e.law(lag,a['states'][ids],48),p,[],f'conditional-old-refine{day}_{lag}')
            for name,s,rate,_,_ in SPECS:
                zT=np.sqrt(-np.expm1(-2*rate*day/365))*eps[ids];mu=np.exp(-rate*lag)*zT;sd=np.sqrt(-np.expm1(-2*rate*lag))
                coarse=k.response(c,mu,sd,s,128);fine=k.response(c,mu,sd,s,256);oldfine=k.response(fine_c,mu,sd,s,256)
                checks.append(dict(kind='new_Z_conditional_128_to_256',day=day,lag=lag,candidate=name,relative_error=float(np.max(abs(coarse/fine-1)))))
                checks.append(dict(kind='old_Gaussian_conditional_32_to_48',day=day,lag=lag,candidate=name,relative_error=float(np.max(abs(fine/oldfine-1)))))
    return checks


def run(out):
    out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();prior=HERE/'runs/m2_maturity_curve_fit_001'
    d=json.loads((prior/'run.json').read_text());p=np.array(d['selected']);cache=Path(d['cache']);controls=Path(d['controls'])
    spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv');vc=json.loads((prior/'vix_controls.json').read_text())
    architecture,_,premium,_=frozen_stage3();assert json.loads(json.dumps(asdict(architecture)))==d['architecture']
    assert json.loads(json.dumps(asdict(premium)))==d['premium'];dated=DatedM2Pricer(architecture,premium);initial=np.array(d['initial_state'])
    np.testing.assert_array_equal(initial,dated.initial_state([1,0,0])[0])
    files=[*sorted(HERE.glob('*.py')),HERE/'stochastic_vov_kernel.cpp',HERE/'cap_tail_kernel.cpp',HERE/'STOCHASTIC_VOV_PLAN.md',prior/'run.json',prior/'selected_model.json',prior/'vix_controls.json',cache/'design.json',cache/'spx_contracts.csv',cache/'vix_contracts.csv',cache/'forward_targets.csv',controls/'coefficients.json',HERE/'runs/m2_vix_targeted_check_001/proposals.json',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(files);save_sources(out,captured,ROOT)
    design=dict(status='running fixed-curve dedicated sensitivity grid',parent=str(prior),baseline=p.tolist(),cache=str(cache),controls=str(controls),caps=CAPS,
         specs=[dict(name=n,sensitivity=s,rate=rate,half_life=half,factor=factor) for n,s,rate,half,factor in SPECS],
         fixed='All original parameters including entire 16-piece xi and floor; independent added Q OU with zero mean and Z0=0',
         source_sha256=source_hashes(captured,ROOT),seed_base=SEED,paths_per_strike=4096,market_out_of_sample=False)
    dump(out/'design.json',design);lib=out/'stochastic_vov.dylib'
    subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'stochastic_vov_kernel.cpp'),'-o',str(lib)],check=True)
    k=LayerKernel(lib,d['ends_years']);e=engine(out,cache,k,dated,initial,spx,vix,f,vc,controls,d)
    def new_engine(vix_only=True):
        ee=LayerEngine(cache,k,dated,initial,spx,vix,f,controls=vc,frozen_controls=controls,ends=d['ends_years'],reference_levels=d['reference_levels'])
        if vix_only:
            ee.spx=spx.iloc[:0];ee.original_quotes=vix.reset_index(drop=True);ee.quotes=canonical_market(ee.original_quotes);ee.loss=vix_only_loss(ee.quotes,f);ee.otm_loss=vix_only_loss(ee.original_quotes,f)
        return ee
    profiles=[];allchecks={};norm_audits={};norms={};exposures={};draws={};tailrows=[];moments=[]
    try:
        for cap,continuation in CAPS:
            k.configure(continuation);e.last_readout=None;e.exposure_readout_key=None
            print('REFERENCE',cap,flush=True);r=e.predict(p,True,False);score=save_result(e,r,p,cap+'_ref',out)
            if continuation==0:
                saved=pd.read_csv(prior/'prices_canonical_selected.csv');saved=saved[saved.instrument.eq('VIX')]
                error=float(np.max(abs(saved.model_price_U.to_numpy()-r['prices'])));assert error<2e-6;dump(out/'baseline_replay.json',dict(maximum_price_difference=error))
            profiles.append(dict(name=cap+'_ref',cap=cap,continuation=continuation,layer=False,shape=shape_metrics(e.quotes,r['prices']),score=score))
            tt,mm=tails(out,cap+'_ref',vix,p,np.array(d['ends_days']));tailrows+=tt;moments+=mm
            checks=[];print('NORMALIZERS',cap,flush=True);nn=make_normalizers(e,k,p,checks,out/f'normalizers_{cap}.npz');norms[cap]=nn;norm_audits[cap]=nn[SPECS[0][0]][1]
            dump(out/'normalizer_audits.json',norm_audits)
            ex,aa,seeds=make_exposures(e,k,p,nn,cache,checks,out/f'exposures_{cap}.npz');exposures[cap]=ex;draws=aa
            dump(out/'training_Z_seeds.json',seeds);allchecks[cap]=checks;dump(out/'interpolation_checks.json',allchecks)
            qe=quadrature_checks(e,k,p,nn,cache);dump(out/f'quadrature_{cap}.json',qe)
            ee=new_engine()
            try:
                for spec in SPECS:
                    name=cap+'_'+spec[0];install(ee,k,spec,nn[spec[0]],ex[spec[0]],aa);r=ee.predict(p,True,False);score=save_result(ee,r,p,name,out)
                    shape=shape_metrics(ee.quotes,r['prices']);profiles.append(dict(name=name,cap=cap,continuation=continuation,layer=True,spec=spec,shape=shape,score=score))
                    tt,mm=tails(out,name,vix,p,np.array(d['ends_days']));tailrows+=tt;moments+=mm
                    print('PROFILE',name,'native VIX+f',score['original_otm_components']['VIX']+score['original_otm_components']['forwards'],'wing',shape['wing_MAE'],'ATM',shape['ATM_MAE'],flush=True)
                    dump(out/'profiles.json',profiles);pd.DataFrame(tailrows).to_csv(out/'training_tails.csv',index=False);pd.DataFrame(moments).to_csv(out/'training_moments.csv',index=False)
            finally:ee.close()
        selected=[]
        for cap,_ in CAPS:
            eligible=[r for r in profiles if r['cap']==cap and r['layer'] and r['shape']['wing_complete']]
            if not eligible:raise RuntimeError(f'No complete wing candidate for {cap}; retain all failures')
            selected.append(min(eligible,key=lambda r:(r['shape']['wing_MAE'],r['score']['original_otm_components']['VIX']+r['score']['original_otm_components']['forwards'])))
        pairs=sorted(set(r['spec'][0] for r in selected));names=[cap+'_ref' for cap,_ in CAPS]+[cap+'_'+name for cap,_ in CAPS for name in pairs]
        selection=dict(training_winners=[r['name'] for r in selected],matched_specs=pairs,full_candidates=names,frozen_before_SPX_and_fresh_sampling=True,rule='Minimum equal-maturity wing IV MAE per cap, complete wing required; no aggregate/ATM/SPX acceptance implied')
        dump(out/'selection.json',selection);print('SELECTION',selection,flush=True)
        full=out/'full';full.mkdir();fullscores=[]
        for cap,continuation in CAPS:
            k.configure(continuation)
            # Each no-layer cap reference receives its own original normalizer.
            k.spec=('reference',0.,np.log(2)*365/30,30,1.)
            ref=engine(out,cache,k,dated,initial,spx,vix,f,vc,controls,d,vix_only=False)
            try:
                print('FULL',cap+'_ref',flush=True);r=ref.predict(p,True,False);fullscores.append(save_result(ref,r,p,cap+'_ref',full))
            finally:ref.close()
            for name in pairs:
                spec=next(s for s in SPECS if s[0]==name);ee=new_engine(False);install(ee,k,spec,norms[cap][name],exposures[cap][name],draws)
                try:
                    print('FULL',cap+'_'+name,flush=True);r=ee.predict(p,True,False);fullscores.append(save_result(ee,r,p,cap+'_'+name,full));dump(full/'scores.json',fullscores);dump(full/'stock_Z_audit.json',k.stock_audit)
                finally:ee.close()
        fresh=fresh_cache(out,dated,initial,cache);freshout=out/'fresh_check';freshout.mkdir();freshscores=[];ft=[];fm=[]
        for cap,continuation in CAPS:
            k.configure(continuation);ref=engine(out,cache,k,dated,initial,spx,vix,f,vc,None,d);ref.cache=fresh;ref.folder=fresh/'production'
            try:r=ref.predict(p,True,False);freshscores.append(save_result(ref,r,p,cap+'_ref',freshout))
            finally:ref.close()
            checks=[];ex,aa,seeds=make_exposures(e,k,p,norms[cap],fresh,checks,freshout/f'exposures_{cap}.npz')
            dump(freshout/f'interpolation_{cap}.json',checks)
            for name in pairs:
                spec=next(s for s in SPECS if s[0]==name);ee=new_engine();install(ee,k,spec,norms[cap][name],ex[name],aa)
                try:r=ee.predict(p,True,False);freshscores.append(save_result(ee,r,p,cap+'_'+name,freshout))
                finally:ee.close()
            for name in [cap+'_ref']+[cap+'_'+n for n in pairs]:
                tt,mm=tails(freshout,name,vix,p,np.array(d['ends_days']));ft+=tt;fm+=mm
            dump(freshout/'scores.json',freshscores);pd.DataFrame(ft).to_csv(freshout/'tails.csv',index=False);pd.DataFrame(fm).to_csv(freshout/'moments.csv',index=False)
        check_sources(captured)
        dump(out/'run.json',design|dict(status='completed dedicated stochastic sensitivity grid, matched full SPX comparisons and new VIX checking',selection=selection,full_scores=fullscores,fresh_scores=freshscores,wall_seconds=time.perf_counter()-begin,baseline_replaced=False,limitations=['Bounded fixed-curve grid; no global calibration.','Full SPX only for matched training selections and both cap references.','No new market OOS, P calibration, Lean proof or population convergence certificate.']))
        print('COMPLETED',time.perf_counter()-begin,flush=True)
    finally:e.close()


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    try:run(args.output.resolve())
    except Exception as error:
        if args.output.exists():dump(args.output/'failure.json',dict(type=type(error).__name__,message=str(error)))
        raise
