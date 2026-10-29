"""Bounded clean cap comparison, actual SPX tradeoffs, genuinely new VIX nets."""
import argparse,ctypes as ct,json,time,subprocess,hashlib
from dataclasses import asdict
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import minimize
from .full_readout_fit import HERE,arr,dump
from .m2_spx import ROOT,M2,frozen_stage3
from .maturity_curve import CurveKernel
from .cap_tail_study import engine,tails
from .cap_controlled_core import window_basis,adjusted,curve_jac,prediction,shape_metrics
from .fit_parity_parameters_v3 import save_result
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture

SHAPES=[(dl,a) for dl in [0.,.375,.75] for a in [0.,.25,.5]]
SEED=2026102511

def configure(k,a):k.lib.configure(0 if a==0 else 2,a)

class ProgressKernel(CurveKernel):
    def stock(self,folder,norm,p,derivatives=False):
        if derivatives:raise ValueError('Experimental readout is value-only')
        if folder.name=='batch000':print('SPX START',folder.parent.name,flush=True);self.begin=time.perf_counter()
        v=super().stock(folder,norm,p,False)
        if folder.name=='batch031':print('SPX COMPLETE',folder.parent.name,'seconds',time.perf_counter()-self.begin,flush=True)
        return v

def fit_shape(e,p,B,limit):
    D=curve_jac(p,B);memo={};history=[]
    def calc(s):
        key=tuple(s)
        if key not in memo:
            memo.clear();x=adjusted(p,s,B);r=prediction(e,x,D);memo[key]=r
            history.append(dict(adjustments=np.asarray(s).tolist(),native_loss=r['native_loss'],common_loss=r['common_loss'],forwards=r['forward'].tolist()))
        return memo[key]
    market=e.forwards.forward_vix_points.to_numpy()
    def slack(s):
        r=calc(s);err=r['forward']-market
        return np.r_[.5-err,.5+err,(limit-r['common_loss'])/max(1.,limit)]
    def cj(s):
        r=calc(s);return np.r_[-r['forward_jac'],r['forward_jac'],-r['common_gradient'][None,:]/max(1.,limit)]
    result=minimize(lambda s:calc(s)['native_loss'],np.zeros(4),jac=lambda s:calc(s)['native_gradient'],bounds=[(-.15,.15)]*4,
                    constraints=[dict(type='ineq',fun=slack,jac=cj)],method='SLSQP',options=dict(maxiter=60,ftol=1e-9))
    s=result.x;r=calc(s);sl=slack(s);x=adjusted(p,s,B)
    return x,dict(adjustments=s.tolist(),success=bool(result.success),termination=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev),history=history,
                  native_loss=float(r['native_loss']),common_loss=float(r['common_loss']),slack=sl.tolist(),feasible=bool(np.min(sl)>=-1e-6),forwards=r['forward'].tolist())

def new_cache(out,dated,initial,design):
    dest=out/'fresh_cache';dest.mkdir(exist_ok=False);production=dest/'production';production.mkdir();old=set()
    for folder in HERE.joinpath('runs').glob('*/production'):
        for path in folder.glob('vix_d*_seeds.json'):
            for row in json.loads(path.read_text()):old.add(int(row['seed']))
    marginal=HERE/'runs/m2_cap_tail_001/marginal_response.csv'
    if marginal.exists():old.update(pd.read_csv(marginal).seed.astype(int).tolist())
    proposals=json.loads((HERE/'runs/m2_vix_targeted_check_001/proposals.json').read_text());allseeds=[];meta=[]
    for day in [27,55,118,153]:
        dec=proposals[str(day)];mixture=DefensiveGaussianMixture(dec['shifts'],dec['probabilities']);sampler=GaussianOuterSampler(dated,day/365,initial);states=[];normals=[];weights=[];components=[];records=[]
        for net in range(8):
            seed=int(np.random.SeedSequence([SEED,2,day,net]).generate_state(1)[0]);assert seed not in old and seed not in allseeds;allseeds.append(seed)
            z,w,c=mixture.draw(9,seed);states.append(sampler.mean+z@sampler.root.T);normals.append(z);weights.append(w);components.append(c);records.append(dict(net=net,seed=seed,power=9))
        src=np.load(Path(design['cache'])/'production'/f'vix_d{day}.npz')
        path=production/f'vix_d{day}.npz';np.savez_compressed(path,states=np.concatenate(states),normals=np.concatenate(normals),weights=np.concatenate(weights),components=np.concatenate(components),nets=np.repeat(np.arange(8),512),times=src['times'],time_weights=src['time_weights'])
        dump(production/f'vix_d{day}_seeds.json',records);meta.append(dict(days=day,seeds=records,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    dump(dest/'design.json',dict(seed_base=SEED,seed_formula='SeedSequence([base,2,days,net]).generate_state(1)[0]',previous_seed_count=len(old),new_seeds_disjoint=True,files=meta,proposal=proposals,created_after_selection=True,paths_per_strike=4096))
    return dest

def run(out):
    out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();prior=HERE/'runs/m2_maturity_curve_fit_001';design=json.loads((prior/'run.json').read_text());base=np.array(design['selected']);B=window_basis(design['ends_days']);assert np.array_equal(np.flatnonzero(B.sum(1)),[1,2,3,4,7,8,10])
    cache=Path(design['cache']);controls=Path(design['controls']);spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv');vc=json.loads((prior/'vix_controls.json').read_text())
    files=[Path(__file__),HERE/'cap_controlled_core.py',HERE/'test_cap_controlled.py',HERE/'CAP_CONTROLLED_PLAN.md',HERE/'cap_tail_kernel.cpp',prior/'run.json',prior/'selected_model.json',prior/'vix_controls.json',cache/'design.json',cache/'spx_contracts.csv',cache/'vix_contracts.csv',cache/'forward_targets.csv',controls/'coefficients.json',HERE/'runs/m2_vix_targeted_check_001/proposals.json',*sorted(HERE.glob('*.py')),*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(files);save_sources(out,captured,ROOT);lib=out/'cap.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'cap_tail_kernel.cpp'),'-o',str(lib)],check=True)
    k=ProgressKernel(lib,lib,np.array(design['ends_years']));k.lib.configure.argtypes=[ct.c_int,ct.c_double];k.lib.direct.argtypes=k.lib.expectation.argtypes
    a,_,premium,_=frozen_stage3();assert json.loads(json.dumps(asdict(a)))==design['architecture'];assert json.loads(json.dumps(asdict(premium)))==design['premium'];dated=DatedM2Pricer(a,premium);initial=np.array(design['initial_state']);np.testing.assert_array_equal(initial,dated.initial_state([1,0,0])[0])
    current=dict(status='training shape comparison',parent=str(prior),baseline=base.tolist(),basis=B.tolist(),shapes=SHAPES,adjustment_bounds=[-.15,.15],paths_per_strike=4096,source_sha256=source_hashes(captured,ROOT),cache=str(cache),controls=str(controls),seed_base=SEED,market_out_of_sample=False)
    dump(out/'design.json',current);records=[];tailrows=[];mom=[];gradchecks=[]
    e=engine(out,cache,k,dated,initial,spx,vix,f,vc,controls,design)
    try:
        configure(k,0);r=e.predict(base,True,False);save_result(e,r,base,'baseline',out);baseline_shape=shape_metrics(e.quotes,r['prices']);baseline_native=r['original_otm_components']['VIX']+r['original_otm_components']['forwards'];limit=r['components']['VIX']+r['components']['forwards']
        old=pd.read_csv(prior/'prices_canonical_selected.csv');old=old[old.instrument.eq('VIX')];err=float(np.max(abs(old.model_price_U.to_numpy()-r['prices'])));assert err<2e-6;dump(out/'baseline_replay.json',dict(maximum_price_difference=err,shape=baseline_shape,native_vix_plus_forward=baseline_native,common_vix_plus_forward=limit))
        for i,(dl,strength) in enumerate(SHAPES):
            name=f'shape{i:02d}';tick=time.perf_counter();configure(k,strength);p=base.copy();p[5]+=dl;e.last_readout=None;e.exposure_readout_key=None
            print('FIT',name,'L shift',dl,'continuation',strength,flush=True);selected,fit=fit_shape(e,p,B,limit)
            # Finite-difference the four actual controlled price/forward objective columns.
            zero=np.zeros(4);D=curve_jac(p,B);rr=prediction(e,p,D);fd=[]
            for j in range(4):
                h=2e-6;ss=zero.copy();ss[j]=h;pp=prediction(e,adjusted(p,ss,B),D);ss[j]=-h;mm=prediction(e,adjusted(p,ss,B),D);fd.append((pp['native_loss']-mm['native_loss'])/(2*h))
            rel=float(np.linalg.norm(np.array(fd)-rr['native_gradient'])/max(1.,np.linalg.norm(rr['native_gradient'])));assert rel<.003;gradchecks.append(dict(candidate=name,relative_error=rel))
            result=e.predict(selected,True,False);score=save_result(e,result,selected,name,out);shape=shape_metrics(e.quotes,result['prices']);t,m=tails(out,name,vix,selected,np.array(design['ends_days']));tailrows+=t;mom+=m
            # Preserve exact exposures to independently replay the four-parameter fit.
            np.savez_compressed(out/f'exposures_{name}.npz',**{f'd{day}':e.readout_exposures[(day,'production')][0] for day in [27,55,118,153]})
            preference=fit['feasible'] and fit['native_loss']<baseline_native and shape['ATM_MAE']<=baseline_shape['ATM_MAE']+5 and shape['wing_complete']
            row=dict(name=name,cap_shift=dl,continuation=strength,parameters=selected.tolist(),fit=fit,shape=shape,preferred=bool(preference),seconds=time.perf_counter()-tick,score=score);records.append(row)
            dump(out/'profiles.json',records);dump(out/'gradient_checks.json',gradchecks);pd.DataFrame(tailrows).to_csv(out/'training_tails.csv',index=False);pd.DataFrame(mom).to_csv(out/'training_moments.csv',index=False)
            print('PROFILE',name,fit['native_loss'],'wingMAE',shape['wing_MAE'],'ATM',shape['ATM_MAE'],'feasible',fit['feasible'],'preferred',preference,flush=True)
    finally:e.close()
    selected=[]
    for family in ['capped','continued']:
        pool=[r for r in records if (r['continuation']==0)==(family=='capped') and r['fit']['feasible']];preferred=[r for r in pool if r['preferred']]
        if pool:
            best=min(preferred or pool,key=lambda r:r['shape']['wing_MAE']);selected.append(dict(name=best['name'],family=family,parameters=best['parameters'],continuation=best['continuation'],preferred=best['preferred'],selection_rule='minimum equal-maturity wing IV MAE among preferred feasible candidates, or feasible fallback'))
    dump(out/'selection.json',dict(candidates=selected,frozen_before_fresh_sampling=True,selection_used_training_only=True));print('SELECTION FROZEN',[(r['family'],r['name']) for r in selected],flush=True)
    full=out/'full';full.mkdir();fullscores=[]
    for row in selected:
        configure(k,row['continuation']);x=np.array(row['parameters']);e=engine(out,cache,k,dated,initial,spx,vix,f,vc,controls,design,False);tick=time.perf_counter()
        try:
            print('FULL',row['name'],flush=True);r=e.predict(x,True,False);s=save_result(e,r,x,row['name'],full);s['wall_seconds']=time.perf_counter()-tick;fullscores.append(s);dump(full/'scores.json',fullscores);print('FULL SCORE',row['name'],r['original_otm_components'],flush=True)
        finally:e.close()
    # Fresh states are generated only after training-based representative selection.
    fresh=new_cache(out,dated,initial,design);check=out/'fresh_check';check.mkdir();checks=[];ft=[];fm=[]
    for row in [dict(name='baseline',parameters=base.tolist(),continuation=0),*selected]:
        configure(k,row['continuation']);x=np.array(row['parameters']);e=engine(out,cache,k,dated,initial,spx,vix,f,vc,None,design);e.cache=fresh;e.folder=fresh/'production';print('FRESH CHECK',row['name'],flush=True)
        try:
            r=e.predict(x,True,False);checks.append(save_result(e,r,x,row['name'],check));t,m=tails(check,row['name'],vix,x,np.array(design['ends_days']));ft+=t;fm+=m
        finally:e.close()
    pd.DataFrame(ft).to_csv(check/'tails.csv',index=False);pd.DataFrame(fm).to_csv(check/'moments.csv',index=False);dump(check/'scores.json',checks)
    check_sources(captured);dump(out/'run.json',current|dict(status='completed controlled cap comparison, actual SPX tradeoffs and fresh VIX numerical checks',selection=selected,full_scores=fullscores,fresh_scores=checks,wall_seconds=time.perf_counter()-begin,baseline_adopted_or_replaced=False,new_layer_implemented=False))
    print('COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);run(ap.parse_args().output.resolve())
