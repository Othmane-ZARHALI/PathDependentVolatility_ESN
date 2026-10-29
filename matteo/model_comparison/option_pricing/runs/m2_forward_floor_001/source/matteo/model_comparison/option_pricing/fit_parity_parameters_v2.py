"""Bounded joint11-parameter fit, common calls and4096 paths everywhere."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from .m2_spx import ROOT,M2,frozen_stage3
from .full_readout_fit import HERE,NAMES,SCALE,Kernel,dump
from .controlled_readout_fast_fit import LOWER,UPPER
from .parity_parameter_core import accept_proposal
from .parity_parameter_constraints import ConstrainedFitEngine as FitEngine,constrained_local,anchor_criteria
from .parity_stock import StockValues
from .audit_continuous_strike_controls import shape
from .run_parity_adaptive import ivs
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.pricing import DatedM2Pricer


def save_result(engine,r,p,name,out):
    q,f=engine.loss.tables(r['prices'],r['errors'],r['forwards'],r['forward_errors'])
    q['actual_paths']=r['actual_paths'];q['effective_paths']=r['effective_paths'];q['native_otm_model_price_U']=r['native_prices'];q['native_otm_mc_se_U']=r['native_errors']
    q['model_iv_percent']=ivs(q,r['prices']);q['iv_unavailable']=~np.isfinite(q.model_iv_percent)
    q.to_csv(out/f'prices_canonical_{name}.csv',index=False);f.to_csv(out/f'forwards_{name}.csv',index=False)
    native,_=engine.otm_loss.tables(r['native_prices'],r['native_errors'],r['forwards'],r['forward_errors']);native.to_csv(out/f'prices_otm_{name}.csv',index=False)
    np.savez_compressed(out/f'paths_{name}.npz',**r['raw'])
    record=dict(candidate=name,parameters=p.tolist(),canonical_loss=r['total'],canonical_components=r['components'],original_otm_loss=r['original_otm_total'],original_otm_components=r['original_otm_components'],
                inside_canonical_bidask=int(q.interval_residual.eq(0).sum()),inside_original_otm_bidask=int(native.interval_residual.eq(0).sum()),canonical_negative_prices=int(q.negative_estimate.sum()),
                native_otm_negative_prices=int(native.negative_estimate.sum()),iv_unavailable=int(q.iv_unavailable.sum()),normalizer=r['normalizer'],allocations=r['allocations'],moments=r['moments'],
                shapes=[dict(instrument=ins,days=int(day),**shape(g,g.model_price_U.to_numpy())) for (ins,day),g in q.groupby(['instrument','expiry_days'])])
    dump(out/f'score_{name}.json',record)
    if r.get('derivatives_computed'):
        np.savez_compressed(out/f'jacobian_{name}.npz',price=r['price_jacobian'],forward=r['forward_jacobian'],residual=r['jacobian'],residuals=r['residuals'])
    return record


def run(out,check_cache,check_controls):
    out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();parent=HERE/'runs/m2_parity_uniform_001'
    cache=HERE/'runs/m2_parity_adaptive_cache_001';controls=HERE/'runs/m2_parity_adaptive_controls_001';prior=json.loads((parent/'run.json').read_text());initial_p=np.array(prior['selected'])
    vix_controls=json.loads((parent/'vix_controls.json').read_text());spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');forwards=pd.read_csv(cache/'forward_targets.csv')
    files=['fit_parity_parameters_v2.py','fit_parity_parameters.py','parity_parameter_constraints.py','test_parity_parameter_constraints.py','parity_parameter_core.py','test_parity_parameter_core.py','parity_adaptive_engine.py','parity_adaptive.py','parity_stock.py','parity_stock_kernel.cpp',
           'parity_hedge_controls.py','controlled_readout_kernel.cpp','controlled_readout_convolution_kernel.cpp','controlled_readout_fast_fit.py','full_readout_fit.py','all_maturity_core.py',
           'joint_curve_core.py','m2_spx.py','PARITY_PARAMETER_FIT_PLAN.md','PARITY_PARAMETER_FIT_PLAN_V2.md','run_parity_adaptive.py']
    sources=[*[HERE/n for n in files],HERE/'runs/m2_parity_parameter_fit_001/jacobian_check.json',HERE/'runs/m2_parity_parameter_fit_001/design.json',parent/'run.json',parent/'prices_canonical_selected.csv',parent/'vix_controls.json',cache/'design.json',controls/'coefficients.json',
             check_cache/'design.json',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(sources);save_sources(out,captured,ROOT)
    lib=out/'convolution.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'controlled_readout_convolution_kernel.cpp'),'-o',str(lib)],check=True)
    slib=out/'stock.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'parity_stock_kernel.cpp'),'-o',str(slib)],check=True)
    a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial,_=dated.initial_state([1,0,0])
    design=dict(status='joint parameter adjustment in progress',source_sha256=source_hashes(captured,ROOT),cache=str(cache),controls=str(controls),check_cache=str(check_cache),check_controls=str(check_controls),
                baseline=initial_p.tolist(),parameter_names=NAMES,bounds_lower=LOWER.tolist(),bounds_upper=UPPER.tolist(),architecture=asdict(a),premium=asdict(premium),initial_state=initial.tolist(),
                paths_per_strike=4096,allocation='uniform',market_otm_inputs=1511,market_out_of_sample=False,fitted_in_this_run=True,
                budget=dict(outer_proposals=3,inner_iterations=12,inner_nfev=20,optimization_seconds=4500,initial_radius=.1,maximum_radius=.3),
                method='VIX-focused SLSQP with both-convention VIX/joint constraints; local SPX/exact VIX guidance and full true4096 acceptance.',
                checking='New independent4096-path/8x512 numerical streams; no check prices used for selection.')
    dump(out/'design.json',design);dump(out/'vix_controls.json',vix_controls)
    e=FitEngine(cache,Kernel(lib),dated,initial,spx,vix,forwards,controls=vix_controls,frozen_controls=controls,stock_values=StockValues(slib))
    print('INITIAL FULL4096 JACOBIAN',flush=True);tick=time.perf_counter();anchor=e.predict(initial_p,True,True);print('INITIAL',anchor['total'],anchor['components'],'seconds',time.perf_counter()-tick,flush=True)
    baseline_record=save_result(e,anchor,initial_p,'baseline',out);baseline_result=anchor;baseline_p=initial_p.copy();anchor_p=initial_p.copy()
    saved=pd.read_csv(parent/'prices_canonical_selected.csv');replay=float(np.max(abs(saved.model_price_U.to_numpy()-anchor['prices'])/saved.price_scale_U.to_numpy()))
    dump(out/'baseline_replay.json',dict(maximum_quote_scale_difference=replay,tolerance=1e-6));assert replay<1e-6,replay
    previous=json.loads((HERE/'runs/m2_parity_parameter_fit_001/jacobian_check.json').read_text())
    assert previous['full_price_relative_error']<.02
    direction=np.array([.01,.002,.03,.004,.02,.04,.02,.001,-.001,.002,.00001]);h=1e-4
    c0=constrained_local(e,initial_p,initial_p,anchor);plus=constrained_local(e,initial_p+h*direction,initial_p,anchor);minus=constrained_local(e,initial_p-h*direction,initial_p,anchor)
    fd=(plus['criteria']-minus['criteria'])/(2*h);analytical=c0['criteria_gradients']@direction
    error=float(np.linalg.norm(fd-analytical)/max(1.,np.linalg.norm(fd)))
    dump(out/'jacobian_check.json',dict(constrained_criteria_relative_error=error,tolerance=.01,full_true_pricer_production_reference=previous,
                                      note='True stock/VIX pricer unchanged; previous full-production directional check retained. New native/canonical constraint derivatives screened here.'))
    assert error<.01,error;print('CONSTRAINED PRODUCTION JACOBIAN',error,flush=True)
    opt_begin=time.perf_counter();radius=.1;inner=[];outer=[];events=[];selected_record=baseline_record
    class InnerBudget(Exception):pass
    for cycle in range(3):
        if time.perf_counter()-opt_begin>4500:break
        lo=np.maximum((LOWER-anchor_p)/SCALE,-radius);hi=np.minimum((UPPER-anchor_p)/SCALE,radius)
        limits=anchor_criteria(anchor);margins=np.maximum(1e-4,np.array([.02,.02,.005,.005])*limits*(radius/.1));limits=limits-margins
        normalization=np.maximum(1.,anchor_criteria(anchor));memo={};best_feasible=[np.inf,None,None];cycle_evaluations=[0]
        def current(z):
            x=anchor_p+SCALE*np.asarray(z);key=tuple(x)
            if key not in memo:
                if cycle_evaluations[0]>=20 or time.perf_counter()-opt_begin>4500:raise InnerBudget('Declared inner/optimization evaluation budget')
                tick=time.perf_counter();r=constrained_local(e,x,anchor_p,anchor);cycle_evaluations[0]+=1;memo.clear();memo[key]=r
                slack=(limits-r['criteria'])/normalization
                if np.min(slack)>=-1e-8 and r['merit']<best_feasible[0]:best_feasible[:]=[r['merit'],x.copy(),r]
                entry=dict(cycle=cycle,index=len(inner),parameters=x.tolist(),surrogate_canonical_loss=r['canonical_joint'],surrogate_native_loss=r['native_joint'],
                           VIX_merit=r['merit'],canonical_VIX=r['canonical_vix'],native_VIX=r['native_vix'],constraint_slack=slack.tolist(),seconds=time.perf_counter()-tick,SPX_prices_are_linearized=True)
                inner.append(entry);dump(out/'inner_evaluations.json',inner);print('INNER CONSTRAINED',cycle,len(inner),r['merit'],r['canonical_joint'],float(np.min(slack)),entry['seconds'],flush=True)
            return memo[key]
        def constraint(z):return (limits-current(z)['criteria'])/normalization
        def constraint_jac(z):return -current(z)['criteria_gradients']*SCALE[None,:]/normalization[:,None]
        fit=None;message=''
        try:
            fit=minimize(lambda z:current(z)['merit'],np.zeros(11),jac=lambda z:current(z)['merit_gradient']*SCALE,
                         bounds=list(zip(lo,hi)),constraints=[dict(type='ineq',fun=constraint,jac=constraint_jac)],method='SLSQP',options=dict(maxiter=12,ftol=1e-6))
            message=str(fit.message)
        except InnerBudget as error:message=str(error)
        if best_feasible[1] is None:
            events.append(dict(cycle=cycle,status='no locally feasible improvement found',radius=radius,inner_nfev=cycle_evaluations[0],message=message));dump(out/'search_events.json',events);radius*=.5
            print('NO FEASIBLE INNER PROPOSAL',cycle,message,flush=True);continue
        candidate=best_feasible[1];local=best_feasible[2];predicted=anchor['total']-local['canonical_joint']
        print('FULL CONSTRAINED PROPOSAL',cycle,'predicted decrease',predicted,flush=True);tick=time.perf_counter();trial=e.predict(candidate,True,True)
        decision=accept_proposal(anchor,trial,predicted);record=save_result(e,trial,candidate,f'proposal{cycle}',out)
        row=dict(cycle=cycle,radius=radius,parameters=candidate.tolist(),canonical_loss=trial['total'],original_otm_loss=trial['original_otm_total'],components=trial['components'],native_components=trial['original_otm_components'],
                 decision=decision,seconds=time.perf_counter()-tick,inner_nfev=cycle_evaluations[0],inner_termination=message,inner_success=bool(fit.success) if fit is not None else False,
                 local_constraint_margins=margins.tolist(),VIX_guidance_merit=best_feasible[0])
        outer.append(row);dump(out/'outer_proposals.json',outer);print('OUTER',json.dumps(row),flush=True)
        if decision['accepted']:
            anchor_p=candidate.copy();anchor=trial;selected_record=record
            if decision['ratio'] is not None and decision['ratio']>.7:radius=min(.3,radius*1.5)
        else:radius*=.5
    # Local guidance may stop at its budget; incumbent remains a true evaluation.
    # Freeze before any checking prediction. Retaining the baseline is allowed.
    selected=anchor_p.copy();assert selected[4]>0
    dump(out/'selected.json',dict(parameters=selected.tolist(),levels=selected[7:10].tolist(),training_loss=anchor['total'],training_components=anchor['components'],checking_not_used=True,
                                 accepted_steps=sum(x['decision']['accepted'] for x in outer),optimizer_seconds=time.perf_counter()-opt_begin))
    selected_record=save_result(e,anchor,selected,'selected',out);e.close();print('SELECTION FROZEN',selected.tolist(),flush=True)
    while not (check_cache/'run.json').exists() or not (check_controls/'run.json').exists():
        if time.perf_counter()-begin>18000:raise TimeoutError('Independent checking preparation exceeded total budget')
        time.sleep(3)
    checking=out/'checking';checking.mkdir();check_design=design|dict(cache=str(check_cache),controls=str(check_controls),selected=selected.tolist(),status='independent checking in progress')
    dump(checking/'design.json',check_design);chk=FitEngine(check_cache,Kernel(lib),dated,initial,spx,vix,forwards,controls=vix_controls,frozen_controls=check_controls,stock_values=StockValues(slib));check_records=[]
    try:
        for name,p in [('baseline',baseline_p),('selected',selected)]:
            print('CHECK',name,flush=True);r=chk.predict(p,True,False);check_records.append(save_result(chk,r,p,name,checking));print('CHECK SCORE',name,r['total'],r['components'],flush=True)
    finally:chk.close()
    check_sources(captured)
    dump(checking/'run.json',check_design|dict(status='complete independent numerical checking',scores=check_records))
    model=json.loads((HERE/'runs/m2_controlled_readout_fit_002/selected_model.json').read_text())
    model['readout']={**model.get('readout',{}),**dict(zip(NAMES[:7],selected[:7].tolist()))};model['variance_floor']=float(selected[10]);model['forward_variance']['levels']=selected[7:10].tolist();model['limitations']=['Bounded local calibration; independent checking and shape/tail diagnostics in run.json.','Price convergence, global optimality and P compatibility unresolved.'];model['status']='new exploratory calibration; see full run and independent checking'
    dump(out/'selected_model.json',model)
    dump(out/'run.json',design|dict(status='completed bounded joint parameter adjustment',selected=selected.tolist(),selected_parameters=selected.tolist(),selected_levels=selected[7:10].tolist(),
                                    training_scores=[baseline_record,selected_record],checking_scores=check_records,outer_proposals=outer,wall_seconds=time.perf_counter()-begin))
    print('FIT COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--check-cache',type=Path,required=True);p.add_argument('--check-controls',type=Path,required=True);a=p.parse_args()
    try:run(a.output.resolve(),a.check_cache.resolve(),a.check_controls.resolve())
    except Exception as error:
        if a.output.exists():dump(a.output/'failure.json',dict(type=type(error).__name__,message=str(error)))
        raise
