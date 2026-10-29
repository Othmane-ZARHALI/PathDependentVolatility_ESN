"""Resume paused run001 from its complete frozen baseline; retain original sources."""
from .parity_checkpoint import load_result
import hashlib
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .full_readout_fit import HERE,NAMES,dump
from .maturity_curve import ENDS_DAYS,CurveKernel,CurveEngine,embed_three_piece,bounds,validate,vix_only_loss
from .maturity_curve_local import CurvePriceLocal
from .maturity_curve_search import propose
from .parity_parameter_core import accept_proposal
from .parity_parameter_constraints import anchor_criteria
from .fit_parity_parameters_v3 import save_result
from .parity_adaptive import canonical_market
from .all_maturity_core import InSampleLoss
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from esn_eq.pricing import DatedM2Pricer


def vix_criteria(engine,p,anchor):
    v=engine.vix_prediction(p);n=len(engine.spx);prices=anchor['prices'].copy();prices[n:]=v['prices']
    native=anchor['native_prices'].copy();native[n:]=v['prices'];index=0
    for j,(_,g) in enumerate(engine.vix.groupby('expiry_days',sort=True)):
        m=len(g);native[n+index:n+index+m]-=g.option_type.eq('put').to_numpy()*(v['forwards'][j]-g.strike.to_numpy());index+=m
    c=engine.loss.evaluate(prices,v['forwards'],p[7:-1]);o=engine.otm_loss.evaluate(native,v['forwards'],p[7:-1])
    return np.array([c['components']['VIX']+c['components']['forwards'],o['components']['VIX']+o['components']['forwards']])


def curve_grid_audit(engine,p):
    records=[];left=np.r_[0,engine.ends[:-1]]
    for day,(folders,_) in engine.stock.items():
        grid=np.load(folders[0]/'grid.npy');seg=np.minimum(np.searchsorted(engine.ends,grid[:-1],side='right'),engine.ncurve-1)
        measured=np.bincount(seg,weights=np.diff(grid),minlength=engine.ncurve)
        exact=np.maximum(0,np.minimum(engine.ends,grid[-1])-left)
        records.append(dict(days=int(day),maximum_step_days=float(np.max(np.diff(grid))*365),
                            left_integrated_variance=float(measured@p[7:-1]),exact_integrated_variance=float(exact@p[7:-1]),
                            difference=float((measured-exact)@p[7:-1]),note='Deterministic curve integration only; not a full stock timestep error bound.'))
    return records


def production_curve_derivatives(engine,p,anchor):
    n=len(engine.spx);base=engine.vix_prediction(p)
    np.testing.assert_allclose(base['prices'],anchor['prices'][n:],rtol=1e-11,atol=1e-11)
    np.testing.assert_allclose(base['jacobian'],anchor['price_jacobian'][n:],rtol=1e-10,atol=1e-9)
    _,_,scale=bounds(engine.ncurve);rows=[]
    for j in range(7,len(p)):
        h=2e-7*scale[j];plus=p.copy();minus=p.copy();plus[j]+=h;minus[j]-=h
        a=engine.vix_prediction(plus);b=engine.vix_prediction(minus)
        fd=np.r_[(a['prices']-b['prices'])/(2*h),(a['forwards']-b['forwards'])/(2*h)]
        jac=np.r_[base['jacobian'][:,j],base['forward_jacobian'][:,j]]
        err=float(np.linalg.norm(fd-jac)/max(1.,np.linalg.norm(jac)));rows.append(dict(column=j,relative_error=err));assert err<.005,(j,err)
    return rows


def run(out):
    assert (out/'paused.json').exists() and not (out/'run.json').exists()
    frozen_design=json.loads((out/'design.json').read_text())
    for name,digest in frozen_design['source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    begin=time.perf_counter();prior=HERE/'runs/m2_joint_local_fit_001'
    old=json.loads((prior/'run.json').read_text());oldp=np.array(old['selected']);p=embed_three_piece(oldp);validate(p)
    cache=HERE/'runs/m2_parity_adaptive_cache_001';controls=HERE/'runs/m2_parity_adaptive_controls_001'
    check_cache=HERE/'runs/m2_parity_check_cache_001';check_controls=HERE/'runs/m2_parity_check_controls_001'
    spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv');vc=json.loads((prior/'vix_controls.json').read_text())
    files=['resume_maturity_curve.py','parity_checkpoint.py','fit_maturity_curve.py','maturity_curve.py','maturity_curve_predict.py','maturity_curve_stock.cpp','maturity_curve_local.py','maturity_curve_search.py','test_maturity_curve.py','MATURITY_CURVE_FIT_PLAN.md',
           'parity_adaptive_engine.py','parity_adaptive.py','parity_hedge_controls.py','parity_parameter_core.py','parity_parameter_constraints.py','calibration_policy.py','forward_variance_constraints.py',
           'full_readout_fit.py','controlled_readout_fast_fit.py','controlled_readout_kernel.cpp','controlled_readout_convolution_kernel.cpp','all_maturity_core.py','joint_curve_core.py','m2_spx.py','fit_parity_parameters_v3.py','run_parity_adaptive.py']
    captured=capture_sources([*[HERE/n for n in files],prior/'run.json',prior/'vix_controls.json',cache/'design.json',controls/'coefficients.json',check_cache/'design.json',check_controls/'coefficients.json',
                              cache/'spx_contracts.csv',cache/'vix_contracts.csv',cache/'forward_targets.csv',*sorted((M2/'src/esn_eq').glob('*.py'))]);save_sources(out/'continuation',captured,ROOT)
    conv=out/'convolution.dylib';stock=out/'curve.dylib'
    assert conv.exists() and stock.exists()
    a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial,_=dated.initial_state([1,0,0]);ends=ENDS_DAYS/365;lo,hi,_=bounds(len(ends))
    names=NAMES[:7]+[f'xi_{i+1:02d}' for i in range(len(ends))]+['variance_floor']
    design=dict(status='nested baseline in progress',baseline=p.tolist(),prior_parameters=oldp.tolist(),parameter_names=names,ends_days=ENDS_DAYS.tolist(),ends_years=ends.tolist(),
                bounds_lower=lo.tolist(),bounds_upper=hi.tolist(),architecture=asdict(a),premium=asdict(premium),initial_state=initial.tolist(),reference_levels=oldp[7:10].tolist(),
                source_sha256=source_hashes(captured,ROOT),cache=str(cache),controls=str(controls),check_cache=str(check_cache),check_controls=str(check_controls),paths_per_strike=4096,
                market_otm_inputs=1511,stock_steps_per_trading_day=512,method='Curve-first then all24 local joint-price search; exact VIX precheck; full SPX/VIX acceptance.',
                budget=dict(proposals=3,seconds_after_baseline=1800),checking='VIX-only separate reused numerical streams after selection; not pristine holdout/market OOS.')
    dump(out/'resume_design.json',design);dump(out/'resume_status.json',dict(status='resumed by user',date='2026-10-25',baseline_reused=True))
    def engine(c,ctrl):return CurveEngine(c,CurveKernel(conv,stock,ends),dated,initial,spx,vix,f,controls=vc,frozen_controls=ctrl,ends=ends,reference_levels=oldp[7:10])
    e=engine(cache,controls);full_calls=0
    try:
        print('RESTORE COMPLETE BASELINE',len(p),'parameters',flush=True)
        anchor=load_result(out,'baseline');base=json.loads((out/'score_baseline.json').read_text())
        np.testing.assert_array_equal(p,np.array(base['parameters']))
        anchor_p=p.copy();selected=base;full_calls=1;baseline_seconds=341.46298262500204
        print('BASELINE RESTORED',anchor['total'],anchor['original_otm_total'],flush=True)
        outer=[];searches=[];screens=[];radius=.1;opt_begin=time.perf_counter()
        for cycle in range(3):
            if time.perf_counter()-opt_begin>1800:break
            curve_only=cycle==0
            if cycle==1:radius=.05
            local=propose(CurvePriceLocal(anchor_p,anchor,e.quotes,e.original_quotes,f),radius=radius,curve_only=curve_only)
            entry={k:v for k,v in local.items() if k not in ['parameters','prediction']};entry|=dict(cycle=cycle,anchor_parameters=anchor_p.tolist(),parameters=None if local['parameters'] is None else local['parameters'].tolist(),
                       approximate_criteria=None if local['prediction'] is None else local['prediction']['criteria'].tolist());
            if cycle==0:
                prior_search=json.loads((out/'local_searches.json').read_text())[0]
                np.testing.assert_allclose(entry['parameters'],prior_search['parameters'],rtol=0,atol=1e-12)
            searches.append(entry);dump(out/'local_searches.json',searches)
            print('LOCAL',cycle,local['seconds'],local['nfev'],entry['approximate_criteria'],flush=True)
            if local['parameters'] is None:radius*=.5;continue
            candidate=validate(local['parameters']);predicted=anchor['total']-local['prediction']['criteria'][0]
            if predicted<=1e-8:radius*=.5;continue
            tick=time.perf_counter();screen=vix_criteria(e,candidate,anchor);eligible=bool(np.all(screen<anchor_criteria(anchor)[2:]-1e-8))
            screens.append(dict(cycle=cycle,parameters=candidate.tolist(),true_VIX_criteria=screen.tolist(),eligible=eligible,seconds=time.perf_counter()-tick));dump(out/'vix_screens.json',screens)
            print('EXACT VIX SCREEN',cycle,screen.tolist(),eligible,flush=True)
            if not eligible:radius*=.5;continue
            print('FULL PROPOSAL',cycle,flush=True);tick=time.perf_counter();trial=e.predict(candidate,True,True);full_calls+=1;seconds=time.perf_counter()-tick
            decision=accept_proposal(anchor,trial,predicted);record=save_result(e,trial,candidate,f'proposal{cycle}',out)
            row=dict(cycle=cycle,parameters=candidate.tolist(),radius=radius,curve_only=curve_only,decision=decision,seconds=seconds,approximate_criteria=local['prediction']['criteria'].tolist(),
                     true_criteria=anchor_criteria(trial).tolist(),canonical_loss=trial['total'],original_otm_loss=trial['original_otm_total'])
            outer.append(row);dump(out/'outer_proposals.json',outer);print('TRUE DECISION',json.dumps(row),flush=True)
            if decision['accepted']:
                anchor=trial;anchor_p=candidate.copy();selected=record
                radius=min(.1,radius*1.4) if decision['ratio']>.75 else (radius*.5 if decision['ratio']<.25 else radius)
            else:radius*=.5
        validate(anchor_p);selected=save_result(e,anchor,anchor_p,'selected',out);dump(out/'grid_selected.json',curve_grid_audit(e,anchor_p))
        dump(out/'selected.json',dict(parameters=anchor_p.tolist(),levels=anchor_p[7:-1].tolist(),ends_days=ENDS_DAYS.tolist(),accepted_steps=sum(x['decision']['accepted'] for x in outer),checking_not_used=True))
        print('SELECTION FROZEN',anchor_p.tolist(),flush=True)
    finally:e.close()
    check=out/'checking_vix';check.mkdir();ce=engine(check_cache,check_controls)
    ce.spx=spx.iloc[:0];ce.original_quotes=vix.reset_index(drop=True);ce.quotes=canonical_market(ce.original_quotes);ce.loss=vix_only_loss(ce.quotes,f);ce.otm_loss=vix_only_loss(ce.original_quotes,f)
    checking=[]
    try:
        for name,x in [('baseline',p),('selected',anchor_p)]:
            print('VIX CHECK',name,flush=True);value=ce.predict(x,True,True);checking.append(save_result(ce,value,x,name,check));print('VIX CHECK SCORE',name,value['components'],value['original_otm_components'],flush=True)
    finally:ce.close()
    config=json.loads((prior/'selected_model.json').read_text());config['readout'].update(dict(zip(NAMES[:7],anchor_p[:7].tolist())))
    config['variance_floor']=float(anchor_p[-1]);config['forward_variance']=dict(ends_years=ends.tolist(),levels=anchor_p[7:-1].tolist())
    config['status']='Exploratory maturity-curve calibration; VIX-only separate numerical checking'
    config['limitations']=['Not global optimum or identified curve.','SPX checking not repeated; price convergence and VIX shape fit unresolved.']
    dump(out/'selected_model.json',config);check_sources(captured)
    dump(out/'run.json',design|dict(status='completed maturity-curve calibration attempt',selected=anchor_p.tolist(),selected_parameters=anchor_p.tolist(),selected_levels=anchor_p[7:-1].tolist(),
                                   training_scores=[base,selected],checking_scores=checking,outer_proposals=outer,full_stock_calls=full_calls,baseline_seconds=baseline_seconds,wall_seconds=time.perf_counter()-begin))
    dump(check/'run.json',dict(baseline=p.tolist(),selected=anchor_p.tolist(),scores=checking,status='VIX-only reused separate-stream numerical checking'))
    print('COMPLETE',time.perf_counter()-begin,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    try:run(args.output.resolve())
    except Exception as error:
        if args.output.exists():dump(args.output/'failure.json',dict(type=type(error).__name__,message=str(error)))
        raise
