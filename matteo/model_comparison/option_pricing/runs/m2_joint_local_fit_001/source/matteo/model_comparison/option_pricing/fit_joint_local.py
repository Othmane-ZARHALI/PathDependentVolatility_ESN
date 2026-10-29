"""Point1 pilot: cheap joint SPX/VIX local search, full-price acceptance."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT, M2, frozen_stage3
from .full_readout_fit import HERE, NAMES, Kernel, dump
from .calibration_policy import LOWER, UPPER, validate
from .joint_price_local import JointPriceLocal
from .joint_local_search import propose
from .parity_parameter_core import FitEngine, accept_proposal
from .parity_parameter_constraints import anchor_criteria
from .parity_stock import StockValues
from .fit_parity_parameters_v3 import save_result
from .provenance import capture_sources, source_hashes, save_sources, check_sources
from esn_eq.pricing import DatedM2Pricer


def run(out, proposals=2):
    out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    prior=HERE/'runs/m2_parity_parameter_fit_003'
    cache=HERE/'runs/m2_parity_adaptive_cache_001';controls=HERE/'runs/m2_parity_adaptive_controls_001'
    old=json.loads((prior/'run.json').read_text());baseline=np.clip(np.array(old['selected']),LOWER,UPPER);validate(baseline)
    spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');forwards=pd.read_csv(cache/'forward_targets.csv')
    native=pd.concat([spx,vix],ignore_index=True)
    vc=json.loads((prior/'vix_controls.json').read_text())
    files=['fit_joint_local.py','joint_price_local.py','joint_local_search.py','calibration_policy.py','test_joint_price_local.py','POINT1_ACCELERATION_PLAN.md',
           'fit_parity_parameters_v3.py','parity_parameter_core.py','parity_parameter_constraints.py','parity_adaptive_engine.py','parity_adaptive.py',
           'parity_stock.py','parity_stock_kernel.cpp','parity_hedge_controls.py','controlled_readout_kernel.cpp','controlled_readout_convolution_kernel.cpp',
           'controlled_readout_fast_fit.py','full_readout_fit.py','all_maturity_core.py','joint_curve_core.py','m2_spx.py','run_parity_adaptive.py']
    sources=[*[HERE/n for n in files],prior/'run.json',prior/'vix_controls.json',cache/'design.json',controls/'coefficients.json',
             cache/'spx_contracts.csv',cache/'vix_contracts.csv',cache/'forward_targets.csv',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(sources);save_sources(out,captured,ROOT)
    lib=out/'convolution.dylib';slib=out/'stock.dylib'
    for source,target in [('controlled_readout_convolution_kernel.cpp',lib),('parity_stock_kernel.cpp',slib)]:
        subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/source),'-o',str(target)],check=True)
    a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial,_=dated.initial_state([1,0,0])
    design=dict(status='full baseline in progress',baseline=baseline.tolist(),parameter_names=NAMES,bounds_lower=LOWER.tolist(),bounds_upper=UPPER.tolist(),
                minimum_volatility=.05,maximum_floor_volatility=.06,prior_parameters=old['selected'],cache=str(cache),controls=str(controls),
                architecture=asdict(a),premium=asdict(premium),initial_state=initial.tolist(),paths_per_strike=4096,allocation='uniform',market_otm_inputs=1511,
                source_sha256=source_hashes(captured,ROOT),checking='Training-ensemble implementation pilot; no new independent checking priced.',
                market_out_of_sample=False,method='Affine joint SPX/VIX price/forward trust model; original loss and full true-price guards.',
                budget=dict(full_proposals=proposals,seconds_after_baseline=1200,initial_radius=.02))
    dump(out/'design.json',design);dump(out/'vix_controls.json',vc)
    e=FitEngine(cache,Kernel(lib),dated,initial,spx,vix,forwards,controls=vc,frozen_controls=controls,stock_values=StockValues(slib))
    try:
        print('FULL BASELINE: exact new5% bound',baseline.tolist(),flush=True);tick=time.perf_counter()
        anchor=e.predict(baseline,True,True);baseline_seconds=time.perf_counter()-tick
        base_record=save_result(e,anchor,baseline,'baseline',out);anchor_p=baseline.copy();selected_record=base_record
        print('BASELINE',anchor['total'],'seconds',baseline_seconds,flush=True)
        outer=[];local_runs=[];events=[];radius=.02;opt_begin=time.perf_counter();full_calls=1
        for cycle in range(proposals):
            if time.perf_counter()-opt_begin>1200:break
            model=JointPriceLocal(anchor_p,anchor,e.quotes,native,forwards)
            tick=time.perf_counter()
            for _ in range(200):model.evaluate(anchor_p)
            latency=(time.perf_counter()-tick)/200
            local=propose(model,radius=radius)
            entry={k:v for k,v in local.items() if k not in ['parameters','prediction']}
            entry|=dict(cycle=cycle,per_evaluation_seconds=latency,parameters=None if local['parameters'] is None else local['parameters'].tolist(),
                        anchor_parameters=anchor_p.tolist(),anchor_criteria=anchor_criteria(anchor).tolist(),
                        approximate_criteria=None if local['prediction'] is None else local['prediction']['criteria'].tolist())
            local_runs.append(entry);dump(out/'local_searches.json',local_runs)
            print('LOCAL',cycle,local['seconds'],'seconds',local['nfev'],'evaluations',entry['approximate_criteria'],flush=True)
            if local['parameters'] is None:
                events.append(dict(cycle=cycle,reason='No local reserve-feasible candidate',radius=radius));radius*=.5;continue
            candidate=validate(local['parameters']);predicted=anchor['total']-local['prediction']['criteria'][0]
            if predicted<=1e-8:
                events.append(dict(cycle=cycle,reason='No predicted canonical decrease',radius=radius));radius*=.5;continue
            print('FULL PROPOSAL',cycle,candidate.tolist(),flush=True);tick=time.perf_counter()
            trial=e.predict(candidate,True,True);seconds=time.perf_counter()-tick;full_calls+=1
            decision=accept_proposal(anchor,trial,predicted)
            rec=save_result(e,trial,candidate,f'proposal{cycle}',out)
            true_criteria=anchor_criteria(trial)
            pred_nv=anchor_criteria(anchor)[3]-local['prediction']['criteria'][3]
            nv_ratio=float((anchor_criteria(anchor)[3]-true_criteria[3])/pred_nv) if pred_nv>0 else None
            row=dict(cycle=cycle,parameters=candidate.tolist(),radius=radius,seconds=seconds,decision=decision,
                     canonical_loss=trial['total'],original_otm_loss=trial['original_otm_total'],components=trial['components'],native_components=trial['original_otm_components'],
                     approximate_criteria=local['prediction']['criteria'].tolist(),true_criteria=true_criteria.tolist(),native_vix_agreement_ratio=nv_ratio,
                     maximum_price_prediction_error_in_quote_scales=float(np.max(abs(trial['prices']-local['prediction']['prices'])/e.quotes.price_scale_U.to_numpy())))
            outer.append(row);dump(out/'outer_proposals.json',outer);print('TRUE DECISION',json.dumps(row),flush=True)
            if decision['accepted']:
                anchor_p=candidate.copy();anchor=trial;selected_record=rec
                agreement=min(decision['ratio'],nv_ratio if nv_ratio is not None else decision['ratio'])
                radius=min(.05,radius*1.4) if agreement>.75 else (radius*.5 if agreement<.25 else radius)
            else:radius*=.5
        validate(anchor_p);selected_record=save_result(e,anchor,anchor_p,'selected',out)
        dump(out/'search_events.json',events)
        dump(out/'selected.json',dict(parameters=anchor_p.tolist(),levels=anchor_p[7:10].tolist(),training_loss=anchor['total'],accepted_steps=sum(x['decision']['accepted'] for x in outer),
                                     minimum_instantaneous_volatility=float(np.sqrt(anchor_p[10])),checking_not_used=True))
        config=json.loads((prior/'selected_model.json').read_text());config['readout'].update(dict(zip(NAMES[:7],anchor_p[:7].tolist())))
        config['variance_floor']=float(anchor_p[10]);config['forward_variance']['levels']=anchor_p[7:10].tolist()
        config['calibration_policy']=dict(minimum_volatility=.05,floor_volatility_bounds=[.05,.06],beta_J_positive=True)
        config['status']='Exploratory joint-local acceleration pilot; training-ensemble full-price verification only'
        config['limitations']=['No new independent checking; no convergence/global optimality claim.','Prior SPX signed-estimator failures and VIX smile misfit remain unresolved.']
        dump(out/'selected_model.json',config);check_sources(captured)
        previous_inner=json.loads((prior/'inner_evaluations.json').read_text())
        timing=dict(full_baseline_seconds=baseline_seconds,full_proposal_seconds=sum(x['seconds'] for x in outer),full_pricer_calls=full_calls,
                    local_search_seconds=sum(x['seconds'] for x in local_runs),local_evaluations=sum(len(x['history']) for x in local_runs),
                    local_evaluation_median_seconds=float(np.median([x['per_evaluation_seconds'] for x in local_runs])),
                    historical_nonlinear_VIX_local_median_seconds=float(np.median([x['seconds'] for x in previous_inner])),
                    caveat='Approximate vs nonlinear computations. Local latency speed ratio is not end-to-end speedup; model/floor and search differ.')
        dump(out/'timing.json',timing)
        dump(out/'run.json',design|dict(status='completed full-price-guarded acceleration pilot',selected=anchor_p.tolist(),selected_parameters=anchor_p.tolist(),selected_levels=anchor_p[7:10].tolist(),
                                       training_scores=[base_record,selected_record],outer_proposals=outer,timing=timing,wall_seconds=time.perf_counter()-begin))
        print('COMPLETE',timing,flush=True)
    finally:e.close()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--proposals',type=int,default=2);a=p.parse_args()
    try:run(a.output.resolve(),a.proposals)
    except Exception as error:
        if a.output.exists():dump(a.output/'failure.json',dict(type=type(error).__name__,message=str(error)))
        raise
