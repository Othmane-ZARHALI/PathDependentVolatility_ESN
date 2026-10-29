"""Five-piece refinement at the searching-selected restricted state.

Select on the original searching draws. Evaluate both curves on the same
three independently generated state-run checking ensembles, explicitly reused
across comparisons rather than claimed to be a new independent experiment.
"""
import argparse
from pathlib import Path
import json
import time
import shutil

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .joint_state_fit import HERE, OLD, dump
from .m2_spx import ROOT, M2
from .spx_conditional_control_study import setup
from .joint_state_core import AdaptiveDatedNormalizer, StateDatedPricer
from .joint_curve_fit import SPXData,VIXData
from .curve_refinement_core import FlexibleCurvePathBatch,FlexibleFitData
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.pricing import ForwardVarianceCurve
from esn_eq.vix_pricing import ConditionalVIXPricer
from esn_eq.vix_curve_profile import prepare_curve_exposures


def run(args):
    resume=Path(args.resume_search).resolve() if args.resume_search else None
    args.state_run=args.state_run.resolve()
    args.output=args.output.resolve()
    args.output.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    selected=json.loads((args.state_run/'selected.json').read_text())
    design_state=json.loads((args.state_run/'design.json').read_text())
    state=selected['state'];j=state['index'];initial=np.array(state['initial']);oldlevels=np.array(selected['levels'])
    baseline=oldlevels[[0,0,1,1,2]];statepenalty=.1*state['prior_squared_norm']/3
    targets,p,parameters,zero,oldcurve,panel,master,zeroM,hashes=setup()
    t57=float(master[np.argmin(abs(master-57/365))]);t151=float(master[np.argmin(abs(master-151/365))])
    ends=(t57,.25,t151,.75,oldcurve.ends[-1]);template=ForwardVarianceCurve(ends,tuple(baseline))
    assert abs(t57*365-57)<1e-10 and abs(t151*365-151)<.01
    scratch=Path(design_state['scratch']);oldscratch=Path(design_state['old_scratch'])
    spx=pd.read_csv(OLD/'spx_contracts.csv');vix=pd.read_csv(OLD/'vix_contracts.csv')
    membership=pd.read_csv(OLD/'spx_membership.csv');primary=membership.priority_selected.to_numpy(bool)
    forwards=targets.vix_forwards.sort_values('expiry').reset_index(drop=True)
    assert primary.sum()==707 and len(vix)==175
    sources=[Path(__file__),HERE/'curve_refinement_core.py',HERE/'joint_curve_fit.py',HERE/'joint_curve_core.py',
             HERE/'joint_state_core.py',HERE/'CURVE_REFINEMENT_PLAN.md',HERE/'active_scope.json',
             args.state_run/'selected.json',args.state_run/'design.json',
             *sorted((M2/'src/esn_eq').glob('*.py'))]
    if resume:
        prior=json.load(open(resume/'design.json'))
        assert prior['ends_years']==list(ends) and prior['state']==state
        sources += [resume/'selected.json',resume/'design.json',resume/'optimizer_ledger.json',resume/'failure.json']
    captured=capture_sources(sources);save_sources(args.output,captured,ROOT)
    design=dict(state=state,three_levels=oldlevels.tolist(),embedded_baseline_levels=baseline.tolist(),ends_years=list(ends),
                added_dates_days=[t57*365,t151*365],state_run=str(args.state_run),target_sha256=targets.artifact_sha256,
                implementation_sha256=source_hashes(captured,ROOT),budgets=dict(starts=2,max_nfev=args.max_nfev,
                    search_wall_seconds=1800,paths=512,vix_nets=8,points_per_net=64,checks=3,conditional_GH=32,normalizer_GH=48,time_order=8),
                penalty='Existing 0.1*sum(diff(xi)^2), now over five adjacent levels; same penalty on embedded old curve.',
                checking='Reuse three independent state-run checking ensembles, never use checking scores in optimization or selection. Paired, correlated comparison with state experiment, not a new independent experiment.',
                scope='Unchanged 707 cutoff-3 SPX +175 VIX options +4 forwards. No held-out market pricing.',
                selection='Only original searching paths, fixed state selected using that same searching sample.',
                approximation='New VIX time splitting at added knots: report embedded-old-curve repricing separately from fitting gain.')
    if resume:design.update(resumed_search_from=str(resume),search_repeated=False,
        checking_note='Frozen five-piece selection is unchanged. Resume only on NEW checking ensembles from resumed state run, after its earlier pre-score implementation failure.')
    dump(args.output/'design.json',design);membership.to_csv(args.output/'spx_membership.csv',index=False)
    vix.to_csv(args.output/'vix_contracts.csv',index=False);forwards.to_csv(args.output/'forward_targets.csv',index=False)
    g=AdaptiveDatedNormalizer(p,parameters,initial,float(master[-1]));dump(args.output/'normalizer_audit.json',g.audit)
    if j:
        dated=StateDatedPricer(p.architecture,p.premium);dated.state_normalizer=g
    else:dated=p
    engine=ConditionalVIXPricer(dated,parameters,template,initial,normalizer_order=48)
    windows={days:engine.gaussian_window(float(vix[vix.expiry_days==days].T_years.iloc[0]),time_order=8) for days in [27,55,118,153]}
    def status(stage,**details):
        dump(args.output/'status.json',dict(stage=stage,elapsed_seconds=time.perf_counter()-begin,**details))
        print(stage,json.dumps(details),flush=True)
    def await_file(path):
        start=time.perf_counter()
        while not path.exists():
            if (args.state_run/'failure.json').exists():raise RuntimeError('Parent state experiment failed')
            if time.perf_counter()-start>3600:raise TimeoutError('Parent checking artifact wait exhausted')
            time.sleep(1)
    embedding=[]
    def build(role,trial):
        # Parent may recreate a memory-mapped selected-state ratio while
        # recording its comparison. A CSV sentinel is written only AFTER all
        # candidate paths/exposures finish, preventing reads during truncation.
        await_file(args.state_run/f'prices_{role}_t{trial}_selected.csv')
        ds={};dv={}
        for days in [29,57,211,302]:
            q=spx[spx.expiry_days==days].reset_index(drop=True)
            grid=master[:np.searchsorted(master,q.T_years.iloc[0])+1]
            assert all(k in grid for k in ends if k<grid[-1])
            f=args.state_run/f'controls_{role}_t{trial}_d{days}.npz';await_file(f)
            with np.load(f) as a:controls,weights,beta,stocks=[a[k].copy() for k in ['controls','weights','beta','anchor_stock']]
            batches=[]
            for batch in range(4):
                stem=scratch/f'{role}_t{trial}_d{days}_b{batch}'
                if j==0:
                    ratiofile=Path(str(oldscratch/f'r1_t0_d{days}_b{batch}')+'_ratio.npy') if role=='search' else Path(str(stem)+'_zero_ratio.npy')
                else:ratiofile=Path(str(stem)+f'_state{j}_ratio.npy')
                await_file(ratiofile)
                # The parent finishes each selected comparison before these are
                # consumed; wait for its VIX archive below before using stocks.
                batches.append(FlexibleCurvePathBatch(np.load(ratiofile,mmap_mode='r'),np.load(str(stem)+'_equity.npy',mmap_mode='r'),
                    np.diff(grid),np.minimum(np.searchsorted(ends,grid[:-1],side='right'),4),5))
            d=SPXData(q,batches,weights,controls,stocks);d.beta=beta;ds[days]=d
        for days in [27,55,118,153]:
            f=args.state_run/f'vix_state{j}_{role}_t{trial}_d{days}.npz';await_file(f)
            with np.load(f) as a:x,w,c,beta,oldA,oldfloor,delta=[a[k].copy() for k in ['states','weights','controls','beta','exposures','floor_integral','delta']]
            status('Five-piece VIX exposure',role=role,trial=trial,days=days)
            if resume and role=='search':
                with np.load(resume/f'vix_search_t0_d{days}.npz') as old5:
                    np.testing.assert_array_equal(x,old5['states']);np.testing.assert_array_equal(w,old5['weights'])
                    coef=old5['exposures'].copy();floor_integral=float(old5['floor_integral']);window_delta=float(old5['delta'])
            else:
                cache=prepare_curve_exposures(engine,x,windows[days],order=32)
                coef=cache.coefficients;floor_integral=cache.floor_integral;window_delta=cache.delta
            old_integral=float(oldfloor)+oldA@(oldlevels-1e-4)
            embedded=floor_integral+coef@(baseline-1e-4)
            oldv=100*np.sqrt(old_integral/float(delta));newv=100*np.sqrt(embedded/float(delta))
            b=(w*(newv-oldv)).reshape(8,-1).mean(axis=1)
            embedding.append(dict(role=role,trial=trial,days=days,max_vix_change=float(abs(newv-oldv).max()),
                weighted_forward_change=float(b.mean()),paired_scramble_se=float(b.std(ddof=1)/np.sqrt(8)),
                max_relative_integral_change=float(np.max(abs(embedded/old_integral-1))),
                interpretation='Same embedded variance curve; differences are quadrature/time-splitting, not fitting.'))
            q=vix[vix.expiry_days==days].reset_index(drop=True)
            dv[days]=VIXData(q,coef,floor_integral,window_delta,w,c,beta,8)
            np.savez_compressed(args.output/f'vix_{role}_t{trial}_d{days}.npz',states=x,weights=w,controls=c,beta=beta,
                exposures=coef,floor_integral=floor_integral,delta=window_delta)
        pd.DataFrame(embedding).to_csv(args.output/'embedding_quadrature_changes.csv',index=False)
        return FlexibleFitData(ds,dv,spx,vix,forwards,primary,5)
    train=build('search',0);search_start=time.perf_counter();memo={};evaluations=[];ledger=[];candidates=[]
    def evaluate(x):
        key=np.asarray(x).tobytes()
        if key not in memo:
            if time.perf_counter()-search_start>1800:raise TimeoutError('Curve search budget exhausted')
            value=train.predict(x);memo.clear();memo[key]=value
            row=dict(evaluation=len(evaluations),levels=np.asarray(x).tolist(),total=value['total']+statepenalty,
                components=value['components']|dict(state_penalty=statepenalty),negative_count=value['negative_count'])
            evaluations.append(row)
            with (args.output/'evaluations.jsonl').open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n')
        return memo[key]
    r=evaluate(baseline);candidates.append((r['total'],baseline))
    for number,start in enumerate([] if resume else [baseline,np.full(5,.06)]):
        fit=least_squares(lambda x:evaluate(x)['residuals'],start,jac=lambda x:evaluate(x)['jacobian'],
            bounds=(np.full(5,.00010001),np.full(5,.25)),max_nfev=args.max_nfev,ftol=1e-6,xtol=1e-6,gtol=1e-5)
        r=evaluate(fit.x);candidates.append((r['total'],fit.x.copy()))
        ledger.append(dict(start=number,levels=fit.x.tolist(),total=r['total']+statepenalty,success=bool(fit.success),nfev=fit.nfev,message=fit.message))
    if resume:
        prior_selection=json.load(open(resume/'selected.json'));levels=np.array(prior_selection['levels'])
        ledger=json.load(open(resume/'optimizer_ledger.json'))
        assert prior_selection['state']==state and prior_selection['baseline_levels']==baseline.tolist()
    else:_,levels=min(candidates,key=lambda x:x[0])
    dump(args.output/'selected.json',dict(levels=levels.tolist(),baseline_levels=baseline.tolist(),state=state,
        selected_using='search only',checking_not_used=True,
        search_seconds=prior_selection['search_seconds'] if resume else time.perf_counter()-search_start,
        resumed_search_from=str(resume) if resume else None,selection_repeated=False if resume else None))
    dump(args.output/'optimizer_ledger.json',ledger)
    status('Five-piece selection frozen',levels=levels.tolist())
    records=[];blocks={};paths={}
    def record(data,role,trial):
        for name,x in [('baseline',baseline),('selected',levels)]:
            r=data.predict(x,derivative=False);q,f=data.loss.tables(r['prices'][data.primary],r['errors'][data.primary],r['forwards'],r['forward_errors'])
            q.to_csv(args.output/f'prices_{role}_t{trial}_{name}.csv',index=False);f.to_csv(args.output/f'forwards_{role}_t{trial}_{name}.csv',index=False)
            full=data.all_quotes.copy();full['model_price_U']=r['prices'];full['mc_se_U']=r['errors'];full['priority_objective_member']=data.primary
            full.to_csv(args.output/f'full_training_{role}_t{trial}_{name}.csv',index=False)
            records.append(dict(role=role,trial=trial,candidate=name,levels=x.tolist(),total=r['total']+statepenalty,
                components=r['components']|dict(state_penalty=statepenalty),negative_count=r['negative_count']))
            blocks[f'{role}_t{trial}_{name}_prices']=r['blocks'];blocks[f'{role}_t{trial}_{name}_forwards']=r['forward_blocks']
            for key,array in r['archive'].items():paths[f'{role}_t{trial}_{name}_{key}']=array
            status('Refinement comparison',role=role,trial=trial,candidate=name,total=r['total']+statepenalty)
    record(train,'search',0)
    for trial in range(3):record(build('fresh',trial),'fresh',trial)
    check_sources(captured);np.savez_compressed(args.output/'comparison_blocks.npz',**blocks)
    np.savez_compressed(args.output/'comparison_paths.npz',**paths);dump(args.output/'scores.json',records)
    dump(args.output/'run.json',design|dict(status='completed exploratory five-piece refinement',selected_levels=levels.tolist(),scores=records,
        wall_seconds=time.perf_counter()-begin,normalizer_audit=g.audit))
    status('completed')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--state-run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--max-nfev',type=int,default=30)
    parser.add_argument('--resume-search',type=Path)
    args=parser.parse_args()
    try:run(args)
    except Exception as e:
        if args.output.exists():dump(args.output/'failure.json',dict(type=type(e).__name__,message=str(e)))
        raise
