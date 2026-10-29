"""Raw stock/mass and controlled VIX-squared checks, separate from fit scores."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .joint_state_fit import OLD,VIX_PARENT,HERE
from .spx_conditional_control_study import setup
from .curve_refinement_core import FlexibleCurvePathBatch
from esn_eq.pricing import ForwardVarianceCurve
from matteo.preprocessing.option_variance import vix_option_moment_bound


def stats(x):return float(np.mean(x)),float(np.std(x,ddof=1)/np.sqrt(len(x)))


def audit(run):
    record=json.load(open(run/'run.json'));isstate='states' in record
    if isstate:
        parent=run;state_indices={'baseline':0,'selected':record['selected']['state']['index']}
        levels={'baseline':record['baseline_levels'],'selected':record['selected']['levels']}
        reference=json.load(open(run/'source'/OLD.relative_to(HERE.parents[2])/'run.json'))
        ends=reference['family']['ends']
    else:
        parent=Path(record['state_run']);state_indices={k:record['state']['index'] for k in ['baseline','selected']}
        levels={'baseline':record['embedded_baseline_levels'],'selected':record['selected_levels']};ends=record['ends_years']
    parent_record=json.load(open(parent/'run.json'));scratch=Path(parent_record['scratch'])
    _,_,_,_,_,_,master,_,_=setup()
    spx=pd.read_csv(run/'spx_membership.csv');vix=pd.read_csv(run/'vix_contracts.csv')
    full_beta=json.load(open(VIX_PARENT/'coefficients.json'))
    stock_rows=[];moment_rows=[];pools={}
    for trial in range(3):
        for days,group in spx.groupby('expiry_days'):
            grid=master[:np.searchsorted(master,group.T_years.iloc[0])+1]
            weights=np.load(parent/f'controls_fresh_t{trial}_d{days}.npz')['weights']
            for name,xi in levels.items():
                j=state_indices[name];ys=[]
                for batch in range(4):
                    stem=scratch/f'fresh_t{trial}_d{days}_b{batch}'
                    ratio=Path(str(stem)+('_zero_ratio.npy' if j==0 else f'_state{j}_ratio.npy'))
                    path=FlexibleCurvePathBatch(np.load(ratio,mmap_mode='r'),np.load(str(stem)+'_equity.npy',mmap_mode='r'),
                        np.diff(grid),np.minimum(np.searchsorted(ends,grid[:-1],side='right'),len(xi)-1),len(xi))
                    ys.append(path.terminal(xi,derivative=False)[0])
                for target,values in [('mass',weights),('stock',weights*np.concatenate(ys))]:
                    mean,se=stats(values);stock_rows.append(dict(trial=trial,days=int(days),candidate=name,target=target,
                        estimate=mean,se=se,expected=1.,warning_2se=bool(abs(mean-1)>2*se)))
        for days,group in vix.groupby('expiry_days'):
            for name,xi in levels.items():
                j=state_indices[name]
                path=run/f'vix_state{j}_fresh_t{trial}_d{days}.npz' if isstate else run/f'vix_fresh_t{trial}_d{days}.npz'
                with np.load(path) as a:
                    squared=10000*(float(a['floor_integral'])+a['exposures']@(np.asarray(xi)-1e-4))/float(a['delta'])
                    beta=np.asarray(full_beta[str(days)]['coefficients']['gaussian'])[:,1]
                    b=(a['weights']*squared-a['controls']@beta).reshape(8,64).mean(axis=1)
                    delta=float(a['delta'])
                curve=ForwardVarianceCurve(tuple(ends),tuple(xi));T=float(group.T_years.iloc[0]);target=10000*curve.integral(T,T+delta)/delta
                mean,se=stats(b);pools.setdefault((int(days),name),[]).append(b)
                moment_rows.append(dict(trial=trial,days=int(days),candidate=name,model_target=target,estimate=mean,se=se,
                    residual=mean-target,warning_2se=bool(abs(mean-target)>2*se)))
    moment_pool=[];compat=[]
    for (days,name),arrays in pools.items():
        mean,se=stats(np.concatenate(arrays));target=next(r['model_target'] for r in moment_rows if r['days']==days and r['candidate']==name)
        moment_pool.append(dict(days=days,candidate=name,model_target=target,estimate=mean,se=se,residual=mean-target,warning_2se=bool(abs(mean-target)>2*se)))
        f=pd.read_csv(run/f'forwards_fresh_pooled_{name}.csv');r=f[f.expiry_days==days].iloc[0]
        bound=vix_option_moment_bound(vix[vix.expiry_days==days],float(r.forward_vix_points),float(r.model_forward))
        compat.append(dict(days=days,candidate=name,model_mean=float(r.model_forward),forward_se=float(r.mc_se),
            model_second_moment=target,market_bid_necessary_bound=bound['moment_lower_bound'],necessary_margin=target-bound['moment_lower_bound'],
            jensen_margin=target-float(r.model_forward)**2))
    for name,data in [('stock_forward_checks',stock_rows),('vix_moment_checks',moment_rows),('vix_moment_pooled',moment_pool),('necessary_moment_compatibility',compat)]:
        pd.DataFrame(data).to_csv(run/(name+'.csv'),index=False)
    report=dict(stock_warning_rows=sum(r['warning_2se'] for r in stock_rows),pooled_vix_moment_warnings=[r for r in moment_pool if r['warning_2se']],
                negative_necessary_margins=[r for r in compat if r['necessary_margin']<0],
                caveat='Dependent numerical diagnostics at estimated model means. Necessary market bounds plus moment/forward uncertainty and unresolved numerical bias; no certified infeasibility or price convergence.')
    (run/'identity_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True)
    audit(parser.parse_args().run)
