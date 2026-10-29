"""Descriptive MLMC efficiency, independent-price consistency and finest bias checks."""
import argparse
import json
from pathlib import Path
from .spx_mlmc import draw_estimators
import numpy as np
import pandas as pd


def audit(output):
    record=json.loads((output/'run.json').read_text());allocation=json.loads((output/'allocation.json').read_text())
    e=pd.read_csv(output/'estimates.csv');levels=pd.read_csv(output/'levels.csv');costs=pd.read_csv(output/'costs.csv');d=pd.read_csv(output/'diagnostics.csv')
    pooled=e[e.estimator=='pooled'];central=[2,3,8,9]
    result=dict(allocation=allocation['counts'],predicted_mlmc_cost_seconds=allocation['predicted_cost_seconds'],
        predicted_direct_budget_seconds=allocation['budget_seconds'],trials=[],sampling=[],finest_corrections=[],pilot_stability=[])
    for trial in range(record['arguments']['trials']):
        m=pooled[(pooled.trial==trial)&(pooled.method=='mlmc')];direct=pooled[(pooled.trial==trial)&(pooled.method=='direct')]
        table=m.merge(direct,on=['trial','contract'],suffixes=('_mlmc','_direct'))
        valid=(table.se_mlmc>0)&(table.se_direct>0);ratios=(table.loc[valid,'se_direct']/table.loc[valid,'se_mlmc'])**2
        cm=costs[(costs.trial==trial)&(costs.role==2)].seconds.sum();cd=costs[(costs.trial==trial)&(costs.role==3)].seconds.sum()
        centered=table[table.contract.isin(central)]
        aggregate_gain=float(np.mean((centered.se_direct/centered.scale_direct)**2)/np.mean((centered.se_mlmc/centered.scale_mlmc)**2))
        z=(table.price_mlmc-table.price_direct)/np.sqrt(table.se_mlmc**2+table.se_direct**2)
        result['trials'].append(dict(trial=trial,mlmc_seconds=float(cm),direct_seconds=float(cd),
            realized_cost_ratio_mlmc_over_direct=float(cm/cd),median_sampling_variance_gain=float(ratios.median()),
            improved=int((ratios>1).sum()),valid_contracts=int(valid.sum()),central_mean_scaled_variance_gain=aggregate_gain,
            central_cost_adjusted_efficiency_gain=float(aggregate_gain*cd/cm),
            price_discrepancies_over_two_se=int((z.abs()>2).sum()),max_abs_price_difference_z=float(z.abs().max())))
    for method,g in pooled.groupby('method'):
        result['sampling'].append(dict(method=method,nominal_se_screens=int((g.se<=g.scale/8).sum()),
            median_se_over_scale=float((g.se/g.scale).median()),negative_estimates=int((g.price<0).sum()),
            note='Nominal screen only; all independent levels and tail-hit counts must be inspected.'))
    for trial,g in levels[(levels.level==3)&(levels.estimator=='pooled')].groupby('trial'):
        result['finest_corrections'].append(dict(trial=int(trial),median_abs_change_over_scale=float((g['mean'].abs()/g.scale).median()),
            above_scale_eighth=int((g['mean'].abs()>g.scale/8).sum()),above_two_se=int((g['mean'].abs()>2*g.se).sum()),
            max_abs_change_over_scale=float((g['mean'].abs()/g.scale).max())))
    for (method,estimator),g in e[e.estimator.isin(['halfA','halfB'])].groupby(['method','estimator']):
        t=g.merge(pooled[pooled.method==method],on=['trial','contract'],suffixes=('_half','_pooled'))
        valid=t.se_pooled>0;ratio=t.loc[valid,'se_half']/t.loc[valid,'se_pooled']
        result['pilot_stability'].append(dict(method=method,fit=estimator,median_se_ratio=float(ratio.median()),maximum_se_ratio=float(ratio.max()),negative_estimates=int((g.price<0).sum())))
    # Signed correction means can be negative; only final negative option estimates are flagged.
    result['zero_hit_level_rows']=levels[(levels.estimator=='pooled')&((levels.fine_hits==0)|(levels.coarse_hits==0))].to_dict(orient='records')
    z=(d['mean']-1)/d.se
    result['expectation_warnings']=d[z.abs()>2].assign(z=z[z.abs()>2]).to_dict(orient='records')
    result['design_seconds']=float(costs[costs.role<2].seconds.sum());result['setup_seconds']=record['setup_seconds'];result['wall_seconds']=record['wall_seconds']
    quotes=pd.read_csv(output/'contracts.csv')
    with np.load(output/'design_draws.npz') as design:
        x=draw_estimators(design['alloc_l3_fine_y'],design['alloc_l3_fine_c'],np.asarray(allocation['direct_coefficients']),quotes.forward.to_numpy())[:,:,1]
        direct_v=float((x[:,central]/quotes.price_scale_U.to_numpy()[central]).var(axis=0,ddof=1).mean()/allocation['direct_paths'])
    multilevel_v=float(np.sum(np.asarray(allocation['variances'])/np.asarray(allocation['counts'])))
    result['design_variance_prediction']=dict(direct=direct_v,mlmc=multilevel_v,direct_over_mlmc=direct_v/multilevel_v,
        note='Descriptive diagnosis; MLMC allocation was frozen without choosing between estimators.')
    result['correction_control_gains']=[]
    for level,g in levels[levels.contract.isin(central)].groupby('level'):
        raw=g[g.estimator=='raw'];controlled=g[g.estimator=='pooled']
        ratio=float((raw.variance/raw.scale**2).mean()/(controlled.variance/controlled.scale**2).mean())
        result['correction_control_gains'].append(dict(level=int(level),raw_over_controlled_variance=ratio,
            note='Aggregate observed central-contract path variance; fixed checking paths and coefficients.'))
    result['limitations']=['Cost match is pilot-predicted; actual times and variance-cost products govern interpretation.',
        'Empirical ratio estimates are conditional on frozen allocation/control fits and have sampling noise.',
        '512-step price remains discretized and GH48 approximate; nominal SE screens do not certify total error.',
        'No checking-based reallocation or estimator selection is applied.']
    (output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    audit(parser.parse_args().output)
