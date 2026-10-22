"""Observed reflected-pair correction variance gains with full measured costs."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def audit(output):
    record=json.loads((output/'run.json').read_text());e=pd.read_csv(output/'estimates.csv');c=pd.read_csv(output/'costs.csv')
    delta=pd.read_csv(output/'paired_differences.csv');d=pd.read_csv(output/'diagnostics.csv');central=[2,3,8,9]
    result=dict(comparisons=[],pilot_stability=[],zero_payoff_rows=[],kernel_endpoint_error=float(c.endpoint_error.max()),
                kernel_equity_error=float(c.equity_error.max()),setup_seconds=record['setup_seconds'],wall_seconds=record['wall_seconds'])
    for (trial,level,estimator),g in e[e.estimator.isin(['raw','pooled'])].groupby(['trial','level','estimator']):
        table=g[g.method=='plain'].merge(g[g.method=='bridge'],on=['trial','level','contract'],suffixes=('_plain','_bridge'))
        valid=(table.variance_plain>0)&(table.variance_bridge>0)
        ratio=table.loc[valid,'variance_plain']/table.loc[valid,'variance_bridge']
        group=table[table.contract.isin(central)]
        var_plain=float((group.variance_plain/group.scale_plain**2).mean());var_bridge=float((group.variance_bridge/group.scale_bridge**2).mean())
        timing=c[(c.role==1)&(c.trial==trial)&(c.level==level)]
        ordinary=float(timing.ordinary_seconds.sum());bridge=float(timing.bridge_seconds.sum())
        result['comparisons'].append(dict(trial=int(trial),level=int(level),estimator=estimator,
            improved_contracts=int((ratio>1).sum()),valid_contracts=int(valid.sum()),median_variance_gain=float(ratio.median()),
            central_variance_gain=var_plain/var_bridge,ordinary_seconds=ordinary,bridge_seconds=bridge,
            cost_ratio_bridge_over_ordinary=bridge/ordinary,central_variance_cost_gain=var_plain/var_bridge*ordinary/bridge))
    for (level,method,fit),g in e[e.estimator.isin(['halfA','halfB'])].groupby(['level','method','estimator']):
        base=e[(e.level==level)&(e.method==method)&(e.estimator=='pooled')]
        paired=g.merge(base,on=['trial','contract'],suffixes=('_half','_pooled'));valid=paired.se_pooled>0
        ratio=paired.loc[valid,'se_half']/paired.loc[valid,'se_pooled']
        result['pilot_stability'].append(dict(level=int(level),method=method,fit=fit,median_se_ratio=float(ratio.median()),max_se_ratio=float(ratio.max())))
    for (level,method),g in e[e.estimator=='pooled'].groupby(['level','method']):
        result.setdefault('correction_sizes',[]).append(dict(level=int(level),method=method,
            median_abs_correction_over_scale=float((g.correction.abs()/g.scale).median()),
            changes_above_scale_eighth=int((g.correction.abs()>g.scale/8).sum()),above_two_se=int((g.correction.abs()>2*g.se).sum())))
    zero=e[(e.estimator=='pooled')&(e.method=='bridge')&(e.fine_hits==0)&(e.reflected_hits==0)&(e.coarse_hits==0)]
    result['zero_payoff_rows']=zero.to_dict(orient='records')
    check=delta[(delta.estimator=='pooled')&(delta.se>0)].copy();check['z']=check.difference/check.se
    result['paired_mean_warnings']=check[check.z.abs()>2].to_dict(orient='records')
    z=(d['mean']-1)/d.se;result['expectation_warnings']=d[z.abs()>2].assign(z=z[z.abs()>2]).to_dict(orient='records')
    # Explanatory replay with the SAME frozen plain coefficients, not refitting.
    frozen=json.loads((output/'frozen.json').read_text());quotes=pd.read_csv(output/'contracts.csv')
    result['same_coefficient_decomposition']=[]
    with np.load(output/'checking_draws.npz') as draws:
        for level in (2,3):
            beta=np.asarray(frozen['coefficients'][str(level)]['plain'])[0]
            for trial in range(record['arguments']['trials']):
                prefix=f't{trial}_l{level}'
                def values(method):
                    y,c=draws[prefix+'_'+method+'_y'],draws[prefix+'_'+method+'_c']
                    return ((y-np.einsum('njk,jk->nj',c,beta))*quotes.forward.to_numpy()/quotes.price_scale_U.to_numpy())[:,central]
                plain=values('plain');average=values('bridge');difference=plain-average
                vp=float(plain.var(axis=0,ddof=1).mean());va=float(average.var(axis=0,ddof=1).mean());vd=float(difference.var(axis=0,ddof=1).mean())
                covariance=float(((average-average.mean(axis=0))*(difference-difference.mean(axis=0))).sum(axis=0).mean()/(len(average)-1))
                np.testing.assert_allclose(vp,va+vd+2*covariance,rtol=1e-12)
                result['same_coefficient_decomposition'].append(dict(level=level,trial=trial,averaged_variance_fraction=va/vp,
                    difference_variance_fraction=vd/vp,twice_covariance_fraction=2*covariance/vp,
                    note='Post-hoc explanation only; no checking-based coefficient change or selection.'))
    result['limitations']=['Methods share original paths; pairs are iid, the two fine paths within a pair are dependent.',
        'Variance-cost gains use actual costs at equal pair counts; they are descriptive, not a matched-runtime allocation experiment.',
        'No checking-based coupling/fit selection; corrections only, no complete pricing clearance.']
    (output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    audit(parser.parse_args().output)
