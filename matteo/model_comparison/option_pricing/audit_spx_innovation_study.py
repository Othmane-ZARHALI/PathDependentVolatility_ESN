"""Observed control-dictionary gains, pilot stability and unrepaired warnings."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_innovation_study import datasets


def audit(output):
    r=json.loads((output/'run.json').read_text());e=pd.read_csv(output/'estimates.csv');c=pd.read_csv(output/'costs.csv');d=pd.read_csv(output/'diagnostics.csv')
    result=dict(comparisons=[],pilot_stability=[],price_screens=[],correction_sizes=[],paired_mean_warnings=[],expectation_warnings=[],
        setup_seconds=r['setup_seconds'],wall_seconds=r['wall_seconds'])
    for (trial,target),g in e.groupby(['trial','target']):
        m=g[g.method=='base'].merge(g[g.method=='augmented'],on=['trial','target','contract'],suffixes=('_base','_aug'))
        valid=(m.variance_base>0)&(m.variance_aug>0);ratio=m.loc[valid,'variance_base']/m.loc[valid,'variance_aug']
        a=m[m.contract.isin([2,3,8,9])];gain=float((a.variance_base/a.scale_base**2).mean()/(a.variance_aug/a.scale_aug**2).mean())
        times=c[(c.role==1)&(c.trial==trial)];ordinary=float(times.ordinary_seconds.sum());augmented=float(times.augmented_seconds.sum())
        result['comparisons'].append(dict(trial=int(trial),target=target,improved=int((ratio>1).sum()),valid_contracts=int(valid.sum()),
            median_variance_gain=float(ratio.median()),central_variance_gain=gain,ordinary_seconds=ordinary,augmented_seconds=augmented,
            joint_pipeline_cost_ratio=augmented/ordinary,central_joint_pipeline_efficiency_gain=gain*ordinary/augmented))
    for (target,method),g in e[e.method.isin(['base_A','base_B','aug_A','aug_B'])].groupby(['target','method']):
        primary='augmented' if method.startswith('aug_') else 'base';base=e[(e.target==target)&(e.method==primary)]
        m=g.merge(base,on=['trial','contract'],suffixes=('_half','_pooled'));valid=m.se_pooled>0;ratio=m.loc[valid,'se_half']/m.loc[valid,'se_pooled']
        result['pilot_stability'].append(dict(target=target,method=method,median_se_ratio=float(ratio.median()),max_se_ratio=float(ratio.max())))
    for (target,method),g in e[e.method.isin(['legacy','base','augmented'])].groupby(['target','method']):
        if target.startswith('price'):
            result['price_screens'].append(dict(target=target,method=method,positive_hit_screens=int(((g.nonzero_raw>0)&(g.se<=g.scale/8)).sum()),
                zero_raw_hits=int((g.nonzero_raw==0).sum()),negative_adjusted_prices=int((g.estimate<0).sum()),median_se_over_scale=float((g.se/g.scale).median())))
        else:
            result['correction_sizes'].append(dict(target=target,method=method,median_abs_change_over_scale=float((g.estimate.abs()/g.scale).median()),
                above_scale_eighth=int((g.estimate.abs()>g.scale/8).sum()),above_two_se=int((g.estimate.abs()>2*g.se).sum())))
    paired=pd.read_csv(output/'paired_differences.csv');valid=paired.se>0;z=paired.loc[valid,'difference']/paired.loc[valid,'se']
    result['paired_mean_warnings']=paired.loc[valid].assign(z=z).loc[z.abs()>2].to_dict(orient='records')
    z=(d['mean']-d.expectation)/d.se;result['expectation_warnings']=d[z.abs()>2].assign(z=z[z.abs()>2]).to_dict(orient='records')
    # Contract-weighted gains can have poorly sampled tails even when bare innovation sums pass.
    controls=[];labels=['stock','reference_payoff','deterministic_delta','spot_delta','mass',
        'ou_linear_C','ou_linear_Y','ou_linear_F','ou_square_C','ou_square_Y','ou_square_F','ou_mixed_C','ou_mixed_Y','ou_mixed_F']
    with np.load(output/'checking_draws.npz') as checking:
        for trial in range(r['arguments']['trials']):
            prefix=f't{trial}_';data={key[len(prefix):]:checking[key] for key in checking.files if key.startswith(prefix) and not key.endswith('_adjusted')}
            for target,(y,base,vol) in datasets(data).items():
                values=np.concatenate([base,vol],axis=2);means=values.mean(axis=0);ses=values.std(axis=0,ddof=1)/np.sqrt(len(values))
                for contract in range(values.shape[1]):
                    for k,label in enumerate(labels):
                        z=float(means[contract,k]/ses[contract,k]) if ses[contract,k]>0 else None
                        controls.append(dict(trial=trial,target=target,contract=contract,control=label,mean=float(means[contract,k]),se=float(ses[contract,k]),z=z))
    pd.DataFrame(controls).to_csv(output/'control_means.csv',index=False)
    result['contract_control_mean_warnings']=[row for row in controls if row['z'] is not None and abs(row['z'])>2]
    result['contract_control_mean_warning_count']=len(result['contract_control_mean_warnings'])
    result['limitations']=['Primary compares ridge14 to ridge5 at the same regularization; legacy OLS5 retained separately.',
        'Costs apply to the shared 128/256/512-grid pipeline, not standalone price-only simulation.',
        'Sample ratios and multiple correlated flags are descriptive; no checking-based selection or clipping.',
        'Variance reduction alone does not remove stock-integration or normalization bias.']
    (output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path);audit(parser.parse_args().output)
