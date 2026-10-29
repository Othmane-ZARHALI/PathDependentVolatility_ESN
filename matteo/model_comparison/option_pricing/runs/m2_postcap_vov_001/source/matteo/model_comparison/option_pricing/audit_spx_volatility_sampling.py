"""Retain design/fresh discrepancies and pilot instability without re-selection."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def gains(e, candidate, central=False):
    data=e[e.contract.isin([2,3,8,9])] if central else e
    table=data.pivot(index=['trial','contract'],columns='method',values='se')
    valid=(table['equity']>0)&(table[candidate]>0)
    ratio=(table.loc[valid,'equity']/table.loc[valid,candidate])**2
    return dict(cases=int(valid.sum()),improved=int((ratio>1+1e-12).sum()),median=float(ratio.median()),
                geometric_mean=float(np.exp(np.log(ratio).mean())),minimum=float(ratio.min()),maximum=float(ratio.max()))


def audit(output):
    record=json.loads((output/'run.json').read_text());selection=json.loads((output/'selection.json').read_text())
    estimates=pd.read_csv(output/'estimates.csv');changes=pd.read_csv(output/'paired_changes.csv')
    fine=estimates[(estimates.steps==256)&(estimates.estimator=='pooled')]
    result=dict(selected=record['selected'],design_objective=selection['objective'],candidates={},pilot_stability=[],refinements=[])
    base_score=next(row['score'] for row in selection['candidates'] if row['method']=='equity')
    for name,g in fine.groupby('method'):
        score=next(row['score'] for row in selection['candidates'] if row['method']==name)
        result['candidates'][name]=dict(design_geometric_variance_gain_vs_equity=float(np.exp(base_score-score)),
            fresh_all_contract_gains=gains(fine,name),fresh_central_contract_gains=gains(fine,name,True),
            positive_hit_sampling_screens=int(((g.positive_payoffs>0)&(g.se<=g.scale/8)).sum()),
            zero_hits=int((g.positive_payoffs==0).sum()),negative_prices=int((g.price<0).sum()),
            median_se_over_scale=float((g.se/g.scale).median()))
        for trial,gt in fine.groupby('trial'):
            result['candidates'][name][f'trial{trial}_gains']=gains(gt,name)
    for (name,estimator),g in estimates[(estimates.steps==256)&estimates.estimator.isin(['halfA','halfB'])].groupby(['method','estimator']):
        paired=g.merge(fine[['method','trial','contract','se']],on=['method','trial','contract'],suffixes=('_half','_pooled'))
        ratio=paired.se_half/paired.se_pooled
        result['pilot_stability'].append(dict(method=name,fit=estimator,median_se_ratio_to_pooled=float(ratio.median()),max_se_ratio_to_pooled=float(ratio.max()),
            negative_prices=int((g.price<0).sum()),positive_hit_screens=int(((g.positive_payoffs>0)&(g.se<=g.scale/8)).sum())))
    for (name,low),g in changes[changes.estimator=='pooled'].groupby(['method','coarse_steps']):
        f=fine[fine.method==name][['trial','contract','se']]
        m=g.merge(f,on=['trial','contract'],suffixes=('_correction','_fine'))
        ratio=(m.se_correction/m.se_fine)**2
        result['refinements'].append(dict(method=name,coarse_steps=int(low),median_abs_change_over_scale=float((g.change.abs()/g.scale).median()),
            above_scale_eighth=int((g.change.abs()>g.scale/8).sum()),above_two_se=int((g.change.abs()>2*g.se).sum()),
            median_correction_variance_over_fine=float(ratio.median())))
    diagnostics=pd.read_csv(output/'diagnostics.csv')
    result['monitors']=diagnostics.to_dict(orient='records')
    result['monitor_max_abs_z']=float(((diagnostics['mean']-1).abs()/diagnostics.se).max())
    result['costs']=pd.read_csv(output/'costs.csv').groupby(['method','role'])[['state_seconds','stock_control_seconds']].sum().reset_index().to_dict(orient='records')
    with np.load(output/'checking_draws.npz') as draws:
        result['paired_half_fit_differences']=[]
        for name in record['proposals']:
            for trial in range(record['arguments']['trials']):
                x=draws[f'{name}_r2_t{trial}_g2_adjusted'];diff=x[:,:,2]-x[:,:,3]
                scale=pd.read_csv(output/'contracts.csv').price_scale_U.to_numpy()
                se=diff.std(axis=0,ddof=1)/np.sqrt(len(diff))
                valid=se>0
                z=np.abs(diff.mean(axis=0)[valid])/se[valid]
                result['paired_half_fit_differences'].append(dict(method=name,trial=trial,zero_se_contracts=int((~valid).sum()),
                    max_abs_price_difference_over_scale=float(np.max(np.abs(diff.mean(axis=0))/scale)),
                    max_abs_z=float(z.max()) if len(z) else 0.))
    (output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('monitors','costs')},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    audit(parser.parse_args().output)
