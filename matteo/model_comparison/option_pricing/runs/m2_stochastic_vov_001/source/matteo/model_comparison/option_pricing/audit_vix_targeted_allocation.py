"""Descriptive fresh-proposal comparisons, with scramble uncertainty kept explicit."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def run(folder):
    meta=json.loads((folder/'run.json').read_text());table=pd.read_csv(folder/'summary.csv')
    changes=pd.read_csv(folder/'prefix_changes.csv');fine=max(meta['arguments']['powers'])
    rows=[];option_rows=[]
    for days in meta['scope']['vix_fit_expiry_days']:
        block=table[(table.days==days)&(table.power==fine)]
        a=block[block.method=='random_targeted'].reset_index(drop=True)
        b=block[block.method=='fixed_targeted'].reset_index(drop=True)
        assert len(a)==len(b)
        for i in range(len(a)):
            x,y=a.iloc[i],b.iloc[i]
            assert x.metric==y.metric
            delta=y.estimate-x.estimate
            independent_se=float(np.hypot(x.replication_se,y.replication_se))
            if x.metric=='option':
                assert x.strike==y.strike and x.option_type==y.option_type
                option_rows.append(dict(days=days,strike=x.strike,option_type=x.option_type,scale=x.scale,
                    random_targeted=x.estimate,fixed_targeted=y.estimate,fixed_minus_random=delta,
                    independent_difference_se=independent_se,
                    random_se_over_scale=x.se_over_scale,fixed_se_over_scale=y.se_over_scale,
                    se_ratio=y.replication_se/x.replication_se if x.replication_se>0 else None))
            else:
                rows.append(dict(days=days,metric=x.metric,random_targeted=x.estimate,random_se=x.replication_se,
                    fixed_targeted=y.estimate,fixed_se=y.replication_se,fixed_minus_random=delta,
                    independent_difference_se=independent_se,
                    model_moment_target=y.model_moment_target if x.metric=='second_moment' else None,
                    random_moment_residual=x.estimate-x.model_moment_target if x.metric=='second_moment' else None,
                    fixed_moment_residual=y.estimate-y.model_moment_target if x.metric=='second_moment' else None))
    option=pd.DataFrame(option_rows);by_method=[]
    for days in meta['scope']['vix_fit_expiry_days']:
        for method in ('random_targeted','fixed_targeted'):
            x=table[(table.days==days)&(table.power==fine)&(table.method==method)&(table.metric=='option')]
            c=changes[(changes.days==days)&(changes.method==method)&(changes.metric=='option')]
            ratio=c.change.abs()/c.scale
            by_method.append(dict(days=days,method=method,quotes=len(x),positive_estimates=int((x.estimate>0).sum()),
                se_screen_passes=int(((x.se_over_scale<=.125)&(x.positive_payoff_count>=2)).sum()),
                median_se_over_scale=float(x.se_over_scale.median()),max_se_over_scale=float(x.se_over_scale.max()),
                median_max_pooled_share=float(x.max_pooled_share.median()),max_max_pooled_share=float(x.max_pooled_share.max()),
                zero_hit_quotes=int((x.positive_payoff_count==0).sum()),quotes_with_a_zero_hit_scramble=int((x.zero_hit_replications>0).sum()),
                median_abs_prefix_change_over_scale=float(ratio.median()),max_abs_prefix_change_over_scale=float(ratio.max())))
    mass=[]
    raw=pd.read_csv(folder/'outer_values.csv')
    for (days,method),x in raw.groupby(['days','method']):
        replicates=x.groupby('replicate').weight.mean()
        mass.append(dict(days=int(days),method=method,estimate=float(replicates.mean()),replication_se=float(replicates.std(ddof=1)/np.sqrt(len(replicates))),minimum_weight=float(x.weight.min()),maximum_weight=float(x.weight.max())))
    splits=[]
    reps=pd.read_csv(folder/'replications.csv');reps=reps[reps.power==fine]
    for keys,group in reps.groupby(['days','method','metric','strike','option_type'],dropna=False):
        x=group.sort_values('replicate').estimate.to_numpy();half=len(x)//2
        a,b=x[:half],x[half:];sa=a.std(ddof=1)/np.sqrt(len(a));sb=b.std(ddof=1)/np.sqrt(len(b))
        splits.append(dict(days=int(keys[0]),method=keys[1],metric=keys[2],strike=keys[3],option_type=keys[4],first_mean=float(a.mean()),second_mean=float(b.mean()),first_se=float(sa),second_se=float(sb),second_minus_first=float(b.mean()-a.mean()),independent_difference_se=float(np.hypot(sa,sb)),se_ratio_second_over_first=float(sb/sa) if sa>0 else None))
    pd.DataFrame(splits).to_csv(folder/'replication_splits.csv',index=False)
    summary=[]
    for days,x in option.groupby('days'):
        summary.append(dict(days=int(days),quotes=len(x),median_option_se_ratio=float(x.se_ratio.median()),
            quotes_with_lower_fixed_se=int((x.se_ratio<1).sum()),
            max_abs_price_difference_over_scale=float((x.fixed_minus_random.abs()/x.scale).max()),
            differences_exceed_two_combined_se=int((x.fixed_minus_random.abs()>2*x.independent_difference_se).sum())))
    record=dict(status='descriptive frozen proposal comparison; no accuracy acceptance',
        inputs={name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in ('run.json','summary.csv','prefix_changes.csv','outer_values.csv','replications.csv')},
        method_summaries=by_method,option_comparisons=summary,weighted_mass_diagnostics=mass,
        audit_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limitations=['Observed SE screens are necessary diagnostics, not price-error guarantees. Half-replication splits reuse checking data and are descriptive, not new validation.',
            'Two combined SE comparisons are descriptive and not simultaneous or calibrated tail tests.',
            'Finite scramble SE can miss rare contributions; prefix stability cannot rule them out.',
            'Identical frozen targeted density; allocation and net construction change, without new pilot fitting.'])
    pd.DataFrame(rows).to_csv(folder/'method_comparison.csv',index=False)
    option.to_csv(folder/'option_comparison.csv',index=False)
    (folder/'comparison.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path);run(parser.parse_args().folder)
