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
        a=block[block.method=='learned_mean'].reset_index(drop=True)
        b=block[block.method=='score_inflated'].reset_index(drop=True)
        assert len(a)==len(b)
        for i in range(len(a)):
            x,y=a.iloc[i],b.iloc[i]
            assert x.metric==y.metric
            delta=y.estimate-x.estimate
            independent_se=float(np.hypot(x.replication_se,y.replication_se))
            if x.metric=='option':
                assert x.strike==y.strike and x.option_type==y.option_type
                option_rows.append(dict(days=days,strike=x.strike,option_type=x.option_type,scale=x.scale,
                    learned_mean=x.estimate,score_inflated=y.estimate,score_minus_learned=delta,
                    independent_difference_se=independent_se,
                    learned_se_over_scale=x.se_over_scale,score_se_over_scale=y.se_over_scale,
                    se_ratio=y.replication_se/x.replication_se if x.replication_se>0 else None))
            else:
                rows.append(dict(days=days,metric=x.metric,learned_mean=x.estimate,learned_se=x.replication_se,
                    score_inflated=y.estimate,score_se=y.replication_se,score_minus_learned=delta,
                    independent_difference_se=independent_se,
                    model_moment_target=y.model_moment_target if x.metric=='second_moment' else None,
                    learned_moment_residual=x.estimate-x.model_moment_target if x.metric=='second_moment' else None,
                    score_moment_residual=y.estimate-y.model_moment_target if x.metric=='second_moment' else None))
    option=pd.DataFrame(option_rows);by_method=[]
    for days in meta['scope']['vix_fit_expiry_days']:
        for method in ('learned_mean','score_inflated'):
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
    summary=[]
    for days,x in option.groupby('days'):
        summary.append(dict(days=int(days),quotes=len(x),median_option_se_ratio=float(x.se_ratio.median()),
            quotes_with_lower_score_se=int((x.se_ratio<1).sum()),
            max_abs_price_difference_over_scale=float((x.score_minus_learned.abs()/x.scale).max()),
            differences_exceed_two_combined_se=int((x.score_minus_learned.abs()>2*x.independent_difference_se).sum())))
    record=dict(status='descriptive frozen proposal comparison; no accuracy acceptance',
        inputs={name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in ('run.json','summary.csv','prefix_changes.csv','outer_values.csv')},
        method_summaries=by_method,option_comparisons=summary,weighted_mass_diagnostics=mass,
        audit_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limitations=['Observed SE screens are necessary diagnostics, not price-error guarantees.',
            'Two combined SE comparisons are descriptive and not simultaneous or calibrated tail tests.',
            'Finite scramble SE can miss rare contributions; prefix stability cannot rule them out.',
            'One frozen prior pilot and one proposed rank/scale design; no full pilot-variation study.'])
    pd.DataFrame(rows).to_csv(folder/'method_comparison.csv',index=False)
    option.to_csv(folder/'option_comparison.csv',index=False)
    (folder/'comparison.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path);run(parser.parse_args().folder)
