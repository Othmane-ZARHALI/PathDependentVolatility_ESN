"""Descriptive paired-control results; no accuracy or superiority certificate."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def run(folder):
    r=json.loads((folder/'run.json').read_text());s=pd.read_csv(folder/'summary.csv');fine=max(r['arguments']['powers']);x=s[s.power==fine];rows=[]
    for (days,method),block in x.groupby(['days','method']):
        q=block[block.metric=='option'];f=block[block.metric=='forward'].iloc[0];m=block[block.metric=='second_moment'].iloc[0]
        rows.append(dict(days=int(days),method=method,option_quotes=len(q),positive_option_estimates=int((q.estimate>0).sum()),negative_option_estimates=int((q.estimate<0).sum()),option_se_screens=int(((q.estimate>0)&(q.positive_payoff_count>=2)&(q.se_over_scale<=.125)).sum()),median_option_se_ratio_to_raw=float(q.se_ratio_to_raw.median()),max_option_se_over_scale=float(q.se_over_scale.max()),forward_estimate=float(f.estimate),forward_se=float(f.replication_se),forward_se_ratio_to_raw=float(f.se_ratio_to_raw),forward_se_screen=bool(f.estimate>0 and f.se_over_scale<=.125),forward_change_from_raw=float(f.change_from_raw),forward_paired_change_se=float(f.paired_change_se),second_moment=float(m.estimate),moment_se=float(m.replication_se),moment_target=float(m.model_moment_target),moment_residual=float(m.estimate-m.model_moment_target)))
    table=pd.DataFrame(rows);table.to_csv(folder/'comparison.csv',index=False)
    warnings=table[table.moment_residual.abs()>2*table.moment_se].copy()
    warnings['residual_over_reported_se']=warnings.moment_residual/warnings.moment_se
    warnings.to_csv(folder/'moment_warnings.csv',index=False)
    changes=pd.read_csv(folder/'prefix_changes.csv');prefix=[]
    for (days,method),q in changes[changes.metric=='option'].groupby(['days','method']):prefix.append(dict(days=int(days),method=method,max_abs_change_over_scale=float((q.change.abs()/q.scale).max()),median_abs_change_over_scale=float((q.change.abs()/q.scale).median())))
    pd.DataFrame(prefix).to_csv(folder/'prefix_summary.csv',index=False)
    # Fits share these checking draws; sensitivity is not new validation.
    sensitivity=x[x.method.isin(['gaussian','gaussian_batch_a','gaussian_batch_b','legacy12'])].copy()
    sensitivity.to_csv(folder/'pilot_fit_sensitivity.csv',index=False)
    stability=[]
    for days in r['scope']['vix_fit_expiry_days']:
        block=x[x.days==days];a=block[block.method=='gaussian'];b=block[block.method=='legacy12']
        fa=a[a.metric=='forward'].iloc[0];fb=b[b.metric=='forward'].iloc[0]
        qa=a[a.metric=='option'].reset_index(drop=True);qb=b[b.metric=='option'].reset_index(drop=True)
        batches=block[(block.metric=='forward')&block.method.isin(['gaussian_batch_a','gaussian_batch_b'])]
        stability.append(dict(days=days,primary_forward_se=float(fa.replication_se),legacy_forward_se=float(fb.replication_se),primary_to_legacy_forward_se_ratio=float(fa.replication_se/fb.replication_se),median_option_se_ratio_primary_to_legacy=float((qa.replication_se/qb.replication_se).median()),batch_forward_se_min=float(batches.replication_se.min()),batch_forward_se_max=float(batches.replication_se.max()),batch_forward_estimate_min=float(batches.estimate.min()),batch_forward_estimate_max=float(batches.estimate.max())))
    pd.DataFrame(stability).to_csv(folder/'stability_summary.csv',index=False)
    record=dict(status='descriptive fresh paired control comparison; no convergence acceptance',inputs={n:hashlib.sha256((folder/n).read_bytes()).hexdigest() for n in ('run.json','summary.csv','prefix_changes.csv','replications.csv','coefficients.json')},primary_method=r['primary_method'],audit_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limitations=['SE screens exclude negative controlled estimates and retain original raw hit counts; they do not certify tail or quadrature error.','Mass-only, balanced pilot and legacy alternatives were frozen before checking, with pooled R28 primary; no favorable alternative selected afterward.','Control-variant differences and count changes are paired; repeated-pilot efficiency and deterministic bias are separate from conditional checking SEs. Moment-warning flags across correlated fits are descriptive, not calibrated multiple tests.'])
    (folder/'comparison.json').write_text(json.dumps(record,indent=2)+'\n');print(table.to_string(index=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path);run(parser.parse_args().folder)
