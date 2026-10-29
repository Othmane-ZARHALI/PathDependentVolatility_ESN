"""Replay positive conditional prices, raw payoffs, statistics and covariance gates."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import ndtr
from .m2_spx import frozen_stage3
from esn_eq.conditional_spx import lognormal_payoffs,partial_transition
from esn_eq.pricing import DatedM2Pricer


def audit(output):
    r=json.loads((output/'run.json').read_text());q=pd.read_csv(output/'contracts.csv');e=pd.read_csv(output/'estimates.csv');paired=pd.read_csv(output/'paired_differences.csv')
    for path,digest in r['implementation_sha256'].items():assert hashlib.sha256((output/'source'/path).read_bytes()).hexdigest()==digest
    seeds=[seed for row in r['seed_records'] for seed in row['seeds']];assert len(seeds)==len(set(seeds))
    architecture,parameters,premium,_=frozen_stage3()
    for name,value in [('architecture',architecture),('parameters',parameters),('premium',premium)]:assert r[name]==json.loads(json.dumps(asdict(value)))
    p=DatedM2Pricer(architecture,premium);transition=json.loads((output/'transition_audit.json').read_text())
    for row in transition:
        _,_,actual=partial_transition(p,row['dt'])
        for key,value in actual.items():np.testing.assert_allclose(value,row[key],rtol=1e-12,atol=1e-20)
    monitors=[]
    with np.load(output/'draws.npz') as data:
        for trial in range(r['arguments']['trials']):
            raw=data[f't{trial}_raw'];conditional=data[f't{trial}_conditional'];stock=data[f't{trial}_raw_stock']
            mean=data[f't{trial}_log_mean'];variance=data[f't{trial}_log_variance'];forward=data[f't{trial}_conditional_stock']
            assert np.all(raw>=0)&np.all(conditional>=0)&np.all(variance>=0)
            np.testing.assert_allclose(forward,np.exp(mean+.5*variance),rtol=1e-14,atol=1e-14)
            np.testing.assert_allclose(raw,np.maximum(np.where(q.option_type=='call',1.,-1.)[None,:]*(stock-q.K_over_F.to_numpy()[None,:]),0.),rtol=1e-14,atol=1e-14)
            np.testing.assert_allclose(conditional,lognormal_payoffs(mean,variance,q.K_over_F.to_numpy()[None,:],q.option_type.to_numpy()[None,:]),rtol=1e-12,atol=1e-14)
            # Independent direct-CDF Black expression (absolute tolerance covers cancellation in tiny tails).
            sd=np.sqrt(variance);d2=(mean-np.log(q.K_over_F.to_numpy()[None,:]))/sd;d1=d2+sd
            plain=np.where(q.option_type.to_numpy()[None,:]=='call',forward*ndtr(d1)-q.K_over_F.to_numpy()[None,:]*ndtr(d2),
                q.K_over_F.to_numpy()[None,:]*ndtr(-d2)-forward*ndtr(-d1))
            np.testing.assert_allclose(conditional,plain,rtol=1e-10,atol=1e-12)
            for row in e[e.trial==trial].itertuples():
                values=(raw if row.method=='raw' else conditional)[:,row.contract]*q.forward.iloc[row.contract]
                np.testing.assert_allclose([values.mean(),values.std(ddof=1)/np.sqrt(len(values)),values.var(ddof=1)],
                    [row.estimate,row.se,row.variance],rtol=1e-11,atol=1e-10)
            for row in paired[paired.trial==trial].itertuples():
                values=(conditional[:,row.contract]-raw[:,row.contract])*q.forward.iloc[row.contract]
                np.testing.assert_allclose([values.mean(),values.std(ddof=1)/np.sqrt(len(values))],[row.difference,row.se],rtol=1e-11,atol=1e-10)
            for days,g in q.groupby('expiry_days'):
                j=int(g.index[0])
                for name,values,target in [('raw_stock',stock[:,j],1.),('conditional_stock',forward[:,j],1.),('paired_stock',forward[:,j]-stock[:,j],0.)]:
                    m=float(values.mean());se=float(values.std(ddof=1)/np.sqrt(len(values)))
                    monitors.append(dict(trial=trial,expiry_days=int(days),statistic=name,mean=m,se=se,expectation=target,z=(m-target)/se if se else None))
    combined=e[e.method=='raw'].merge(e[e.method=='conditional'],on=['trial','contract'],suffixes=('_raw','_conditional'))
    combined['variance_gain']=combined.variance_raw/combined.variance_conditional.replace(0,np.nan)
    combined.to_csv(output/'comparisons.csv',index=False);pd.DataFrame(monitors).to_csv(output/'monitors.csv',index=False)
    summary=[]
    for (days,method),g in e.groupby(['expiry_days','method']):
        summary.append(dict(expiry_days=int(days),method=method,negative=int((g.estimate<0).sum()),zero_raw_hit_cases=int((g.raw_positive_hits==0).sum()),
            conservative_sampling_passes=int(((g.raw_positive_hits>0)&(g.se<=g.scale/8)&(g.estimate>=0)).sum()),
            median_se_over_scale=float((g.se/g.scale).median()),min_raw_hits=int(g.raw_positive_hits.min())))
    paired['z']=paired.difference/paired.se.replace(0,np.nan);monitors=pd.DataFrame(monitors)
    result=dict(auditor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        verification='Snapshot hashes, model, transition audit, raw and positive conditional payoffs, independent Black formula and all means/SEs/variances replay.',
        estimates=len(e),paired=len(paired),independent_equity_variance_fraction=r['independent_equity_variance_fraction'],summary=summary,
        variance_improved=int((combined.variance_gain>1).sum()),variance_comparisons=len(combined),median_variance_gain=float(combined.variance_gain.median()),
        paired_mean_warnings=paired[paired.z.abs()>2].to_dict(orient='records'),stock_mean_warnings=monitors[monitors.z.abs()>2].to_dict(orient='records'),
        max_standardized_covariance_error=max(row['max_standardized_covariance_error'] for row in transition),
        max_absolute_covariance_error=max(row['max_abs_legacy_covariance_error'] for row in transition),
        limitations=['Partial conditional vs raw target-Q only; no comparison with controlled importance-sampling baseline or separate standalone cost.',
            'Sample variance gains can fluctuate around the population conditional-variance inequality.',
            'Zero raw payoff hits cannot establish tail accuracy even with positive/tiny conditional payoffs.',
            'Finite precision covariance audits are not formal nonlinear payoff-error bounds; time integration/GH refinement remain pending.'])
    (output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path);audit(parser.parse_args().output)
