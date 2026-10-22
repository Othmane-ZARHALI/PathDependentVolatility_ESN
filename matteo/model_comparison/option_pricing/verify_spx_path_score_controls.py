"""Replay frozen Gaussian score coefficients, estimates and descriptive gains."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_path_score_controls import ROOT, dictionary, ridge_fit


def verify(parent,output):
    frozen=json.loads((output/'frozen.json').read_text())
    quotes=pd.read_csv(parent/'contracts.csv')
    for relative,digest in frozen['input_and_source_sha256'].items():
        assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==digest
    with np.load(parent/'design_draws.npz') as design:
        for name,fits in frozen['coefficients'].items():
            prefix=f'{name}_r0_t0';y,c,w,z=[design[prefix+'_'+key] for key in ('y','c','w','z')]
            for k,indices in enumerate((np.arange(len(w)),np.arange(len(w)//2),np.arange(len(w)//2,len(w)))):
                for j in range(y.shape[1]):
                    beta=ridge_fit(w[indices]*y[indices,j],dictionary(c[indices,j],w[indices],z[indices]))
                    np.testing.assert_allclose(beta,fits[k][j],rtol=1e-11,atol=1e-11)
    estimates=pd.read_csv(output/'estimates.csv')
    with np.load(parent/'checking_draws.npz') as checking, np.load(output/'checking_draws.npz') as adjusted:
        for name,fits in frozen['coefficients'].items():
            fits=np.asarray(fits)
            for trial in estimates.trial.unique():
                prefix=f'{name}_r2_t{trial}';w=checking[prefix+'_w'];z=checking[prefix+'_z']
                for level in range(3):
                    y,c=checking[prefix+f'_g{level}_y'],checking[prefix+f'_g{level}_c']
                    for j,row in enumerate(quotes.itertuples()):
                        for k in range(3):
                            values=(w*y[:,j]-dictionary(c[:,j],w,z)@fits[k,j])*row.forward
                            np.testing.assert_allclose(values,adjusted[prefix+f'_g{level}'][:,j,k],rtol=1e-11,atol=1e-10)
        for row in estimates.itertuples():
            level={64:0,128:1,256:2}[row.steps];k=['pooled','halfA','halfB'].index(row.fit)
            values=adjusted[f'{row.method}_r2_t{row.trial}_g{level}'][:,row.contract,k]
            np.testing.assert_allclose([values.mean(),values.std(ddof=1)/np.sqrt(len(values))],[row.price,row.se],rtol=1e-11,atol=1e-10)
    baseline=pd.read_csv(parent/'estimates.csv');baseline=baseline[(baseline.steps==256)&(baseline.estimator=='pooled')]
    summaries=[]
    for (name,fit),group in estimates[estimates.steps==256].groupby(['method','fit']):
        merged=group.merge(baseline[baseline.method==name],on=['trial','method','contract'],suffixes=('_new','_old'))
        ratio=(merged.se_old/merged.se_new)**2
        summaries.append(dict(method=name,fit=fit,median_variance_gain_vs_original_pooled=float(ratio.median()),
            improved=int((ratio>1).sum()),contract_trials=len(merged),negative_prices=int((group.price<0).sum()),
            median_se_over_scale=float((group.se/group.scale).median()),
            positive_hit_screens=int(((group.se<=group.scale/8)&(group.positive_payoffs>0)).sum())))
    result=dict(selected=frozen['selected'],comparisons=summaries,
        limitations=['Half fits and pooled fits share checking paths.','Variance ratios are descriptive, not certified population gains.',
                     'Ridge plus dictionary differs from original OLS fit; effect is combined.'])
    (output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f'Validated coefficients and {len(estimates)} score-control estimates; primary proposal remains {frozen["selected"]}.')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('parent',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();verify(args.parent,args.output)
