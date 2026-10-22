"""Replay proposal likelihoods, design selection, pilot fits and fresh SPX checks."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_volatility_sampling import fit_variants, adjusted
from esn_eq.importance_sampling import DefensiveGaussianMixture


def verify(output):
    record=json.loads((output/'run.json').read_text())
    selection=json.loads((output/'selection.json').read_text())
    quotes=pd.read_csv(output/'contracts.csv')
    for relative,digest in record['implementation_sha256'].items():
        assert hashlib.sha256((output/'source'/relative).read_bytes()).hexdigest()==digest
    seeds=[s for row in record['seed_records'] for s in row['seeds']]
    assert len(seeds)==len(set(seeds)), 'Pilot/selection/check seed collision'
    with np.load(output/'design_draws.npz') as design, np.load(output/'checking_draws.npz') as checking:
        scores={}
        for name,definition in record['proposals'].items():
            mixture=DefensiveGaussianMixture(definition['shifts'],definition['probabilities'])
            for role,trial,archive in [(0,0,design),(1,0,design)]+[(2,t,checking) for t in range(record['arguments']['trials'])]:
                prefix=f'{name}_r{role}_t{trial}'
                w,z=archive[prefix+'_w'],archive[prefix+'_z']
                np.testing.assert_allclose(mixture.weights(z),w,rtol=1e-12,atol=1e-14)
                assert np.all((w>0)&(w<=2))
                assert np.all((archive[prefix+'_components']>=0)&(archive[prefix+'_components']<len(mixture.shifts)))
            prefix=f'{name}_r0_t0'
            fits=fit_variants(design[prefix+'_y'],design[prefix+'_c'],design[prefix+'_w'])
            np.testing.assert_allclose(fits,record['fitted_controls'][name],rtol=1e-11,atol=1e-11)
            prefix=f'{name}_r1_t0'
            x=adjusted(design[prefix+'_y'],design[prefix+'_c'],design[prefix+'_w'],fits,quotes.forward.to_numpy())[:,:,1]
            variance=x.var(axis=0,ddof=1)
            scores[name]=float(np.log(np.maximum(variance[[2,3,8,9]]/quotes.price_scale_U.to_numpy()[[2,3,8,9]]**2,1e-24)).mean())
            saved=next(row for row in selection['candidates'] if row['method']==name)
            np.testing.assert_allclose(scores[name],saved['score'],rtol=1e-11,atol=1e-12)
            np.testing.assert_allclose(variance,saved['variance'],rtol=1e-11,atol=1e-12)
            for trial in range(record['arguments']['trials']):
                prefix=f'{name}_r2_t{trial}'; w=checking[prefix+'_w']
                for level in range(3):
                    y,c=checking[prefix+f'_g{level}_y'],checking[prefix+f'_g{level}_c']
                    x=adjusted(y,c,w,fits,quotes.forward.to_numpy())
                    np.testing.assert_allclose(x,checking[prefix+f'_g{level}_adjusted'],rtol=1e-11,atol=1e-10)
        assert min(scores,key=scores.get)==selection['selected']==record['selected']
        estimates=pd.read_csv(output/'estimates.csv')
        changes=pd.read_csv(output/'paired_changes.csv')
        labels=['raw','pooled','halfA','halfB']
        for row in estimates.itertuples():
            level={64:0,128:1,256:2}[row.steps]; prefix=f'{row.method}_r2_t{row.trial}'
            x=checking[prefix+f'_g{level}_adjusted'][:,row.contract,labels.index(row.estimator)]
            np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x))],[row.price,row.se],rtol=1e-11,atol=1e-10)
            assert row.positive_payoffs==int((checking[prefix+f'_g{level}_y'][:,row.contract]>0).sum())
        for row in changes.itertuples():
            low={64:0,128:1}[row.coarse_steps]; prefix=f'{row.method}_r2_t{row.trial}'
            delta=(checking[prefix+f'_g{low+1}_adjusted']-checking[prefix+f'_g{low}_adjusted'])[:,row.contract,labels.index(row.estimator)]
            np.testing.assert_allclose([delta.mean(),delta.std(ddof=1)/np.sqrt(len(delta)),delta.var(ddof=1)],
                                      [row.change,row.se,row.variance],rtol=1e-11,atol=1e-10)
        for row in pd.read_csv(output/'diagnostics.csv').itertuples():
            index=['mass','stock29','stock57','reference29','reference57'].index(row.statistic)
            x=checking[f'{row.method}_r2_t{row.trial}_monitors'][:,index]
            np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x))],[row.mean,row.se],rtol=1e-11,atol=1e-12)
    assert len(estimates)==record['arguments']['trials']*3*3*12*4
    assert len(changes)==record['arguments']['trials']*3*2*12*4
    print(f'Validated likelihoods, pilot fits, design selection, {len(estimates)} estimates and {len(changes)} paired corrections; numerical accuracy remains open.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    verify(parser.parse_args().output)
