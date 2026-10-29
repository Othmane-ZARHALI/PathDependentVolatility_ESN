"""Replay frozen allocation, level likelihoods, controls and MLMC uncertainty."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_mlmc import fit_controls,draw_estimators
from .spx_volatility_sampling import proposal
from esn_eq.mlmc import allocate_levels,independent_level_statistics


def verify(output):
    record=json.loads((output/'run.json').read_text());allocation=json.loads((output/'allocation.json').read_text())
    quotes=pd.read_csv(output/'contracts.csv');costs=pd.read_csv(output/'costs.csv')
    for relative,digest in record['implementation_sha256'].items():
        assert hashlib.sha256((output/'source'/relative).read_bytes()).hexdigest()==digest
    seeds=[s for row in record['seed_records'] for s in row['seeds']]
    assert len(seeds)==len(set(seeds)), 'Seed collision'
    coefficients={int(k):np.asarray(v) for k,v in allocation['coefficients'].items()}
    direct_coefficients=np.asarray(allocation['direct_coefficients']);forwards=quotes.forward.to_numpy()
    def check_weights(archive,prefix):
        w,z=archive[prefix+'_w'],archive[prefix+'_z']
        np.testing.assert_allclose(proposal(1).weights(z),w,rtol=1e-12,atol=1e-14)
        assert np.all((w>0)&(w<=2))
    with np.load(output/'design_draws.npz') as design:
        for level in range(4):
            prefix=f'train_l{level}';check_weights(design,prefix)
            y,c=design[prefix+'_y'],design[prefix+'_c'];n=len(y)
            splits=[np.arange(n),np.arange(n//2),np.arange(n//2,n)]
            for k,i in enumerate(splits):
                np.testing.assert_allclose(fit_controls(y[i],c[i]),coefficients[level][k],rtol=1e-10,atol=1e-10)
                if level==3:
                    np.testing.assert_allclose(fit_controls(design[prefix+'_fine_y'][i],design[prefix+'_fine_c'][i]),direct_coefficients[k],rtol=1e-10,atol=1e-10)
            prefix=f'alloc_l{level}';check_weights(design,prefix)
            x=draw_estimators(design[prefix+'_y'],design[prefix+'_c'],coefficients[level],forwards)[:,:,1]
            v=(x[:,[2,3,8,9]]/quotes.price_scale_U.to_numpy()[[2,3,8,9]]).var(axis=0,ddof=1).mean()
            np.testing.assert_allclose(v,allocation['variances'][level],rtol=1e-10)
            rows=costs[(costs.role<2)&(costs.level==level)]
            np.testing.assert_allclose(rows.seconds.sum()/rows.paths.sum(),allocation['cost_per_path'][level],rtol=1e-10)
            if level==3:
                np.testing.assert_allclose(rows.fine_seconds.sum()/rows.paths.sum(),allocation['direct_cost_per_path'],rtol=1e-10)
    counts=allocate_levels(allocation['variances'],allocation['cost_per_path'],allocation['budget_seconds'],allocation['minimum'],allocation['quantum'])
    np.testing.assert_array_equal(counts,allocation['counts'])
    assert counts@allocation['cost_per_path']<=allocation['budget_seconds']*(1+1e-11)
    estimates=pd.read_csv(output/'estimates.csv');levels=pd.read_csv(output/'levels.csv');diagnostics=pd.read_csv(output/'diagnostics.csv')
    labels=['raw','pooled','halfA','halfB'];summary={}
    with np.load(output/'checking_draws.npz') as checking:
        for trial in range(record['arguments']['trials']):
            draws=[]
            for level in range(4):
                prefix=f'mlmc_t{trial}_l{level}';check_weights(checking,prefix)
                x=draw_estimators(checking[prefix+'_y'],checking[prefix+'_c'],coefficients[level],forwards)
                np.testing.assert_allclose(x,checking[prefix+'_adjusted'],rtol=1e-11,atol=1e-10)
                assert len(x)==counts[level];draws.append(x)
            summary[('mlmc',trial)]=independent_level_statistics(draws)
            prefix=f'direct_t{trial}';check_weights(checking,prefix)
            x=draw_estimators(checking[prefix+'_y'],checking[prefix+'_c'],direct_coefficients,forwards)
            np.testing.assert_allclose(x,checking[prefix+'_adjusted'],rtol=1e-11,atol=1e-10)
            summary[('direct',trial)]=(x.mean(axis=0),x.std(axis=0,ddof=1)/np.sqrt(len(x)))
        for row in estimates.itertuples():
            mean,se=summary[(row.method,row.trial)];k=labels.index(row.estimator)
            np.testing.assert_allclose([mean[row.contract,k],se[row.contract,k]],[row.price,row.se],rtol=1e-11,atol=1e-10)
        for row in levels.itertuples():
            prefix=f'mlmc_t{row.trial}_l{row.level}';x=checking[prefix+'_adjusted'][:,row.contract,labels.index(row.estimator)]
            np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x)),x.var(ddof=1)],[row.mean,row.se,row.variance],rtol=1e-11,atol=1e-10)
            hit=checking[prefix+'_hits'][:,row.contract]
            np.testing.assert_array_equal(hit.sum(axis=0),[row.fine_hits,row.coarse_hits,row.nonzero_corrections])
        for row in diagnostics.itertuples():
            prefix=f'mlmc_t{row.trial}_l{row.level}' if row.method=='mlmc' else f'direct_t{row.trial}'
            index=['mass','stock29','stock57','reference29','reference57'].index(row.statistic)
            x=checking[prefix+'_monitors'][:,index]
            np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x))],[row.mean,row.se],rtol=1e-11,atol=1e-12)
    assert len(estimates)==record['arguments']['trials']*12*4*2
    assert len(levels)==record['arguments']['trials']*4*12*4
    print(f'Validated allocation {counts.tolist()}, exact path weights, fits, {len(estimates)} estimates, {len(levels)} level rows and {len(diagnostics)} monitors; bias remains unresolved.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    verify(parser.parse_args().output)
