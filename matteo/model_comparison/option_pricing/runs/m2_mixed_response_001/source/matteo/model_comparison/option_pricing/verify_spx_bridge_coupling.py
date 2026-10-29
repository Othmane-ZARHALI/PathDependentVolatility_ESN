"""Replay reflected-pair likelihood bookkeeping, pilot fits and paired uncertainty."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_bridge_coupling import fit_variants
from .spx_mlmc import draw_estimators
from .spx_path_importance import weighted_controls
from .spx_volatility_sampling import proposal


def verify(output):
    record=json.loads((output/'run.json').read_text());frozen=json.loads((output/'frozen.json').read_text())
    quotes=pd.read_csv(output/'contracts.csv');costs=pd.read_csv(output/'costs.csv')
    for relative,digest in record['implementation_sha256'].items():
        assert hashlib.sha256((output/'source'/relative).read_bytes()).hexdigest()==digest
    seeds=[s for row in record['seed_records'] for s in row['seeds']]
    assert len(seeds)==len(set(seeds)), 'Pilot/checking seed collision'
    assert costs.endpoint_error.max()<5e-9 and costs.equity_error.max()<5e-10
    for diagnostics in record['kernel_diagnostics'].values():
        for row in diagnostics:
            assert max(row[k] for k in ('endpoint_residual','orthogonal_residual','involution_residual'))<1e-10
    def check_likelihood(archive,prefix):
        w,z=archive[prefix+'_w'],archive[prefix+'_z']
        np.testing.assert_allclose(proposal(1).weights(z),w,rtol=1e-12,atol=1e-14)
        assert np.all((w>0)&(w<=2))
        f,a,c=[archive[prefix+'_'+k] for k in ('fine_payoffs','anti_payoffs','coarse_payoffs')]
        np.testing.assert_allclose(archive[prefix+'_plain_y'],w[:,None]*(f-c),rtol=1e-12,atol=1e-14)
        np.testing.assert_allclose(archive[prefix+'_bridge_y'],w[:,None]*(.5*(f+a)-c),rtol=1e-12,atol=1e-14)
        # This deliberately uses the ORIGINAL weight, not f/q at anti_z.
    with np.load(output/'pilot_draws.npz') as pilot:
        for level in (2,3):
            prefix=f'pilot_l{level}';check_likelihood(pilot,prefix)
            data={key:pilot[prefix+'_'+key] for key in ('plain_y','plain_c','bridge_y','bridge_c')}
            for method in ('plain','bridge'):
                np.testing.assert_allclose(fit_variants(data,method),frozen['coefficients'][str(level)][method],rtol=1e-10,atol=1e-10)
    labels=['raw','pooled','halfA','halfB'];rows=pd.read_csv(output/'estimates.csv');differences=pd.read_csv(output/'paired_differences.csv')
    with np.load(output/'checking_draws.npz') as checking:
        for trial in range(record['arguments']['trials']):
            for level in (2,3):
                prefix=f't{trial}_l{level}';check_likelihood(checking,prefix)
                for method in ('plain','bridge'):
                    x=draw_estimators(checking[prefix+'_'+method+'_y'],checking[prefix+'_'+method+'_c'],np.asarray(frozen['coefficients'][str(level)][method]),quotes.forward.to_numpy())
                    np.testing.assert_allclose(x,checking[prefix+'_'+method+'_adjusted'],rtol=1e-11,atol=1e-10)
        for row in rows.itertuples():
            prefix=f't{row.trial}_l{row.level}';x=checking[prefix+'_'+row.method+'_adjusted'][:,row.contract,labels.index(row.estimator)]
            assert len(x)==row.pairs
            np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x)),x.var(ddof=1)],[row.correction,row.se,row.variance],rtol=1e-11,atol=1e-10)
            for label,hit in (('fine_payoffs',row.fine_hits),('anti_payoffs',row.reflected_hits),('coarse_payoffs',row.coarse_hits)):
                assert int((checking[prefix+'_'+label][:,row.contract]>0).sum())==hit
        for row in differences.itertuples():
            prefix=f't{row.trial}_l{row.level}'
            x=(checking[prefix+'_bridge_adjusted']-checking[prefix+'_plain_adjusted'])[:,row.contract,labels.index(row.estimator)]
            np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x))],[row.difference,row.se],rtol=1e-11,atol=1e-10)
        for row in pd.read_csv(output/'diagnostics.csv').itertuples():
            index=['mass','stock29','stock57','reflected_stock29','reflected_stock57','reference29','reference57'].index(row.statistic)
            x=checking[f't{row.trial}_l{row.level}_monitors'][:,index]
            np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x))],[row.mean,row.se],rtol=1e-11,atol=1e-12)
    assert len(rows)==record['arguments']['trials']*2*2*12*4
    assert len(differences)==record['arguments']['trials']*2*12*4
    print(f'Validated reflection kernels, original weights on both payoffs, frozen fits, {len(rows)} correction estimates and {len(differences)} paired differences; no accuracy certificate.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    verify(parser.parse_args().output)
