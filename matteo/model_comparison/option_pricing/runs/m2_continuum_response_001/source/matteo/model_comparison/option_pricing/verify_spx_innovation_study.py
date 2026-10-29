"""Replay innovation-control dictionaries, pilot fits and fresh paired estimates."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_innovation_study import METHODS,datasets,train,adjusted
from .spx_volatility_sampling import proposal


def verify(output):
    record=json.loads((output/'run.json').read_text());frozen=json.loads((output/'frozen.json').read_text());quotes=pd.read_csv(output/'contracts.csv')
    for relative,digest in record['implementation_sha256'].items():
        assert hashlib.sha256((output/'source'/relative).read_bytes()).hexdigest()==digest
    seeds=[s for row in record['seed_records'] for s in row['seeds']];assert len(seeds)==len(set(seeds))
    def check(data):
        np.testing.assert_allclose(proposal(1).weights(data['z']),data['w'],rtol=1e-12,atol=1e-14)
        assert np.all((data['w']>0)&(data['w']<=2))
        for name,(y,base,vol) in datasets(data).items():
            assert base.shape==(*y.shape,5) and vol.shape==(*y.shape,3 if record['arguments'].get('linear_only') else 9)
            np.testing.assert_array_equal(base[:,:,-1],np.broadcast_to(data['w'][:,None]-1.,y.shape))
    with np.load(output/'pilot_draws.npz') as pilot:
        data={key:pilot[key] for key in pilot.files};check(data)
        for name,triple in datasets(data).items():
            fits=train(*triple)
            for method,beta in fits.items():np.testing.assert_allclose(beta,frozen['coefficients'][name][method],rtol=1e-10,atol=1e-10)
    e=pd.read_csv(output/'estimates.csv');paired=pd.read_csv(output/'paired_differences.csv');diagnostics=pd.read_csv(output/'diagnostics.csv')
    controls_by_target={}
    with np.load(output/'checking_draws.npz') as checking:
        for trial in range(record['arguments']['trials']):
            prefix=f't{trial}_';data={key[len(prefix):]:checking[key] for key in checking.files if key.startswith(prefix) and not key.endswith('_adjusted')};check(data)
            for name,triple in datasets(data).items():
                controls_by_target[(trial,name)]=np.concatenate([triple[1],triple[2]],axis=2)
                x=adjusted(*triple,frozen['coefficients'][name],quotes.forward.to_numpy())
                np.testing.assert_allclose(x,checking[prefix+name+'_adjusted'],rtol=1e-11,atol=1e-10)
        for row in e.itertuples():
            x=checking[f't{row.trial}_{row.target}_adjusted'][:,row.contract,METHODS.index(row.method)]
            np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x)),x.var(ddof=1)],[row.estimate,row.se,row.variance],rtol=1e-11,atol=1e-10)
        for row in paired.itertuples():
            x=checking[f't{row.trial}_{row.target}_adjusted'][:,row.contract]
            d=x[:,3]-x[:,2]
            np.testing.assert_allclose([d.mean(),d.std(ddof=1)/np.sqrt(len(d))],[row.difference,row.se],rtol=1e-11,atol=1e-10)
        for row in diagnostics.itertuples():
            if row.target=='forward':
                index=['mass','stock29','stock57','reference29','reference57'].index(row.statistic)
                x=checking[f't{row.trial}_monitors'][:,index]
            else:
                level={128:1,256:2,512:3}[int(row.target.removeprefix('innovations'))]
                x=checking[f't{row.trial}_g{level}_innovation_moments'][:,int(row.statistic)]
            np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x))],[row.mean,row.se],rtol=1e-11,atol=1e-12)
    control_count=0
    if (output/'control_means.csv').exists():
        control_labels=['stock','reference_payoff','deterministic_delta','spot_delta','mass','ou_linear_C','ou_linear_Y','ou_linear_F',
                        'ou_square_C','ou_square_Y','ou_square_F','ou_mixed_C','ou_mixed_Y','ou_mixed_F']
        control_rows=pd.read_csv(output/'control_means.csv');control_count=len(control_rows)
        for row in control_rows.itertuples():
            x=controls_by_target[(row.trial,row.target)][:,row.contract,control_labels.index(row.control)]
            np.testing.assert_allclose([x.mean(),x.std(ddof=1)/np.sqrt(len(x))],[row.mean,row.se],rtol=1e-11,atol=1e-12)
    if (output/'audit.json').exists():
        audit=json.loads((output/'audit.json').read_text())
        with np.load(output/'checking_draws.npz') as checking:
            for row in audit.get('paired_tail_concentration',[]):
                x=checking[f"t{row['trial']}_{row['target']}_adjusted"][:,row['contract']]
                difference=x[:,3]-x[:,2];squares=difference*difference;index=int(np.argmax(squares))
                assert row['largest_path']==index
                share=float(squares[index]/squares.sum()) if squares.sum() else 0.
                np.testing.assert_allclose([difference[index],difference[index]/len(difference),share],
                    [row['largest_difference'],row['contribution_to_mean'],row['largest_share_uncentered_squares']],rtol=1e-12,atol=1e-12)
    targets=1 if record['arguments'].get('direct_only') else 5
    assert len(e)==record['arguments']['trials']*targets*12*len(METHODS)
    assert len(paired)==record['arguments']['trials']*targets*12
    print(f'Validated exact-weight bookkeeping, correction mass columns, frozen fits, {len(e)} estimates, {len(paired)} paired differences and {len(diagnostics)} bare monitors plus {control_count} contract-control means; bias/tails remain open.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path);verify(parser.parse_args().output)
