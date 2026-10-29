"""Recompute all stored SPX pilot estimates and paired corrections from draws."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def verify(output):
    record = json.loads((output/'run.json').read_text())
    for relative, digest in record['implementation_sha256'].items():
        assert hashlib.sha256((output/'source'/relative).read_bytes()).hexdigest() == digest
    seeds = [seed for row in record['seed_records'] for role in ('pilot','evaluation') for seed in row[role]]
    assert len(seeds) == len(set(seeds)), 'Pilot/evaluation seed collision'
    estimates = pd.read_csv(output/'estimates.csv')
    changes = pd.read_csv(output/'paired_changes.csv')
    diagnostics = pd.read_csv(output/'diagnostics.csv')
    with np.load(output/'checking_draws.npz') as draws:
        for row in estimates.itertuples():
            level = {64:0,128:1,256:2}[row.steps]
            values = draws[f't{row.trial}_{row.method}_g{level}_q{row.normalizer_order}'][:,row.contract, int(row.estimator == 'hedging')]
            np.testing.assert_allclose([values.mean(), values.std(ddof=1)/np.sqrt(len(values))], [row.price,row.se], rtol=1e-11, atol=1e-12)
            assert 0 <= row.positive_payoffs <= row.paths
        levels = {'64_to_128':((0,48),(1,48)), '128_to_256':((1,48),(2,48)), 'GH48_to_64':((2,48),(2,64))}
        for row in changes.itertuples():
            low, high = levels[row.change_type]
            def get(v):
                return draws[f't{row.trial}_{row.method}_g{v[0]}_q{v[1]}'][:,row.contract,int(row.estimator == 'hedging')]
            diff = get(high)-get(low)
            np.testing.assert_allclose([diff.mean(), diff.std(ddof=1)/np.sqrt(len(diff)),diff.var(ddof=1)], [row.change,row.se,row.correction_variance], rtol=1e-11, atol=1e-12)
        labels = ['mass','stock29','stock57','reference29','reference57']
        for row in diagnostics.itertuples():
            values = draws[f't{row.trial}_{row.method}_monitors'][:, labels.index(row.statistic)]
            np.testing.assert_allclose([values.mean(),values.std(ddof=1)/np.sqrt(len(values))], [row.estimate,row.se], rtol=1e-11, atol=1e-12)
            if row.statistic == 'mass':
                assert np.all((values > 0) & (values <= 2))
                if row.method == 'target':
                    np.testing.assert_array_equal(values, 1.)
    assert len(estimates) == record['arguments']['trials']*2*4*12*2
    assert len(changes) == record['arguments']['trials']*2*3*12*2
    print(f'Validated {len(estimates)} estimates, {len(changes)} paired changes, {len(diagnostics)} monitors and source hashes; accuracy remains unqualified.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    verify(parser.parse_args().output)
