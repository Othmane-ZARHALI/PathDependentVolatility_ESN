"""Reconstruct searching predictions from preserved fixed-shock disk caches.

This is a same-sample correspondence check, never a new independent check.
The temporary caches are optional acceleration artifacts; seeds/source snapshots
remain available for regeneration if they are later removed.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .spx_conditional_control_study import setup
from .joint_curve_core import CurvePathBatch
from .joint_curve_fit import SPXData, VIXData, FitData


def load_search(run):
    record = json.loads((run/'design.json').read_text())
    _, _, _, _, curve, _, master, _, _ = setup()
    scratch = Path(record['scratch'])
    spx = pd.read_csv(run/'spx_contracts.csv'); vix = pd.read_csv(run/'vix_contracts.csv')
    primary = pd.read_csv(run/'spx_membership.csv').priority_selected.to_numpy(bool)
    forwards = pd.read_csv(run/'forwards_search_t0_anchor.csv')
    spx_data, vix_data = {}, {}
    for days, q in spx.groupby('expiry_days', sort=True):
        grid = master[:np.searchsorted(master, q.T_years.iloc[0])+1]
        paths = []
        for batch in range(4):
            stem = scratch/f'r1_t0_d{days}_b{batch}'
            paths.append(CurvePathBatch(np.load(str(stem)+'_ratio.npy', mmap_mode='r'),
                np.load(str(stem)+'_equity.npy', mmap_mode='r'), np.diff(grid),
                np.minimum(np.searchsorted(curve.ends, grid[:-1], side='right'), 2)))
        with np.load(run/f'spx_r1_t0_d{days}.npz') as archive:
            data = SPXData(q.reset_index(drop=True), paths, archive['weights'].copy(), archive['controls'].copy(), archive['anchor_stock'].copy())
            data.beta = archive['beta'].copy()
            spx_data[int(days)] = data
    for days, q in vix.groupby('expiry_days', sort=True):
        with np.load(run/f'vix_r1_t0_d{days}.npz') as a:
            vix_data[int(days)] = VIXData(q.reset_index(drop=True), a['exposures'].copy(), float(a['floor_integral']),
                float(a['delta']), a['weights'].copy(), a['controls'].copy(), a['beta'].copy(), 8)
    return FitData(spx_data, vix_data, spx, vix, forwards, primary)


def replay(run):
    data = load_search(run)
    selected = np.asarray(json.loads((run/'selected.json').read_text())['selected_levels'])
    result = data.predict(selected)
    q = pd.read_csv(run/'prices_search_t0_selected.csv')
    np.testing.assert_allclose(result['prices'][data.primary], q.model_price_U, rtol=1e-10, atol=1e-10)
    singular = np.linalg.svd(result['jacobian'], compute_uv=False)
    relative = np.linalg.svd(result['jacobian']*selected[None, :], compute_uv=False)
    negatives = []
    rows = [json.loads(line) for line in (run/'evaluations.jsonl').read_text().splitlines()]
    for row in rows:
        if row['negative_count']:
            candidate = data.predict(np.asarray(row['levels']), derivative=False)
            assert candidate['negative_count'] == row['negative_count']
            failed = data.all_quotes.copy(); failed['model_price_U'] = candidate['prices']; failed['mc_se_U'] = candidate['errors']
            failed['objective_member'] = data.primary
            failed.to_csv(run/f'search_signed_evaluation_{row["evaluation"]:03d}.csv', index=False)
            negatives.append(dict(evaluation=row['evaluation'], levels=row['levels'], negative_count=row['negative_count']))
    gradient_errors = []
    for days, spx in data.spx.items():
        batch = spx.batches[0]
        _, derivative = batch.terminal(selected)
        for j in range(3):
            d = np.eye(3)[j]*1e-6
            fd = (batch.terminal(selected+d, derivative=False)[0]-batch.terminal(selected-d, derivative=False)[0])/2e-6
            error = float(np.max(abs(fd-derivative[:, j])))
            np.testing.assert_allclose(fd, derivative[:, j], rtol=2e-5, atol=2e-7)
            gradient_errors.append(dict(days=days, coordinate=j, max_absolute_stock_derivative_error=error))
    audit = dict(selected_prices_replayed=len(q), negative_search_points_retained=negatives,
        selected_residual_jacobian_singular_values=singular.tolist(), relative_parameter_jacobian_singular_values=relative.tolist(),
        stock_derivative_checks=gradient_errors,
        scope='Same frozen searching paths, 128 paths per expiry for derivative checks. Singular values describe this numerical curve-only problem, not population or full-model identification.')
    (run/'cache_replay.json').write_text(json.dumps(audit, indent=2, allow_nan=False)+'\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--run', type=Path, required=True)
    replay(parser.parse_args().run)
