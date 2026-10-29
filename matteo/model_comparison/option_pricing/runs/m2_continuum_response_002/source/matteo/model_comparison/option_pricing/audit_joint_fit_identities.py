"""Numerical counterparts to Lean's stock-step check and VIX moment diagnostics."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .spx_conditional_control_study import setup
from .joint_curve_core import CurvePathBatch
from .joint_curve_fit import VIX_PARENT
from esn_eq.pricing import ForwardVarianceCurve
from matteo.preprocessing.option_variance import vix_option_moment_bound


def stats(values):
    return float(np.mean(values)), float(np.std(values, ddof=1)/np.sqrt(len(values)))


def audit(run):
    record = json.loads((run/'run.json').read_text())
    _, _, _, _, anchor_curve, _, master, _, _ = setup()
    scratch = Path(record['scratch'])
    spx = pd.read_csv(run/'spx_contracts.csv'); vix = pd.read_csv(run/'vix_contracts.csv')
    levels = dict(anchor=np.asarray(record['anchor_levels']), selected=np.asarray(record['selected_levels']))
    stock_rows, moment_rows, moment_pool = [], [], {}
    full_beta = json.loads((VIX_PARENT/'coefficients.json').read_text())
    for trial in range(3):
        for days, group in spx.groupby('expiry_days', sort=True):
            grid = master[:np.searchsorted(master, group.T_years.iloc[0])+1]
            with np.load(run/f'spx_r2_t{trial}_d{days}.npz') as a:
                weights, anchor_stock = a['weights'], a['anchor_stock']
            selected_stock = []
            for batch in range(4):
                stem = scratch/f'r2_t{trial}_d{days}_b{batch}'
                path = CurvePathBatch(np.load(str(stem)+'_ratio.npy', mmap_mode='r'), np.load(str(stem)+'_equity.npy', mmap_mode='r'),
                    np.diff(grid), np.minimum(np.searchsorted(anchor_curve.ends, grid[:-1], side='right'), 2))
                selected_stock.append(path.terminal(levels['selected'], derivative=False)[0])
            for candidate, stock in [('anchor', anchor_stock), ('selected', np.concatenate(selected_stock))]:
                for target, value in [('mass', weights), ('stock', weights*stock)]:
                    mean, se = stats(value)
                    stock_rows.append(dict(trial=trial, expiry_days=int(days), candidate=candidate, target=target,
                        estimate=mean, se=se, expected=1., warning_2se=bool(abs(mean-1)>2*se)))
        for days, group in vix.groupby('expiry_days', sort=True):
            with np.load(run/f'vix_r2_t{trial}_d{days}.npz') as a:
                beta = np.asarray(full_beta[str(days)]['coefficients']['gaussian'])[:, 1]
                for candidate, xi in levels.items():
                    integrated = float(a['floor_integral'])+a['exposures']@(xi-1e-4)
                    squared = 10000*integrated/float(a['delta'])
                    adj = a['weights']*squared-a['controls']@beta
                    repetitions = adj.reshape(8, 64).mean(axis=1)
                    moment_pool.setdefault((int(days), candidate), []).append(repetitions)
                    curve = ForwardVarianceCurve(tuple(record['family']['ends']), tuple(xi))
                    t = float(group.T_years.iloc[0]); target = 10000/float(a['delta'])*curve.integral(t, t+float(a['delta']))
                    mean, se = stats(repetitions)
                    moment_rows.append(dict(trial=trial, expiry_days=int(days), candidate=candidate,
                        model_target=target, estimate=mean, replication_se=se, residual=mean-target,
                        warning_2se=bool(abs(mean-target)>2*se)))
    pooled_rows, compatibility = [], []
    for (days, candidate), arrays in moment_pool.items():
        a = np.concatenate(arrays); mean, se = stats(a)
        target = next(row['model_target'] for row in moment_rows if row['expiry_days'] == days and row['candidate'] == candidate)
        pooled_rows.append(dict(expiry_days=days, candidate=candidate, model_target=target, estimate=mean, replication_se=se,
            residual=mean-target, warning_2se=bool(abs(mean-target)>2*se)))
        table = pd.read_csv(run/f'forwards_fresh_pooled_{candidate}.csv')
        row = table[table.expiry_days == days].iloc[0]
        group = vix[vix.expiry_days == days]
        bound = vix_option_moment_bound(group, float(row.forward_vix_points), float(row.model_forward))
        compatibility.append(dict(expiry_days=days, candidate=candidate, model_mean=float(row.model_forward),
            forward_sampling_se=float(row.mc_se), model_second_moment=target,
            necessary_option_second_moment_bound=bound['moment_lower_bound'],
            necessary_margin=target-bound['moment_lower_bound'], jensen_margin=target-float(row.model_forward)**2,
            caveat='Necessary market-bid bound evaluated at estimated model mean. Negative margin is retained with mean uncertainty; not a certified incompatibility verdict.'))
    for name, rows in [('stock_forward_checks', stock_rows), ('vix_moment_checks', moment_rows),
            ('vix_moment_pooled', pooled_rows), ('necessary_moment_compatibility', compatibility)]:
        pd.DataFrame(rows).to_csv(run/(name+'.csv'), index=False)
    result = dict(stock_monitor_warning_count=sum(row['warning_2se'] for row in stock_rows),
        pooled_vix_moment_warnings=[row for row in pooled_rows if row['warning_2se']],
        necessary_negative_margins=[row for row in compatibility if row['necessary_margin'] < 0],
        scope='New fixed selected/anchor curves on fresh saved draws; all four maturities retained. Statistical warnings are dependent diagnostics, not proofs of bias or price convergence.')
    (run/'identity_audit.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--run', type=Path, required=True)
    audit(parser.parse_args().run)
