"""Replay frozen conditional-control artifacts and report bounded decision metrics."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .m2_spx import ROOT
from .spx_conditional_control_study import setup, fit_pilots, adjusted, METHODS
from .spx_volatility_sampling import proposal
from esn_eq.conditional_spx import partial_transition


def audit(run):
    record = json.loads((run/'run.json').read_text())
    for name, digest in record['implementation_sha256'].items():
        snapshot = run/'source'/name
        assert hashlib.sha256(snapshot.read_bytes()).hexdigest() == digest, name
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    seeds = [seed for expiry in record['expiry_records'] for stream in expiry['seed_records'] for seed in stream['seeds']]
    assert len(seeds) == len(set(seeds)), 'Pilot/check or cross-method seed collision'
    estimates, costs, diagnostics = (pd.read_csv(run/(name+'.csv')) for name in ['estimates', 'costs', 'diagnostics'])
    comparisons, sensitivity, differences, pilots, budget = [], [], [], [], []
    replayed = 0
    for expiry in record['expiry_records']:
        days = expiry['expiry_days']
        out = run/f'expiry_{days}'
        quotes = pd.read_csv(out/'contracts.csv')
        frozen = json.loads((out/'frozen.json').read_text())
        with np.load(out/'pilot_draws.npz') as archive:
            for method in METHODS:
                coefficients = fit_pilots(archive[method+'_y'], archive[method+'_c'])
                np.testing.assert_allclose(coefficients, frozen['coefficients'][method], rtol=1e-12, atol=1e-12)
                np.testing.assert_allclose(archive[method+'_w'], proposal(1).weights(archive[method+'_scores']), atol=1e-14)
                for j, row in enumerate(quotes.itertuples()):
                    pilots.append(dict(expiry_days=days, method=method, contract=j, band=row.band,
                        positive_paths=int(archive[method+'_positive'][:, j].sum()),
                        max_weighted_pilot_payoff_U=float(archive[method+'_y'][:, j].max()*row.forward),
                        mean_weighted_pilot_payoff_U=float(archive[method+'_y'][:, j].mean()*row.forward)))
        data_by_method = {}
        with np.load(out/'checking_draws.npz') as archive:
            for trial in range(2):
                data_by_method.clear()
                for method in METHODS:
                    prefix = f'{method}_t{trial}_'
                    data = {key: archive[prefix+key] for key in ['y', 'c', 'w', 'scores', 'components', 'positive', 'monitors']}
                    np.testing.assert_allclose(data['w'], proposal(1).weights(data['scores']), atol=1e-14)
                    assert data['w'].max() <= 2.+1e-14
                    np.testing.assert_allclose(data['c'][:, :, -1], np.broadcast_to((data['w']-1.)[:, None], data['c'].shape[:2]), atol=1e-14)
                    x = adjusted(data, np.asarray(frozen['coefficients'][method]), quotes.forward.to_numpy())
                    np.testing.assert_allclose(x, archive[prefix+'adjusted'], rtol=1e-12, atol=1e-12)
                    rows = estimates[(estimates.expiry_days == days)&(estimates.trial == trial)&(estimates.method == method)]
                    for j, quote in enumerate(quotes.itertuples()):
                        for v, variant in enumerate(['full', 'half_A', 'half_B']):
                            a = x[:, j, v]
                            row = rows[(rows.contract == j)&(rows.variant == variant)].iloc[0]
                            squares = (a-a.mean())**2
                            share = squares.max()/squares.sum() if squares.sum() else 0.
                            np.testing.assert_allclose([a.mean(), a.std(ddof=1)/np.sqrt(len(a)), a.var(ddof=1), share],
                                row[['estimate', 'se', 'variance', 'largest_residual_share']].to_numpy(float), rtol=1e-11, atol=1e-12)
                            replayed += 1
                        base_se = x[:, j, 0].std(ddof=1)/np.sqrt(len(x))
                        for v, variant in [(1, 'half_A'), (2, 'half_B')]:
                            delta = x[:, j, v]-x[:, j, 0]
                            sensitivity.append(dict(expiry_days=days, trial=trial, method=method, contract=j,
                                strike=quote.strike, band=quote.band, variant=variant,
                                se_ratio=float(x[:, j, v].std(ddof=1)/x[:, j, 0].std(ddof=1)),
                                paired_mean=float(delta.mean()), paired_se=float(delta.std(ddof=1)/np.sqrt(len(delta)))))
                        budget.append(dict(expiry_days=days, trial=trial, method=method, contract=j, band=quote.band,
                            observed_paths=len(x), plugin_required_paths=float(len(x)*(base_se/(quote.price_scale_U/8))**2),
                            caveat='Fixed observed variance only; excludes missed tails, pilot changes and deterministic bias.'))
                    # Replay all known-mean monitors/control rows.
                    for row in diagnostics[(diagnostics.expiry_days == days)&(diagnostics.trial == trial)&(diagnostics.method == method)].itertuples():
                        if row.target.startswith('contract'):
                            contract, control = row.target.removeprefix('contract').split('_control')
                            a = data['c'][:, int(contract), int(control)]
                        else:
                            a = data['monitors'][:, ['mass', 'stock', 'reference_stock'].index(row.target)]
                        np.testing.assert_allclose([a.mean(), a.std(ddof=1)/np.sqrt(len(a))], [row.mean, row.se], rtol=1e-11, atol=1e-12)
                    data_by_method[method] = x
                # Independently sampled estimators: never use a paired cross-method SE.
                a, b = data_by_method['reference'][:, :, 0], data_by_method['conditional'][:, :, 0]
                for j, q in enumerate(quotes.itertuples()):
                    differences.append(dict(expiry_days=days, trial=trial, contract=j, strike=q.strike,
                        difference=float(b[:, j].mean()-a[:, j].mean()),
                        independent_se=float(np.sqrt(a[:, j].var(ddof=1)/len(a)+b[:, j].var(ddof=1)/len(b)))))
                for band_name, selected in [('all', np.ones(len(quotes), bool)), ('inner', quotes.band.to_numpy() == 0),
                        ('middle', quotes.band.to_numpy() == 1), ('outer', quotes.band.to_numpy() == 2)]:
                    scale = quotes.price_scale_U.to_numpy()[selected]
                    va = np.mean(a[:, selected].var(axis=0, ddof=1)/scale**2)
                    vb = np.mean(b[:, selected].var(axis=0, ddof=1)/scale**2)
                    seconds = {method: costs[(costs.expiry_days == days)&(costs.method == method)&(costs.role == 1)&(costs.trial == trial)].seconds.sum() for method in METHODS}
                    comparisons.append(dict(expiry_days=days, trial=trial, band=band_name,
                        reference_scaled_variance=float(va), conditional_scaled_variance=float(vb),
                        variance_gain=float(va/vb), reference_seconds=float(seconds['reference']),
                        conditional_seconds=float(seconds['conditional']),
                        efficiency_gain=float(va*seconds['reference']/(vb*seconds['conditional']))))
    result = estimates.copy()
    result['nominal_se_screen'] = (result.estimate > 0)&(result.se <= result.scale/8)
    # A conditional probability above floating-point underflow is not a raw
    # exercise hit. The first audit mistakenly called the count-based proxy a
    # sampling screen; preserve its outputs and use accurate labels here.
    result['positive_path_count_at_least_two'] = result.positive_paths >= 2
    result['positive_path_count_is_raw_exercise_hits'] = result.method == 'reference'
    result['negative_estimate'] = result.estimate < 0
    result.to_csv(run/'screens.csv', index=False)
    for name, rows in [('efficiency', comparisons), ('pilot_sensitivity', sensitivity), ('independent_differences', differences),
            ('pilot_hits', pilots), ('plugin_budgets', budget)]:
        pd.DataFrame(rows).to_csv(run/(name+'.csv'), index=False)
    diagnostics['z_warning'] = abs(diagnostics['mean']-diagnostics.expectation) > 2*diagnostics.se
    diagnostics.to_csv(run/'control_mean_audit.csv', index=False)
    # Recreate covariance audits without repeating any priced path experiments.
    _, p, _, _, _, _, master, _, _ = setup()
    transitions = [partial_transition(p, float(dt))[2] for dt in np.unique(np.diff(master))]
    (run/'transition_audit.json').write_text(json.dumps(transitions, indent=2, allow_nan=False)+'\n')
    efficiency = pd.DataFrame(comparisons)
    global_efficiency = []
    for trial in range(2):
        rows = efficiency[(efficiency.trial == trial)&(efficiency.band == 'all')]
        gain = (rows.reference_scaled_variance.mean()*rows.reference_seconds.sum()/
            (rows.conditional_scaled_variance.mean()*rows.conditional_seconds.sum()))
        global_efficiency.append(dict(trial=trial, efficiency_gain=float(gain),
            reference_seconds=float(rows.reference_seconds.sum()), conditional_seconds=float(rows.conditional_seconds.sum())))
    summary = dict(replayed_prices=replayed, captured_sources=len(record['implementation_sha256']),
        unique_pilot_check_batch_seeds=len(seeds),
        primary=result[result.variant == 'full'].groupby('method').agg(prices=('estimate', 'size'),
            negative=('negative_estimate', 'sum'), nominal_se_screens=('nominal_se_screen', 'sum'),
            median_se_scale=('se', lambda x: float(np.median(x/result.loc[x.index, 'scale'])))).reset_index().to_dict('records'),
        numerical_quote_clearance=False,
        screen_caveat='SE threshold only; conditional positive-path counts are not exercise hits. Both nominal passes have tiny pilot targets and >99% single-path residual concentration; tail precision remains unresolved.',
        max_standardized_covariance_error=max(t['max_standardized_covariance_error'] for t in transitions),
        known_mean_warning_count=int(diagnostics.z_warning.sum()),
        all_expiry_efficiency=efficiency.query('band == "all"').to_dict('records'),
        global_efficiency=global_efficiency,
        wall_seconds=record['wall_seconds'],
        conclusion='Saved-artifact and covariance correspondence only; not a population error certificate.')
    (run/'audit.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    audit(parser.parse_args().run)
