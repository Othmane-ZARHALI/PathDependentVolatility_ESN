"""Independent score reconstruction, paired fresh checks and frozen-fit diagnostics."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .m2_spx import ROOT


def score(quotes, targets, prices, forwards, levels):
    """Independent scalar implementation; signed numerical estimates retained."""
    residual = (np.minimum(prices-quotes.price_bid_U.to_numpy(), 0.)+
        np.maximum(prices-quotes.price_ask_U.to_numpy(), 0.))/quotes.price_scale_U.to_numpy()
    components = {}
    for family in ['SPX', 'VIX']:
        expiry_means = []
        for _, group in quotes[quotes.instrument == family].groupby('expiry'):
            expiry_means.append(float(np.mean(residual[group.index]**2)))
        components[family] = float(np.mean(expiry_means))
    components['forwards'] = float(np.mean(((forwards-targets.forward_vix_points.to_numpy())/.5)**2))
    components['curve_penalty'] = float(.1*np.sum(np.diff(levels)**2))
    return float(sum(components.values())), components, residual


def paired_loss(quotes, targets, a, b, fa, fb, anchor, selected):
    n = len(a)
    am, bm, fam, fbm = a.mean(axis=0), b.mean(axis=0), fa.mean(axis=0), fb.mean(axis=0)
    delta = score(quotes, targets, bm, fbm, selected)[0]-score(quotes, targets, am, fam, anchor)[0]
    leave = []
    for j in range(n):
        aa, bb = (a.sum(axis=0)-a[j])/(n-1), (b.sum(axis=0)-b[j])/(n-1)
        af, bf = (fa.sum(axis=0)-fa[j])/(n-1), (fb.sum(axis=0)-fb[j])/(n-1)
        leave.append(score(quotes, targets, bb, bf, selected)[0]-score(quotes, targets, aa, af, anchor)[0])
    se = float(np.sqrt((n-1)/n*np.sum((leave-np.mean(leave))**2)))
    return dict(selected_minus_anchor=float(delta), descriptive_block_jackknife_se=se,
        leave_one_block_range=[float(min(leave)), float(max(leave))], blocks=n,
        caveat='Nonlinear loss jackknife with few blocks and unresolved tails; not a calibrated confidence certificate.')


def audit(run):
    record = json.loads((run/'run.json').read_text())
    for name, digest in record['implementation_sha256'].items():
        assert hashlib.sha256((run/'source'/name).read_bytes()).hexdigest() == digest, name
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    seeds = [r['seed'] for r in record['seed_records']]
    assert len(seeds) == len(set(seeds)), 'Pilot/search/fresh seed collision'
    selected = np.asarray(record['selected_levels']); anchor = np.asarray(record['anchor_levels'])
    assert json.loads((run/'selected.json').read_text())['fresh_checking_not_used']
    spx = pd.read_csv(run/'spx_membership.csv'); vix = pd.read_csv(run/'vix_contracts.csv')
    all_quotes = pd.concat([spx, vix], ignore_index=True)
    mask = np.r_[spx.priority_selected.to_numpy(bool), np.ones(len(vix), bool)]
    q = all_quotes[mask].reset_index(drop=True)
    targets = pd.read_csv(run/'forwards_search_t0_anchor.csv')
    assert len(q) == 882 and mask[:len(spx)].sum() == 707 and len(vix) == 175
    rows, summaries, changes, blocks_report = [], [], [], []
    pooled = {name: [] for name in ['a', 'b', 'fa', 'fb']}
    raw_sources = {}
    with np.load(run/'comparison_blocks.npz') as blocks, np.load(run/'comparison_paths.npz') as paths:
        for role, trials in [('search', [0]), ('fresh', range(3))]:
            for trial in trials:
                result = {}
                for candidate, levels in [('anchor', anchor), ('selected', selected)]:
                    p = blocks[f'{role}_t{trial}_{candidate}_prices'][:, mask]
                    f = blocks[f'{role}_t{trial}_{candidate}_forwards']
                    prices, forwards = p.mean(axis=0), f.mean(axis=0)
                    saved = pd.read_csv(run/f'prices_{role}_t{trial}_{candidate}.csv')
                    savedf = pd.read_csv(run/f'forwards_{role}_t{trial}_{candidate}.csv')
                    np.testing.assert_allclose(prices, saved.model_price_U, rtol=1e-10, atol=1e-10)
                    np.testing.assert_allclose(forwards, savedf.model_forward, rtol=1e-12, atol=1e-12)
                    total, components, residual = score(q, targets, prices, forwards, levels)
                    entry = next(r for r in record['scores'] if (r['role'], r['trial'], r['candidate']) == (role, trial, candidate))
                    np.testing.assert_allclose(total, entry['total'], rtol=1e-11, atol=1e-10)
                    for key, value in components.items():
                        np.testing.assert_allclose(value, entry['components'][key], rtol=1e-11, atol=1e-10)
                    results = dict(role=role, trial=trial, candidate=candidate, total=total, components=components,
                        negatives=int(np.sum(prices < 0)), nominal_se_passes=int(saved.nominal_se_screen.sum()),
                        inside_count=int(np.sum(residual == 0)), median_se_over_scale=float(np.median(saved.mc_se_U/saved.price_scale_U)))
                    summaries.append(results)
                    for family in ['SPX', 'VIX']:
                        for days, group in q[q.instrument == family].groupby('expiry_days'):
                            ix = group.index.to_numpy()
                            values = paths[f'{role}_t{trial}_{candidate}_{family.lower()}_{days}']
                            if family == 'SPX':
                                fullgroup = spx[spx.expiry_days == days].reset_index(drop=True)
                                subset = fullgroup.priority_selected.to_numpy(bool)
                                np.testing.assert_allclose(values.mean(axis=0)[subset], prices[ix], rtol=1e-10, atol=1e-10)
                                np.testing.assert_allclose(values.std(axis=0, ddof=1)[subset]/np.sqrt(len(values)), saved.mc_se_U.to_numpy()[ix], rtol=1e-10, atol=1e-10)
                            else:
                                archive = np.load(run/f'vix_r{1 if role == "search" else 2}_t{trial}_d{days}.npz')
                                integral = float(archive['floor_integral'])+archive['exposures']@(levels-1e-4)
                                vv = 100*np.sqrt(integral/float(archive['delta']))
                                np.testing.assert_allclose(vv, values, rtol=1e-13, atol=1e-12)
                                sign = np.where(group.option_type.to_numpy() == 'call', 1., -1.)
                                h = np.c_[vv, np.maximum(sign*(vv[:, None]-group.strike.to_numpy()), 0.)]
                                adj = archive['weights'][:, None]*h-archive['controls']@archive['beta']
                                means = adj.reshape(8, 64, -1).mean(axis=1)
                                np.testing.assert_allclose(means[:, 1:], p[:, ix], rtol=1e-10, atol=1e-10)
                                np.testing.assert_allclose(means[:, 1:].std(axis=0, ddof=1)/np.sqrt(8), saved.mc_se_U.to_numpy()[ix], rtol=1e-10, atol=1e-10)
                                archive.close()
                            rows.append(dict(role=role, trial=trial, candidate=candidate, instrument=family, expiry_days=int(days),
                                quotes=len(ix), interval_mse=float(np.mean(residual[ix]**2)),
                                inside_count=int(np.sum(residual[ix] == 0)), negative_count=int(np.sum(prices[ix] < 0)),
                                max_abs_scaled_distance=float(np.max(abs(residual[ix]))),
                                median_se_scale=float(np.median(saved.mc_se_U.to_numpy()[ix]/group.price_scale_U.to_numpy()))))
                    result[candidate] = (p, f)
                a, fa = result['anchor']; b, fb = result['selected']
                diagnostic = paired_loss(q, targets, a, b, fa, fb, anchor, selected)
                changes.append(dict(role=role, trial=trial, **diagnostic))
                if role == 'fresh':
                    for key, array in [('a', a), ('b', b), ('fa', fa), ('fb', fb)]:
                        pooled[key].append(array)
        pooled = {key: np.concatenate(value) for key, value in pooled.items()}
        pooled_loss = paired_loss(q, targets, pooled['a'], pooled['b'], pooled['fa'], pooled['fb'], anchor, selected)
        for candidate, pname, fname, levels in [('anchor', 'a', 'fa', anchor), ('selected', 'b', 'fb', selected)]:
            prices, forwards = pooled[pname].mean(axis=0), pooled[fname].mean(axis=0)
            total, components, residual = score(q, targets, prices, forwards, levels)
            # Pool iid SPX paths, but use independent scramble means for VIX SE.
            errors = np.empty(len(q))
            for days, group in q[q.instrument == 'SPX'].groupby('expiry_days'):
                fullgroup = spx[spx.expiry_days == days].reset_index(drop=True)
                values = np.concatenate([paths[f'fresh_t{t}_{candidate}_spx_{days}'] for t in range(3)])
                errors[group.index] = values[:, fullgroup.priority_selected.to_numpy(bool)].std(axis=0, ddof=1)/np.sqrt(len(values))
            vi = q.instrument.to_numpy() == 'VIX'
            errors[vi] = pooled[pname][:, vi].std(axis=0, ddof=1)/np.sqrt(24)
            report = q.copy(); report['model_price_U'] = prices; report['mc_se_U'] = errors; report['interval_residual'] = residual
            report['negative_estimate'] = prices < 0; report['nominal_se_screen'] = (prices > 0)&(errors <= report.price_scale_U.to_numpy()/8)
            report.to_csv(run/f'prices_fresh_pooled_{candidate}.csv', index=False)
            ftable = targets[['expiry', 'expiry_days', 'T_years', 'forward_vix_points']].copy()
            ftable['model_forward'] = forwards; ftable['mc_se'] = pooled[fname].std(axis=0, ddof=1)/np.sqrt(24)
            ftable['error_points'] = forwards-ftable.forward_vix_points
            ftable.to_csv(run/f'forwards_fresh_pooled_{candidate}.csv', index=False)
            summaries.append(dict(role='fresh_pooled', candidate=candidate, total=total, components=components,
                negatives=int(sum(prices < 0)), nominal_se_passes=int(report.nominal_se_screen.sum()),
                inside_count=int(sum(residual == 0)), median_se_over_scale=float(np.median(errors/report.price_scale_U.to_numpy()))))
    pd.DataFrame(rows).to_csv(run/'fit_blocks.csv', index=False)
    pd.DataFrame(summaries).to_json(run/'score_audit.json', orient='records', indent=2)
    (run/'paired_loss_audit.json').write_text(json.dumps(dict(per_sample=changes, pooled=pooled_loss), indent=2, allow_nan=False)+'\n')
    result = dict(source_hashes_replayed=len(record['implementation_sha256']), unique_pilot_search_check_seeds=len(seeds),
        searching_and_fresh_quote_scores_replayed=8*882, paired_fresh_loss=pooled_loss,
        scope=dict(spx_priority=707, spx_wing_diagnostics=273, vix_options=175, vix_maturities=[27, 55, 118, 153], vix_forwards=4),
        score_audit=summaries, final_status='Exploratory fit only; stock-grid, normalizer, conditional/time quadrature and rare-tail accuracy unresolved.')
    (run/'audit.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--run', type=Path, required=True)
    audit(parser.parse_args().run)
