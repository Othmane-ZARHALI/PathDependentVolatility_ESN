"""Bounded retained-shock conditional-control versus original vol7/ridge5 check."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd

from .m2_spx import ROOT, M2, frozen_stage3
from .fine_grid import CachedNormalizerPricer
from .spx_hedging_study import refinement, controls
from .spx_path_importance import weighted_controls
from .spx_path_score_controls import ridge_fit
from .spx_volatility_sampling import proposal
from .provenance import capture_sources, source_hashes, save_sources, check_sources
from esn_eq.conditional_controls import retained_path_basis, simulate_conditional_controls
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import dated_grid
from esn_eq.variance_reduction import deterministic_stock_controls
from esn_eq.volatility_path_importance import volatility_path_basis, prepare_volatility_paths


FOLDER = Path(__file__).resolve().parent
PARENT = FOLDER/'runs/m2_spx_central_panel_001'
METHODS = ('reference', 'conditional')


def setup():
    targets = MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    architecture, parameters, premium, hashes = frozen_stage3()
    p = CachedNormalizerPricer(architecture, premium)
    initial, _ = p.initial_state([0., 0., 0.])
    prior = json.loads((PARENT/'run.json').read_text())
    fit = json.loads((FOLDER/'runs/m2_vix_curve_fit_001/run.json').read_text())
    for record in (prior, fit):
        if record['target_sha256'] != targets.artifact_sha256 or record['initial_state'] != initial.tolist():
            raise ValueError('Frozen target/state changed')
        for key, value in dict(architecture=architecture, parameters=parameters, premium=premium).items():
            if record[key] != json.loads(json.dumps(asdict(value))):
                raise ValueError('Frozen model changed')
    prefix = str((M2/'src/esn_eq').relative_to(ROOT))+'/'
    for name, digest in prior['implementation_sha256'].items():
        if name.startswith(prefix) and hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Original normalizer source changed')
    provenance = json.loads((PARENT/'normalizer_prefix_check.json').read_text())
    if (hashlib.sha256((PARENT/'normalizers.npz').read_bytes()).hexdigest() != provenance['new_sha256']
            or prior['normalizer_order'] != 48):
        raise ValueError('Normalizer cache hash/order changed')
    curve = NuisanceCurveFamily.from_targets(targets).build(fit['selected_levels'])
    if list(curve.levels) != prior['curve_levels']:
        raise ValueError('Frozen curve changed')
    scope = json.loads((FOLDER/'active_scope.json').read_text())
    if (scope['spx_research_priority']['criterion']['cutoff'] != 3
            or scope['vix_fit_expiry_days'] != [27, 55, 118, 153]):
        raise ValueError('Active scope changed')
    panel = pd.read_csv(PARENT/'contracts.csv')
    current = targets.select('SPX').set_index(['expiry', 'strike', 'option_type'])
    for row in panel.itertuples():
        actual = current.loc[(row.expiry, row.strike, row.option_type)]
        np.testing.assert_allclose(actual[['T_years', 'forward', 'price_scale_U']].to_numpy(float),
            [row.T_years, row.forward, row.price_scale_U], rtol=0, atol=1e-12)
    events = np.r_[panel.T_years.unique(), 7/365, 21/365, [t for t in curve.ends if t < panel.T_years.max()]]
    times = dated_grid(events, 1/(252*64))
    for _ in range(3):
        times = refinement(times)
    with np.load(PARENT/'normalizers.npz') as archive:
        np.testing.assert_array_equal(times[:-1], archive['times'])
        values = archive['values'].copy()
    values.setflags(write=False)
    return targets, p, parameters, initial, curve, panel, times, values, hashes


def sample(method, p, parameters, initial, curve, quotes, times, basis, n, seed):
    maturity = quotes.T_years.to_numpy()
    strikes, kinds = quotes.K_over_F.to_numpy(), quotes.option_type.to_numpy()
    if method == 'conditional':
        r = simulate_conditional_controls(p, times, initial, parameters, curve,
            maturity, strikes, kinds, n, seed, basis=basis, mixture=proposal(1))
        y, c, w, z, component = r.payoff, r.controls, r.weights, r.scores, r.components
        stock, reference = r.stock, r.reference_stock
    else:
        r = prepare_volatility_paths(p, times, initial, n, seed, basis, proposal(1))
        h = spx_hedging_controls(p, r.cache, parameters, curve, maturity, strikes, kinds)
        stock = np.stack([h.terminal_ratios[t] for t in maturity], axis=1)
        deterministic = deterministic_stock_controls(r.cache, curve, np.unique(maturity))
        reference = np.stack([deterministic[t] for t in maturity], axis=1)
        pairs = controls(dict(stock=stock, hedges=h.hedges), reference, quotes, curve)
        y, c = np.stack([a for a, b in pairs], axis=1), np.stack([b for a, b in pairs], axis=1)
        w, z, component = r.weights, r.score, r.components
    weighted = [weighted_controls(y[:, j], c[:, j], w) for j in range(y.shape[1])]
    return dict(y=np.stack([a for a, b in weighted], axis=1),
        c=np.stack([b for a, b in weighted], axis=1), w=w, scores=z, components=component,
        positive=(y > 0), monitors=np.c_[w, w*stock[:, 0], w*reference[:, 0]])


def fit_pilots(y, c):
    indices = (np.arange(len(y)), np.arange(len(y)//2), np.arange(len(y)//2, len(y)))
    return np.stack([np.stack([ridge_fit(y[ix, j], c[ix, j]) for j in range(y.shape[1])]) for ix in indices])


def adjusted(data, coefficients, forwards):
    return (data['y'][:, :, None]-np.einsum('njk,vjk->njv', data['c'], coefficients))*forwards[None, :, None]


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    begin = time.perf_counter()
    helpers = ['m2_spx', 'fine_grid', 'spx_hedging_study', 'spx_path_importance',
        'spx_path_score_controls', 'spx_volatility_sampling']
    sources = [Path(__file__).resolve(), FOLDER/'active_scope.json',
        PARENT/'run.json', PARENT/'normalizers.npz', PARENT/'normalizer_prefix_check.json', PARENT/'contracts.csv',
        FOLDER/'runs/m2_vix_curve_fit_001/run.json',
        *[FOLDER/(name+'.py') for name in helpers], *sorted((M2/'src/esn_eq').glob('*.py'))]
    captured = capture_sources(sources)
    save_sources(args.output, captured, ROOT)
    targets, p, parameters, initial, curve, panel, master, normalizers, hashes = setup()
    panel.to_csv(args.output/'contracts.csv', index=False)
    design = dict(arguments=vars(args)|dict(output=str(args.output)), architecture=asdict(p.architecture),
        parameters=asdict(parameters), premium=asdict(p.premium), initial_state=initial.tolist(),
        curve_levels=list(curve.levels), target_sha256=targets.artifact_sha256, frozen_model_sha256=hashes,
        implementation_sha256=source_hashes(captured, ROOT), normalizer_order=48, steps_per_trading_day=512,
        software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__),
        primary='Independent pilot ridge5, penalty=1; same five controls, analytically conditioned for new method.',
        proposal='Original full-path vol7 versus NEW retained-coordinate four-score defensive mixture with same shift sizes/masses. Different proposal laws; compare whole estimators, not pure conditioning gain.',
        uncertainty_unit='iid checking paths conditional on independent frozen pilot; methods use independent streams; pilot-half comparisons paired within method.',
        fixed_budget='256 pilot + two 256-path fresh checks per method/expiry by default; no adaptive selection or extension.',
        adoption='Require both fresh trials to improve scale-weighted whole-panel variance times measured runtime, without worse tail/pilot warnings. Passing is still not price convergence.',
        limitations=['24 representatives of 707 SPX cutoff-3 priority quotes; no scope reduction; all four VIX maturities remain active.',
            'No timestep or normalizer refinement in this run; exact OU transitions do not make stock update exact.',
            'Negative controlled estimates retained; no clipping or weight self-normalization.',
            'Conditional positive-tail counts include unrepresentably tiny values and cannot certify adequate tail exposure.',
            'Different proposal laws mean cross-method differences use independent combined SE, not paired SE.',
            'Primary costs are standalone checking pipelines with warmed common GH48; setup and pilot costs also retained.'])
    (args.output/'design.json').write_text(json.dumps(design, indent=2, allow_nan=False)+'\n')
    all_rows, costs, diagnostics, records = [], [], [], []
    for days, quotes in panel.groupby('expiry_days', sort=True):
        quotes = quotes.reset_index(drop=True)
        out = args.output/f'expiry_{int(days)}'
        out.mkdir()
        quotes.to_csv(out/'contracts.csv', index=False)
        times = master[:np.searchsorted(master, quotes.T_years.iloc[0])+1]
        values = normalizers[:len(times)-1]
        values.setflags(write=False)
        p.normalizers[(times[:-1].tobytes(), initial.tobytes(), parameters, 48)] = values
        bases, basis_seconds = {}, {}
        for method, builder in [('reference', volatility_path_basis), ('conditional', retained_path_basis)]:
            tick = time.perf_counter()
            bases[method] = builder(p, times, 7/365)
            basis_seconds[method] = time.perf_counter()-tick
        archive, seed_records, coefficients = {}, [], {}

        def collect(method, role, trial):
            count = args.pilot if role == 0 else args.paths
            pieces, seeds = {}, []
            for batch, start in enumerate(range(0, count, args.batch_size)):
                n = min(args.batch_size, count-start)
                seed = int(np.random.SeedSequence([args.seed, int(days), METHODS.index(method), role, trial, batch]).generate_state(1)[0])
                seeds.append(seed)
                (args.output/'status.json').write_text(json.dumps(dict(expiry_days=int(days), method=method, role=role, trial=trial, paths=start))+'\n')
                tick = time.perf_counter()
                result = sample(method, p, parameters, initial, curve, quotes, times, bases[method], n, seed)
                costs.append(dict(expiry_days=int(days), method=method, role=role, trial=trial, batch=batch, paths=n, seconds=time.perf_counter()-tick))
                for key, value in result.items():
                    pieces.setdefault(key, []).append(value)
                print(f'{days}d {method} role={role} trial={trial} paths={start+n}/{count}', flush=True)
            seed_records.append(dict(method=method, role=role, trial=trial, seeds=seeds))
            return {key: np.concatenate(value) for key, value in pieces.items()}

        # Complete and freeze both independent pilots before any checking draw.
        pilot_archive = {}
        for method in METHODS:
            data = collect(method, 0, 0)
            coefficients[method] = fit_pilots(data['y'], data['c'])
            for key, value in data.items():
                pilot_archive[f'{method}_{key}'] = value
        frozen = dict(coefficients={k: v.tolist() for k, v in coefficients.items()},
            variants=['full', 'half_A', 'half_B'], checking_not_used=True, pilot_seeds=seed_records.copy())
        (out/'frozen.json').write_text(json.dumps(frozen, indent=2, allow_nan=False)+'\n')
        np.savez_compressed(out/'pilot_draws.npz', **pilot_archive)
        for trial in range(args.trials):
            # Alternate ordering to reduce systematic warmup/thermal timing bias.
            for method in (METHODS if trial % 2 == 0 else METHODS[::-1]):
                data = collect(method, 1, trial)
                x = adjusted(data, coefficients[method], quotes.forward.to_numpy())
                for key, value in data.items():
                    archive[f'{method}_t{trial}_{key}'] = value
                archive[f'{method}_t{trial}_adjusted'] = x
                for j, row in enumerate(quotes.itertuples()):
                    for variant, name in enumerate(['full', 'half_A', 'half_B']):
                        a = x[:, j, variant]
                        squares = (a-a.mean())**2
                        all_rows.append(dict(expiry_days=int(days), trial=trial, method=method,
                            variant=name, contract=j, strike=row.strike, option_type=row.option_type, band=row.band,
                            estimate=float(a.mean()), se=float(a.std(ddof=1)/np.sqrt(len(a))),
                            variance=float(a.var(ddof=1)), scale=row.price_scale_U,
                            positive_paths=int(data['positive'][:, j].sum()),
                            largest_residual_share=float(squares.max()/squares.sum()) if squares.sum() else 0.))
                    for control in range(5):
                        a = data['c'][:, j, control]
                        diagnostics.append(dict(expiry_days=int(days), trial=trial, method=method,
                            target=f'contract{j}_control{control}', mean=float(a.mean()),
                            se=float(a.std(ddof=1)/np.sqrt(len(a))), expectation=0.))
                for j, label in enumerate(['mass', 'stock', 'reference_stock']):
                    a = data['monitors'][:, j]
                    diagnostics.append(dict(expiry_days=int(days), trial=trial, method=method,
                        target=label, mean=float(a.mean()), se=float(a.std(ddof=1)/np.sqrt(len(a))), expectation=1.))
        np.savez_compressed(out/'checking_draws.npz', **archive)
        record = dict(expiry_days=int(days), grid_intervals=len(times)-1, basis_seconds=basis_seconds, seed_records=seed_records)
        (out/'run.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
        records.append(record)
    check_sources(captured)
    for name, rows in [('estimates', all_rows), ('costs', costs), ('diagnostics', diagnostics)]:
        pd.DataFrame(rows).to_csv(args.output/(name+'.csv'), index=False)
    result = design|dict(status='completed bounded conditional-control check; price accuracy remains unqualified',
        expiry_records=records, wall_seconds=time.perf_counter()-begin)
    (args.output/'run.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    (args.output/'status.json').write_text('{"status":"completed"}\n')
    print(f'Completed in {result["wall_seconds"]:.1f}s', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--pilot', type=int, default=256)
    parser.add_argument('--paths', type=int, default=256)
    parser.add_argument('--batch-size', type=int, default=128)
    parser.add_argument('--trials', type=int, default=2)
    parser.add_argument('--seed', type=int, default=2026102301)
    args = parser.parse_args()
    if (args.output.exists() or min(args.pilot, args.paths, args.batch_size) < 8
            or args.pilot % 2 or args.trials != 2
            or any(n % args.batch_size == 1 for n in [args.pilot, args.paths])):
        parser.error('Fresh output, even pilot, valid batch sizes and exactly two trials required')
    run(args)
