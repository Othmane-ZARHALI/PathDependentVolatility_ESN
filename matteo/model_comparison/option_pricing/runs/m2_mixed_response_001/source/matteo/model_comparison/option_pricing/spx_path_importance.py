"""Bounded frozen-SPX path importance pilot; separate sampling/grid/GH changes."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd

from .m2_spx import ROOT, M2, frozen_stage3
from .fine_grid import CachedNormalizerPricer
from .spx_hedging_study import choose_quotes, refinement, controls, beta_fit
from .provenance import capture_sources, source_hashes, check_sources, save_sources
from esn_eq.path_importance import prepare_importance_paths
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import dated_grid, coarsen_paths
from esn_eq.variance_reduction import deterministic_stock_controls


def weighted_controls(payoff, control, weights):
    """Target-zero controls must be weighted too; add known-zero mass control."""
    return weights*payoff, np.c_[weights[:, None]*control, weights-1.]


def collect(p, parameters, curve, initial, quotes, grids, args, trial, method, role):
    count = args.pilot if role == 0 else args.paths
    variants = [(2, 48)] if role == 0 else [(0, 48), (1, 48), (2, 48), (2, 64)]
    collected = {v: [] for v in variants}
    monitors, seeds = [], []
    maturities = quotes.T_years.to_numpy()
    for batch, start in enumerate(range(0, count, args.batch_size)):
        n = min(args.batch_size, count-start)
        seed = int(np.random.SeedSequence([args.seed, trial, method, role, batch]).generate_state(1)[0])
        seeds.append(seed)
        if method:
            proposal = prepare_importance_paths(p, grids[-1], initial, n, seed, args.magnitude)
            cache, weights = proposal.cache, proposal.weights
        else:
            cache, weights = p.prepare(grids[-1], initial, n, seed), np.ones(n)
        ref = deterministic_stock_controls(cache, curve, np.unique(maturities))
        reference = np.stack([ref[t] for t in maturities], axis=1)
        monitor = [weights]
        for level, order in variants:
            c = cache if level == 2 else coarsen_paths(cache, grids[level])
            hedged = spx_hedging_controls(p, c, parameters, curve, maturities,
                        quotes.K_over_F.to_numpy(), quotes.option_type.to_numpy(), order=order)
            data = dict(stock=np.stack([hedged.terminal_ratios[t] for t in maturities], axis=1), hedges=hedged.hedges)
            pairs = controls(data, reference, quotes, curve)
            collected[(level, order)].append((np.stack([x[0] for x in pairs], axis=1),
                                             np.stack([x[1] for x in pairs], axis=1), weights))
            if level == 2 and order == 48:
                monitor += [weights*data['stock'][:, j] for j in (0, len(quotes)-1)]
                monitor += [weights*reference[:, j] for j in (0, len(quotes)-1)]
        monitors.append(np.stack(monitor, axis=1))
        print(f'trial={trial} method={method} role={role} paths={start+n}/{count}', flush=True)
    return {v: tuple(np.concatenate([part[j] for part in parts]) for j in range(3))
            for v, parts in collected.items()}, np.concatenate(monitors), seeds


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    fit_path = Path(__file__).resolve().with_name('runs')/'m2_vix_curve_fit_001/run.json'
    sources = capture_sources([Path(__file__).resolve(), fit_path,
        *[Path(__file__).resolve().with_name(n+'.py') for n in ('m2_spx','fine_grid','spx_hedging_study')],
        *sorted((M2/'src/esn_eq').glob('*.py'))])
    save_sources(args.output, sources, ROOT)
    (args.output/'status.json').write_text(json.dumps(dict(status='running', arguments=vars(args)|dict(output=str(args.output))), indent=2))
    targets = MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    quotes = choose_quotes(targets)
    architecture, parameters, premium, hashes = frozen_stage3()
    p = CachedNormalizerPricer(architecture, premium)
    initial, _ = p.initial_state([0.,0.,0.])
    fitted = json.loads(fit_path.read_text())
    if fitted['target_sha256'] != targets.artifact_sha256 or fitted['initial_state'] != initial.tolist():
        raise ValueError('Frozen inputs changed')
    for k, v in dict(architecture=architecture, parameters=parameters, premium=premium).items():
        if fitted[k] != json.loads(json.dumps(asdict(v))):
            raise ValueError('Frozen model changed')
    curve = NuisanceCurveFamily.from_targets(targets).build(fitted['selected_levels'])
    events = np.r_[quotes.T_years.unique(), [t for t in curve.ends if t < quotes.T_years.max()]]
    g = dated_grid(events, 1/(252*64))
    grids = [g, refinement(g), refinement(refinement(g))]
    rows, changes, diagnostics, coefficients, seed_records, draws = [], [], [], [], [], {}
    for trial in range(args.trials):
        for method in (0, 1):
            name = ('target','importance')[method]
            pilot, _, ps = collect(p, parameters, curve, initial, quotes, grids, args, trial, method, 0)
            py, pc, pw = pilot[(2,48)]
            betas = []
            for j in range(len(quotes)):
                y, c = weighted_controls(py[:, j], pc[:, j], pw)
                b = beta_fit(y, c)
                betas.append(b)
                coefficients.append(dict(trial=trial, method=name, contract=j, beta=b.tolist()))
            del pilot
            evaluation, monitor, es = collect(p, parameters, curve, initial, quotes, grids, args, trial, method, 1)
            seed_records.append(dict(trial=trial, method=name, pilot=ps, evaluation=es))
            draws[f't{trial}_{name}_monitors'] = monitor
            for j, label in enumerate(('mass','stock29','stock57','reference29','reference57')):
                diagnostics.append(dict(trial=trial, method=name, statistic=label, estimate=float(monitor[:,j].mean()),
                    se=float(monitor[:,j].std(ddof=1)/np.sqrt(args.paths)), target=1.))
            adjusted = {}
            for (level, order), (payoffs, cs, w) in evaluation.items():
                vals = []
                for j, row in enumerate(quotes.itertuples()):
                    y, c = weighted_controls(payoffs[:,j], cs[:,j], w)
                    x = np.c_[y, y-c@betas[j]]*row.forward
                    vals.append(x)
                    for k, estimator in enumerate(('raw','hedging')):
                        rows.append(dict(trial=trial, method=name, estimator=estimator, contract=j,
                            expiry_days=row.expiry_days, strike=row.strike, option_type=row.option_type,
                            steps=64*2**level, normalizer_order=order, price=float(x[:,k].mean()),
                            se=float(x[:,k].std(ddof=1)/np.sqrt(args.paths)), scale=row.price_scale_U,
                            positive_payoffs=int(np.count_nonzero(payoffs[:,j])), paths=args.paths))
                adjusted[(level,order)] = np.stack(vals, axis=1)
                draws[f't{trial}_{name}_g{level}_q{order}'] = adjusted[(level,order)]
            for low, high, label in (((0,48),(1,48),'64_to_128'),((1,48),(2,48),'128_to_256'),((2,48),(2,64),'GH48_to_64')):
                diff = adjusted[high]-adjusted[low]
                for j, row in enumerate(quotes.itertuples()):
                    for k, estimator in enumerate(('raw','hedging')):
                        changes.append(dict(trial=trial, method=name, estimator=estimator, contract=j, change_type=label,
                            change=float(diff[:,j,k].mean()), se=float(diff[:,j,k].std(ddof=1)/np.sqrt(args.paths)),
                            correction_variance=float(diff[:,j,k].var(ddof=1)), scale=row.price_scale_U))
    check_sources(sources)
    record = dict(status='completed bounded diagnostic; no convergence or calibration clearance', arguments=vars(args)|dict(output=str(args.output)),
        architecture=asdict(architecture), parameters=asdict(parameters), premium=asdict(premium), initial_state=initial.tolist(),
        curve_levels=list(curve.levels), target_sha256=targets.artifact_sha256, frozen_model_sha256=hashes,
        implementation_sha256=source_hashes(sources, ROOT), seed_records=seed_records, coefficients=coefficients,
        grid_intervals=[len(g)-1 for g in grids], max_lambda_dt=[float(p.layout.rates.max()*np.diff(g).max()) for g in grids],
        uncertainty_unit='iid checking paths conditional on independent pilot coefficients; paired refinements within method',
        proposal='0.5 target + 0.25 Brownian drift +m/sqrt(H) + 0.25 drift -m/sqrt(H); correlated OU shifts; full-H exact nonnormalized weights',
        limitations=['12 training contracts only; all maturity/split scopes unchanged.',
            'Independent target/importance draws: variance comparisons are descriptive, not paired price tests.',
            'One pilot per method/trial; no checking-set coefficient selection.',
            'GH48/64 and 128/256 changes are diagnostics, not population bias bounds.',
            'Zero hits or zero empirical SE do not establish tail price accuracy.',
            'No MLMC allocation/cost optimization; correction variances only.'], wall_seconds=time.perf_counter()-start)
    for name, data in (('estimates',rows),('paired_changes',changes),('diagnostics',diagnostics)):
        pd.DataFrame(data).to_csv(args.output/(name+'.csv'), index=False)
    quotes.to_csv(args.output/'contracts.csv', index=False)
    np.savez_compressed(args.output/'checking_draws.npz', **draws)
    (args.output/'run.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    (args.output/'status.json').write_text(json.dumps(dict(status='completed'))+'\n')
    print(f'Completed in {record["wall_seconds"]:.1f}s', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--paths', type=int, default=1024)
    parser.add_argument('--pilot', type=int, default=512)
    parser.add_argument('--batch-size', type=int, default=256)
    parser.add_argument('--trials', type=int, default=2)
    parser.add_argument('--seed', type=int, default=20261025)
    parser.add_argument('--magnitude', type=float, default=2.)
    args = parser.parse_args()
    if args.output.exists() or min(args.paths,args.pilot,args.batch_size) < 8 or args.trials < 1 or not np.isfinite(args.magnitude) or args.magnitude <= 0:
        parser.error('Require new output, positive magnitude, valid budgets')
    if any(n % args.batch_size == 1 for n in (args.paths,args.pilot)):
        parser.error('Avoid one-path final batch')
    run(args)
