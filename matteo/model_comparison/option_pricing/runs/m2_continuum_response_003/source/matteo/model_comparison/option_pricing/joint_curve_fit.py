"""Bounded actual joint SPX/VIX nuisance-curve optimization, then fresh checks."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import tempfile
import time

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import least_squares

from .m2_spx import ROOT, M2
from .spx_conditional_control_study import setup
from .spx_volatility_sampling import proposal
from .spx_path_score_controls import ridge_fit
from .joint_curve_core import CurvePathBatch, path_payoffs, donor_indices, ExploratoryJointLoss
from .provenance import capture_sources, source_hashes, save_sources, check_sources
from esn_eq.volatility_path_importance import volatility_path_basis, prepare_volatility_paths
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.variance_reduction import deterministic_stock_controls
from esn_eq.conditional_spx import lognormal_payoffs
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.gaussian_controls import gaussian_control_rows
from esn_eq.vix_curve_profile import prepare_curve_exposures
from esn_eq.vix_pricing import ConditionalVIXPricer
from esn_eq.market_calibration import NuisanceCurveFamily


FOLDER = Path(__file__).resolve().parent
VIX_PARENT = FOLDER/'runs/m2_vix_gaussian_stability_001'
MEMBERSHIP = FOLDER/'runs/spx_scaled_moneyness_003/training_membership.csv'


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


class SPXData:
    def __init__(self, quotes, batches, weights, controls, anchor_stock):
        self.quotes, self.batches, self.weights, self.controls = quotes, batches, weights, controls
        self.anchor_stock = anchor_stock
        self.beta = None

    def price(self, levels, derivative=True, terminal=None):
        stocks, gradients = [], []
        for batch in self.batches:
            stock, gradient = batch.terminal(levels, derivative)
            stocks.append(stock); gradients.append(gradient)
        stock, gradient = np.concatenate(stocks), np.concatenate(gradients)
        # Comparison against the original pricer was made when each cache built.
        h, dh = path_payoffs(stock, gradient, self.quotes)
        values = self.weights[:, None]*h-np.einsum('njk,jk->nj', self.controls, self.beta)
        jac = np.mean(self.weights[:, None, None]*dh, axis=0)
        if terminal is not None:
            terminal.update(stock=stock, values=values)
        return values.mean(axis=0), jac, values


class VIXData:
    def __init__(self, quotes, coefficients, floor_integral, delta, weights, controls, beta, replications):
        self.quotes, self.exposures = quotes, coefficients
        self.floor_integral, self.delta, self.weights = floor_integral, delta, weights
        self.controls, self.beta, self.replications = controls, beta, replications

    def price(self, levels):
        integral = self.floor_integral+self.exposures@(np.asarray(levels)-1e-4)
        vix = 100*np.sqrt(integral/self.delta)
        dv = 5000*self.exposures/(self.delta*vix[:, None])
        sign = np.where(self.quotes.option_type.to_numpy() == 'call', 1., -1.)
        distance = sign*(vix[:, None]-self.quotes.strike.to_numpy())
        raw = np.c_[vix, np.maximum(distance, 0.)]
        derivative = np.concatenate([dv[:, None, :], (distance > 0)[:, :, None]*sign[None, :, None]*dv[:, None, :]], axis=1)
        values = self.weights[:, None]*raw-self.controls@self.beta
        jac = np.mean(self.weights[:, None, None]*derivative, axis=0)
        blocks = values.reshape(self.replications, -1, values.shape[1]).mean(axis=1)
        return values.mean(axis=0), jac, blocks, vix


class FitData:
    def __init__(self, spx, vix, spx_quotes, vix_quotes, forwards, primary):
        self.spx, self.vix = spx, vix
        self.spx_quotes, self.vix_quotes, self.forwards = spx_quotes, vix_quotes, forwards
        self.all_quotes = pd.concat([spx_quotes, vix_quotes], ignore_index=True)
        self.primary = np.r_[primary, np.ones(len(vix_quotes), bool)]
        self.loss = ExploratoryJointLoss(self.all_quotes[self.primary], forwards)

    def predict(self, levels, derivative=True):
        prices, jac = np.empty(len(self.all_quotes)), np.empty((len(self.all_quotes), 3))
        errors = np.empty(len(prices)); f = np.empty(len(self.forwards)); fj = np.empty((len(f), 3)); fe = np.empty(len(f))
        blocks = np.empty((8, len(prices))); forward_blocks = np.empty((8, len(f)))
        archive, tail = {}, []
        for days, data in self.spx.items():
            mean, derivative_values, values = data.price(levels, derivative)
            ix = self.spx_quotes.index[self.spx_quotes.expiry_days == days].to_numpy()
            prices[ix], jac[ix], errors[ix] = mean, derivative_values, values.std(axis=0, ddof=1)/np.sqrt(len(values))
            blocks[:, ix] = values.reshape(8, -1, len(ix)).mean(axis=1)
            archive[f'spx_{days}'] = values
            for j, q in enumerate(data.quotes.itertuples()):
                squares = (values[:, j]-mean[j])**2
                tail.append(dict(instrument='SPX', expiry_days=days, strike=q.strike, option_type=q.option_type,
                    largest_centered_square_share=float(squares.max()/squares.sum()) if squares.sum() else 0.))
        for days, data in self.vix.items():
            mean, derivative_values, b, v = data.price(levels)
            ix = len(self.spx_quotes)+self.vix_quotes.index[self.vix_quotes.expiry_days == days].to_numpy()
            k = int(np.flatnonzero(self.forwards.expiry_days.to_numpy() == days)[0])
            prices[ix], jac[ix], errors[ix] = mean[1:], derivative_values[1:], b[:, 1:].std(axis=0, ddof=1)/np.sqrt(8)
            f[k], fj[k], fe[k] = mean[0], derivative_values[0], b[:, 0].std(ddof=1)/np.sqrt(8)
            blocks[:, ix], forward_blocks[:, k] = b[:, 1:], b[:, 0]
            archive[f'vix_{days}'] = v
            # Concentration across independent complete scrambles, not iid points.
            for j, q in enumerate(data.quotes.itertuples()):
                squares = (b[:, j+1]-mean[j+1])**2
                tail.append(dict(instrument='VIX', expiry_days=days, strike=q.strike, option_type=q.option_type,
                    largest_centered_square_share=float(squares.max()/squares.sum()) if squares.sum() else 0.))
        result = self.loss.evaluate(prices[self.primary], f, levels, jac[self.primary], fj)
        return result|dict(prices=prices, price_jac=jac, forwards=f, forward_jac=fj, errors=errors,
            forward_errors=fe, blocks=blocks, forward_blocks=forward_blocks, archive=archive, tail=tail)


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    begin = time.perf_counter()
    helpers = ['m2_spx', 'fine_grid', 'spx_conditional_control_study', 'spx_hedging_study',
        'spx_path_importance', 'spx_path_score_controls', 'spx_volatility_sampling', 'joint_curve_core']
    sources = [Path(__file__).resolve(), MEMBERSHIP, FOLDER/'active_scope.json',
        FOLDER/'runs/m2_spx_central_panel_001/run.json', FOLDER/'runs/m2_spx_central_panel_001/normalizers.npz',
        FOLDER/'runs/m2_spx_central_panel_001/normalizer_prefix_check.json', FOLDER/'runs/m2_spx_central_panel_001/contracts.csv',
        FOLDER/'runs/m2_vix_curve_fit_001/run.json', *[FOLDER/(n+'.py') for n in helpers],
        *[VIX_PARENT/n for n in ['run.json', 'proposals.json', 'coefficients.json']], *sorted((M2/'src/esn_eq').glob('*.py'))]
    captured = capture_sources(sources); save_sources(args.output, captured, ROOT)
    targets, p, parameters, initial, anchor_curve, panel, master, normalizers, model_hashes = setup()
    family = NuisanceCurveFamily.from_targets(targets); anchor = np.array(anchor_curve.levels)
    spx = targets.select('SPX').reset_index(drop=True)
    vix = targets.select('VIX').reset_index(drop=True)
    forwards = targets.vix_forwards.sort_values('expiry').reset_index(drop=True)
    membership = pd.read_csv(MEMBERSHIP).set_index(['expiry', 'strike', 'option_type'])
    keys = pd.MultiIndex.from_frame(spx[['expiry', 'strike', 'option_type']])
    primary = membership.loc[keys, 'priority_selected'].to_numpy(bool)
    assert primary.sum() == 707 and len(spx) == 980 and len(vix) == 175
    for quotes, name in [(spx, 'spx'), (vix, 'vix')]:
        quotes.to_csv(args.output/(name+'_contracts.csv'), index=False)
    spx.assign(priority_selected=primary).to_csv(args.output/'spx_membership.csv', index=False)
    vix_record = json.loads((VIX_PARENT/'run.json').read_text())
    if vix_record['target_sha256'] != targets.artifact_sha256 or vix_record['curve_levels'] != anchor.tolist():
        raise ValueError('VIX reference differs from fixed anchor')
    for key, value in [('architecture', asdict(p.architecture)), ('parameters', asdict(parameters)),
            ('premium', asdict(p.premium)), ('initial_state', initial.tolist())]:
        if vix_record[key] != json.loads(json.dumps(value)):
            raise ValueError('VIX reference model/state mismatch')
    for name, key in [('coefficients.json', 'coefficient_sha256'), ('proposals.json', 'proposal_sha256')]:
        if hashlib.sha256((VIX_PARENT/name).read_bytes()).hexdigest() != vix_record[key]:
            raise ValueError('VIX frozen controls/proposal changed')
    proposals = json.loads((VIX_PARENT/'proposals.json').read_text())
    vix_betas = json.loads((VIX_PARENT/'coefficients.json').read_text())
    scratch = Path(tempfile.mkdtemp(prefix='m2-joint-curve-', dir='/private/tmp'))
    design = dict(arguments=vars(args)|dict(output=str(args.output)), architecture=asdict(p.architecture),
        parameters=asdict(parameters), premium=asdict(p.premium), initial_state=initial.tolist(),
        family=asdict(family), anchor_levels=anchor.tolist(), free_curve_levels=[0, 1, 2],
        target_sha256=targets.artifact_sha256, frozen_model_sha256=model_hashes,
        implementation_sha256=source_hashes(captured, ROOT), scratch=str(scratch),
        objective='707 cutoff-3 SPX training options +175 VIX training options +4 forwards. Other 273 SPX training wings are diagnostics only.',
        spx_estimator='Original per-expiry vol7 law; frozen anchor stock/reference/two hedge controls plus mass, ridge5 penalty1. Nearest same-side representative hedge amortized across all strikes. Not a claim of same variance as per-strike production controls.',
        vix_estimator='Unchanged targeted mixture and pooled28 Gaussian control coefficients; fresh exposures; anchor coefficients frozen across curves.',
        signed_estimates='Same interval-distance formula extended to finite signed Monte Carlo estimates. Negatives retained and flagged. Production JointPriceObjective remains strict and unchanged.',
        budgets=dict(spx_pilot=256, spx_search=512, spx_fresh=512, vix_replications=8, vix_points_per_replication=64,
            stock_steps_per_trading_day=512, conditional_GH=32, normalizer_GH=48, time_order=8,
            starts=8, max_nfev_per_start=args.max_nfev, searching_wall_minutes=60, fresh_trials=3),
        limitations=['Coarse exploratory fit below4096 screening budget; no qualified market fit, physical feasibility or model ranking.',
            'Fixed state/readout/Q and piecewise curve restrict the family; failure is not M2 infeasibility.',
            'SE excludes stock discretization, normalizer and VIX quadrature/time bias, which remain separately unresolved.',
            'Fresh seeds check numerics on training contracts, not held-out market performance.',
            'Analytical cache derivatives ignore measure changes because state/readout/Q are fixed; this cache cannot be reused after changing them.',
            'Signed-score extension and fixed anchor shared controls are explicit exploratory implementation choices; negative prices never become qualified quotes.'])
    dump(args.output/'design.json', design)
    seed_records, costs, beta_spx, bases = [], [], {}, {}
    beta_vix = {days: np.asarray(vix_betas[str(days)]['coefficients']['gaussian'])[:, np.r_[0, np.arange(len(vix[vix.expiry_days == days]))+2]]
        for days in [27, 55, 118, 153]}
    for days in [29, 57, 211, 302]:
        t = float(spx[spx.expiry_days == days].T_years.iloc[0])
        grid = master[:np.searchsorted(master, t)+1]
        value = normalizers[:len(grid)-1]; value.setflags(write=False)
        p.normalizers[(grid[:-1].tobytes(), initial.tobytes(), parameters, 48)] = value
        bases[days] = volatility_path_basis(p, grid, 7/365)

    def build_spx(role, trial, keep_paths=True):
        result = {}
        count = 256 if role == 0 else 512
        for days in [29, 57, 211, 302]:
            q = spx[spx.expiry_days == days].reset_index(drop=True)
            representatives = panel[panel.expiry_days == days].reset_index(drop=True)
            donor = donor_indices(q, representatives)
            grid = master[:np.searchsorted(master, q.T_years.iloc[0])+1]
            ratios, weights, cs, ys = [], [], [], []
            for batch, start in enumerate(range(0, count, 128)):
                dump(args.output/'status.json', dict(stage='SPX paths/anchor controls', role=role, trial=trial, days=days, paths=start))
                tick = time.perf_counter(); seed = int(np.random.SeedSequence([args.seed, 1, role, trial, days, batch]).generate_state(1)[0])
                r = prepare_volatility_paths(p, grid, initial, 128, seed, bases[days], proposal(1))
                cache, w = r.cache, r.weights
                h = spx_hedging_controls(p, cache, parameters, anchor_curve, representatives.T_years.to_numpy(),
                    representatives.K_over_F.to_numpy(), representatives.option_type.to_numpy())
                stock = h.terminal_ratios[float(q.T_years.iloc[0])]
                ref = deterministic_stock_controls(cache, anchor_curve, [float(q.T_years.iloc[0])])[float(q.T_years.iloc[0])]
                total = anchor_curve.integral(0., float(q.T_years.iloc[0]))
                known = lognormal_payoffs(-.5*total, total, q.K_over_F.to_numpy(), q.option_type.to_numpy())
                sign = np.where(q.option_type.to_numpy() == 'call', 1., -1.)
                reference_control = np.maximum(sign*(ref[:, None]-q.K_over_F.to_numpy()), 0.)-known
                control = np.concatenate([np.broadcast_to((stock-1)[:, None, None], (128, len(q), 1)),
                    reference_control[:, :, None], h.hedges[:, donor]], axis=2)
                control = np.concatenate([w[:, None, None]*control, np.broadcast_to((w-1)[:, None, None], (128, len(q), 1))], axis=2)
                if keep_paths:
                    stem = scratch/f'r{role}_t{trial}_d{days}_b{batch}'
                    a = np.lib.format.open_memmap(str(stem)+'_ratio.npy', mode='w+', dtype='float64', shape=cache.equity_increments.shape)
                    normalizer = normalizers[:len(grid)-1]
                    for k in range(0, len(a), 2048):
                        a[k:k+2048] = p.multiplier(cache.projections[k:k+2048], parameters)/normalizer[k:k+2048, None]
                    a.flush(); np.save(str(stem)+'_equity.npy', cache.equity_increments)
                    path = CurvePathBatch(np.load(str(stem)+'_ratio.npy', mmap_mode='r'), np.load(str(stem)+'_equity.npy', mmap_mode='r'),
                        np.diff(grid), np.minimum(np.searchsorted(anchor_curve.ends, grid[:-1], side='right'), 2))
                    oracle, _ = path.terminal(anchor)
                    np.testing.assert_allclose(oracle, stock, rtol=2e-11, atol=2e-12)
                    ratios.append(path)
                weights.append(w); cs.append(control); ys.append(stock)
                seed_records.append(dict(instrument='SPX', role=role, trial=trial, days=days, batch=batch, seed=seed))
                costs.append(dict(instrument='SPX', role=role, trial=trial, days=days, batch=batch, seconds=time.perf_counter()-tick))
                print(f'SPX role={role} trial={trial} {days}d {start+128}/{count}', flush=True)
                del cache, r, h
            data = SPXData(q, ratios, np.concatenate(weights), np.concatenate(cs), np.concatenate(ys))
            if role == 0:
                sign = np.where(q.option_type.to_numpy() == 'call', 1., -1.)
                target = data.weights[:, None]*np.maximum(sign*(data.anchor_stock[:, None]-q.K_over_F.to_numpy()), 0.)
                beta_spx[days] = np.stack([ridge_fit(target[:, j], data.controls[:, j]) for j in range(len(q))])*q.forward.to_numpy()[:, None]
            data.beta = beta_spx[days]
            np.savez_compressed(args.output/f'spx_r{role}_t{trial}_d{days}.npz', controls=data.controls, weights=data.weights,
                beta=data.beta, anchor_stock=data.anchor_stock, donor=donor)
            result[days] = data
        return result

    engine = ConditionalVIXPricer(p, parameters, anchor_curve, initial, normalizer_order=48)
    windows = {days: engine.gaussian_window(float(vix[vix.expiry_days == days].T_years.iloc[0]), time_order=8) for days in [27, 55, 118, 153]}

    def build_vix(role, trial):
        result = {}
        for days in [27, 55, 118, 153]:
            q = vix[vix.expiry_days == days].reset_index(drop=True)
            sampler = GaussianOuterSampler(p, float(q.T_years.iloc[0]), initial)
            mix = DefensiveGaussianMixture(proposals[str(days)]['shifts'], proposals[str(days)]['probabilities'])
            U = np.asarray(vix_betas[str(days)]['basis'])
            states, weights, controls = [], [], []
            for rep in range(8):
                seed = int(np.random.SeedSequence([args.seed, 2, role, trial, days, rep]).generate_state(1)[0])
                z, w, _ = mix.draw(6, seed)
                states.append(sampler.mean+z@sampler.root.T); weights.append(w); controls.append(gaussian_control_rows(z, w, U))
                seed_records.append(dict(instrument='VIX', role=role, trial=trial, days=days, replicate=rep, seed=seed))
            states = np.concatenate(states); tick = time.perf_counter()
            dump(args.output/'status.json', dict(stage='VIX exposures', role=role, trial=trial, days=days))
            cache = prepare_curve_exposures(engine, states, windows[days], order=32)
            coefficients = cache.coefficients
            direct = engine.evaluate(states[:3], windows[days], order=32).vix
            np.testing.assert_allclose(cache.evaluate(anchor).vix[:3], direct, rtol=1e-12, atol=1e-11)
            data = VIXData(q, coefficients, cache.floor_integral, cache.delta, np.concatenate(weights), np.concatenate(controls), beta_vix[days], 8)
            np.savez_compressed(args.output/f'vix_r{role}_t{trial}_d{days}.npz', exposures=coefficients, floor_integral=cache.floor_integral,
                delta=cache.delta, weights=data.weights, controls=data.controls, beta=data.beta, states=states)
            costs.append(dict(instrument='VIX', role=role, trial=trial, days=days, seconds=time.perf_counter()-tick))
            print(f'VIX role={role} trial={trial} {days}d exposures512', flush=True)
            result[days] = data
        return result

    build_spx(0, 0, keep_paths=False)
    dump(args.output/'frozen_spx_controls.json', {str(k): v.tolist() for k, v in beta_spx.items()})
    train = FitData(build_spx(1, 0), build_vix(1, 0), spx, vix, forwards, primary)
    train.loss.quotes.assign(loss_weight=train.loss.weights).to_csv(args.output/'objective_membership.csv', index=False)
    evaluations, ledger, candidates = [], [], []
    memo, search_begin = {}, time.perf_counter()

    def evaluate(x):
        key = np.asarray(x).tobytes()
        if key not in memo:
            if time.perf_counter()-search_begin > 3600:
                raise TimeoutError('Declared searching wall budget exhausted')
            value = train.predict(x)
            memo.clear(); memo[key] = value
            row = dict(evaluation=len(evaluations), levels=np.asarray(x).tolist(), total=value['total'],
                components=value['components'], negative_count=value['negative_count'])
            evaluations.append(row)
            with (args.output/'evaluations.jsonl').open('a') as file:
                file.write(json.dumps(row, allow_nan=False)+'\n')
            if len(evaluations) % 5 == 1:
                print(f'fit evaluation={len(evaluations)} loss={value["total"]:.6f} levels={np.asarray(x)} negatives={value["negative_count"]}', flush=True)
        return memo[key]

    lo = np.asarray(family.lower); hi = np.asarray(family.upper)
    starts = [anchor, [.03, .04, .06], [.04, .06, .08], [.06, .05, .06],
        [.08, .08, .08], [.10, .10, .10], [.02, .10, .04], [.10, .02, .10]]
    for number, start in enumerate(starts):
        x = np.clip(start, lo+1e-8, hi-1e-8)
        try:
            fit = least_squares(lambda y: evaluate(y)['residuals'], x, jac=lambda y: evaluate(y)['jacobian'],
                bounds=(lo, hi), max_nfev=args.max_nfev, ftol=1e-6, xtol=1e-6, gtol=1e-5)
            value = evaluate(fit.x)
            ledger.append(dict(start=number, initial=x.tolist(), levels=fit.x.tolist(), total=value['total'],
                success=bool(fit.success), status=int(fit.status), nfev=int(fit.nfev), message=fit.message,
                negative_count=value['negative_count']))
            candidates.append((value['total'], fit.x.copy()))
            print(f'Completed start{number}: {value["total"]:.6f}, success={fit.success}', flush=True)
        except TimeoutError as error:
            ledger.append(dict(start=number, initial=x.tolist(), status='interrupted', message=str(error)))
            break
        dump(args.output/'optimizer_ledger.json', ledger)
    if not candidates:
        raise RuntimeError('No optimizer candidate completed; retained sources/status/evaluation log')
    # Include anchor among choices; selection cannot silently make train loss worse.
    anchor_result = train.predict(anchor)
    candidates.append((anchor_result['total'], anchor))
    _, selected = min(candidates, key=lambda item: item[0])
    fit_seconds = time.perf_counter()-search_begin
    dump(args.output/'selected.json', dict(selected_levels=selected.tolist(), anchor_levels=anchor.tolist(),
        selected_using='searching simulations only', fresh_checking_not_used=True))
    dump(args.output/'optimizer_ledger.json', ledger)
    records, block_archive, price_archive = [], {}, {}

    def record(role, trial, data):
        results = {}
        for name, levels in [('anchor', anchor), ('selected', selected)]:
            r = data.predict(levels, derivative=False); results[name] = r
            q, f = data.loss.tables(r['prices'][data.primary], r['errors'][data.primary], r['forwards'], r['forward_errors'])
            q.to_csv(args.output/f'prices_{role}_t{trial}_{name}.csv', index=False)
            f.to_csv(args.output/f'forwards_{role}_t{trial}_{name}.csv', index=False)
            fullq = data.all_quotes.copy(); fullq['model_price_U'] = r['prices']; fullq['mc_se_U'] = r['errors']
            fullq['priority_objective_member'] = data.primary
            fullq.to_csv(args.output/f'full_training_diagnostics_{role}_t{trial}_{name}.csv', index=False)
            pd.DataFrame(r['tail']).to_csv(args.output/f'tails_{role}_t{trial}_{name}.csv', index=False)
            block_archive[f'{role}_t{trial}_{name}_prices'] = r['blocks']
            block_archive[f'{role}_t{trial}_{name}_forwards'] = r['forward_blocks']
            for key, array in r['archive'].items():
                price_archive[f'{role}_t{trial}_{name}_{key}'] = array
            records.append(dict(role=role, trial=trial, candidate=name, levels=levels.tolist(), total=r['total'],
                components=r['components'], negative_count=r['negative_count'], invalid_forward_count=r['invalid_forward_count']))
        return results

    record('search', 0, train)
    # Fresh numerical checks are generated only after selected.json exists.
    for trial in range(3):
        fresh = FitData(build_spx(2, trial), build_vix(2, trial), spx, vix, forwards, primary)
        record('fresh', trial, fresh)
        del fresh
    check_sources(captured)
    np.savez_compressed(args.output/'comparison_blocks.npz', **block_archive)
    np.savez_compressed(args.output/'comparison_paths.npz', **price_archive)
    pd.DataFrame(costs).to_csv(args.output/'costs.csv', index=False)
    dump(args.output/'scores.json', records)
    result = design|dict(status='completed coarse joint curve optimization and three fresh numerical checks; unresolved price accuracy',
        selected_levels=selected.tolist(), optimizer_ledger=ledger, evaluations=len(evaluations),
        optimizer_seconds=fit_seconds, seed_records=seed_records, scores=records,
        wall_seconds=time.perf_counter()-begin,
        software=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__))
    dump(args.output/'run.json', result); dump(args.output/'status.json', dict(status='completed'))
    print(json.dumps(dict(selected_levels=selected.tolist(), scores=records, wall_seconds=result['wall_seconds'])), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-nfev', type=int, default=25)
    parser.add_argument('--seed', type=int, default=2026102302)
    args = parser.parse_args()
    if args.output.exists() or args.max_nfev < 2:
        parser.error('Fresh output and valid optimizer budget required')
    run(args)
