"""Full-scope VIX Sobol versus defensive-mixture importance-sampling pilot."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from .m2_spx import ROOT,M2,frozen_stage3
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.outer_sampling import GaussianOuterSampler,replication_statistics
from esn_eq.pricing import DatedM2Pricer
from esn_eq.vix_pricing import ConditionalVIXPricer


def run(args):
    start = time.perf_counter()
    scope_path = Path(__file__).with_name('active_scope.json')
    old_path = Path(__file__).with_name('runs')/'m2_vix_curve_fit_001/run.json'
    scope = json.loads(scope_path.read_text())
    if scope['vix_fit_expiry_days'] != [27,55,118,153] or scope['vix_deferred_expiry_days']:
        raise ValueError('This analysis requires the explicitly restored four-maturity scope')
    sources = capture_sources([Path(__file__),scope_path,old_path,Path(__file__).with_name('m2_spx.py'),
        *sorted((M2/'src/esn_eq').glob('*.py'))])
    targets = MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    architecture,parameters,premium,frozen_hashes = frozen_stage3()
    dated = DatedM2Pricer(architecture,premium)
    initial,_ = dated.initial_state([0.,0.,0.])
    old = json.loads(old_path.read_text())
    expected = json.loads(json.dumps(dict(architecture=asdict(architecture),parameters=asdict(parameters),
        premium=asdict(premium),initial_state=initial.tolist(),target_sha256=targets.artifact_sha256)))
    if any(old[k] != value for k,value in expected.items()):
        raise ValueError('Frozen candidate model/state/targets changed')
    curve = NuisanceCurveFamily.from_targets(targets).build(old['selected_levels'])
    engine = ConditionalVIXPricer(dated,parameters,curve,initial,normalizer_order=48)
    rows,raw,diagnostics,seeds = [],[],[],[]
    for days in scope['vix_fit_expiry_days']:
        quotes = targets.select('VIX',expiry_days=[days]).reset_index(drop=True)
        expiry = float(quotes.T_years.iloc[0])
        sampler = GaussianOuterSampler(dated,expiry,initial)
        mixture = DefensiveGaussianMixture.symmetric_leading(len(initial))
        window = engine.gaussian_window(expiry,time_order=8)
        strikes = quotes.strike.to_numpy()
        sign = np.where(quotes.option_type.to_numpy() == 'call',1.,-1.)
        descriptors = [dict(metric='forward',strike=None,option_type=None,scale=.5),
            dict(metric='second_moment',strike=None,option_type=None,scale=None)]
        descriptors += [dict(metric='option',strike=float(q.strike),option_type=q.option_type,
            scale=float(q.price_scale_U)) for q in quotes.itertuples()]
        for method_id,method in enumerate(('sobol','mixture')):
            estimates = {power:[] for power in args.powers}
            counts = {power:[] for power in args.powers}
            concentration = {power:[] for power in args.powers}
            for replicate in range(args.replications):
                tick = time.perf_counter()
                seed = int(np.random.SeedSequence([args.seed,days,method_id,replicate]).generate_state(1)[0])
                if method == 'mixture':
                    z,w,component = mixture.draw(max(args.powers),seed)
                    states = sampler.mean+z @ sampler.root.T
                else:
                    states = sampler.draw(max(args.powers),seed)
                    w,component = np.ones(len(states)),np.zeros(len(states),int)
                values = engine.evaluate(states,window,order=32).vix
                payoffs = np.maximum(sign*(values[:,None]-strikes),0.)
                contributions = w[:,None]*np.column_stack((values,values**2,payoffs))
                for power in args.powers:
                    n = 2**power
                    estimates[power].append(contributions[:n].mean(axis=0))
                    counts[power].append((payoffs[:n]>0).sum(axis=0))
                    sums = contributions[:n].sum(axis=0)
                    concentration[power].append(np.divide(contributions[:n].max(axis=0),sums,
                        out=np.full(len(sums),np.nan),where=sums>0))
                    diagnostics.append(dict(days=days,method=method,replicate=replicate,power=power,
                        weight_mean=float(w[:n].mean()),weight_min=float(w[:n].min()),weight_max=float(w[:n].max()),
                        weight_ess=float(w[:n].sum()**2/np.sum(w[:n]**2))))
                raw.append(pd.DataFrame(dict(days=days,method=method,replicate=replicate,
                    point=np.arange(len(values)),vix=values,weight=w,component=component)))
                seeds.append(dict(days=days,method=method,replicate=replicate,seed=seed))
                print(json.dumps(dict(days=days,method=method,replicate=replicate+1,
                    seconds=time.perf_counter()-tick)),flush=True)
            for power in args.powers:
                matrix = np.array(estimates[power])
                mean,se = replication_statistics(matrix)
                for j,desc in enumerate(descriptors):
                    c = np.array(counts[power])[:,j-2] if j>=2 else None
                    shares = np.array(concentration[power])[:,j]
                    rows.append(dict(days=days,method=method,power=power,**desc,
                        estimate=float(mean[j]),replication_se=float(se[j]),
                        se_over_scale=float(se[j]/desc['scale']) if desc['scale'] else None,
                        positive_payoff_count=int(c.sum()) if c is not None else None,
                        zero_hit_replications=int(np.sum(c==0)) if c is not None else None,
                        max_single_contribution_share=float(np.nanmax(shares)) if np.isfinite(shares).any() else None))
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    record = dict(status='fixed-candidate sampling diagnostic, no refit or accuracy acceptance',
        arguments=vars(args)|dict(output=str(args.output)),scope=scope,
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),
        initial_state=initial.tolist(),curve_levels=old['selected_levels'],
        target_sha256=targets.artifact_sha256,frozen_model_sha256=frozen_hashes,
        implementation_sha256=source_hashes(sources,ROOT),seed_records=seeds,
        proposal='0.5 unshifted + 0.125 each +/-2 in first two full PCA coordinates; no tuning',
        likelihood='exact phi(z)/sum pi_k phi(z-m_k); no sample normalization or forward rescaling',
        weight_upper_bound=2.,conditional_order=32,normalizer_order=48,time_order=8,
        uncertainty_unit='independent complete Sobol scrambles, including proposal-component selection',
        limitations=['Eight-scramble standard errors and variance comparisons are exploratory.',
            'Zero hits or zero observed SE do not certify rare-event prices.',
            'Weight ESS is diagnostic, not a payoff accuracy certificate.',
            'Outer error excludes fixed conditional, time and normalizer quadrature biases.',
            'The two leading PCA directions need not cover the important nonlinear payoff tails.',
            'Existing market input conflicts remain; no held-out quotes used.'],
        wall_seconds=time.perf_counter()-start)
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(args.output/'weight_diagnostics.csv',index=False)
    pd.concat(raw,ignore_index=True).to_csv(args.output/'outer_values.csv',index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--powers',type=int,nargs=2,default=[5,6])
    parser.add_argument('--replications',type=int,default=8)
    parser.add_argument('--seed',type=int,default=20261026)
    args = parser.parse_args()
    if args.output.exists() or args.replications<4 or not 3<=args.powers[0]<args.powers[1]:
        parser.error('Require new output, >=4 scrambles and increasing powers >=3')
    run(args)
