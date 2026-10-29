"""Matched-cost iid versus independently scrambled Sobol VIX sampling pilot."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd
import scipy

from .m2_spx import ROOT, M2, frozen_stage3
from .provenance import capture_sources, source_hashes, check_sources, save_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.outer_sampling import GaussianOuterSampler, replication_statistics
from esn_eq.pricing import DatedM2Pricer
from esn_eq.vix_pricing import ConditionalVIXPricer


def metrics(values, strikes, kinds):
    sign = np.where(kinds == "call", 1., -1.)
    payoffs = np.maximum(sign*(values[:, None]-strikes), 0.)
    return np.r_[values.mean(), np.mean(values**2), payoffs.mean(axis=0)], (payoffs > 0).sum(axis=0)


def run(args):
    sources = capture_sources([Path(__file__), Path(__file__).with_name("m2_spx.py"),
                               *sorted((M2/"src/esn_eq").glob("*.py"))])
    start = time.perf_counter()
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    quotes = targets.select("VIX", expiry_days=[args.expiry_days]).reset_index(drop=True)
    if quotes.empty:
        raise ValueError("No selected training maturity")
    architecture, parameters, premium, hashes = frozen_stage3()
    dated = DatedM2Pricer(architecture, premium)
    initial, _ = dated.initial_state([0., 0., 0.])
    curve = targets.curve("training_mid_piecewise_v1")
    maturity = float(quotes.T_years.iloc[0])
    sampler = GaussianOuterSampler(dated, maturity, initial)
    engine = ConditionalVIXPricer(dated, parameters, curve, initial, normalizer_order=48)
    window = engine.gaussian_window(maturity, time_order=8)
    strikes, kinds = quotes.strike.to_numpy(), quotes.option_type.to_numpy()
    descriptors = [dict(metric="forward",strike=None,option_type=None,scale=.5),
                   dict(metric="second_moment",strike=None,option_type=None,scale=None)]
    descriptors += [dict(metric="option",strike=float(row.strike),option_type=row.option_type,
                         scale=float(row.price_scale_U)) for row in quotes.itertuples()]
    estimates, hits, timings, seed_records, raw = {}, {}, [], [], []
    refinement_states = None
    for method_id, method in enumerate(("iid", "sobol")):
        for power in args.powers:
            estimates[method,power], hits[method,power] = [], []
        for replicate in range(args.replications):
            tick = time.perf_counter()
            seed = int(np.random.SeedSequence([args.seed,args.expiry_days,method_id,replicate]).generate_state(1)[0])
            states = sampler.draw(max(args.powers), seed, method)
            values = engine.evaluate(states, window, order=32).vix
            for power in args.powers:
                means, counts = metrics(values[:2**power], strikes, kinds)
                estimates[method,power].append(means)
                hits[method,power].append(counts)
            raw.append(pd.DataFrame(dict(method=method,replicate=replicate,seed=seed,
                point=np.arange(len(values)),vix=values)))
            if method == "sobol" and replicate == 0:
                order = np.argsort(values)
                selected = order[np.linspace(0,len(order)-1,8).astype(int)]
                refinement_states = states[selected]
            elapsed = time.perf_counter()-tick
            timings.append(dict(method=method,replicate=replicate,wall_seconds=elapsed))
            seed_records.append(dict(method=method,replicate=replicate,seed=seed))
            print(json.dumps(dict(method=method,replicate=replicate+1,
                replications=args.replications,points=len(values),wall_seconds=elapsed)), flush=True)
    rows, rep_rows, changes = [], [], []
    for method in ("iid", "sobol"):
        for power in args.powers:
            matrix = np.array(estimates[method,power])
            mean, se = replication_statistics(matrix)
            counts = np.array(hits[method,power])
            for metric, desc in enumerate(descriptors):
                row = dict(method=method,power=power,points_per_replication=2**power,
                    replications=args.replications,**desc,estimate=float(mean[metric]),
                    replication_se=float(se[metric]),
                    zero_replication_estimates=int(np.sum(matrix[:,metric] == 0.)))
                if metric >= 2:
                    row["positive_payoff_total"] = int(counts[:,metric-2].sum())
                    row["replications_with_positive_payoff"] = int(np.sum(counts[:,metric-2] > 0))
                rows.append(row)
                for replicate in range(args.replications):
                    rep_rows.append(dict(method=method,power=power,replicate=replicate,**desc,
                        estimate=float(matrix[replicate,metric]),
                        positive_payoff_count=int(counts[replicate,metric-2]) if metric >= 2 else None))
        coarse, fine = min(args.powers), max(args.powers)
        difference = np.array(estimates[method,fine])-np.array(estimates[method,coarse])
        mean, se = replication_statistics(difference)
        for metric, desc in enumerate(descriptors):
            changes.append(dict(method=method,coarse_power=coarse,fine_power=fine,**desc,
                paired_replication_change=float(mean[metric]),paired_replication_se=float(se[metric])))
    summary = pd.DataFrame(rows)
    ratios = []
    for power in args.powers:
        a, b = np.array(estimates["iid",power]), np.array(estimates["sobol",power])
        for metric, desc in enumerate(descriptors):
            va, vb = float(np.var(a[:,metric],ddof=1)), float(np.var(b[:,metric],ddof=1))
            ratios.append(dict(power=power,**desc,iid_replication_variance=va,
                sobol_replication_variance=vb,iid_over_sobol_variance=va/vb if vb > 0 else None))
    # Isolate conditional GH and global normalizer changes at eight selected states.
    base = engine.evaluate(refinement_states, window, order=32).vix
    conditional = engine.evaluate(refinement_states, window, order=48).vix
    refined_engine = ConditionalVIXPricer(dated, parameters, curve, initial, normalizer_order=64)
    normalizer = refined_engine.evaluate(refinement_states,
        refined_engine.gaussian_window(maturity,time_order=8),order=48).vix
    refinement = dict(states="eight ordered VIX quantile positions from first Sobol replication; diagnostic selection",
        base_vix=base.tolist(),conditional_order_48_vix=conditional.tolist(),
        normalizer_order_64_vix=normalizer.tolist(),
        max_conditional_change=float(np.max(np.abs(conditional-base))),
        max_normalizer_change=float(np.max(np.abs(normalizer-conditional))))
    record = dict(status="unfitted matched-cost outer sampling pilot; no price acceptance or fit ranking",
        arguments=vars(args)|dict(output=str(args.output)),architecture=asdict(architecture),
        parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),
        curve_id="training_mid_piecewise_v1",expiry_days=args.expiry_days,maturity=maturity,
        clock="ACT/365, delta=30/365; exact terminal OU distribution",
        input_compatibility=targets.audit["variance_audit"],target_sha256=targets.artifact_sha256,
        frozen_model_sha256=hashes,implementation_sha256=source_hashes(sources, ROOT),
        covariance_factor="all PCA columns retained, ordered by decayed projection variance contribution",
        direction_lags=sampler.direction_lags.tolist(),coordinate_scores=sampler.coordinate_scores.tolist(),
        normalizer_order=48,conditional_order=32,time_order=8,time_nodes=len(window.times),
        uncertainty_unit="independent complete iid batches / independent Sobol scrambles",
        points_within_scramble_are_iid=False,seed_records=seed_records,timings=timings,
        second_moment_target=engine.moment_target(window),quadrature_refinement=refinement,
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        wall_seconds=time.perf_counter()-start,
        limitations=["Eight-replication variance ratios are exploratory and noisy, especially in tails.",
            "Zero payoff estimates/SEs do not certify rare-event expectations.",
            "Outer uncertainty excludes conditional GH, time integration and normalizer approximation.",
            "Quadrature checks cover selected states, not the full outer law or all parameter candidates.",
            "PCA ordering may perform differently at other parameters and expiries; no general advantage is assumed."])
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    (args.output/"run.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    summary.to_csv(args.output/"summary.csv",index=False)
    pd.DataFrame(rep_rows).to_csv(args.output/"replications.csv",index=False)
    pd.DataFrame(ratios).to_csv(args.output/"variance_ratios.csv",index=False)
    pd.DataFrame(changes).to_csv(args.output/"refinement.csv",index=False)
    pd.concat(raw,ignore_index=True).to_csv(args.output/"outer_values.csv",index=False)
    print(summary[summary.metric == "forward"].to_json(orient="records",indent=2),flush=True)
    print(json.dumps(refinement),flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--expiry-days",type=int,default=55)
    parser.add_argument("--powers",type=int,nargs=2,default=[6,7])
    parser.add_argument("--replications",type=int,default=8)
    parser.add_argument("--seed",type=int,default=20261022)
    args = parser.parse_args()
    if args.output.exists() or args.replications < 4 or not 3 <= args.powers[0] < args.powers[1]:
        parser.error("Require a new output, >=4 replications and increasing powers >=3")
    run(args)
