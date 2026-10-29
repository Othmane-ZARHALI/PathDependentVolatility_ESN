"""Matched-path SPX static versus adapted hedging controls on three grids."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from .m2_spx import ROOT,M2,frozen_stage3
from .fine_grid import CachedNormalizerPricer
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.market_targets import MarketTargets
from esn_eq.pricing import dated_grid,coarsen_paths
from esn_eq.variance_reduction import deterministic_stock_controls,_control_data


def choose_quotes(targets):
    rows = []
    for _,group in targets.select("SPX",expiry_days=[29,57]).groupby("expiry"):
        for target,kind in ((.5,"put"),(.8,"put"),(.98,"put"),(1.02,"call"),(1.2,"call"),(1.5,"call")):
            used = {(int(row.expiry),float(row.strike),row.option_type) for row in rows}
            eligible = group[(group.option_type == kind) & ~group.apply(
                lambda row:(int(row.expiry),float(row.strike),row.option_type) in used,axis=1)]
            rows.append(eligible.loc[(eligible.K_over_F-target).abs().idxmin()])
    frame = pd.DataFrame(rows).sort_values(["T_years","strike"]).reset_index(drop=True)
    if frame.duplicated(["expiry","strike","option_type"]).any():
        raise ValueError("Representative targets collided")
    return frame


def refinement(grid):
    return np.r_[grid[0],np.concatenate([np.linspace(a,b,3)[1:] for a,b in zip(grid[:-1],grid[1:])])]


def samples(pricer,parameters,curve,initial,quotes,grids,args,trial,role):
    paths = args.pilot if role == 0 else args.paths
    used = [len(grids)-1] if role == 0 else list(range(len(grids)))
    pieces = {level:dict(stock=[],hedges=[]) for level in used}
    references,seeds = [],[]
    maturities = quotes.T_years.to_numpy()
    strikes,kinds = quotes.K_over_F.to_numpy(),quotes.option_type.to_numpy()
    for batch,start in enumerate(range(0,paths,args.batch_size)):
        count = min(args.batch_size,paths-start)
        seed = int(np.random.SeedSequence([args.seed,trial,role,batch]).generate_state(1)[0])
        seeds.append(seed)
        cache = pricer.prepare(grids[-1],initial,count,seed)
        reference = deterministic_stock_controls(cache,curve,np.unique(maturities))
        references.append(np.stack([reference[t] for t in maturities],axis=1))
        for level in used:
            paths_cache = cache if level == len(grids)-1 else coarsen_paths(cache,grids[level])
            result = spx_hedging_controls(pricer,paths_cache,parameters,curve,maturities,strikes,kinds,order=48)
            pieces[level]["stock"].append(np.stack([result.terminal_ratios[t] for t in maturities],axis=1))
            pieces[level]["hedges"].append(result.hedges)
            del paths_cache,result
        del cache
        print(f"trial={trial} role={role} paths={start+count}/{paths}",flush=True)
    return {level:{key:np.concatenate(value) for key,value in data.items()} for level,data in pieces.items()},np.concatenate(references),seeds


def controls(data,reference,quotes,curve):
    output = []
    for j,row in enumerate(quotes.itertuples()):
        payoff,base = next(_control_data(data["stock"][:,j],reference[:,j],[row.K_over_F],
                                        [row.option_type],curve.integral(0.,row.T_years)))
        output.append((payoff,np.c_[base,data["hedges"][:,j,:]]))
    return output


def beta_fit(payoff,control):
    return np.linalg.lstsq(control-control.mean(axis=0),payoff-payoff.mean(),rcond=1e-10)[0]


def run(args):
    sources = capture_sources([Path(__file__),Path(__file__).with_name("m2_spx.py"),Path(__file__).with_name("fine_grid.py"),
                               *sorted((M2/"src/esn_eq").glob("*.py"))])
    start = time.perf_counter()
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    quotes = choose_quotes(targets)
    architecture,parameters,premium,hashes = frozen_stage3()
    pricer = CachedNormalizerPricer(architecture,premium)
    initial,_ = pricer.initial_state([0.,0.,0.])
    fit_path = Path(__file__).with_name("runs")/"m2_vix_curve_fit_001/run.json"
    fitted = json.loads(fit_path.read_text())
    if fitted["target_sha256"] != targets.artifact_sha256:
        raise ValueError("Target snapshot changed")
    expected = json.loads(json.dumps(dict(architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium))))
    if any(fitted[k] != v for k,v in expected.items()) or fitted["initial_state"] != initial.tolist():
        raise ValueError("Fixed model changed")
    curve = NuisanceCurveFamily.from_targets(targets).build(fitted["selected_levels"])
    events = np.r_[quotes.T_years.unique(),[t for t in curve.ends if t < quotes.T_years.max()]]
    coarse = dated_grid(events,1/(252*64))
    grids = [coarse,refinement(coarse),refinement(refinement(coarse))]
    estimates,paired,coefficients,seeds = [],[],[],[]
    for trial in range(args.trials):
        pilot,reference,pilot_seeds = samples(pricer,parameters,curve,initial,quotes,grids,args,trial,0)
        fitted_controls = []
        for j,(payoff,control) in enumerate(controls(pilot[2],reference,quotes,curve)):
            static,augmented = beta_fit(payoff,control[:,:2]),beta_fit(payoff,control)
            fitted_controls.append((static,augmented))
            coefficients.append(dict(trial=trial,contract=j,static=static.tolist(),augmented=augmented.tolist()))
        del pilot,reference
        evaluation,reference,evaluation_seeds = samples(pricer,parameters,curve,initial,quotes,grids,args,trial,1)
        seeds.append(dict(trial=trial,pilot=pilot_seeds,evaluation=evaluation_seeds))
        adjusted = {}
        for level,data in evaluation.items():
            adjusted[level] = []
            for j,((payoff,control),(static,augmented),row) in enumerate(zip(controls(data,reference,quotes,curve),fitted_controls,quotes.itertuples())):
                draws = np.stack([payoff,payoff-control[:,:2] @ static,payoff-control @ augmented],axis=1)*row.forward
                adjusted[level].append(draws)
                for count in (args.paths//4,args.paths//2,args.paths):
                    x = draws[:count]
                    for method,name in enumerate(("plain","static","hedging")):
                        estimates.append(dict(trial=trial,contract=j,expiry_days=row.expiry_days,strike=row.strike,
                            option_type=row.option_type,K_over_F=row.K_over_F,grid_steps=64*2**level,paths=count,method=name,
                            price=float(x[:,method].mean()),se=float(x[:,method].std(ddof=1)/np.sqrt(count)),
                            scale=row.price_scale_U,positive_payoff_count=int(np.sum(payoff[:count] > 0))))
        for low,high in ((0,1),(1,2)):
            for j,row in enumerate(quotes.itertuples()):
                difference = adjusted[high][j]-adjusted[low][j]
                for method,name in enumerate(("plain","static","hedging")):
                    paired.append(dict(trial=trial,contract=j,expiry_days=row.expiry_days,strike=row.strike,option_type=row.option_type,
                        coarse_steps=64*2**low,fine_steps=64*2**high,method=name,
                        change=float(difference[:,method].mean()),se=float(difference[:,method].std(ddof=1)/np.sqrt(args.paths)),
                        scale=row.price_scale_U))
    record = dict(status="representative-strike matched-path variance-reduction diagnostic; no full-surface accuracy claim",
        arguments=vars(args)|dict(output=str(args.output)),architecture=asdict(architecture),parameters=asdict(parameters),
        premium=asdict(premium),initial_state=initial.tolist(),curve_levels=list(curve.levels),
        target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,implementation_sha256=source_hashes(sources,ROOT),
        grid_intervals=[len(g)-1 for g in grids],max_lambda_dt=[float(pricer.layout.rates.max()*np.diff(g).max()) for g in grids],
        seed_records=seeds,coefficients=coefficients,normalizer_order=48,
        uncertainty_unit="iid evaluation paths conditional on independent pilot coefficients",
        wall_seconds=time.perf_counter()-start,
        limitations=["Twelve preselected actual training contracts; this does not qualify all 648 near SPX quotes.",
            "Three time grids share one finest state/equity path sample in each trial; grid differences are paired.",
            "All regression coefficients come from independent pilot paths at the finest grid and stay fixed across grids.",
            "Control means are zero for the discrete stock scheme; continuous-model time bias and normalizer error remain.",
            "Nested sample-count prefixes are dependent; SE changes are descriptive, not a convergence-order proof.",
            "Zero-hit tail expectations remain unresolved even with variance-reduced SEs."])
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    (args.output/"run.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    pd.DataFrame(estimates).to_csv(args.output/"estimates.csv",index=False)
    pd.DataFrame(paired).to_csv(args.output/"paired_grid_changes.csv",index=False)
    quotes.to_csv(args.output/"contracts.csv",index=False)
    print(f"Completed in {record['wall_seconds']:.1f}s",flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--paths",type=int,default=2048)
    parser.add_argument("--pilot",type=int,default=1024)
    parser.add_argument("--batch-size",type=int,default=256)
    parser.add_argument("--trials",type=int,default=2)
    parser.add_argument("--seed",type=int,default=20261024)
    args = parser.parse_args()
    if args.output.exists() or min(args.paths,args.pilot,args.batch_size) < 8 or args.paths % 4 or args.trials < 1:
        parser.error("Require a new output and valid budgets")
    if any(n % args.batch_size == 1 for n in (args.paths,args.pilot)):
        parser.error("Avoid a final batch with one path")
    run(args)
