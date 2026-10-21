"""Active 27/55-day VIX pilot with one identifiable near variance level."""
import argparse
from dataclasses import asdict,replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import least_squares

from .m2_spx import ROOT,M2,frozen_stage3
from .vix_curve_fit import CurveFitData
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.calibration_loss import JointPriceObjective,LossConfig
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.market_targets import MarketTargets
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.pricing import DatedM2Pricer
from esn_eq.vix_curve_profile import VIXCurveExposures,prepare_curve_exposures
from esn_eq.vix_pricing import ConditionalVIXPricer


def active_targets(targets,scope):
    included,deferred = scope["vix_fit_expiry_days"],scope["vix_deferred_expiry_days"]
    if set(included) != {27,55} or set(deferred) != {118,153} or set(included)&set(deferred):
        raise ValueError("This one-level pilot requires the declared 27/55-day scope")
    quotes = targets.select("VIX",expiry_days=included).reset_index(drop=True)
    forwards = targets.vix_forwards[targets.vix_forwards.expiry_days.isin(included)].sort_values("expiry").reset_index(drop=True)
    if set(quotes.expiry_days) != set(included) or set(forwards.expiry_days) != set(included):
        raise ValueError("Missing active expiry/forward")
    return quotes,forwards


def prefix_cache(cache,replications,points,prefix):
    if not 1 <= prefix <= points:
        raise ValueError("Invalid scramble prefix")
    coefficients = cache.coefficients.reshape(replications,points,-1)[:,:prefix,:].reshape(-1,len(cache.ends)).copy()
    coefficients.setflags(write=False)
    return replace(cache,coefficients=coefficients)


def run(args):
    start = time.perf_counter()
    scope_path = Path(__file__).with_name("active_scope.json")
    scope = json.loads(scope_path.read_text())
    sources = capture_sources([Path(__file__),scope_path,Path(__file__).with_name("vix_curve_fit.py"),
        Path(__file__).with_name("m2_spx.py"),*sorted((M2/"src/esn_eq").glob("*.py"))])
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    quotes,forwards = active_targets(targets,scope)
    family = NuisanceCurveFamily.from_targets(targets)
    architecture,parameters,premium,frozen_hashes = frozen_stage3()
    dated = DatedM2Pricer(architecture,premium)
    initial,_ = dated.initial_state([0.,0.,0.])
    engine = ConditionalVIXPricer(dated,parameters,family.build([.04]*3),initial,normalizer_order=48)
    old_folder = Path(__file__).with_name("runs")/"m2_vix_curve_fit_001"
    old = json.loads((old_folder/"run.json").read_text())
    expected = json.loads(json.dumps(dict(architecture=asdict(architecture),parameters=asdict(parameters),
        premium=asdict(premium),initial_state=initial.tolist(),target_sha256=targets.artifact_sha256)))
    if any(old[key] != value for key,value in expected.items()):
        raise ValueError("Training cache model/state/targets changed")
    if (old["normalizer_order"],old["conditional_order"],old["time_order"]) != (48,32,8):
        raise ValueError("Training cache quadrature convention changed")
    arrays = np.load(old_folder/"exposures.npz")
    train_caches,windows,stored = {},{},{}
    for target in forwards.itertuples():
        expiry = int(target.expiry)
        window = engine.gaussian_window(float(target.T_years),time_order=8)
        if target.T_years+window.delta >= family.ends[0]:
            raise ValueError("More than one curve level affects an active VIX window")
        coefficients = arrays[f"role0_expiry{expiry}"].copy()
        if np.any(coefficients[:,1:] != 0):
            raise ValueError("Cached inactive segments are nonzero")
        coefficients.setflags(write=False)
        train_caches[expiry] = VIXCurveExposures(coefficients,old["exposure_manifest"][str(expiry)]["floor_integral"],
            family.ends,window.delta,32,engine,engine._configuration())
        windows[expiry] = window
        stored[f"train_expiry{expiry}"] = coefficients
    train = CurveFitData(quotes,forwards,train_caches,old["arguments"]["replications"],2**old["arguments"]["power"])
    objective = JointPriceObjective(quotes,forwards,LossConfig(spx_weight=0.))
    def score(data,x):
        levels = [float(x)]*3  # Inactive placeholders; not identified long-run coefficients.
        prices,means = data.predictions(levels)
        return objective.evaluate(prices,means,curve_penalty=family.penalty(levels))
    ledger = []
    for initial_x in (.01,.02,.03,.04,.05,.06,.08,.10):
        result = least_squares(lambda x:score(train,x[0])["residuals"],[initial_x],
            bounds=([family.lower[0]],[family.upper[0]]),max_nfev=300,ftol=1e-10,xtol=1e-10,gtol=1e-10)
        ledger.append(dict(start=initial_x,x=float(result.x[0]),success=bool(result.success),
            total=score(train,result.x[0])["total"],nfev=int(result.nfev)))
    converged = [row for row in ledger if row["success"]]
    if not converged:
        raise RuntimeError("No starting point converged")
    selected = min(converged,key=lambda row:row["total"])["x"]
    candidates = dict(neutral_20pct_vol=.04,previous_start=float(old["starting_levels"][0]),
                      previous_four_maturity_fit=float(old["selected_levels"][0]),near_fit=selected)
    print(json.dumps(dict(near_variance=selected,training_loss=score(train,selected)["total"])),flush=True)
    # Independent states generated only after the one-dimensional selection is frozen.
    fresh_caches,fresh_states,seeds,timings = {},{},{},[]
    points = 2**args.power
    for target in forwards.itertuples():
        tick = time.perf_counter()
        expiry = int(target.expiry)
        sampler = GaussianOuterSampler(dated,float(target.T_years),initial)
        seeds[expiry] = [int(np.random.SeedSequence([args.seed,expiry,replicate,2]).generate_state(1)[0])
                         for replicate in range(args.replications)]
        states = np.concatenate([sampler.draw(args.power,seed) for seed in seeds[expiry]])
        cache = prepare_curve_exposures(engine,states,windows[expiry],order=32)
        fresh_caches[expiry],fresh_states[expiry] = cache,states
        stored[f"fresh_expiry{expiry}"] = cache.coefficients
        elapsed = time.perf_counter()-tick
        timings.append(dict(expiry=expiry,points=len(states),seconds=elapsed))
        print(json.dumps(timings[-1]),flush=True)
    datasets = {"train":train}
    for prefix in (points//2,points):
        caches = {expiry:prefix_cache(cache,args.replications,points,prefix) for expiry,cache in fresh_caches.items()}
        datasets[f"fresh_{prefix}"] = CurveFitData(quotes,forwards,caches,args.replications,prefix)
    scores,tables = [],{}
    for label,data in datasets.items():
        for candidate,x in candidates.items():
            result = score(data,x)
            scores.append(dict(sample=label,candidate=candidate,variance=x,total=result["total"],
                inside_fraction=float(np.mean(result["quotes"].interval_residual == 0)),**result["components"]))
        result = score(data,selected)
        diagnostic = data.predictions([selected]*3,details=True)
        q,f = result["quotes"],result["forwards"]
        q["replication_se_U"] = diagnostic["price_se"]
        q["positive_payoff_count"] = diagnostic["positive_payoff_count"].astype(int)
        q["positive_payoff_replications"] = diagnostic["positive_payoff_replications"].astype(int)
        f["replication_se"] = diagnostic["forward_se"]
        tables[f"prices_{label}.csv"],tables[f"forwards_{label}.csv"] = q,f
    refinements = []
    direct = ConditionalVIXPricer(dated,parameters,family.build([selected]*3),initial,normalizer_order=48)
    high = ConditionalVIXPricer(dated,parameters,family.build([selected]*3),initial,normalizer_order=64)
    for target in forwards.itertuples():
        expiry = int(target.expiry)
        base = fresh_caches[expiry].evaluate([selected]*3).vix
        indices = np.argsort(base)[np.linspace(0,len(base)-1,8).astype(int)]
        states = fresh_states[expiry][indices]
        conditional = direct.evaluate(states,direct.gaussian_window(target.T_years),order=48).vix
        normalizer = high.evaluate(states,high.gaussian_window(target.T_years),order=48).vix
        time_refined = high.evaluate(states,high.gaussian_window(target.T_years,time_order=16),order=48).vix
        refinements.append(dict(expiry_days=int(target.expiry_days),states="eight ordered fresh-sample VIX positions",
            max_conditional_32_48_change=float(np.max(np.abs(conditional-base[indices]))),
            max_normalizer_48_64_change=float(np.max(np.abs(normalizer-conditional))),
            max_time_8_16_change=float(np.max(np.abs(time_refined-normalizer)))))
    record = dict(status="active reduced-scope exploratory VIX-only fit; numerical gates unresolved",
        scope=scope,scope_sha256=hashlib.sha256(scope_path.read_bytes()).hexdigest(),
        arguments=vars(args)|dict(output=str(args.output)),architecture=asdict(architecture),parameters=asdict(parameters),
        premium=asdict(premium),initial_state=initial.tolist(),active_variance_level=selected,
        free_parameter="single near variance level, affecting both active VIX windows",
        inactive_segments="Equal-level placeholders only; no post-0.25-year variance is identified, fitted or forecast here.",
        curve_penalty=0.,training_cache_role=0,fresh_seeds=seeds,fresh_timings=timings,
        quadrature_refinement=refinements,normalizer_order=48,conditional_order=32,time_order=8,
        target_sha256=targets.artifact_sha256,frozen_model_sha256=frozen_hashes,
        input_compatibility=targets.audit["variance_audit"],
        training_manifest_sha256=hashlib.sha256((old_folder/"run.json").read_bytes()).hexdigest(),
        training_cache_sha256=hashlib.sha256((old_folder/"exposures.npz").read_bytes()).hexdigest(),
        implementation_sha256=source_hashes(sources,ROOT),wall_seconds=time.perf_counter()-start,
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__),
        limitations=["Readout/state/risk prices fixed; no SPX objective or physical feasibility constraint is fitted.",
            "Only active training quotes and forwards enter the objective. Historical candidates are fixed comparisons, not reselected.",
            "Training and fresh budgets are below 4096 outer paths per maturity; no final accuracy classification.",
            "Fresh prefix estimates are dependent; selected-state quadrature changes are not full-distribution bias bounds.",
            "No deferred-maturity or held-out market performance claim follows from a near-only result."])
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    (args.output/"run.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    np.savez_compressed(args.output/"exposures.npz",**stored)
    pd.DataFrame(ledger).to_csv(args.output/"candidate_ledger.csv",index=False)
    pd.DataFrame(scores).to_csv(args.output/"scores.csv",index=False)
    for name,frame in tables.items():
        frame.to_csv(args.output/name,index=False)
    print(pd.DataFrame(scores).to_string(index=False),flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--power",type=int,default=7)
    parser.add_argument("--replications",type=int,default=8)
    parser.add_argument("--seed",type=int,default=20261025)
    args = parser.parse_args()
    if args.output.exists() or args.power < 4 or args.replications < 4:
        parser.error("Require a new output, >=16 points per scramble and >=4 scrambles")
    run(args)
