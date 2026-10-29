"""Exploratory two-coefficient VIX curve fit using cached Gaussian exposures.

Readout, architecture, Q dynamics and dated state are fixed. This is a low-budget
VIX-only profile; SPX compatibility and final numerical gates remain open.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import least_squares

from .m2_spx import ROOT, M2, frozen_stage3
from .provenance import capture_sources, source_hashes, check_sources, save_sources
from esn_eq.calibration_loss import JointPriceObjective, LossConfig
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.market_targets import MarketTargets
from esn_eq.outer_sampling import GaussianOuterSampler, replication_statistics
from esn_eq.vix_curve_profile import prepare_curve_exposures
from esn_eq.vix_pricing import ConditionalVIXPricer
from matteo.preprocessing.option_variance import vix_option_moment_bound


class CurveFitData:
    def __init__(self, quotes, forwards, caches, replications, points):
        self.quotes, self.forwards, self.caches = quotes, forwards, caches
        self.replications, self.points = replications, points

    def predictions(self, levels, details=False):
        prices = np.empty(len(self.quotes))
        price_errors, price_hits, repeat_hits = prices.copy(),prices.copy(),prices.copy()
        forwards, forward_errors, moments, moment_errors = {},{},{},{}
        for expiry, group in self.quotes.groupby("expiry"):
            values = self.caches[int(expiry)].evaluate(levels).vix.reshape(self.replications,self.points)
            forward, forward_se = replication_statistics(values.mean(axis=1))
            moment, moment_se = replication_statistics((values**2).mean(axis=1))
            forwards[int(expiry)], forward_errors[int(expiry)] = float(forward),float(forward_se)
            moments[int(expiry)], moment_errors[int(expiry)] = float(moment),float(moment_se)
            sign = np.where(group.option_type.to_numpy() == "call",1.,-1.)
            payoffs = np.maximum(sign*(values[:,:,None]-group.strike.to_numpy()),0.)
            mean,se = replication_statistics(payoffs.mean(axis=1))
            prices[group.index], price_errors[group.index] = mean,se
            price_hits[group.index] = (payoffs > 0).sum(axis=(0,1))
            repeat_hits[group.index] = ((payoffs > 0).sum(axis=1) > 0).sum(axis=0)
        ordered_forwards = np.array([forwards[int(expiry)] for expiry in self.forwards.expiry])
        if not details:
            return prices,ordered_forwards
        return dict(prices=prices,forwards=ordered_forwards,price_se=price_errors,
            positive_payoff_count=price_hits,positive_payoff_replications=repeat_hits,
            forward_se=[forward_errors[int(expiry)] for expiry in self.forwards.expiry],
            second_moment=[moments[int(expiry)] for expiry in self.forwards.expiry],
            second_moment_se=[moment_errors[int(expiry)] for expiry in self.forwards.expiry])


def run(args):
    start = time.perf_counter()
    sources = capture_sources([Path(__file__),Path(__file__).with_name("m2_spx.py"),
        ROOT/"matteo/preprocessing/option_variance.py",*sorted((M2/"src/esn_eq").glob("*.py"))])
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    family = NuisanceCurveFamily.from_targets(targets)
    quotes = targets.select("VIX").reset_index(drop=True)
    forwards = targets.vix_forwards.sort_values("expiry").reset_index(drop=True)
    architecture,parameters,premium,frozen_hashes = frozen_stage3()
    from esn_eq.pricing import DatedM2Pricer
    dated = DatedM2Pricer(architecture,premium)
    initial,_ = dated.initial_state([0.,0.,0.])
    fixed_third = targets.curve("training_mid_piecewise_v1").levels[-1]
    template = family.build([.04,.05,fixed_third])
    engine = ConditionalVIXPricer(dated,parameters,template,initial,normalizer_order=48)
    windows = {int(f.expiry):engine.gaussian_window(float(f.T_years),time_order=8) for f in forwards.itertuples()}
    if any(np.any(w.times >= family.ends[1]) for w in windows.values()):
        raise ValueError("This pilot assumes the third segment is outside all VIX windows")
    points = 2**args.power
    seed_records, exposure_arrays, check_states, timings = [],{},{},[]

    def build(role):
        caches = {}
        for target in forwards.itertuples():
            tick = time.perf_counter()
            expiry = int(target.expiry)
            sampler = GaussianOuterSampler(dated,float(target.T_years),initial)
            states = []
            for replicate in range(args.replications):
                seed = int(np.random.SeedSequence([args.seed,role,expiry,replicate]).generate_state(1)[0])
                states.append(sampler.draw(args.power,seed,"sobol"))
                seed_records.append(dict(role=role,expiry=expiry,replicate=replicate,seed=seed))
            states = np.concatenate(states)
            cache = prepare_curve_exposures(engine,states,windows[expiry],order=32)
            if np.any(cache.coefficients[:,2] != 0):
                raise ValueError("Unidentified third coefficient unexpectedly affects VIX")
            caches[expiry] = cache
            exposure_arrays[f"role{role}_expiry{expiry}"] = cache.coefficients
            if role == 1:
                check_states[expiry] = states
            elapsed = time.perf_counter()-tick
            timings.append(dict(role=role,expiry=expiry,wall_seconds=elapsed))
            print(json.dumps(dict(role=role,expiry_days=int(target.expiry_days),
                points=points*args.replications,exposure_seconds=elapsed)),flush=True)
        return CurveFitData(quotes,forwards,caches,args.replications,points)

    train = build(0)
    config = LossConfig(spx_weight=0.)
    objective = JointPriceObjective(quotes,forwards,config)
    initializer_path = Path(__file__).with_name("runs")/"m2_curve_initializer_001/run.json"
    initializer = json.loads(initializer_path.read_text())
    if initializer["target_sha256"] != targets.artifact_sha256 or initializer["family"] != json.loads(json.dumps(asdict(family))):
        raise ValueError("Curve initializer inputs do not match")
    starting = [*initializer["levels"][:2],fixed_third]
    evaluations = 0
    def levels_from(x):
        levels = np.r_[x,fixed_third]
        family.build(levels)
        return levels
    def score(data,levels):
        prices,means = data.predictions(levels)
        return objective.evaluate(prices,means,curve_penalty=family.penalty(levels))
    def residual(x):
        nonlocal evaluations
        evaluations += 1
        return score(train,levels_from(x))["residuals"]
    fit_start = time.perf_counter()
    starts = [starting[:2],[.03,.04],[.04,.06],[.06,.05],[.08,.08],[.10,.10],[.02,.10],[.10,.02]]
    ledger = []
    for index,x in enumerate(starts):
        fit = least_squares(residual,x,bounds=(family.lower[:2],family.upper[:2]),
            max_nfev=300,ftol=1e-9,xtol=1e-9,gtol=1e-9)
        levels = levels_from(fit.x)
        result = score(train,levels)
        ledger.append(dict(start_id=index,start_xi0=x[0],start_xi1=x[1],xi0=levels[0],xi1=levels[1],
            fixed_xi2=fixed_third,total=result["total"],vix=result["components"]["VIX"],
            forwards=result["components"]["forwards"],curve_penalty=result["components"]["curve_penalty"],
            success=bool(fit.success),nfev=int(fit.nfev),message=fit.message))
    converged = [row for row in ledger if row["success"]]
    if not converged:
        raise RuntimeError("No optimizer start converged")
    selected = min(converged,key=lambda row:row["total"])
    fitted = np.array([selected["xi0"],selected["xi1"],fixed_third])
    fit_seconds = time.perf_counter()-fit_start
    print(json.dumps(dict(selected_levels=fitted.tolist(),training_loss=selected["total"],
        optimizer_seconds=fit_seconds,evaluations=evaluations)),flush=True)
    # Selection is complete before any independently scrambled checking data.
    check = build(1)
    output_tables, moment_checks, scores = {},[],{}
    for role,data in (("train",train),("check",check)):
        result = score(data,fitted)
        diagnostics = data.predictions(fitted,details=True)
        q,f = result["quotes"],result["forwards"]
        q["replication_se_U"] = diagnostics["price_se"]
        q["positive_payoff_count"] = diagnostics["positive_payoff_count"].astype(int)
        q["positive_payoff_replications"] = diagnostics["positive_payoff_replications"].astype(int)
        f["replication_se"] = diagnostics["forward_se"]
        f["second_moment"] = diagnostics["second_moment"]
        f["second_moment_se"] = diagnostics["second_moment_se"]
        output_tables[f"prices_{role}.csv"],output_tables[f"forwards_{role}.csv"] = q,f
        output_tables[f"blocks_{role}.csv"] = result["blocks"]
        scores[role] = dict(total=result["total"],components=result["components"],
            initializer_total=score(data,starting)["total"],
            inside_fraction=float(np.mean(q.interval_residual == 0)),zero_hit_quotes=int(np.sum(q.positive_payoff_count == 0)))
        for target,mean in zip(forwards.itertuples(),diagnostics["forwards"]):
            group = quotes[quotes.expiry == target.expiry]
            bound = vix_option_moment_bound(group,float(target.forward_vix_points),float(mean))
            moment = 1e4/(30/365)*family.build(fitted).integral(target.T_years,target.T_years+30/365)
            moment_checks.append(dict(role=role,expiry=int(target.expiry),model_mean=float(mean),
                curve_second_moment=moment,necessary_bound=bound["moment_lower_bound"],
                moment_margin=moment-bound["moment_lower_bound"],jensen_margin=moment-float(mean)**2))
    # Separate refinement of conditional GH and the global normalizer at selected check states.
    refinement = []
    curve = family.build(fitted)
    direct = ConditionalVIXPricer(dated,parameters,curve,initial,normalizer_order=48)
    high = ConditionalVIXPricer(dated,parameters,curve,initial,normalizer_order=64)
    for target in forwards.itertuples():
        expiry = int(target.expiry)
        base_all = check.caches[expiry].evaluate(fitted).vix
        indices = np.argsort(base_all)[np.linspace(0,len(base_all)-1,8).astype(int)]
        states = check_states[expiry][indices]
        conditional = direct.evaluate(states,direct.gaussian_window(target.T_years),order=48).vix
        normalizer = high.evaluate(states,high.gaussian_window(target.T_years),order=48).vix
        refinement.append(dict(expiry=expiry,states="eight ordered check-sample VIX positions",
            base=base_all[indices].tolist(),conditional48=conditional.tolist(),normalizer64=normalizer.tolist(),
            max_conditional_change=float(np.max(np.abs(conditional-base_all[indices]))),
            max_normalizer_change=float(np.max(np.abs(normalizer-conditional)))))
    record = dict(status="exploratory VIX-only curve fit; numerical and SPX/physical gates unresolved",
        calibration_scope="two nuisance-curve coefficients; all structural/readout/state/risk-price parameters fixed",
        arguments=vars(args)|dict(output=str(args.output)),architecture=asdict(architecture),
        parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),family=asdict(family),
        selected_levels=fitted.tolist(),starting_levels=starting,free_coefficients=[0,1],
        fixed_coefficients={"2":fixed_third},third_coefficient_identified_by_vix=False,
        quote_split="train_only; numerical check uses same training contracts with fresh scrambles",
        loss_config=asdict(config),scores=scores,seed_records=seed_records,exposure_timings=timings,
        optimizer_seconds=fit_seconds,residual_evaluations=evaluations,normalizer_order=48,conditional_order=32,time_order=8,
        quadrature_refinement=refinement,input_compatibility=targets.audit["variance_audit"],
        target_sha256=targets.artifact_sha256,frozen_model_sha256=frozen_hashes,
        initializer_sha256=hashlib.sha256(initializer_path.read_bytes()).hexdigest(),
        implementation_sha256=source_hashes(sources,ROOT),
        exposure_manifest={str(e):dict(delta=w.delta,floor_integral=engine.variance_floor*float(w.weights.sum()))
                           for e,w in windows.items()},
        wall_seconds=time.perf_counter()-start,
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__),
        limitations=["Preliminary curve-only fit uses fewer than the plan's 4096 outer screening paths per maturity.",
            "Independent checking scrambles assess sampling sensitivity, not held-out market performance.",
            "Outer SE excludes conditional quadrature, normalizer and time error; zero-hit tails remain unresolved.",
            "Only curve levels are optimized; no SPX price objective or physical constraints are enforced.",
            "Necessary moment checks at sampled model means are diagnostics, not feasibility certificates."])
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    (args.output/"run.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    np.savez_compressed(args.output/"exposures.npz",**exposure_arrays)
    pd.DataFrame(ledger).to_csv(args.output/"candidate_ledger.csv",index=False)
    pd.DataFrame(moment_checks).to_csv(args.output/"necessary_moment_checks.csv",index=False)
    for name,frame in output_tables.items():
        frame.to_csv(args.output/name,index=False)
    print(json.dumps(dict(scores=scores,selected_levels=fitted.tolist(),wall_seconds=record["wall_seconds"]),indent=2),flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--power",type=int,default=6)
    parser.add_argument("--replications",type=int,default=4)
    parser.add_argument("--seed",type=int,default=20261023)
    args = parser.parse_args()
    if args.output.exists() or args.power < 3 or args.replications < 4:
        parser.error("Require a new output, >=8 points per scramble and >=4 scrambles")
    run(args)
