"""Batched, coupled 64/128-step M2 SPX and fixed-state VIX integration checks."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import platform

import numpy as np
import pandas as pd

from .m2_spx import ROOT, M2, frozen_stage3
from .provenance import capture_sources, source_hashes, check_sources, save_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.pricing import DatedM2Pricer, dated_grid, coarsen_paths, payoff_statistics
from esn_eq.variance_reduction import (deterministic_stock_controls,
    fit_payoff_controls, controlled_payoff_statistics)
from esn_eq.vix_pricing import ConditionalVIXPricer


class CachedNormalizerPricer(DatedM2Pricer):
    """Cache only deterministic normalizers within this fixed-model experiment."""
    def __init__(self, architecture, premium):
        super().__init__(architecture, premium)
        self.normalizers = {}

    def normalizer(self, times, initial, parameters, order):
        key = (np.asarray(times).tobytes(), np.asarray(initial).tobytes(), parameters, order)
        if key not in self.normalizers:
            values = super().normalizer(times, initial, parameters, order)
            values.setflags(write=False)
            self.normalizers[key] = values
        return self.normalizers[key]


def summarize(values):
    return float(np.mean(values)), float(np.std(values, ddof=1)/np.sqrt(len(values)))


def bundle(pricer, fine, coarse, initial, parameters, curve, maturities,
           paths, batch_size, seed, role, order, include_coarse):
    collected = {name: {float(t): [] for t in maturities}
                 for name in ("fine", "coarse", "reference")}
    seeds = []
    for batch, start in enumerate(range(0, paths, batch_size)):
        count = min(batch_size, paths-start)
        if count < 2:
            raise ValueError("Each batch needs at least two paths")
        actual_seed = int(np.random.SeedSequence([seed, role, batch]).generate_state(1)[0])
        seeds.append(actual_seed)
        cache = pricer.prepare(fine, initial, count, actual_seed)
        result = pricer.spx(cache, parameters, curve, maturities, order=order)
        reference = deterministic_stock_controls(cache, curve, maturities)
        for t in maturities:
            collected["fine"][float(t)].append(result.terminal_ratios[float(t)])
            collected["reference"][float(t)].append(reference[float(t)])
        if include_coarse:
            low = pricer.spx(coarsen_paths(cache, coarse), parameters, curve, maturities, order=order)
            for t in maturities:
                collected["coarse"][float(t)].append(low.terminal_ratios[float(t)])
            del low
        del cache, result, reference
        print(f"role={role} batch={batch+1} completed_paths={start+count}/{paths}", flush=True)
    return {name: {t: np.concatenate(parts) for t, parts in values.items() if parts}
            for name, values in collected.items()}, seeds


def spx_check(args, pricer, initial, parameters, curve, targets):
    quotes = targets.select("SPX", "train", [29,57])
    maturities = np.sort(quotes.T_years.unique())
    events = np.r_[maturities, [k for k in curve.ends if k < maturities[-1]]]
    coarse = dated_grid(events, 1/(252*args.coarse_steps))
    fine = np.r_[0., np.concatenate([np.linspace(a,b,3)[1:] for a,b in zip(coarse[:-1],coarse[1:])])]
    pilot, pilot_seeds = bundle(pricer,fine,coarse,initial,parameters,curve,maturities,
        args.pilot_paths,args.batch_size,args.seed,0,args.order,False)
    evaluation, seeds = bundle(pricer,fine,coarse,initial,parameters,curve,maturities,
        args.paths,args.batch_size,args.seed,1,args.order,True)
    rows, reports = [], []
    for t, group in quotes.groupby("T_years"):
        forward = float(group.forward.iloc[0])
        strikes, kinds = group.strike.to_numpy()/forward, group.option_type.to_numpy()
        total_variance = curve.integral(0.,float(t))
        beta = fit_payoff_controls(pilot["fine"][t],pilot["reference"][t],strikes,kinds,total_variance)
        yf,yc,z = evaluation["fine"][t],evaluation["coarse"][t],evaluation["reference"][t]
        pf,sf = payoff_statistics(yf,strikes,kinds,forward)
        pc,sc = payoff_statistics(yc,strikes,kinds,forward)
        cf,ef = controlled_payoff_statistics(yf,z,strikes,kinds,total_variance,beta,forward)
        cc,ec = controlled_payoff_statistics(yc,z,strikes,kinds,total_variance,beta,forward)
        frame = group.copy()
        frame["fine_price_U"],frame["fine_se_U"],frame["coarse_price_U"],frame["coarse_se_U"] = pf,sf,pc,sc
        frame["fine_control_price_U"],frame["fine_control_se_U"] = cf,ef
        frame["coarse_control_price_U"],frame["coarse_control_se_U"] = cc,ec
        diffs,errors,control_diffs,control_errors = [],[],[],[]
        for k,kind,coefficient in zip(strikes,kinds,beta):
            sign = 1 if kind == "call" else -1
            difference = forward*(np.maximum(sign*(yf-k),0)-np.maximum(sign*(yc-k),0))
            mean,se = summarize(difference); diffs.append(mean); errors.append(se)
            mean,se = summarize(difference-forward*coefficient[0]*(yf-yc))
            control_diffs.append(mean); control_errors.append(se)
        frame["paired_change_U"],frame["paired_change_se_U"] = diffs,errors
        frame["control_paired_change_U"],frame["control_paired_change_se_U"] = control_diffs,control_errors
        frame["fine_positive_payoff_count"] = [int(np.sum(yf>k)) if kind=="call" else int(np.sum(yf<k))
                                              for k,kind in zip(strikes,kinds)]
        frame["coarse_positive_payoff_count"] = [int(np.sum(yc>k)) if kind=="call" else int(np.sum(yc<k))
                                                for k,kind in zip(strikes,kinds)]
        rows.append(frame)
        imax = int(np.argmax(np.abs(control_diffs)))
        fy,fse = summarize(yf); cy,cse = summarize(yc)
        reports.append(dict(expiry_days=int(group.expiry_days.iloc[0]),quotes=len(group),
            coarse_forward_ratio=cy,coarse_forward_se=cse,fine_forward_ratio=fy,fine_forward_se=fse,
            max_abs_paired_change_U=float(np.max(np.abs(diffs))),
            max_abs_control_paired_change_U=float(abs(control_diffs[imax])),
            se_at_max_control_change_U=float(control_errors[imax]),
            strike_at_max_control_change=float(group.strike.iloc[imax]),
            max_control_fine_se_U=float(ef.max()),
            max_change_in_price_scales=float(np.max(np.abs(control_diffs)/group.price_scale_U)),
            negative_control_prices=int(np.sum((cf<0)|(cc<0))),
            numerical_status="refinement diagnostic; no converged-price claim"))
    frame = pd.concat(rows)
    return dict(coarse_steps=args.coarse_steps,fine_steps=2*args.coarse_steps,
        coarse_intervals=len(coarse)-1,fine_intervals=len(fine)-1,
        max_lambda_dt_coarse=float(max(pricer.layout.rates)*np.diff(coarse).max()),
        max_lambda_dt_fine=float(max(pricer.layout.rates)*np.diff(fine).max()),
        pilot_batch_seeds=pilot_seeds,evaluation_batch_seeds=seeds,expiries=reports), frame


def vix_check(args,pricer,initial,parameters,curve):
    engine = ConditionalVIXPricer(pricer,parameters,curve,initial,normalizer_order=args.order)
    states = np.array([pricer.initial_state(u)[0] for u in ([0.,0.,0.],[0.,-1.,-1.],[1.,-2.,-2.])])
    reports = []
    for steps in (args.coarse_steps,2*args.coarse_steps):
        rule = engine.window(55/365,max_step=1/(252*steps))
        gh = engine.evaluate(states,rule,order=args.conditional_order)
        nested = engine.evaluate(states,rule,engine="nested_mc",inner_paths=args.vix_inner,seed=args.seed+17)
        reports.append(dict(steps_per_day=steps,nodes=len(rule.times),
            state_projections_C_Y_F=[[0.,0.,0.],[0.,-1.,-1.],[1.,-2.,-2.]],
            quadrature_vix=gh.vix.tolist(),nested_vix=nested.vix.tolist(),
            nested_integral_se=nested.integrated_variance_se.tolist()))
    return dict(expiry_days=55,states="three fixed diagnostic states, not an outer sample",
        conditional_order=args.conditional_order,inner=args.vix_inner,levels=reports,
        quadrature_fine_minus_coarse=(np.array(reports[1]["quadrature_vix"])-reports[0]["quadrature_vix"]).tolist(),
        note="Nested random streams differ across integration grids; no paired nested grid-error claim.")


def run(args):
    sources = capture_sources([Path(__file__),Path(__file__).with_name("m2_spx.py"),*sorted((M2/"src/esn_eq").glob("*.py"))])
    targets = MarketTargets.load(args.snapshot)
    architecture,parameters,premium,frozen_hashes = frozen_stage3()
    pricer = CachedNormalizerPricer(architecture,premium)
    initial,_ = pricer.initial_state([0.,0.,0.])
    curve = targets.curve("training_mid_piecewise_v1")
    start = time.perf_counter()
    record = dict(status="fixed-model fine-grid diagnostic; no calibration",arguments=vars(args).copy(),
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),
        initial_state=initial.tolist(),input_compatibility=targets.audit["variance_audit"],
        target_sha256=targets.artifact_sha256,frozen_model_sha256=frozen_hashes,
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__),
        implementation_sha256=source_hashes(sources,ROOT))
    record["arguments"] = {k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}
    if args.mode in ("spx","both"):
        record["spx"],prices = spx_check(args,pricer,initial,parameters,curve,targets)
    if args.mode in ("vix","both"):
        record["vix"] = vix_check(args,pricer,initial,parameters,curve)
    record["wall_seconds"] = time.perf_counter()-start
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    (args.output/"run.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    if "spx" in record:
        prices.to_csv(args.output/"prices.csv",index=False)
    print(json.dumps({k:v for k,v in record.items() if k in ("spx","vix","wall_seconds")},indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--snapshot",type=Path,default=ROOT/"data/processed/options/2026-02-19")
    parser.add_argument("--mode",choices=["spx","vix","both"],default="both")
    parser.add_argument("--paths",type=int,default=4096)
    parser.add_argument("--pilot-paths",type=int,default=2048)
    parser.add_argument("--batch-size",type=int,default=512)
    parser.add_argument("--coarse-steps",type=int,default=64)
    parser.add_argument("--order",type=int,default=48)
    parser.add_argument("--conditional-order",type=int,default=32)
    parser.add_argument("--vix-inner",type=int,default=512)
    parser.add_argument("--seed",type=int,default=20261021)
    args = parser.parse_args()
    if min(args.paths,args.pilot_paths,args.batch_size,args.order,args.conditional_order,args.vix_inner)<2 or args.coarse_steps<1:
        parser.error("Invalid numerical budget")
    if any(n % args.batch_size == 1 for n in (args.paths,args.pilot_paths)):
        parser.error("Avoid a final batch with one path")
    if args.output.exists():
        parser.error("Output exists")
    run(args)
