"""High-resolution SPX check of a VIX-only profiled curve, without refitting."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from .m2_spx import ROOT, M2, frozen_stage3
from .fine_grid import CachedNormalizerPricer, spx_check
from .provenance import capture_sources, source_hashes, check_sources, save_sources
from esn_eq.calibration_loss import JointPriceObjective, LossConfig
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.market_targets import MarketTargets


def run(args):
    sources = capture_sources([Path(__file__),Path(__file__).with_name("fine_grid.py"),
        Path(__file__).with_name("m2_spx.py"),*sorted((M2/"src/esn_eq").glob("*.py"))])
    fitted = json.loads(args.fit.read_text())
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    architecture,parameters,premium,hashes = frozen_stage3()
    family = NuisanceCurveFamily.from_targets(targets)
    expected = json.loads(json.dumps(dict(architecture=asdict(architecture),parameters=asdict(parameters),
        premium=asdict(premium),family=asdict(family),target_sha256=targets.artifact_sha256)))
    if any(fitted[key] != value for key,value in expected.items()):
        raise ValueError("Fitted curve has different fixed model or input snapshot")
    pricer = CachedNormalizerPricer(architecture,premium)
    initial,_ = pricer.initial_state([0.,0.,0.])
    if fitted["initial_state"] != initial.tolist():
        raise ValueError("Dated state changed")
    curve = family.build(fitted["selected_levels"])
    report,prices = spx_check(args,pricer,initial,parameters,curve,targets)
    quotes = targets.select("SPX",expiry_days=[29,57])
    objective = JointPriceObjective(quotes,targets.vix_forwards.iloc[:0],
        LossConfig(vix_weight=0.,forward_weight=0.))
    # Align predictions to the objective's frozen contract order.
    keys = ["instrument","expiry","strike","option_type"]
    import pandas as pd
    ordered = prices.set_index(keys).loc[pd.MultiIndex.from_frame(quotes[keys])]
    try:
        result = objective.evaluate(ordered.fine_control_price_U.to_numpy(),np.array([]))
        score = dict(total=result["total"],components=result["components"],
            inside_fraction=float(np.mean(result["quotes"].interval_residual == 0)))
        evaluation_failure = None
    except ValueError as error:
        result,score = None,None
        evaluation_failure = str(error)
    record = dict(status="SPX consequence of exploratory VIX-only curve fit; no SPX optimization or feasibility claim",
        arguments=vars(args)|dict(output=str(args.output),fit=str(args.fit)),
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),
        initial_state=initial.tolist(),family=asdict(family),levels=list(curve.levels),
        fit_sha256=hashlib.sha256(args.fit.read_bytes()).hexdigest(),
        target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,
        implementation_sha256=source_hashes(sources,ROOT),input_compatibility=targets.audit["variance_audit"],
        spx=report,score=score,evaluation_failure=evaluation_failure,
        limitations=["Only 29/57-day SPX training expiries are evaluated; long-maturity and physical checks remain open.",
            "SPX quote-relative accuracy remains unresolved; SE excludes normalizer approximation.",
            "Negative controlled prices fail the entire loss evaluation rather than being clipped or removed."])
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    (args.output/"run.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    prices.to_csv(args.output/"prices.csv",index=False)
    if result is not None:
        result["quotes"].to_csv(args.output/"price_residuals.csv",index=False)
        result["blocks"].to_csv(args.output/"blocks.csv",index=False)
    print(json.dumps(dict(score=score,evaluation_failure=evaluation_failure,spx=report),indent=2),flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit",type=Path,default=Path(__file__).with_name("runs")/"m2_vix_curve_fit_001/run.json")
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--paths",type=int,default=4096)
    parser.add_argument("--pilot-paths",type=int,default=2048)
    parser.add_argument("--batch-size",type=int,default=512)
    parser.add_argument("--coarse-steps",type=int,default=128)
    parser.add_argument("--order",type=int,default=48)
    parser.add_argument("--seed",type=int,default=20261021)
    args = parser.parse_args()
    if args.output.exists() or min(args.paths,args.pilot_paths,args.batch_size,args.order) < 2 or args.coarse_steps < 1:
        parser.error("Require a new output and valid numerical budgets")
    if any(n % args.batch_size == 1 for n in (args.paths,args.pilot_paths)):
        parser.error("Avoid a final batch with one path")
    run(args)
