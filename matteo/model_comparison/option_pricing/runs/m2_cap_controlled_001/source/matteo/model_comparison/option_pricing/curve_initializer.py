"""Constrained nuisance-curve initializer, not an option-model calibration."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, minimize

from .m2_spx import ROOT, M2
from .provenance import capture_sources, source_hashes, check_sources, save_sources
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.market_targets import MarketTargets
from matteo.preprocessing.option_variance import vix_option_moment_bound


def overlap(ends, start, end):
    return np.maximum(0.,np.minimum(end,ends)-np.maximum(start,np.r_[0.,ends[:-1]]))


def run(args):
    sources = capture_sources([Path(__file__),Path(__file__).with_name("m2_spx.py"),
        ROOT/"matteo/preprocessing/option_variance.py",*sorted((M2/"src/esn_eq").glob("*.py"))])
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    family = NuisanceCurveFamily.from_targets(targets)
    strip = targets.curve("training_mid_piecewise_v1")
    spx = targets.select("SPX",expiry_days=[29,57,211,302]).drop_duplicates("expiry")
    A = np.array([overlap(family.ends,0.,t) for t in spx.T_years])
    W = np.array([strip.integral(0.,t) for t in spx.T_years])
    constraints, descriptions = [], []
    for target in targets.vix_forwards.itertuples():
        training = targets.select("VIX").query("expiry == @target.expiry")
        bound = vix_option_moment_bound(training,target.forward_vix_points,target.forward_vix_points)
        constraints.append(1e4/(30/365)*overlap(family.ends,target.T_years,target.T_years+30/365))
        descriptions.append(dict(expiry=int(target.expiry),expiry_days=int(target.expiry_days),
            assumed_model_mean=float(target.forward_vix_points),
            necessary_lower_bound=max(bound["moment_lower_bound"],target.forward_vix_points**2)))
    B = np.array(constraints)
    lower = np.array([row["necessary_lower_bound"] for row in descriptions])
    scaled = A/W[:,None]
    difference = np.diff(np.eye(len(family.ends)),axis=0)
    def objective(levels):
        residual = scaled @ levels-1
        return float(np.mean(residual**2)+family.penalty(levels))
    def jacobian(levels):
        return 2*scaled.T @ (scaled @ levels-1)/len(W)+2*family.smoothness_weight*difference.T @ difference @ levels
    initial = np.array([strip.integral(a,b)/(b-a) for a,b in zip(np.r_[0.,family.ends[:-1]],family.ends)])
    # Keep objective evaluations defined outside bounds: penalties are algebraic here.
    def solver_objective(levels):
        return float(np.mean((scaled @ levels-1)**2)+family.smoothness_weight*np.sum(np.diff(levels)**2))
    result = minimize(solver_objective,initial,jac=jacobian,method="SLSQP",
        bounds=Bounds(family.lower,family.upper),constraints=[LinearConstraint(B,lower,np.inf)],
        options=dict(ftol=1e-12,maxiter=1000))
    if not result.success or np.min(B @ result.x-lower) < -1e-7:
        raise RuntimeError(f"Initializer failed: {result.message}")
    candidate = family.build(result.x)
    checks = pd.DataFrame(descriptions)
    checks["candidate_second_moment"] = B @ result.x
    checks["moment_margin"] = checks.candidate_second_moment-checks.necessary_lower_bound
    checks["strip_second_moment"] = [1e4/(30/365)*strip.integral(t.T_years,t.T_years+30/365)
                                        for t in targets.vix_forwards.itertuples()]
    comparison = spx[["expiry","expiry_days","T_years"]].copy()
    comparison["strip_total_variance"] = W
    comparison["initializer_total_variance"] = A @ result.x
    comparison["relative_change"] = (A @ result.x-W)/W
    record = dict(status="nuisance-curve starting point only; not an SPX/VIX option fit",
        family=asdict(family),levels=result.x.tolist(),initial_levels=initial.tolist(),
        objective="mean squared relative deviation from diagnostic stripped SPX W + declared curve penalty",
        objective_value=objective(result.x),curve_penalty=family.penalty(result.x),
        solver=dict(method="SLSQP",success=bool(result.success),iterations=int(result.nit),message=result.message),
        assumed_vix_mean="quoted market forwards; candidate checks must be recomputed at actual model means",
        minimum_necessary_moment_margin=float(np.min(B @ result.x-lower)),
        input_compatibility=targets.audit["variance_audit"],target_sha256=targets.artifact_sha256,
        implementation_sha256=source_hashes(sources,ROOT),
        limitations=["Necessary moment inequalities are not sufficient for matching VIX option intervals.",
            "Stripped SPX variance is a finite-strike diagnostic; this auxiliary initializer is not an independent market loss.",
            "No model option prices or physical compatibility are assessed; original audit flags remain unchanged.",
            "All candidate ranking must use the common quote loss, not this starting-point objective."])
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    (args.output/"run.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    checks.to_csv(args.output/"necessary_moment_checks.csv",index=False)
    comparison.to_csv(args.output/"strip_comparison.csv",index=False)
    print(json.dumps(dict(levels=record["levels"],objective=record["objective_value"]),indent=2))
    print(checks.to_string(index=False))
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output exists")
    run(args)
