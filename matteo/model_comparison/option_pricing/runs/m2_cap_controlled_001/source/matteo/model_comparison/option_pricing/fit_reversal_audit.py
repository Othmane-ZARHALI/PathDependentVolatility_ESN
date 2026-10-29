"""Paired, fixed-candidate diagnosis of the archived curve-fit reversal.

No refitting or candidate reselection. Complete Sobol scrambles are the units;
leave-one-scramble-out sensitivity is descriptive with only four scrambles.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .m2_spx import ROOT, M2
from .provenance import capture_sources, source_hashes, check_sources, save_sources
from esn_eq.calibration_loss import JointPriceObjective, LossConfig
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.market_targets import MarketTargets


def run(args):
    sources = capture_sources([Path(__file__),Path(__file__).with_name("m2_spx.py"),
                               *sorted((M2/"src/esn_eq").glob("*.py"))])
    inputs = {name:(args.fit/name).read_bytes() for name in
              ("run.json","exposures.npz","prices_train.csv","prices_check.csv")}
    record = json.loads(inputs["run.json"])
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    if record["target_sha256"] != targets.artifact_sha256:
        raise ValueError("Target snapshot changed")
    quotes = targets.select("VIX").reset_index(drop=True)
    forwards = targets.vix_forwards.sort_values("expiry").reset_index(drop=True)
    family = NuisanceCurveFamily.from_targets(targets)
    objective = JointPriceObjective(quotes,forwards,LossConfig(**record["loss_config"]))
    candidates = {"starting":np.array(record["starting_levels"]),"selected":np.array(record["selected_levels"])}
    R,N = record["arguments"]["replications"],2**record["arguments"]["power"]
    archive = np.load(args.fit/"exposures.npz")
    floor = 1e-4  # Original run's declared fixed variance floor; checked by repricing below.
    scores,loo,contributions,prices,summaries,profiles = [],[],[],[],[],[]
    for role_id,role in enumerate(("train","check")):
        data = {}
        for candidate,levels in candidates.items():
            p,f = np.zeros((R,len(quotes))),np.zeros((R,len(forwards)))
            hits = np.zeros(len(quotes),int)
            for fi,target in enumerate(forwards.itertuples()):
                expiry = int(target.expiry)
                manifest = record["exposure_manifest"][str(expiry)]
                coeff = archive[f"role{role_id}_expiry{expiry}"].reshape(R,N,-1)
                values = 100*np.sqrt((manifest["floor_integral"]+coeff @ (levels-floor))/manifest["delta"])
                group = quotes[quotes.expiry == expiry]
                sign = np.where(group.option_type.to_numpy() == "call",1.,-1.)
                payoff = np.maximum(sign*(values[:,:,None]-group.strike.to_numpy()),0.)
                p[:,group.index],f[:,fi] = payoff.mean(axis=1),values.mean(axis=1)
                hits[group.index] = (payoff > 0).sum(axis=(0,1))
            result = objective.evaluate(p.mean(axis=0),f.mean(axis=0),curve_penalty=family.penalty(levels))
            data[candidate] = dict(p=p,f=f,result=result)
            scores.append(dict(role=role,candidate=candidate,total=result["total"],**result["components"]))
            if candidate == "selected":
                np.testing.assert_allclose(result["total"],record["scores"][role]["total"],rtol=1e-12)
                old = pd.read_csv(args.fit/f"prices_{role}.csv")
                np.testing.assert_allclose(p.mean(axis=0),old.model_price_U,rtol=1e-11,atol=1e-12)
            frame = result["quotes"].copy()
            frame["role"],frame["candidate"] = role,candidate
            frame["positive_payoff_count"] = hits
            frame["scramble_se_U"] = p.std(axis=0,ddof=1)/np.sqrt(R)
            frame["se_in_frozen_scales"] = frame.scramble_se_U/frame.price_scale_U
            prices.append(frame)
        a,b = data["selected"],data["starting"]
        delta = a["result"]["total"]-b["result"]["total"]
        quote_delta = a["result"]["quotes"].copy()
        quote_delta["role"] = role
        quote_delta["weighted_loss_change"] = quote_delta.loss_weight*(
            a["result"]["quotes"].interval_residual**2-b["result"]["quotes"].interval_residual**2)
        contributions.append(quote_delta)
        leave_differences = []
        for omitted in range(R):
            keep = np.arange(R) != omitted
            values = {}
            for candidate,d in data.items():
                values[candidate] = objective.evaluate(d["p"][keep].mean(axis=0),d["f"][keep].mean(axis=0),
                    curve_penalty=family.penalty(candidates[candidate]))["total"]
            difference = values["selected"]-values["starting"]
            loo.append(dict(role=role,omitted_scramble=omitted,starting=values["starting"],
                            selected=values["selected"],selected_minus_starting=difference))
            leave_differences.append(difference)
        jackknife_se = np.sqrt((R-1)/R*np.sum((leave_differences-np.mean(leave_differences))**2))
        summaries.append(dict(role=role,selected_minus_starting=delta,leave_one_out_changes=leave_differences,
            descriptive_jackknife_se=float(jackknife_se),
            note="Nonlinear interval loss, four scrambles and unresolved tails: no confidence certificate; train includes selection."))
    check_sources(sources)
    if any((args.fit/name).read_bytes() != raw for name,raw in inputs.items()):
        raise RuntimeError("Input artifacts changed during audit")
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    report = dict(status="fixed-candidate paired audit; no refitting or reselection",
        fit=str(args.fit),input_sha256={k:hashlib.sha256(v).hexdigest() for k,v in inputs.items()},
        target_sha256=targets.artifact_sha256,implementation_sha256=source_hashes(sources,ROOT),
        candidates={k:v.tolist() for k,v in candidates.items()},floor=floor,replications=R,points=N,
        summaries=summaries,loss_config=asdict(objective.config))
    (args.output/"run.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    pd.DataFrame(scores).to_csv(args.output/"scores.csv",index=False)
    pd.DataFrame(loo).to_csv(args.output/"leave_one_scramble_out.csv",index=False)
    pd.concat(contributions,ignore_index=True).to_csv(args.output/"quote_contributions.csv",index=False)
    pd.concat(prices,ignore_index=True).to_csv(args.output/"paired_prices.csv",index=False)
    print(json.dumps(summaries,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit",type=Path,default=Path(__file__).with_name("runs")/"m2_vix_curve_fit_001")
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output exists")
    run(args)
