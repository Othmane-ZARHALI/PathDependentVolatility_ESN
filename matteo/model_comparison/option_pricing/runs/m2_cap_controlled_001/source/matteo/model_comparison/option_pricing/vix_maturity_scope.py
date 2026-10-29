"""Fixed-candidate maturity-exclusion diagnostic; no refitting or selection."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .m2_spx import ROOT,M2
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.calibration_loss import JointPriceObjective,LossConfig
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.market_targets import MarketTargets


def run(args):
    sources = capture_sources([Path(__file__),Path(__file__).with_name("m2_spx.py"),
                               *sorted((M2/"src/esn_eq").glob("*.py"))])
    folder = Path(__file__).with_name("runs")/"m2_vix_curve_fit_001"
    fit = json.loads((folder/"run.json").read_text())
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    if targets.artifact_sha256 != fit["target_sha256"]:
        raise ValueError("Target snapshot changed")
    family = NuisanceCurveFamily.from_targets(targets)
    arrays = np.load(folder/"exposures.npz")
    rows = []
    for scope,days in (("all_four",[27,55,118,153]),("without_153",[27,55,118]),("near_two",[27,55])):
        quotes = targets.select("VIX",expiry_days=days).reset_index(drop=True)
        forwards = targets.vix_forwards[targets.vix_forwards.expiry_days.isin(days)].sort_values("expiry").reset_index(drop=True)
        objective = JointPriceObjective(quotes,forwards,LossConfig(**fit["loss_config"]))
        for role,label in ((0,"train"),(1,"check")):
            for candidate in ("starting","selected"):
                levels = np.array(fit[f"{candidate}_levels"])
                prices = np.zeros(len(quotes)); means = []
                for target in forwards.itertuples():
                    expiry = int(target.expiry)
                    manifest = fit["exposure_manifest"][str(expiry)]
                    values = 100*np.sqrt((manifest["floor_integral"]+
                        arrays[f"role{role}_expiry{expiry}"] @ (levels-1e-4))/manifest["delta"])
                    means.append(values.mean())
                    group = quotes[quotes.expiry == expiry]
                    sign = np.where(group.option_type.to_numpy() == "call",1.,-1.)
                    prices[group.index] = np.maximum(sign*(values[:,None]-group.strike.to_numpy()),0.).mean(axis=0)
                result = objective.evaluate(prices,means,curve_penalty=family.penalty(levels))
                if scope == "all_four" and candidate == "selected":
                    np.testing.assert_allclose(result["total"],fit["scores"][label]["total"],rtol=1e-12)
                rows.append(dict(scope=scope,expiry_days="|".join(map(str,days)),role=label,candidate=candidate,
                    quotes=len(quotes),forward_targets=len(forwards),total=result["total"],**result["components"]))
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    record = dict(status="fixed-candidate scope sensitivity; primary four-maturity protocol unchanged",
        excluded_rows="both VIX options and associated forward for each excluded expiry",
        weights="expiry averages recomputed over included maturities; market quote scales unchanged",
        calibration_performed=False,candidate_reselection=False,
        fit_sha256=hashlib.sha256((folder/"run.json").read_bytes()).hexdigest(),
        cache_sha256=hashlib.sha256((folder/"exposures.npz").read_bytes()).hexdigest(),
        target_sha256=targets.artifact_sha256,implementation_sha256=source_hashes(sources,ROOT))
    (args.output/"run.json").write_text(json.dumps(record,indent=2)+"\n")
    table = pd.DataFrame(rows)
    table.to_csv(args.output/"scores.csv",index=False)
    print(table[["scope","role","candidate","total"]].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output exists")
    run(args)
