"""Regenerate saved independent pilot streams to audit rare payoff coverage."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .m2_spx import ROOT,M2,frozen_stage3
from .fine_grid import CachedNormalizerPricer
from .spx_hedging_study import choose_quotes,refinement
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.pricing import dated_grid


def run(args):
    sources = capture_sources([Path(__file__),*[Path(__file__).with_name(n) for n in
        ("m2_spx.py","fine_grid.py","spx_hedging_study.py")],*sorted((M2/"src/esn_eq").glob("*.py"))])
    original = Path(__file__).with_name("runs")/"m2_spx_hedging_001"
    record = json.loads((original/"run.json").read_text())
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    if targets.artifact_sha256 != record["target_sha256"]:
        raise ValueError("Target snapshot changed")
    quotes = choose_quotes(targets)
    pd.testing.assert_frame_equal(quotes,pd.read_csv(original/"contracts.csv"),check_dtype=False,check_exact=False)
    architecture,parameters,premium,_ = frozen_stage3()
    pricer = CachedNormalizerPricer(architecture,premium)
    initial,_ = pricer.initial_state([0.,0.,0.])
    curve = NuisanceCurveFamily.from_targets(targets).build(record["curve_levels"])
    events = np.r_[quotes.T_years.unique(),[t for t in curve.ends if t < quotes.T_years.max()]]
    fine = refinement(refinement(dated_grid(events,1/(252*64))))
    if len(fine)-1 != record["grid_intervals"][-1]:
        raise ValueError("Grid changed")
    output = []
    for item in record["seed_records"]:
        counts = np.zeros(len(quotes),int)
        margins = np.full(len(quotes),np.inf)
        for batch,seed in enumerate(item["pilot"]):
            count = min(record["arguments"]["batch_size"],record["arguments"]["pilot"]-batch*record["arguments"]["batch_size"])
            cache = pricer.prepare(fine,initial,count,seed)
            result = pricer.spx(cache,parameters,curve,quotes.T_years.unique(),order=48)
            for j,row in enumerate(quotes.itertuples()):
                difference = result.terminal_ratios[row.T_years]-row.K_over_F
                counts[j] += np.sum(difference > 0) if row.option_type == "call" else np.sum(difference < 0)
                margins[j] = min(margins[j],float(np.min(np.abs(difference))))
            print(f"pilot audit trial={item['trial']} batch={batch+1}",flush=True)
        for j,row in enumerate(quotes.itertuples()):
            output.append(dict(trial=item["trial"],contract=j,expiry_days=row.expiry_days,strike=row.strike,
                option_type=row.option_type,positive_pilot_payoffs=int(counts[j]),min_distance_to_strike=float(margins[j])))
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    (args.output/"run.json").write_text(json.dumps(dict(status="independent pilot coverage audit",
        source_run=str(original.relative_to(ROOT)),target_sha256=targets.artifact_sha256,
        implementation_sha256=source_hashes(sources,ROOT)),indent=2)+"\n")
    table = pd.DataFrame(output)
    table.to_csv(args.output/"pilot_hits.csv",index=False)
    print(table[table.positive_pilot_payoffs == 0].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output exists")
    run(args)
