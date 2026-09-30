"""Coarse historical M2 readout screen with an early/late date split.

Select among the 32 Stage-3 M2 readouts on each asset's pre-2015 return
facts. Evaluate the selected readout and synthetic anchor on later historical
targets with independent simulation seeds. ESN's archived parameters used the
full historical sample, so its late-period row is descriptive, not held out.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy

from compare import ROOT, ArchivedESN, Stage3M2, historical_target, load_legacy

sys.path.insert(0, str(ROOT))
from matteo.preprocessing import kurtosis, market_data  # noqa: E402


HERE = Path(__file__).resolve().parent
FACTS = ("abs_acf", "sq_acf", "kurtosis")
SPLIT = "2014-12-31"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def return_target(legacy, frame):
    returns, _ = market_data.returns_and_variance(frame)
    return {
        "abs_acf": legacy.acf_at_lags(np.abs(returns), legacy.TAYLOR_LAGS),
        "sq_acf": legacy.acf_at_lags(returns ** 2, legacy.TAYLOR_LAGS),
        "kurtosis": kurtosis.corrected_excess_kurtosis(returns),
    }


def return_path_stats(legacy, paths):
    values = {fact: [] for fact in FACTS}
    for returns, _, _ in paths:
        values["abs_acf"].append(legacy.acf_at_lags(np.abs(returns), legacy.TAYLOR_LAGS))
        values["sq_acf"].append(legacy.acf_at_lags(returns ** 2, legacy.TAYLOR_LAGS))
        values["kurtosis"].append(legacy.excess_kurtosis(returns))
    return {fact: np.asarray(values[fact]) for fact in FACTS}


def return_stats(legacy, paths):
    values = return_path_stats(legacy, paths)
    return {fact: np.mean(values[fact], axis=0) for fact in FACTS}


def return_score(legacy, stats, target):
    tolerances = legacy.TOLERANCES
    return float(
        np.sum(((stats["abs_acf"] - target["abs_acf"]) /
                (tolerances["abs_acf"] * np.sqrt(6))) ** 2)
        + np.sum(((stats["sq_acf"] - target["sq_acf"]) /
                  (tolerances["sq_acf"] * np.sqrt(6))) ** 2)
        + ((stats["kurtosis"] - target["kurtosis"]) /
           tolerances["kurtosis"]) ** 2
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen-days", type=int, default=2000)
    parser.add_argument("--screen-paths", type=int, default=8)
    parser.add_argument("--screen-seed", type=int, default=850001)
    parser.add_argument("--eval-days", type=int, default=4000)
    parser.add_argument("--eval-paths", type=int, default=20)
    parser.add_argument("--eval-seed", type=int, default=720001)
    parser.add_argument("--output", type=Path, default=HERE / "m2_readout_screen_001.json")
    args = parser.parse_args()
    legacy, _ = load_legacy()

    frames = market_data.load_bench9()
    targets = {}
    counts = {}
    for asset, frame in frames.items():
        early = frame.loc[frame.index <= SPLIT]
        late = frame.loc[frame.index > SPLIT]
        targets[asset] = {
            "early": return_target(legacy, early),
            "late": return_target(legacy, late),
        }
        counts[asset] = {"early_ohlc_rows": len(early), "late_ohlc_rows": len(late),
                         "early_first": str(early.index.min().date()),
                         "early_last": str(early.index.max().date()),
                         "late_first": str(late.index.min().date()),
                         "late_last": str(late.index.max().date())}
        historical_target(legacy, asset)  # Ensure saved full-sample target files are available.

    anchor = Stage3M2()
    with anchor.candidates_path.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["variant"] == "M2"]
    candidates = {}
    for row in rows:
        if row["structural_scenario_id"] == "0":
            candidates[int(row["candidate_id"])] = row
    if len(candidates) != 32:
        raise AssertionError(f"Expected 32 M2 readouts, got {len(candidates)}")

    screened = {}
    screen_stats = {}
    for candidate_id, row in sorted(candidates.items()):
        adapter = Stage3M2()
        adapter.candidate = row
        paths = list(adapter.paths(args.screen_days, 252, args.screen_paths, args.screen_seed))
        stats = return_stats(legacy, paths)
        screen_stats[candidate_id] = stats
        screened[candidate_id] = {
            asset: return_score(legacy, stats, targets[asset]["early"])
            for asset in legacy.NAMES9
        }
        print(f"screened M2 readout {candidate_id:2d}", flush=True)

    selected_ids = {asset: min(screened, key=lambda i: screened[i][asset])
                    for asset in legacy.NAMES9}
    selected_unique = sorted(set(selected_ids.values()) | {int(anchor.candidate["candidate_id"])})
    eval_stats = {}
    eval_path_stats = {}
    for candidate_id in selected_unique:
        adapter = Stage3M2()
        adapter.candidate = candidates[candidate_id]
        paths = list(adapter.paths(args.eval_days, 504, args.eval_paths, args.eval_seed))
        eval_path_stats[candidate_id] = return_path_stats(legacy, paths)
        eval_stats[candidate_id] = {
            fact: np.mean(values, axis=0)
            for fact, values in eval_path_stats[candidate_id].items()
        }
        print(f"evaluated M2 readout {candidate_id:2d}", flush=True)

    comparison = {}
    anchor_id = int(anchor.candidate["candidate_id"])
    for asset in legacy.NAMES9:
        selected = selected_ids[asset]
        late = targets[asset]["late"]
        early = targets[asset]["early"]
        comparison[asset] = {
            "selected_candidate_id": selected,
            "early_screen_score_selected": screened[selected][asset],
            "early_screen_score_anchor": screened[anchor_id][asset],
            "late_score_selected_fresh_seeds": return_score(legacy, eval_stats[selected], late),
            "late_score_anchor_fresh_seeds": return_score(legacy, eval_stats[anchor_id], late),
            "late_score_selected_vs_early_target": return_score(legacy, eval_stats[selected], early),
            "late_score_anchor_vs_early_target": return_score(legacy, eval_stats[anchor_id], early),
        }
        print(f"{asset:12s} candidate {selected:2d}: late score "
              f"{comparison[asset]['late_score_selected_fresh_seeds']:.3f} "
              f"vs anchor {comparison[asset]['late_score_anchor_fresh_seeds']:.3f}", flush=True)

    output = {
        "status": "coarse_readout_screen_return_only_variance_proxy_free",
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__},
        "sha256": {"stock_ohlc": sha256(market_data.STOCK_FILE),
                   "index_ohlc": sha256(market_data.INDEX_FILE),
                   "candidate_table": sha256(anchor.candidates_path),
                   "m2_manifest": sha256(anchor.manifest_path),
                   "screen_code": sha256(__file__)},
        "split": {"early_through": SPLIT, "late_after": SPLIT,
                  "counts_and_dates": counts},
        "screen_simulation": {"days": args.screen_days, "burn_days": 252,
                              "paths": args.screen_paths, "seed_first": args.screen_seed},
        "evaluation_simulation": {"days": args.eval_days, "burn_days": 504,
                                  "paths": args.eval_paths, "seed_first": args.eval_seed},
        "candidate_rule": "all 32 M2 Stage-3 readouts, structural scenario 0; choose minimum early-period return-only SSE separately for each asset",
        "anchor_candidate_id": anchor_id,
        "candidate_parameters": {str(i): candidates[i] for i in candidates},
        "targets": {asset: {period: {fact: np.asarray(value).tolist()
                                      for fact, value in period_target.items()}
                            for period, period_target in by_period.items()}
                    for asset, by_period in targets.items()},
        "screen_scores": {str(i): scores for i, scores in screened.items()},
        "screen_statistics": {str(i): {fact: np.asarray(value).tolist()
                                       for fact, value in stats.items()}
                              for i, stats in screen_stats.items()},
        "evaluation_statistics": {str(i): {fact: np.asarray(value).tolist()
                                           for fact, value in stats.items()}
                                  for i, stats in eval_stats.items()},
        "evaluation_path_statistics": {str(i): {fact: np.asarray(value).tolist()
                                                for fact, value in stats.items()}
                                       for i, stats in eval_path_stats.items()},
        "comparison": comparison,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
