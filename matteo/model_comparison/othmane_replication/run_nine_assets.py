"""Fresh-seed nine-asset diagnostic for the archived ESN fit and M2 anchor.

This is an exploratory comparison of saved parameter sets, not a new fit.
The historical variance is Garman--Klass; the simulators expose raw ESN
variance and M2 average instantaneous variance respectively.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy

from compare import ArchivedESN, Stage3M2, historical_target, load_legacy


HERE = Path(__file__).resolve().parent
FACTS = ("Hhat", "leverage", "abs_acf", "sq_acf", "zumbach", "kurtosis")
BLOCK_SIZES = (1, 40, 6, 6, 5, 1)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarise(legacy, target, paths, common_h_paths):
    stats = [legacy.all_stats_one_path(x, v, raw) for x, v, raw in paths]
    arrays = {fact: np.asarray([stat[fact] for stat in stats], float) for fact in FACTS}
    if not all(np.isfinite(values).all() for values in arrays.values()):
        raise ValueError("Nonfinite statistic on a simulated path")
    average = {fact: np.mean(values, axis=0) for fact, values in arrays.items()}
    residual = legacy.joint_residual(average, target)
    if residual.size != 59 or not np.isfinite(residual).all():
        raise ValueError("Unexpected residual length or nonfinite residual")
    contributions = {}
    cursor = 0
    for fact, length in zip(FACTS, BLOCK_SIZES):
        contributions[fact] = float(np.sum(residual[cursor:cursor + length] ** 2))
        cursor += length
    path_scores = [float(np.sum(legacy.joint_residual(stat, target) ** 2)) for stat in stats]
    common_h = np.asarray(common_h_paths, float)
    return {
        "statistics": {fact: average[fact].tolist() for fact in FACTS},
        "path_statistics": {fact: arrays[fact].tolist() for fact in FACTS},
        "path_band_p5": {fact: np.percentile(values, 5, axis=0).tolist()
                         for fact, values in arrays.items()},
        "path_band_p95": {fact: np.percentile(values, 95, axis=0).tolist()
                          for fact, values in arrays.items()},
        "weighted_squared_error_by_block": contributions,
        "weighted_squared_error_total": float(np.sum(residual ** 2)),
        "mean_pathwise_weighted_squared_error": float(np.mean(path_scores)),
        "pathwise_weighted_squared_errors": path_scores,
        "Hhat_common_lags_2_40": float(np.mean(common_h)),
        "Hhat_common_lags_2_40_paths": common_h.tolist(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--days", type=int, default=4000)
    parser.add_argument("--burn-days", type=int, default=504)
    parser.add_argument("--paths", type=int, default=20)
    parser.add_argument("--seed", type=int, default=720001)
    parser.add_argument("--output", type=Path, default=HERE / "nine_assets_fresh_001.json")
    args = parser.parse_args()
    if args.days < 200 or args.burn_days < 0 or args.paths < 2:
        parser.error("Require days >= 200, burn-days >= 0 and paths >= 2")

    legacy, engine = load_legacy()
    targets = {}
    hashes = None
    for asset in legacy.NAMES9:
        targets[asset], current_hashes = historical_target(legacy, asset)
        if hashes is None:
            hashes = current_hashes
        elif hashes != current_hashes:
            raise AssertionError("Target file hashes changed during the run")

    m2_adapter = Stage3M2()
    m2_paths = list(m2_adapter.paths(args.days, args.burn_days, args.paths, args.seed))
    if any(len(x) != args.days or len(v) != args.days for x, v, _ in m2_paths):
        raise ValueError("M2 observation count differs from requested retained days")
    m2_h_common = [legacy.estimate_H(0.5 * np.log(np.maximum(v, 1e-30)),
                                     np.arange(2, 41)) for _, v, _ in m2_paths]

    assets = {}
    archive_ci = json.loads((HERE.parents[2] / "SimultanuousAllFeatures" /
                             "joint_ci.json").read_text())
    for asset in legacy.NAMES9:
        esn_adapter = ArchivedESN(legacy, engine, asset)
        esn_paths = list(esn_adapter.paths(args.days, args.burn_days, args.paths, args.seed))
        esn_h_common = [legacy.estimate_H(0.5 * np.log(np.maximum(raw, 1e-30)),
                                         np.arange(2, 41)) for _, _, raw in esn_paths]
        esn = summarise(legacy, targets[asset], esn_paths, esn_h_common)
        m2 = summarise(legacy, targets[asset], m2_paths, m2_h_common)
        esn["parameters"] = esn_adapter.metadata()
        m2["parameters"] = m2_adapter.metadata()
        esn["archived_band_mean_absolute_difference"] = {
            fact: float(np.mean(np.abs(np.asarray(esn["statistics"][fact]) -
                                       np.asarray(archive_ci[asset][fact]["mean"]))))
            for fact in FACTS
        }
        assets[asset] = {
            "target": {fact: np.asarray(targets[asset][fact]).tolist() for fact in FACTS},
            "models": {"othmane_archived_fit": esn, "m2_stage3_anchor": m2},
        }
        print(f"{asset:12s} ESN {esn['weighted_squared_error_total']:8.3f} "
              f"M2 {m2['weighted_squared_error_total']:8.3f}", flush=True)

    result = {
        "status": "exploratory_unmatched_variance_proxy_and_in_sample_targets",
        "simulation": {"days": args.days, "burn_days": args.burn_days,
                       "paths": args.paths, "seed_first": args.seed,
                       "seed_last": args.seed + args.paths - 1},
        "objective": "Othmane 59-term residual, squared after averaging path statistics",
        "variance_observation": {
            "historical": "same-day Garman--Klass OHLC variance",
            "othmane": "raw model variance for H/leverage/Zumbach; EWMA variance drives returns",
            "m2": "daily mean instantaneous variance; no OHLC proxy",
        },
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__},
        "sha256": {
            "historical_targets": hashes,
            "legacy_calibration": sha256(HERE.parents[2] / "SimultanuousAllFeatures" /
                                         "joint_calibration.py"),
            "legacy_simulator": sha256(HERE.parents[2] / "SimultanuousAllFeatures" /
                                       "esn_base.py"),
            "archived_fit": sha256(HERE.parents[2] / "SimultanuousAllFeatures" /
                                   "joint_calibration_final.json"),
            "archived_bands": sha256(HERE.parents[2] / "SimultanuousAllFeatures" /
                                     "joint_ci.json"),
            "run_code": sha256(Path(__file__)),
        },
        "assets": assets,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
