"""Evaluate archived ESN A2 and M2 parameters on Othmane's historical facts.

Run from any working directory. This is a baseline evaluator, not a fitter.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
LEGACY = ROOT / "SimultanuousAllFeatures"
M2 = ROOT / "matteo" / "exogenous_reservoir_eq_m2"
TARGETS = ROOT / "data" / "processed"


def load_legacy():
    """Import the archived estimator without modifying its code or cwd."""
    sys.path.insert(0, str(LEGACY))
    return importlib.import_module("joint_calibration"), importlib.import_module("esn_base")


def load_m2():
    """Import M2 source in place, without an editable installation."""
    sys.path.insert(0, str(M2 / "src"))
    return (
        importlib.import_module("esn_eq.config"),
        importlib.import_module("esn_eq.model"),
        importlib.import_module("esn_eq.states"),
    )


def _read_json(path: Path):
    return json.loads(path.read_text())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def historical_target(legacy, asset: str):
    """Use the archived target construction and record exact input hashes."""
    stems = (
        "bench9_real_hurst", "bench9_real_leverage_wide",
        "bench9_real_abs_acf_40", "bench9_real_sq_acf_40",
        "bench9_real_zumbach_profile", "bench9_real_kurtosis",
    )
    paths = [TARGETS / (stem + ".json") for stem in stems]
    values = [_read_json(path) for path in paths]
    target = legacy.build_target_vector(
        values[0][asset], values[1][asset], values[2][asset],
        values[3][asset], values[4], values[5][asset], asset,
    )
    return target, {path.name: _sha256(path) for path in paths}


class ArchivedESN:
    """Othmane's saved per-asset fit and its original simulator."""

    def __init__(self, legacy, engine, asset: str):
        self.legacy = legacy
        self.engine = engine
        self.asset = asset
        self.fit = _read_json(LEGACY / "joint_calibration_final.json")[asset]

    def paths(self, days: int, burn_days: int, count: int, seed: int):
        fit = self.fit
        params = dict(
            self.legacy.P_FIXED,
            H=fit["H_fit"], rough_scale=fit["rough_scale_fit"],
            lam_lo=fit["lam_lo_fit"], lam_hi=fit["lam_hi_fit"],
            m1=fit["m1_fit"],
        )
        matrices = self.engine.build_esn_matrices(self.legacy.OPT_ARCH, params)
        if matrices["kappa0"] >= 0:
            raise ValueError("Archived ESN configuration violates kappa0 < 0")
        for index in range(count):
            returns, variance, raw_variance = self.engine._sim_esn_with_params(
                seed + index, days + burn_days, self.legacy.OPT_ARCH, matrices,
                alpha=fit["alpha_fit"],
            )
            yield returns[burn_days:], variance[burn_days:], raw_variance[burn_days:]

    def metadata(self):
        return {"fit": self.fit, "source": "SimultanuousAllFeatures/joint_calibration_final.json"}


class Stage3M2:
    """The first all-seed Stage-3 M2 configuration, evaluated as an anchor."""

    def __init__(self):
        self.manifest_path = M2 / "experiments" / "stage3_001" / "manifest.json"
        self.candidates_path = M2 / "experiments" / "stage3_001" / "top_configurations.csv"
        self.manifest = _read_json(self.manifest_path)
        with self.candidates_path.open(newline="") as handle:
            candidates = list(csv.DictReader(handle))
        self.candidate = next(
            row for row in candidates
            if row["variant"] == "M2" and row["passes_all_seeds"] == "True"
        )

    def paths(self, days: int, burn_days: int, count: int, seed: int):
        config, model_module, states = load_m2()
        row = self.candidate
        architecture = config.ArchitectureConfig.from_mapping(self.manifest["architecture"])
        architecture = replace(
            architecture,
            feedback=replace(
                architecture.feedback,
                asset_correlation=float(row["feedback_asset_correlation"]),
            ),
            spike=replace(
                architecture.spike,
                asset_correlation=float(row["spike_asset_correlation"]),
            ),
        )
        parameters = config.ReadoutParameters(
            cluster_loading=float(row["cluster_loading"]),
            feedback_curvature=float(row["feedback_curvature"]),
            feedback_shift=float(row["feedback_shift"]),
            spike_curvature=float(row["spike_curvature"]),
            orthogonal_curvature=0.0,
            echo_loading=0.0,
            cap_level=float(row["cap_level"]),
            cap_sharpness=float(row["cap_sharpness"]),
        )
        premium = config.RiskPremium(
            equity=float(row["equity_risk_price"]),
            feedback_idiosyncratic=float(row["feedback_idiosyncratic_risk_price"]),
            spike_idiosyncratic=float(row["spike_idiosyncratic_risk_price"]),
        )
        simulation = config.SimulationConfig(
            years=(days + burn_days) / 252,
            burn_years=burn_days / 252,
            paths=count,
            observations_per_year=252,
            steps_per_observation=2,
        )
        cache = states.ReservoirFactory(architecture, simulation).build(seed)
        model = model_module.ExogenousReservoirVolatilityModel(architecture)
        result = model.simulate(model.prepare(cache, premium), parameters,
                                config.MarketEnvironment(), measure="P")
        for index in range(count):
            yield result.log_returns[index], result.realized_variance[index], None

    def metadata(self):
        return {
            "candidate_id": int(self.candidate["candidate_id"]),
            "structural_scenario_id": int(self.candidate["structural_scenario_id"]),
            "source": "matteo/exogenous_reservoir_eq_m2/experiments/stage3_001",
            "manifest_sha256": _sha256(self.manifest_path),
            "candidate_table_sha256": _sha256(self.candidates_path),
        }


def score_paths(legacy, adapter, target, days: int, burn_days: int, count: int, seed: int):
    """Apply the same archived statistic extractor and objective to each path."""
    samples = [legacy.all_stats_one_path(x, v, v_raw)
               for x, v, v_raw in adapter.paths(days, burn_days, count, seed)]
    names = ("Hhat", "leverage", "abs_acf", "sq_acf", "zumbach", "kurtosis")
    averaged = {name: np.mean([sample[name] for sample in samples], axis=0)
                for name in names}
    if any(not np.isfinite(value).all() for value in averaged.values()):
        raise ValueError("Non-finite diagnostics; increase sample length or inspect simulation")
    residual = legacy.joint_residual(averaged, target)
    if residual.size != 59:
        raise AssertionError(f"Expected 59 objective entries, got {residual.size}")
    block_sizes = (1, 40, 6, 6, 5, 1)
    offset = 0
    contributions = {}
    for name, length in zip(names, block_sizes):
        contributions[name] = float(np.sum(residual[offset:offset + length] ** 2))
        offset += length
    return {
        "statistics": {name: np.asarray(averaged[name]).tolist() for name in names},
        "weighted_squared_error_by_block": contributions,
        "weighted_squared_error_total": float(np.sum(residual ** 2)),
        "parameters": adapter.metadata(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset", default="SP500")
    parser.add_argument("--models", choices=("othmane", "m2", "both"), default="both")
    parser.add_argument("--days", type=int, default=1000)
    parser.add_argument("--burn-days", type=int, default=252)
    parser.add_argument("--paths", type=int, default=8)
    parser.add_argument("--seed", type=int, default=93271)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.days < 200 or args.burn_days < 0 or args.paths < 2:
        parser.error("Require days >= 200, burn-days >= 0 and paths >= 2")
    legacy, engine = load_legacy()
    if args.asset not in legacy.NAMES9:
        parser.error(f"Unknown asset {args.asset!r}; choose from {legacy.NAMES9}")
    target, hashes = historical_target(legacy, args.asset)
    adapters = {}
    if args.models in ("othmane", "both"):
        adapters["othmane_archived_fit"] = ArchivedESN(legacy, engine, args.asset)
    if args.models in ("m2", "both"):
        adapters["m2_stage3_anchor"] = Stage3M2()
    output = {
        "status": "exploratory_unmatched_variance_proxy",
        "asset": args.asset,
        "simulation": {"days": args.days, "burn_days": args.burn_days,
                       "paths": args.paths, "seed": args.seed},
        "target_sha256": hashes,
        "target": {key: np.asarray(value).tolist() for key, value in target.items()},
        "models": {name: score_paths(legacy, adapter, target, args.days,
                                     args.burn_days, args.paths, args.seed)
                   for name, adapter in adapters.items()},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(args.output)
    for name, result in output["models"].items():
        print(f"{name}: weighted squared error {result['weighted_squared_error_total']:.4f}")


if __name__ == "__main__":
    main()
