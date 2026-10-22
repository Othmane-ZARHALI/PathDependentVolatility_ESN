"""Run a checkpointed M2/M3/M2+O/M4 physical-measure pilot."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, replace
from math import isclose
from pathlib import Path

import pandas as pd

from .calibration import (
    ArchitectureSearchBounds,
    LatentArchitectureCalibrator,
    LatentCalibrationEnsemble,
    LatentCalibrationResult,
    LatentTargets,
)
from .config import (
    ArchitectureConfig,
    BankConfig,
    MarketEnvironment,
    SimulationConfig,
    geometric_rates,
)
from .experiment import AcceptanceCriteria
from .factorial import (
    FactorialDesign,
    FactorialExperimentConfig,
    FactorialExperimentRecorder,
    MultiSeedFactorialExperiment,
)
from .reporting import FactorialMarkdownReport
from .sampling import QmcStructuralSampler, StructuralParameters


def build_parser() -> argparse.ArgumentParser:
    """Expose design sizes, independent seeds, and simulation resolution."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("factorial_pilot"))
    parser.add_argument("--candidate-samples", type=int, default=32)
    parser.add_argument("--structural-samples", type=int, default=4)
    parser.add_argument("--seeds", type=int, nargs="+", default=[104729, 130363, 155921])
    parser.add_argument("--design-seed", type=int, default=20260901)
    parser.add_argument(
        "--calibration-seeds",
        type=int,
        nargs="+",
        default=[20260902, 20260903, 20260904],
    )
    parser.add_argument("--calibration-iterations", type=int, default=60)
    parser.add_argument("--years", type=float, default=6.0)
    parser.add_argument("--burn-years", type=float, default=1.0)
    parser.add_argument("--paths", type=int, default=8)
    parser.add_argument("--steps-per-observation", type=int, default=2)
    parser.add_argument("--sampler", choices=("latin_hypercube", "sobol"), default="sobol")
    parser.add_argument("--skip-latent-calibration", action="store_true")
    return parser


def main() -> None:
    """Calibrate latent shapes, run checkpoints, and render the comparison."""
    args = build_parser().parse_args()
    base, latent = _architecture(args)
    architecture = _maximal_architecture(base)
    simulation = SimulationConfig(
        years=args.years,
        burn_years=args.burn_years,
        paths=args.paths,
        steps_per_observation=args.steps_per_observation,
    )
    structures = _structural_design(args.structural_samples, args.design_seed + 1)
    design = FactorialDesign.sample(args.candidate_samples, args.design_seed, args.sampler)
    config = FactorialExperimentConfig(tuple(args.seeds))
    runner = MultiSeedFactorialExperiment(
        architecture,
        simulation,
        MarketEnvironment(),
        AcceptanceCriteria.illustrative(),
    )
    recorder = FactorialExperimentRecorder(args.output)
    result = runner.run(structures, design, config, recorder)
    latent_record = _latent_record(latent, base)
    _save_json(args.output / "latent_calibration.json", latent_record)
    _save_latent_replicates(
        args.output / "latent_calibration_replicates.csv", latent_record
    )
    report = FactorialMarkdownReport().save(
        args.output / "RESULTS.md", result, latent_record
    )
    FactorialMarkdownReport.stress_summary(result.runs).to_csv(
        args.output / "posthoc_threshold_sensitivity.csv", index=False
    )
    FactorialMarkdownReport.feedback_curvature_screen(result.runs).to_csv(
        args.output / "feedback_curvature_screen.csv", index=False
    )
    print(f"Completed {len(config.seeds)} seeds and {len(result.runs)} paired rows.")
    print(report)


def _architecture(
    args: argparse.Namespace,
) -> tuple[ArchitectureConfig, LatentCalibrationEnsemble | None]:
    """Fit only the deterministic latent targets unless explicitly skipped."""
    base = ArchitectureConfig()
    if args.skip_latent_calibration:
        return base, None
    ensemble = LatentArchitectureCalibrator().fit_replicated(
        base,
        LatentTargets(),
        seeds=tuple(args.calibration_seeds),
        max_iterations=args.calibration_iterations,
    )
    return ensemble.best.architecture, ensemble


def _maximal_architecture(base: ArchitectureConfig) -> ArchitectureConfig:
    """Enable both optional feature engines while keeping their loadings cheap."""
    return replace(
        base,
        orthogonal=BankConfig(geometric_rates(0.5, 100.0, 12)),
        echo=replace(base.echo, enabled=True),
    )


def _structural_design(count: int, seed: int) -> tuple[StructuralParameters, ...]:
    """Use the default scenario exactly when no outer exploration is requested."""
    if count < 1:
        raise ValueError("structural-samples must be positive.")
    if count == 1:
        return (StructuralParameters(),)
    return QmcStructuralSampler().sample(count, seed)


def _latent_record(
    ensemble: LatentCalibrationEnsemble | None,
    architecture: ArchitectureConfig,
) -> dict[str, object]:
    """Keep latent-fit evidence separate from the stochastic experiment tables."""
    if ensemble is None:
        return {
            "objective": None,
            "optimizer_message": "latent calibration skipped",
            "feedback_rate_range": [
                architecture.feedback.rates[0],
                architecture.feedback.rates[-1],
            ],
            "cluster_rate_range": [
                architecture.cluster.rates[0],
                architecture.cluster.rates[-1],
            ],
            "boundary_hits": [],
            "replicates": [],
        }
    result = ensemble.best
    return {
        "optimizer_seeds": list(ensemble.seeds),
        "selected_seed": ensemble.best_seed,
        "objective": result.objective,
        "objective_range": [
            min(item.objective for item in ensemble.results),
            max(item.objective for item in ensemble.results),
        ],
        "optimizer_message": result.optimizer_message,
        "rough_variogram": result.rough_variogram.tolist(),
        "cluster_acf": result.cluster_acf.tolist(),
        "feedback_rate_range": [
            result.architecture.feedback.rates[0],
            result.architecture.feedback.rates[-1],
        ],
        "cluster_rate_range": [
            result.architecture.cluster.rates[0],
            result.architecture.cluster.rates[-1],
        ],
        "boundary_hits": _calibration_boundary_hits(result.architecture),
        "replicates": [
            _calibration_replicate(seed, item)
            for seed, item in zip(ensemble.seeds, ensemble.results, strict=True)
        ],
        "architecture": asdict(result.architecture),
    }


def _calibration_replicate(
    seed: int,
    result: LatentCalibrationResult,
) -> dict[str, object]:
    """Serialise one optimizer replication without discarding fitted coordinates."""
    architecture = result.architecture
    return {
        "seed": seed,
        "objective": result.objective,
        "optimizer_message": result.optimizer_message,
        "feedback_low": architecture.feedback.rates[0],
        "feedback_high": architecture.feedback.rates[-1],
        "feedback_power": architecture.feedback_weight_power,
        "cluster_low": architecture.cluster.rates[0],
        "cluster_high": architecture.cluster.rates[-1],
        "cluster_power": architecture.cluster_weight_power,
        "boundary_hits": _calibration_boundary_hits(architecture),
    }


def _calibration_boundary_hits(architecture: ArchitectureConfig) -> list[str]:
    """Flag endpoint solutions that need a wider-box or regularisation check."""
    bounds = ArchitectureSearchBounds()
    values = {
        "feedback_low": architecture.feedback.rates[0],
        "feedback_high": architecture.feedback.rates[-1],
        "feedback_power": architecture.feedback_weight_power,
        "cluster_low": architecture.cluster.rates[0],
        "cluster_high": architecture.cluster.rates[-1],
        "cluster_power": architecture.cluster_weight_power,
    }
    hits = []
    for name, value in values.items():
        lower, upper = getattr(bounds, name)
        if isclose(value, lower, rel_tol=1e-6, abs_tol=1e-12):
            hits.append(f"{name}:lower")
        if isclose(value, upper, rel_tol=1e-6, abs_tol=1e-12):
            hits.append(f"{name}:upper")
    return hits


def _save_json(path: Path, payload: dict[str, object]) -> None:
    """Write one deterministic, human-readable calibration record."""
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _save_latent_replicates(path: Path, payload: dict[str, object]) -> None:
    """Expose optimizer-seed stability in a directly comparable table."""
    pd.DataFrame(payload.get("replicates", [])).to_csv(path, index=False)


if __name__ == "__main__":
    main()
