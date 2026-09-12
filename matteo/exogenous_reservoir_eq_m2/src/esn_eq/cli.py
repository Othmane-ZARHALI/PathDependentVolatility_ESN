"""Command-line entry point for a reproducible parameter-space experiment."""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import ArchitectureConfig, MarketEnvironment, SimulationConfig
from .experiment import AcceptanceCriteria, HierarchicalExplorer
from .sampling import (
    ParameterSpace,
    QmcParameterSampler,
    QmcStructuralSampler,
    StructuralParameters,
)


def build_parser() -> argparse.ArgumentParser:
    """Expose only reproducibility and experiment-size choices."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=64)
    parser.add_argument("--structural-samples", type=int, default=1)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--years", type=float, default=8.0)
    parser.add_argument("--burn-years", type=float, default=2.0)
    parser.add_argument("--paths", type=int, default=8)
    parser.add_argument("--steps-per-observation", type=int, default=4)
    parser.add_argument(
        "--sampler", choices=("latin_hypercube", "sobol"), default="latin_hypercube"
    )
    parser.add_argument("--output", type=Path, default=Path("exploration_results"))
    return parser


def main() -> None:
    """Prepare states once, explore readouts, and persist the full experiment."""
    args = build_parser().parse_args()
    architecture = ArchitectureConfig()
    simulation = SimulationConfig(
        years=args.years,
        burn_years=args.burn_years,
        paths=args.paths,
        steps_per_observation=args.steps_per_observation,
    )
    sampler = QmcParameterSampler(ParameterSpace.plausible_m2())
    parameters = sampler.sample(args.samples, args.seed + 1, args.sampler)
    structural = (
        (StructuralParameters(),)
        if args.structural_samples == 1
        else QmcStructuralSampler().sample(args.structural_samples, args.seed + 2)
    )
    explorer = HierarchicalExplorer(
        architecture,
        simulation,
        MarketEnvironment(),
        AcceptanceCriteria.illustrative(),
    )
    result = explorer.run(structural, parameters, args.seed)
    paths = result.save(args.output)
    print(f"Acceptance rate under illustrative bands: {result.acceptance_rate:.1%}")
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
