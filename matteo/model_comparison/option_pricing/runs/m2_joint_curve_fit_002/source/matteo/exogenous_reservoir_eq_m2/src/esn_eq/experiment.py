"""Joint acceptance, sensitivity, and persistence of parameter exploration."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .config import ArchitectureConfig, MarketEnvironment, ReadoutParameters, SimulationConfig
from .diagnostics import StylisedFactDiagnostics
from .model import ExogenousReservoirVolatilityModel, PreparedPathCache
from .sampling import QmcParameterSampler, StructuralParameters
from .states import ReservoirFactory


@dataclass(frozen=True, slots=True)
class MetricBand:
    """An empirical acceptance interval and its composite-score weight."""

    lower: float
    upper: float
    weight: float = 1.0

    def __post_init__(self) -> None:
        """Require a non-empty interval and non-negative weight."""
        if self.upper <= self.lower or self.weight < 0.0:
            raise ValueError("Require lower < upper and weight >= 0.")

    def distance(self, value: float) -> float:
        """Return zero inside the band and scaled distance outside it."""
        width = self.upper - self.lower
        if value < self.lower:
            return (self.lower - value) / width
        if value > self.upper:
            return (value - self.upper) / width
        return 0.0


@dataclass(frozen=True, slots=True)
class AcceptanceCriteria:
    """User-supplied joint bands; no single metric can hide another failure."""

    bands: Mapping[str, MetricBand]

    @classmethod
    def illustrative(cls) -> "AcceptanceCriteria":
        """Provide placeholders that must be replaced by empirical bands."""
        return cls(
            {
                "max_abs_return_acf": MetricBand(0.0, 0.10),
                "mean_log_variance_acf": MetricBand(0.10, 0.99),
                "rough_hurst": MetricBand(0.02, 0.30),
                "leverage_rank": MetricBand(-0.60, -0.005),
                "zumbach_rank": MetricBand(0.005, 0.60),
                "excess_kurtosis": MetricBand(0.5, 30.0),
                "taylor_gap": MetricBand(-0.05, 0.30, 0.25),
            }
        )

    @classmethod
    def stage2_proxy(cls) -> "AcceptanceCriteria":
        """Declare more demanding synthetic gates for independent validation."""
        return cls(
            {
                "max_abs_return_acf": MetricBand(0.0, 0.10),
                "mean_log_variance_acf": MetricBand(0.15, 0.99),
                "rough_hurst": MetricBand(0.07, 0.16),
                "leverage_rank": MetricBand(-0.60, -0.008),
                "zumbach_rank": MetricBand(0.015, 0.60),
                "excess_kurtosis": MetricBand(3.0, 30.0),
                "taylor_gap": MetricBand(-0.05, 0.30, 0.25),
            }
        )

    def evaluate(self, metrics: Mapping[str, float]) -> tuple[bool, float]:
        """Return strict joint acceptance and a continuous failure distance."""
        distances = []
        for name, band in self.bands.items():
            value = float(metrics[name])
            if not np.isfinite(value):
                return False, float("inf")
            distances.append(band.weight * band.distance(value) ** 2)
        return all(distance == 0.0 for distance in distances), float(np.sum(distances))


@dataclass(frozen=True, slots=True)
class ExplorationResult:
    """Full candidate table, rank sensitivities, and acceptance rate."""

    candidates: pd.DataFrame
    sensitivities: pd.DataFrame
    acceptance_rate: float

    def save(self, directory: str | Path) -> tuple[Path, Path, Path]:
        """Write machine-readable results without discarding failed candidates."""
        output = Path(directory)
        output.mkdir(parents=True, exist_ok=True)
        candidates_path = output / "candidates.csv"
        sensitivity_path = output / "rank_sensitivities.csv"
        summary_path = output / "summary.json"
        self.candidates.to_csv(candidates_path, index=False)
        self.sensitivities.to_csv(sensitivity_path, index=False)
        summary_path.write_text(
            json.dumps(
                {"acceptance_rate": self.acceptance_rate, "candidate_count": len(self.candidates)},
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return candidates_path, sensitivity_path, summary_path


class ParameterExplorer:
    """Evaluate a space-filling readout design on one prepared path cache."""

    def __init__(
        self,
        model: ExogenousReservoirVolatilityModel,
        prepared: PreparedPathCache,
        market: MarketEnvironment,
        criteria: AcceptanceCriteria,
    ) -> None:
        """Keep the feature engine fixed across all candidate evaluations."""
        self.model = model
        self.prepared = prepared
        self.market = market
        self.criteria = criteria
        self.diagnostics = StylisedFactDiagnostics()

    def run(self, parameters: tuple[ReadoutParameters, ...]) -> ExplorationResult:
        """Retain every result so general capacity is not confused with a best fit."""
        rows = [self._evaluate(index, item) for index, item in enumerate(parameters)]
        candidates = pd.DataFrame(rows).sort_values("distance", ignore_index=True)
        sensitivities = self._sensitivities(candidates)
        acceptance_rate = float(candidates["accepted"].mean())
        return ExplorationResult(candidates, sensitivities, acceptance_rate)

    def _evaluate(self, index: int, parameters: ReadoutParameters) -> dict[str, object]:
        """Simulate one cheap readout and attach estimates plus uncertainty."""
        result = self.model.simulate(self.prepared, parameters, self.market, measure="P")
        summary = self.diagnostics.evaluate(result)
        accepted, distance = self.criteria.evaluate(summary.mean)
        row: dict[str, object] = {"candidate_id": index, "accepted": accepted, "distance": distance}
        row.update(QmcParameterSampler.to_record(parameters))
        row.update(summary.mean)
        row.update({f"se_{name}": value for name, value in summary.standard_error.items()})
        row.update(asdict(result.metadata))
        return row

    def _sensitivities(self, candidates: pd.DataFrame) -> pd.DataFrame:
        """Report marginal rank effects as screening evidence, not causality."""
        parameter_names = tuple(QmcParameterSampler.to_record(ReadoutParameters()))
        metric_names = tuple(self.criteria.bands)
        return rank_sensitivities(candidates, parameter_names, metric_names)


def rank_sensitivities(
    candidates: pd.DataFrame,
    parameter_names: tuple[str, ...],
    metric_names: tuple[str, ...],
) -> pd.DataFrame:
    """Compute marginal Spearman screens while retaining p-values."""
    rows = []
    for parameter in parameter_names:
        if candidates[parameter].nunique() < 2:
            continue
        for metric in metric_names:
            coefficient, p_value = spearmanr(candidates[parameter], candidates[metric])
            rows.append(
                {
                    "parameter": parameter,
                    "metric": metric,
                    "spearman_rho": float(coefficient),
                    "p_value": float(p_value),
                }
            )
    return pd.DataFrame(rows)


class HierarchicalExplorer:
    """Nest cheap readouts inside structural scenarios on common random numbers."""

    def __init__(
        self,
        base_architecture: ArchitectureConfig,
        simulation: SimulationConfig,
        market: MarketEnvironment,
        criteria: AcceptanceCriteria,
    ) -> None:
        """Fix dimensions so every outer scenario can share primitive normals."""
        self.base_architecture = base_architecture
        self.simulation = simulation
        self.market = market
        self.criteria = criteria

    def run(
        self,
        structural: tuple[StructuralParameters, ...],
        readouts: tuple[ReadoutParameters, ...],
        seed: int,
    ) -> ExplorationResult:
        """Report the joint capacity region across both hierarchy levels."""
        base_factory = ReservoirFactory(self.base_architecture, self.simulation)
        randomness = base_factory.draw(seed)
        frames = []
        for scenario_id, scenario in enumerate(structural):
            architecture = scenario.architecture(self.base_architecture)
            factory = ReservoirFactory(architecture, self.simulation)
            reservoir = factory.from_randomness(randomness)
            model = ExogenousReservoirVolatilityModel(architecture)
            prepared = model.prepare(reservoir, scenario.risk_premium())
            frame = ParameterExplorer(
                model, prepared, self.market, self.criteria
            ).run(readouts).candidates
            frame.insert(0, "structural_scenario_id", scenario_id)
            for name, value in reversed(tuple(scenario.to_record().items())):
                frame.insert(1, name, value)
            frames.append(frame)
        candidates = pd.concat(frames, ignore_index=True).sort_values(
            "distance", ignore_index=True
        )
        parameter_names = self._parameter_names()
        sensitivity = rank_sensitivities(
            candidates, parameter_names, tuple(self.criteria.bands)
        )
        return ExplorationResult(candidates, sensitivity, float(candidates.accepted.mean()))

    @staticmethod
    def _parameter_names() -> tuple[str, ...]:
        """Return all inner and outer parameter columns used in screening."""
        readout = tuple(QmcParameterSampler.to_record(ReadoutParameters()))
        structural = tuple(StructuralParameters().to_record())
        return structural + readout
