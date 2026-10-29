"""Paired multi-seed ablations for the optional M3 and M4 mechanisms."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, replace
from enum import Enum
from pathlib import Path
from typing import Literal, Mapping, cast

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from .config import ArchitectureConfig, MarketEnvironment, ReadoutParameters, SimulationConfig
from .diagnostics import DiagnosticConfig, StylisedFactDiagnostics
from .experiment import AcceptanceCriteria
from .model import ExogenousReservoirVolatilityModel, PreparedPathCache
from .sampling import (
    ParameterSpace,
    QmcParameterSampler,
    StructuralParameters,
)
from .states import ReservoirFactory


COMPARISON_METRICS = (
    "max_abs_return_acf",
    "mean_log_variance_acf",
    "rough_hurst",
    "leverage",
    "leverage_rank",
    "zumbach",
    "zumbach_rank",
    "excess_kurtosis",
    "absolute_return_tail_ratio",
    "hill_tail_index",
    "hill_threshold_instability",
    "taylor_gap",
    "taylor_fraction",
    "mean_volatility",
    "q995_volatility",
    "q_forward_mean_max_error",
    "p_cap_exceedance_fraction",
    "q_cap_exceedance_fraction",
)

SamplerMethod = Literal["latin_hypercube", "sobol"]


class ModelVariant(str, Enum):
    """Four cells of the echo-by-orthogonal factorial experiment."""

    M2 = "M2"
    M3 = "M3"
    M2_ORTHOGONAL = "M2+O"
    M4 = "M4"

    @property
    def uses_echo(self) -> bool:
        """Return whether the cell activates the nonlinear echo loading."""
        return self in {ModelVariant.M3, ModelVariant.M4}

    @property
    def uses_orthogonal(self) -> bool:
        """Return whether the cell activates the orthogonal quadratic loading."""
        return self in {ModelVariant.M2_ORTHOGONAL, ModelVariant.M4}


VARIANT_ORDER = tuple(ModelVariant)


@dataclass(frozen=True, slots=True)
class FactorialCell:
    """One paired model cell for one common full-parameter draw."""

    candidate_id: int
    variant: ModelVariant
    parameters: ReadoutParameters


@dataclass(frozen=True, slots=True)
class FactorialDesign:
    """A fixed M4 design whose zeroed loadings generate all nested cells."""

    full_parameters: tuple[ReadoutParameters, ...]
    design_seed: int
    sampler: SamplerMethod
    parameter_space: ParameterSpace
    space_name: str

    @classmethod
    def sample(
        cls,
        count: int,
        seed: int,
        method: SamplerMethod = "sobol",
        space: ParameterSpace | None = None,
        space_name: str | None = None,
    ) -> "FactorialDesign":
        """Draw full parameters once so every model cell shares its core values."""
        selected_space = space or ParameterSpace.plausible_m4()
        sampler = QmcParameterSampler(selected_space)
        parameters = sampler.sample(count, seed, method)
        label = space_name or ("plausible_m4" if space is None else "custom")
        return cls(parameters, seed, method, selected_space, label)

    def cells(self) -> tuple[FactorialCell, ...]:
        """Expand each draw into exact-zero controls and active optional cells."""
        return tuple(
            FactorialCell(index, variant, self._parameters(item, variant))
            for index, item in enumerate(self.full_parameters)
            for variant in VARIANT_ORDER
        )

    @staticmethod
    def _parameters(
        parameters: ReadoutParameters,
        variant: ModelVariant,
    ) -> ReadoutParameters:
        """Apply only the two declared factorial switches."""
        return replace(
            parameters,
            echo_loading=parameters.echo_loading if variant.uses_echo else 0.0,
            orthogonal_curvature=(
                parameters.orthogonal_curvature if variant.uses_orthogonal else 0.0
            ),
        )


@dataclass(frozen=True, slots=True)
class FactorialExperimentConfig:
    """Independent cache seeds and between-seed confidence convention."""

    seeds: tuple[int, ...] = (104729, 130363, 155921)
    confidence_level: float = 0.95

    def __post_init__(self) -> None:
        """Require genuine replications rather than repeated labels."""
        if len(self.seeds) < 2 or len(set(self.seeds)) != len(self.seeds):
            raise ValueError("At least two distinct cache seeds are required.")
        if not 0.0 < self.confidence_level < 1.0:
            raise ValueError("confidence_level must lie in (0, 1).")


@dataclass(frozen=True, slots=True)
class FactorialExperimentResult:
    """Raw replications, seed summaries, and paired factorial contrasts."""

    runs: pd.DataFrame
    seed_summary: pd.DataFrame
    scenario_seed_summary: pd.DataFrame
    scenario_summary: pd.DataFrame
    criterion_seed_summary: pd.DataFrame
    criterion_summary: pd.DataFrame
    metric_summary: pd.DataFrame
    paired_effects: pd.DataFrame
    effect_summary: pd.DataFrame
    manifest: Mapping[str, object]

    def save(self, directory: str | Path) -> tuple[Path, ...]:
        """Persist all raw and aggregated evidence in auditable formats."""
        output = Path(directory)
        output.mkdir(parents=True, exist_ok=True)
        frames = {
            "runs.csv": self.runs,
            "seed_summary.csv": self.seed_summary,
            "scenario_seed_summary.csv": self.scenario_seed_summary,
            "scenario_summary.csv": self.scenario_summary,
            "criterion_seed_summary.csv": self.criterion_seed_summary,
            "criterion_summary.csv": self.criterion_summary,
            "metric_summary.csv": self.metric_summary,
            "paired_effects.csv": self.paired_effects,
            "effect_summary.csv": self.effect_summary,
        }
        paths = tuple(_atomic_csv(output / name, frame) for name, frame in frames.items())
        summary_path = _atomic_json(output / "manifest.json", dict(self.manifest))
        return (*paths, summary_path)


class FactorialExperimentRecorder:
    """Write a durable checkpoint after every completed simulation seed."""

    def __init__(self, directory: str | Path) -> None:
        """Reserve one output directory for a single experiment."""
        self.directory = Path(directory)
        self.checkpoints = self.directory / "checkpoints"
        self.completed_seeds: list[int] = []
        self.manifest: dict[str, object] = {}

    def start(self, manifest: Mapping[str, object]) -> None:
        """Create the manifest before the first expensive replication."""
        self.checkpoints.mkdir(parents=True, exist_ok=True)
        self.manifest = dict(manifest)
        progress = {**self.manifest, "status": "running", "completed_seeds": []}
        _atomic_json(self.directory / "progress.json", progress)

    def record_seed(self, seed: int, frame: pd.DataFrame) -> None:
        """Commit one complete seed and update the progress ledger."""
        _atomic_csv(self.checkpoints / f"seed_{seed}.csv", frame)
        self.completed_seeds.append(seed)
        progress = {
            **self.manifest,
            "status": "running",
            "completed_seeds": self.completed_seeds,
            "completed_seed_count": len(self.completed_seeds),
        }
        _atomic_json(self.directory / "progress.json", progress)

    def complete(self, result: FactorialExperimentResult) -> tuple[Path, ...]:
        """Write aggregate tables only after all declared seeds succeed."""
        paths = result.save(self.directory)
        progress = {
            **self.manifest,
            "status": "complete",
            "completed_seeds": self.completed_seeds,
            "completed_seed_count": len(self.completed_seeds),
        }
        _atomic_json(self.directory / "progress.json", progress)
        return paths


class MultiSeedFactorialExperiment:
    """Evaluate all four cells on common states, candidates, and seed replications."""

    def __init__(
        self,
        architecture: ArchitectureConfig,
        simulation: SimulationConfig,
        market: MarketEnvironment,
        criteria: AcceptanceCriteria,
        diagnostics: DiagnosticConfig | None = None,
        study_name: str = "M2/M3/M2+O/M4 factorial pilot",
        criteria_label: str = "illustrative proxy",
    ) -> None:
        """Fix the maximal feature engine and numerical conventions."""
        self._validate_architecture(architecture)
        self.architecture = architecture
        self.simulation = simulation
        self.market = market
        self.criteria = criteria
        self.diagnostics = StylisedFactDiagnostics(diagnostics)
        self.study_name = study_name
        self.criteria_label = criteria_label

    def run(
        self,
        structural: tuple[StructuralParameters, ...],
        design: FactorialDesign,
        config: FactorialExperimentConfig,
        recorder: FactorialExperimentRecorder | None = None,
    ) -> FactorialExperimentResult:
        """Run independent caches while retaining exact within-seed pairing."""
        if not structural:
            raise ValueError("At least one structural scenario is required.")
        manifest = self._manifest(structural, design, config)
        if recorder is not None:
            recorder.start(manifest)
        frames = []
        for seed in config.seeds:
            frame = self._run_seed(seed, structural, design)
            frames.append(frame)
            if recorder is not None:
                recorder.record_seed(seed, frame)
        runs = pd.concat(frames, ignore_index=True)
        result = FactorialAggregator(config.confidence_level).build(runs, manifest)
        if recorder is not None:
            recorder.complete(result)
        return result

    def _run_seed(
        self,
        seed: int,
        structural: tuple[StructuralParameters, ...],
        design: FactorialDesign,
    ) -> pd.DataFrame:
        """Draw one maximal cache and reuse it across every paired comparison."""
        base_factory = ReservoirFactory(self.architecture, self.simulation)
        randomness = base_factory.draw(seed)
        frames = []
        for scenario_id, scenario in enumerate(structural):
            architecture = scenario.architecture(self.architecture)
            reservoir = ReservoirFactory(architecture, self.simulation).from_randomness(
                randomness
            )
            model = ExogenousReservoirVolatilityModel(architecture)
            prepared = model.prepare(reservoir, scenario.risk_premium())
            frames.append(self._run_scenario(seed, scenario_id, scenario, model, prepared, design))
        return pd.concat(frames, ignore_index=True)

    def _run_scenario(
        self,
        seed: int,
        scenario_id: int,
        scenario: StructuralParameters,
        model: ExogenousReservoirVolatilityModel,
        prepared: PreparedPathCache,
        design: FactorialDesign,
    ) -> pd.DataFrame:
        """Evaluate all readouts without rebuilding any feature path."""
        rows: list[dict[str, object]] = []
        evaluated: dict[tuple[int, ReadoutParameters], dict[str, object]] = {}
        for cell in design.cells():
            key = (cell.candidate_id, cell.parameters)
            if key not in evaluated:
                evaluated[key] = self._evaluate_cell(
                    seed, scenario_id, scenario, cell, model, prepared
                )
            row = dict(evaluated[key])
            row["variant"] = cell.variant.value
            row["variant_order"] = VARIANT_ORDER.index(cell.variant)
            rows.append(row)
        return pd.DataFrame(rows)

    def _evaluate_cell(
        self,
        seed: int,
        scenario_id: int,
        scenario: StructuralParameters,
        cell: FactorialCell,
        model: ExogenousReservoirVolatilityModel,
        prepared: PreparedPathCache,
    ) -> dict[str, object]:
        """Simulate one cell and retain both estimates and pathwise uncertainty."""
        result = model.simulate(prepared, cell.parameters, self.market, measure="P")
        summary = self.diagnostics.evaluate(result)
        accepted, distance = self.criteria.evaluate(summary.mean)
        row: dict[str, object] = {
            "seed": seed,
            "structural_scenario_id": scenario_id,
            "candidate_id": cell.candidate_id,
            "variant": cell.variant.value,
            "variant_order": VARIANT_ORDER.index(cell.variant),
            "accepted": accepted,
            "distance": distance,
        }
        row.update(scenario.to_record())
        row.update(QmcParameterSampler.to_record(cell.parameters))
        row.update(summary.mean)
        row.update({f"se_{name}": value for name, value in summary.standard_error.items()})
        row.update(asdict(result.metadata))
        return row

    def _manifest(
        self,
        structural: tuple[StructuralParameters, ...],
        design: FactorialDesign,
        config: FactorialExperimentConfig,
    ) -> dict[str, object]:
        """Capture every choice needed to reproduce the result tables."""
        return {
            "schema_version": 3,
            "model_package_version": "0.4.0",
            "study_name": self.study_name,
            "seeds": list(config.seeds),
            "confidence_level": config.confidence_level,
            "design_seed": design.design_seed,
            "sampler": design.sampler,
            "candidate_count": len(design.full_parameters),
            "parameter_space": {
                "name": design.space_name,
                "bounds": {
                    name: asdict(bound)
                    for name, bound in design.parameter_space.bounds.items()
                },
                "fixed": dict(design.parameter_space.fixed),
            },
            "structural_scenarios": [item.to_record() for item in structural],
            "architecture": asdict(self.architecture),
            "simulation": asdict(self.simulation),
            "market": asdict(self.market),
            "acceptance_criteria": asdict(self.criteria),
            "acceptance_criteria_label": self.criteria_label,
            "acceptance_bands_are_illustrative": True,
        }

    @staticmethod
    def _validate_architecture(architecture: ArchitectureConfig) -> None:
        """Protect the interpretation of the two-by-two factorial design."""
        if not architecture.echo.enabled or not architecture.orthogonal.rates:
            raise ValueError("The maximal architecture needs echo and orthogonal features.")
        if "orthogonal" in architecture.echo.input_banks:
            raise ValueError("Echo inputs must exclude orthogonal factors for this ablation.")


class FactorialAggregator:
    """Aggregate fixed-design results across genuinely independent cache seeds."""

    def __init__(self, confidence_level: float) -> None:
        """Use a Student interval because the initial seed count is small."""
        self.confidence_level = confidence_level

    def build(
        self,
        runs: pd.DataFrame,
        manifest: Mapping[str, object],
    ) -> FactorialExperimentResult:
        """Create seed-level summaries before estimating Monte Carlo uncertainty."""
        seed_summary = self._seed_summary(runs)
        scenario_seed_summary, scenario_summary = self._scenario_summaries(runs)
        criterion_seed_summary, criterion_summary = self._criterion_summaries(
            runs, manifest
        )
        metric_summary = self._metric_summary(runs, seed_summary)
        paired_effects = self._paired_effects(runs)
        effect_summary = self._effect_summary(paired_effects)
        return FactorialExperimentResult(
            runs.sort_values(
                ["seed", "structural_scenario_id", "candidate_id", "variant_order"]
            ).reset_index(drop=True),
            seed_summary,
            scenario_seed_summary,
            scenario_summary,
            criterion_seed_summary,
            criterion_summary,
            metric_summary,
            paired_effects,
            effect_summary,
            manifest,
        )

    @staticmethod
    def _seed_summary(runs: pd.DataFrame) -> pd.DataFrame:
        """Average the fixed parameter design separately within each cache seed."""
        rows = []
        for (seed, variant), frame in runs.groupby(["seed", "variant"], sort=False):
            row: dict[str, object] = {
                "seed": int(seed),
                "variant": variant,
                "variant_order": int(frame["variant_order"].iloc[0]),
                "candidate_evaluations": len(frame),
                "acceptance_rate": float(frame["accepted"].mean()),
                "mean_distance": float(frame["distance"].mean()),
            }
            row.update({name: float(frame[name].mean()) for name in COMPARISON_METRICS})
            rows.append(row)
        return pd.DataFrame(rows).sort_values(["seed", "variant_order"]).reset_index(drop=True)

    def _scenario_summaries(
        self,
        runs: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Separate structural-scenario variation from cache-seed variation."""
        metrics = {"acceptance_rate": "accepted", "distance": "distance"}
        metrics.update({name: name for name in COMPARISON_METRICS})
        seed_rows = []
        groups = runs.groupby(["seed", "structural_scenario_id", "variant"], sort=False)
        for (seed, scenario, variant), frame in groups:
            row: dict[str, object] = {
                "seed": int(seed),
                "structural_scenario_id": int(scenario),
                "variant": variant,
            }
            row.update({name: float(frame[column].mean()) for name, column in metrics.items()})
            seed_rows.append(row)
        seed_summary = pd.DataFrame(seed_rows)
        aggregate_rows = []
        for (scenario, variant), frame in seed_summary.groupby(
            ["structural_scenario_id", "variant"], sort=False
        ):
            for metric in metrics:
                aggregate_rows.append(
                    {
                        "structural_scenario_id": int(scenario),
                        "variant": variant,
                        "metric": metric,
                        **self._statistics(frame[metric].to_numpy(dtype=float)),
                    }
                )
        return seed_summary, pd.DataFrame(aggregate_rows)

    def _criterion_summaries(
        self,
        runs: pd.DataFrame,
        manifest: Mapping[str, object],
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Expose individual gate pass rates behind the joint acceptance result."""
        criteria = cast(Mapping[str, object], manifest["acceptance_criteria"])
        bands = cast(Mapping[str, Mapping[str, float]], criteria["bands"])
        seed_rows = []
        for (seed, variant), frame in runs.groupby(["seed", "variant"], sort=False):
            for metric, band in bands.items():
                passed = frame[metric].between(
                    float(band["lower"]), float(band["upper"]), inclusive="both"
                )
                seed_rows.append(
                    {
                        "seed": int(seed),
                        "variant": variant,
                        "criterion": metric,
                        "pass_rate": float(passed.mean()),
                    }
                )
        seed_summary = pd.DataFrame(seed_rows)
        aggregate_rows = []
        for (variant, criterion), frame in seed_summary.groupby(
            ["variant", "criterion"], sort=False
        ):
            aggregate_rows.append(
                {
                    "variant": variant,
                    "criterion": criterion,
                    **self._statistics(frame.pass_rate.to_numpy(dtype=float)),
                }
            )
        return seed_summary, pd.DataFrame(aggregate_rows)

    def _metric_summary(
        self,
        runs: pd.DataFrame,
        seed_summary: pd.DataFrame,
    ) -> pd.DataFrame:
        """Report seed uncertainty and pooled candidate-distribution quantiles."""
        mappings = {
            "acceptance_rate": "accepted",
            "mean_distance": "distance",
            **{name: name for name in COMPARISON_METRICS},
        }
        rows = []
        for variant in (item.value for item in VARIANT_ORDER):
            seed_frame = seed_summary[seed_summary.variant == variant]
            pooled = runs[runs.variant == variant]
            for summary_name, raw_name in mappings.items():
                stats = self._statistics(seed_frame[summary_name].to_numpy(dtype=float))
                raw = pooled[raw_name].to_numpy(dtype=float)
                rows.append(
                    {
                        "variant": variant,
                        "metric": summary_name,
                        **stats,
                        "pooled_q10": _finite_quantile(raw, 0.10),
                        "pooled_median": _finite_quantile(raw, 0.50),
                        "pooled_q90": _finite_quantile(raw, 0.90),
                    }
                )
        return pd.DataFrame(rows)

    @staticmethod
    def _paired_effects(runs: pd.DataFrame) -> pd.DataFrame:
        """Compute both conditional main effects and the factorial interaction."""
        index = ["seed", "structural_scenario_id", "candidate_id"]
        metrics = ("accepted", "distance", *COMPARISON_METRICS)
        effects = {
            "echo_at_O0": ("M3", "M2", None, None),
            "echo_at_O1": ("M4", "M2+O", None, None),
            "orthogonal_at_echo0": ("M2+O", "M2", None, None),
            "orthogonal_at_echo1": ("M4", "M3", None, None),
            "echo_orthogonal_interaction": ("M4", "M3", "M2+O", "M2"),
        }
        rows = []
        for metric in metrics:
            pivot = runs.pivot(index=index, columns="variant", values=metric).reset_index()
            for effect, labels in effects.items():
                values = pivot[labels[0]].astype(float) - pivot[labels[1]].astype(float)
                if labels[2] is not None and labels[3] is not None:
                    values = (
                        values
                        - pivot[labels[2]].astype(float)
                        + pivot[labels[3]].astype(float)
                    )
                name = "acceptance_rate" if metric == "accepted" else metric
                rows.extend(
                    {
                        **{key: pivot.iloc[row][key] for key in index},
                        "effect": effect,
                        "metric": name,
                        "value": float(values.iloc[row]),
                    }
                    for row in range(len(pivot))
                )
        return pd.DataFrame(rows)

    def _effect_summary(self, effects: pd.DataFrame) -> pd.DataFrame:
        """Estimate paired effects from seed-level averages, not pseudo-replicates."""
        per_seed = (
            effects.groupby(["seed", "effect", "metric"], as_index=False)["value"].mean()
        )
        rows = []
        for (effect, metric), frame in per_seed.groupby(["effect", "metric"], sort=False):
            values = frame["value"].to_numpy(dtype=float)
            pooled = effects[
                (effects.effect == effect) & (effects.metric == metric)
            ].value.to_numpy(dtype=float)
            rows.append(
                {
                    "effect": effect,
                    "metric": metric,
                    **self._statistics(values),
                    "pooled_q10": _finite_quantile(pooled, 0.10),
                    "pooled_median": _finite_quantile(pooled, 0.50),
                    "pooled_q90": _finite_quantile(pooled, 0.90),
                }
            )
        return pd.DataFrame(rows)

    def _statistics(self, values: np.ndarray) -> dict[str, float | int]:
        """Summarise independent seed means with a small-sample t interval."""
        finite = values[np.isfinite(values)]
        count = int(finite.size)
        mean = float(np.mean(finite)) if count else float("nan")
        standard_deviation = float(np.std(finite, ddof=1)) if count > 1 else float("nan")
        standard_error = standard_deviation / np.sqrt(count) if count > 1 else float("nan")
        critical = float(
            student_t.ppf((1.0 + self.confidence_level) / 2.0, count - 1)
        ) if count > 1 else float("nan")
        radius = critical * standard_error
        return {
            "seed_count": count,
            "mean": mean,
            "between_seed_sd": standard_deviation,
            "between_seed_se": standard_error,
            "ci_lower": mean - radius,
            "ci_upper": mean + radius,
            "seed_min": float(np.min(finite)) if count else float("nan"),
            "seed_max": float(np.max(finite)) if count else float("nan"),
        }


def _finite_quantile(values: np.ndarray, probability: float) -> float:
    """Return a quantile after excluding failed numerical estimates."""
    finite = values[np.isfinite(values)]
    return float(np.quantile(finite, probability)) if finite.size else float("nan")


def _atomic_csv(path: Path, frame: pd.DataFrame) -> Path:
    """Replace a CSV only after its complete temporary representation exists."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)
    return path


def _atomic_json(path: Path, payload: Mapping[str, object]) -> Path:
    """Replace a JSON ledger atomically to keep interrupted runs readable."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)
    return path
