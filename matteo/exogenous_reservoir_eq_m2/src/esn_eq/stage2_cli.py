"""Run the independent targeted Stage-2 physical-measure validation study."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, replace
from pathlib import Path
from typing import Mapping, cast

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

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
    FactorialExperimentResult,
    MultiSeedFactorialExperiment,
)
from .reporting import FactorialMarkdownReport
from .sampling import (
    ParameterBound,
    ParameterSpace,
    QmcStructuralSampler,
    StructuralParameters,
)


KEY_METRICS = (
    "rough_hurst",
    "leverage_rank",
    "zumbach_rank",
    "excess_kurtosis",
    "max_abs_return_acf",
    "mean_log_variance_acf",
    "p_cap_exceedance_fraction",
)

READOUT_PARAMETERS = (
    "cluster_loading",
    "feedback_curvature",
    "feedback_shift",
    "spike_curvature",
    "orthogonal_curvature",
    "echo_loading",
    "cap_level",
    "cap_sharpness",
)


def build_parser() -> argparse.ArgumentParser:
    """Expose an independent validation design with conservative defaults."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("experiments/stage2_001"))
    parser.add_argument("--pilot", type=Path, default=Path("experiments/pilot_001"))
    parser.add_argument("--candidate-samples", type=int, default=32)
    parser.add_argument("--structural-samples", type=int, default=4)
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=[196613, 221251, 262147, 294001, 324503],
    )
    parser.add_argument("--design-seed", type=int, default=20261021)
    parser.add_argument("--structural-seed", type=int, default=20261022)
    parser.add_argument("--years", type=float, default=5.0)
    parser.add_argument("--burn-years", type=float, default=1.0)
    parser.add_argument("--paths", type=int, default=16)
    parser.add_argument("--steps-per-observation", type=int, default=2)
    return parser


def main() -> None:
    """Run checkpoints, comparisons, sensitivities, and a self-contained report."""
    args = build_parser().parse_args()
    latent_record = _load_json(args.pilot / "latent_calibration.json")
    architecture = _stage2_architecture(latent_record)
    criteria = AcceptanceCriteria.stage2_proxy()
    structures = _targeted_structures(args.structural_samples, args.structural_seed)
    design = FactorialDesign.sample(
        args.candidate_samples,
        args.design_seed,
        method="sobol",
        space=ParameterSpace.targeted_stage2(),
        space_name="targeted_stage2",
    )
    result = _run(args, architecture, criteria, structures, design)
    _save_supplements(args, result, design, structures, criteria, latent_record)
    print(f"Completed {len(args.seeds)} seeds and {len(result.runs)} paired rows.")
    print(args.output / "RESULTS.md")


def _run(
    args: argparse.Namespace,
    architecture: ArchitectureConfig,
    criteria: AcceptanceCriteria,
    structures: tuple[StructuralParameters, ...],
    design: FactorialDesign,
) -> FactorialExperimentResult:
    """Construct and execute the maximal paired feature engine."""
    simulation = SimulationConfig(
        years=args.years,
        burn_years=args.burn_years,
        paths=args.paths,
        steps_per_observation=args.steps_per_observation,
    )
    runner = MultiSeedFactorialExperiment(
        architecture,
        simulation,
        MarketEnvironment(),
        criteria,
        study_name="Targeted Stage-2 M2/M3/M2+O/M4 validation",
        criteria_label="stage2 synthetic proxy",
    )
    recorder = FactorialExperimentRecorder(args.output)
    return runner.run(
        structures,
        design,
        FactorialExperimentConfig(tuple(args.seeds)),
        recorder,
    )


def _stage2_architecture(record: Mapping[str, object]) -> ArchitectureConfig:
    """Freeze the Pilot-1 lift and change only the declared O-feature scaling."""
    raw = cast(Mapping[str, object], record["architecture"])
    base = ArchitectureConfig.from_mapping(raw)
    return replace(
        base,
        orthogonal=BankConfig(geometric_rates(0.5, 100.0, 12)),
        echo=replace(base.echo, enabled=True),
        orthogonal_energy_mode="standardized_centered",
    )


def _targeted_structures(count: int, seed: int) -> tuple[StructuralParameters, ...]:
    """Vary only P-measure couplings; keep all risk-premium entries fixed."""
    bounds = {
        "feedback_asset_correlation": ParameterBound(0.75, 0.99),
        "spike_asset_correlation": ParameterBound(0.55, 0.95),
    }
    return QmcStructuralSampler(bounds).sample(count, seed)


def _save_supplements(
    args: argparse.Namespace,
    result: FactorialExperimentResult,
    design: FactorialDesign,
    structures: tuple[StructuralParameters, ...],
    criteria: AcceptanceCriteria,
    latent_record: Mapping[str, object],
) -> None:
    """Persist validation-only tables and append their interpretation."""
    output = args.output
    pilot_runs = pd.read_csv(args.pilot / "runs.csv")
    comparison = _baseline_comparison(pilot_runs, result.runs, criteria)
    stability = _candidate_stability(result.runs)
    sensitivities = _rank_sensitivities(result.runs)
    _atomic_csv(output / "pilot_stage2_comparison.csv", comparison)
    _atomic_csv(output / "candidate_stability.csv", stability)
    _atomic_csv(output / "rank_sensitivities.csv", sensitivities)
    _atomic_csv(output / "readout_design.csv", _design_frame(design))
    _atomic_csv(output / "structural_design.csv", _structure_frame(structures))
    _atomic_json(output / "frozen_latent_source.json", _latent_source(latent_record))
    core = FactorialMarkdownReport().build(result, latent_record)
    extra = _validation_sections(result, comparison, stability, args.pilot)
    _atomic_text(output / "RESULTS.md", core.rstrip() + "\n\n" + extra)


def _baseline_comparison(
    pilot: pd.DataFrame,
    stage2: pd.DataFrame,
    criteria: AcceptanceCriteria,
) -> pd.DataFrame:
    """Re-evaluate both studies under exactly the same Stage-2 proxy gates."""
    rows = []
    for study, source in (("Pilot 1", pilot), ("Stage 2", stage2)):
        frame = source.copy()
        frame["proxy_accepted"] = _criteria_mask(frame, criteria)
        per_seed = _per_seed_metrics(frame)
        for (variant, metric), group in per_seed.groupby(["variant", "metric"]):
            values = group["value"].to_numpy(dtype=float)
            rows.append(
                {
                    "study": study,
                    "variant": variant,
                    "metric": metric,
                    "seed_count": len(values),
                    "mean": float(np.mean(values)),
                    "between_seed_sd": float(np.std(values, ddof=1)),
                    "seed_min": float(np.min(values)),
                    "seed_max": float(np.max(values)),
                }
            )
    return pd.DataFrame(rows)


def _criteria_mask(frame: pd.DataFrame, criteria: AcceptanceCriteria) -> np.ndarray:
    """Vectorise the strict conjunction of all declared proxy bands."""
    mask = np.ones(len(frame), dtype=bool)
    for name, band in criteria.bands.items():
        values = frame[name].to_numpy(dtype=float)
        mask &= np.isfinite(values) & (values >= band.lower) & (values <= band.upper)
    return mask


def _per_seed_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    """Average a fixed design before treating cache seeds as replications."""
    rows = []
    for (seed, variant), group in frame.groupby(["seed", "variant"], sort=False):
        values = {"stage2_proxy_acceptance": float(group.proxy_accepted.mean())}
        values.update({name: float(group[name].mean()) for name in KEY_METRICS})
        rows.extend(
            {"seed": int(seed), "variant": variant, "metric": name, "value": value}
            for name, value in values.items()
        )
    return pd.DataFrame(rows)


def _candidate_stability(runs: pd.DataFrame) -> pd.DataFrame:
    """Measure whether the same configuration survives independent cache seeds."""
    keys = ["structural_scenario_id", "candidate_id", "variant"]
    columns = ["accepted", "distance", *KEY_METRICS]
    stable = runs.groupby(keys, as_index=False)[columns].mean()
    stable = stable.rename(
        columns={"accepted": "seed_pass_fraction", "distance": "mean_distance"}
    )
    stable["passes_majority"] = stable.seed_pass_fraction >= 0.60
    stable["passes_80pct"] = stable.seed_pass_fraction >= 0.80
    stable["passes_all_seeds"] = stable.seed_pass_fraction >= 1.0
    return stable.sort_values(
        ["variant", "seed_pass_fraction", "mean_distance"],
        ascending=[True, False, True],
    ).reset_index(drop=True)


def _rank_sensitivities(runs: pd.DataFrame) -> pd.DataFrame:
    """Screen readout and structural directions after averaging replications."""
    metrics = ("accepted", "distance", *KEY_METRICS[:4])
    rows = _readout_sensitivities(runs, metrics)
    rows.extend(_structural_sensitivities(runs, metrics))
    return pd.DataFrame(rows).sort_values(
        ["parameter_class", "variant", "metric", "parameter"]
    ).reset_index(drop=True)


def _readout_sensitivities(
    runs: pd.DataFrame,
    metrics: tuple[str, ...],
) -> list[dict[str, object]]:
    """Estimate candidate-level rank screens without pseudo-replicating seeds."""
    aggregate = runs.groupby(["variant", "candidate_id"], as_index=False)[
        [*READOUT_PARAMETERS, *metrics]
    ].mean()
    rows = []
    for variant, frame in aggregate.groupby("variant"):
        active = [name for name in READOUT_PARAMETERS if frame[name].nunique() > 1]
        rows.extend(_correlation_rows(frame, variant, "readout", active, metrics))
    return rows


def _structural_sensitivities(
    runs: pd.DataFrame,
    metrics: tuple[str, ...],
) -> list[dict[str, object]]:
    """Estimate scenario-level rank screens for the two varied correlations."""
    parameters = ("feedback_asset_correlation", "spike_asset_correlation")
    aggregate = runs.groupby(["variant", "structural_scenario_id"], as_index=False)[
        [*parameters, *metrics]
    ].mean()
    rows = []
    for variant, frame in aggregate.groupby("variant"):
        rows.extend(_correlation_rows(frame, variant, "structural", parameters, metrics))
    return rows


def _correlation_rows(
    frame: pd.DataFrame,
    variant: str,
    parameter_class: str,
    parameters: tuple[str, ...] | list[str],
    metrics: tuple[str, ...],
) -> list[dict[str, object]]:
    """Return Spearman coefficients and p-values for one aggregation level."""
    rows = []
    for parameter in parameters:
        for metric in metrics:
            if frame[parameter].nunique() < 2 or frame[metric].nunique() < 2:
                statistic, pvalue = float("nan"), float("nan")
            else:
                estimate = spearmanr(frame[parameter], frame[metric])
                statistic, pvalue = float(estimate.statistic), float(estimate.pvalue)
            rows.append(
                {
                    "parameter_class": parameter_class,
                    "variant": variant,
                    "parameter": parameter,
                    "metric": metric,
                    "spearman": statistic,
                    "pvalue": pvalue,
                    "sample_count": len(frame),
                }
            )
    return rows


def _validation_sections(
    result: FactorialExperimentResult,
    comparison: pd.DataFrame,
    stability: pd.DataFrame,
    pilot: Path,
) -> str:
    """Explain the independent validation comparison and seed robustness."""
    lines = [
        "## Pilot-to-Stage-2 validation comparison",
        "",
        "Pilot 1 was used only to choose the targeted box. Stage 2 uses a new scrambled "
        "design, new structural design, and disjoint cache seeds. Both studies below are "
        "re-scored under the same Stage-2 proxy gates. Because the sampling boxes differ, "
        "this compares targeted capacity, not global parameter-space volume.",
        "",
        _comparison_table(comparison),
        "",
        "## Configuration stability across Stage-2 seeds",
        "",
        "A configuration is one structural scenario, readout candidate, and model cell. "
        "The rates below show how often the same configuration passes across seeds; this "
        "is stricter than pooling all rows.",
        "",
        _stability_table(stability),
        "",
        "## Stage-2 reading",
        "",
        _decision_text(result, stability),
        "",
        "## Selection boundary",
        "",
        f"The frozen latent architecture came from `{pilot}` and was not refitted. Stage-2 "
        "results may guide a later experiment, but any newly selected subregion or best "
        "candidate requires another untouched design before it can be called validated. "
        "Option-surface calibration remains untested because no SPX/VIX option data or "
        "pricing loss was supplied.",
    ]
    return "\n".join(lines).rstrip() + "\n"


def _comparison_table(comparison: pd.DataFrame) -> str:
    """Display same-gate acceptance and four principal diagnostics."""
    metrics = (
        "stage2_proxy_acceptance",
        "rough_hurst",
        "leverage_rank",
        "zumbach_rank",
        "excess_kurtosis",
    )
    headers = ["model", *metrics]
    rows = []
    for variant in ("M2", "M3", "M2+O", "M4"):
        row = [variant]
        for metric in metrics:
            values = comparison[
                (comparison.variant == variant) & (comparison.metric == metric)
            ].set_index("study")
            pilot = _compact(values.loc["Pilot 1", "mean"])
            stage2 = _compact(values.loc["Stage 2", "mean"])
            row.append(f"{pilot} → {stage2}")
        rows.append(row)
    return _markdown_table(headers, rows)


def _stability_table(stability: pd.DataFrame) -> str:
    """Summarise stable configuration fractions for each factorial cell."""
    rows = []
    for variant in ("M2", "M3", "M2+O", "M4"):
        frame = stability[stability.variant == variant]
        rows.append(
            [
                variant,
                str(len(frame)),
                _percent(frame.passes_majority.mean()),
                _percent(frame.passes_80pct.mean()),
                _percent(frame.passes_all_seeds.mean()),
                _percent(frame.seed_pass_fraction.max()),
            ]
        )
    headers = ["model", "configs", "pass ≥3/5", "pass ≥4/5", "pass 5/5", "best rate"]
    return _markdown_table(headers, rows)


def _decision_text(
    result: FactorialExperimentResult,
    stability: pd.DataFrame,
) -> str:
    """Apply the preregistered direction-of-improvement rules conservatively."""
    metrics = result.metric_summary.set_index(["variant", "metric"])
    effects = result.effect_summary.set_index(["effect", "metric"])
    acceptance = {
        variant: float(metrics.loc[(variant, "acceptance_rate"), "mean"])
        for variant in ("M2", "M3", "M2+O", "M4")
    }
    best = max(acceptance, key=acceptance.get)
    orthogonal = effects.loc[("orthogonal_at_echo0", "acceptance_rate")]
    echo = effects.loc[("echo_at_O0", "acceptance_rate")]
    stable = stability.groupby("variant").passes_80pct.mean()
    return "\n".join(
        (
            f"- `{best}` has the largest mean Stage-2 proxy pass rate at "
            f"`{_percent(acceptance[best])}`; the corresponding ≥4/5-seed stable-configuration "
            f"rate is `{_percent(stable[best])}`.",
            f"- Adding standardized orthogonal energy to M2 changes acceptance by "
            f"`{_compact(orthogonal['mean'])}` with a seed interval "
            f"`[{_compact(orthogonal['ci_lower'])}, {_compact(orthogonal['ci_upper'])}]`.",
            f"- Adding the echo to M2 changes acceptance by `{_compact(echo['mean'])}` "
            f"with a seed interval `[{_compact(echo['ci_lower'])}, "
            f"{_compact(echo['ci_upper'])}]`.",
            "- Retain an optional block only if its paired benefit is directionally stable "
            "and its leverage, Zumbach, roughness, tail, return-ACF, and cap diagnostics do "
            "not reveal a compensating deterioration.",
        )
    )


def _design_frame(design: FactorialDesign) -> pd.DataFrame:
    """Expose the exact full-M4 Sobol points before factorial zeroing."""
    frame = pd.DataFrame([asdict(item) for item in design.full_parameters])
    frame.insert(0, "candidate_id", np.arange(len(frame)))
    return frame


def _structure_frame(structures: tuple[StructuralParameters, ...]) -> pd.DataFrame:
    """Expose targeted correlation scenarios and fixed risk-premium entries."""
    frame = pd.DataFrame([item.to_record() for item in structures])
    frame.insert(0, "structural_scenario_id", np.arange(len(frame)))
    return frame


def _latent_source(record: Mapping[str, object]) -> dict[str, object]:
    """Record what was frozen without pretending to recalibrate the lift."""
    return {
        "source": "Pilot 1 latent calibration",
        "selected_seed": record.get("selected_seed"),
        "objective": record.get("objective"),
        "boundary_hits": record.get("boundary_hits", []),
        "architecture": record.get("architecture"),
    }


def _load_json(path: Path) -> dict[str, object]:
    """Load a required JSON object with an explicit shape check."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}.")
    return payload


def _markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    """Render a compact CommonMark table."""
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def _compact(value: object) -> str:
    """Format rates and diagnostics without concealing small magnitudes."""
    number = float(value)
    return f"{number:.4g}" if np.isfinite(number) else "NA"


def _percent(value: object) -> str:
    """Format one fraction in percent units."""
    number = float(value)
    return f"{100.0 * number:.1f}%" if np.isfinite(number) else "NA"


def _atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    """Write a table completely before replacing its public path."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def _atomic_json(path: Path, payload: Mapping[str, object]) -> None:
    """Write a JSON record completely before replacing its public path."""
    _atomic_text(path, json.dumps(payload, indent=2) + "\n")


def _atomic_text(path: Path, content: str) -> None:
    """Write UTF-8 text completely before replacing its public path."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


if __name__ == "__main__":
    main()
