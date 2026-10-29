"""Confirm the slower-spike/high-shift M2 ridge on untouched simulations."""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path
from typing import Mapping, cast

import numpy as np
import pandas as pd

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
from .stage2_cli import (
    KEY_METRICS,
    _atomic_csv,
    _atomic_json,
    _atomic_text,
    _candidate_stability,
    _compact,
    _criteria_mask,
    _design_frame,
    _load_json,
    _markdown_table,
    _per_seed_metrics,
    _percent,
    _rank_sensitivities,
    _structure_frame,
)


def build_parser() -> argparse.ArgumentParser:
    """Expose fresh validation seeds and a more precise path ensemble."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("experiments/stage3_001"))
    parser.add_argument("--pilot", type=Path, default=Path("experiments/pilot_001"))
    parser.add_argument("--training", type=Path, default=Path("experiments/stage2_001"))
    parser.add_argument("--candidate-samples", type=int, default=32)
    parser.add_argument("--structural-samples", type=int, default=4)
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=[360007, 390001, 420001, 450007, 480013],
    )
    parser.add_argument("--design-seed", type=int, default=20261201)
    parser.add_argument("--structural-seed", type=int, default=20261202)
    parser.add_argument("--years", type=float, default=6.0)
    parser.add_argument("--burn-years", type=float, default=1.0)
    parser.add_argument("--paths", type=int, default=24)
    parser.add_argument("--steps-per-observation", type=int, default=2)
    return parser


def main() -> None:
    """Run the frozen confirmation study and its robustness summaries."""
    args = build_parser().parse_args()
    latent_record = _load_json(args.pilot / "latent_calibration.json")
    architecture = _stage3_architecture(latent_record)
    structures = _targeted_structures(args.structural_samples, args.structural_seed)
    design = FactorialDesign.sample(
        args.candidate_samples,
        args.design_seed,
        method="sobol",
        space=ParameterSpace.stage3_local(),
        space_name="stage3_local_slow_spike_ridge",
    )
    criteria = AcceptanceCriteria.stage2_proxy()
    result = _run(args, architecture, structures, design, criteria)
    _save_supplements(args, result, design, structures, criteria, latent_record)
    print(f"Completed {len(args.seeds)} seeds and {len(result.runs)} paired rows.")
    print(args.output / "RESULTS.md")


def _run(
    args: argparse.Namespace,
    architecture: ArchitectureConfig,
    structures: tuple[StructuralParameters, ...],
    design: FactorialDesign,
    criteria: AcceptanceCriteria,
) -> FactorialExperimentResult:
    """Execute a checkpointed study with duplicate zero-echo cells memoized."""
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
        study_name="Stage-3 confirmation of the slower-spike/high-shift ridge",
        criteria_label="stage2 synthetic proxy retained unchanged",
    )
    return runner.run(
        structures,
        design,
        FactorialExperimentConfig(tuple(args.seeds)),
        FactorialExperimentRecorder(args.output),
    )


def _stage3_architecture(record: Mapping[str, object]) -> ArchitectureConfig:
    """Apply the selected slow spike lift to the otherwise frozen Pilot-1 lift."""
    raw = cast(Mapping[str, object], record["architecture"])
    base = ArchitectureConfig.from_mapping(raw)
    return replace(
        base,
        spike=replace(base.spike, rates=geometric_rates(10.0, 500.0, 4)),
        spike_weight_power=-0.5,
        orthogonal=BankConfig(geometric_rates(0.5, 100.0, 12)),
        echo=replace(base.echo, enabled=True),
        orthogonal_energy_mode="standardized_centered",
    )


def _targeted_structures(count: int, seed: int) -> tuple[StructuralParameters, ...]:
    """Validate a small high-correlation neighborhood with fixed risk premia."""
    bounds = {
        "feedback_asset_correlation": ParameterBound(0.970, 0.999),
        "spike_asset_correlation": ParameterBound(0.940, 0.995),
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
    """Persist stability, strong-Zumbach, comparison, and selection records."""
    stability = _candidate_stability(result.runs)
    top = _top_configurations(stability, design, structures)
    strong = _strong_zumbach_summary(result.runs, criteria)
    comparison = _study_comparison(args.training / "runs.csv", result.runs, criteria)
    _atomic_csv(args.output / "candidate_stability.csv", stability)
    _atomic_csv(args.output / "top_configurations.csv", top)
    _atomic_csv(args.output / "strong_zumbach_sensitivity.csv", strong)
    _atomic_csv(args.output / "stage2_stage3_comparison.csv", comparison)
    _atomic_csv(args.output / "rank_sensitivities.csv", _rank_sensitivities(result.runs))
    _atomic_csv(args.output / "readout_design.csv", _design_frame(design))
    _atomic_csv(args.output / "structural_design.csv", _structure_frame(structures))
    _atomic_json(args.output / "selection_record.json", _selection_record(latent_record))
    core = FactorialMarkdownReport().build(result, latent_record)
    supplement = _report_supplement(result, stability, top, strong, comparison)
    _atomic_text(args.output / "RESULTS.md", core.rstrip() + "\n\n" + supplement)


def _top_configurations(
    stability: pd.DataFrame,
    design: FactorialDesign,
    structures: tuple[StructuralParameters, ...],
) -> pd.DataFrame:
    """Attach exact readout and structural values to seed-stability rankings."""
    design_frame = _design_frame(design)
    structure_frame = _structure_frame(structures)
    merged = stability.merge(design_frame, on="candidate_id", how="left")
    merged = merged.merge(structure_frame, on="structural_scenario_id", how="left")
    return merged.sort_values(
        ["variant", "seed_pass_fraction", "mean_distance"],
        ascending=[True, False, True],
    ).reset_index(drop=True)


def _strong_zumbach_summary(
    runs: pd.DataFrame,
    criteria: AcceptanceCriteria,
) -> pd.DataFrame:
    """Report stronger asymmetry checks without redefining the primary gate."""
    rows = []
    for (seed, variant), frame in runs.groupby(["seed", "variant"], sort=False):
        other = np.ones(len(frame), dtype=bool)
        for name, band in criteria.bands.items():
            if name == "zumbach_rank":
                continue
            values = frame[name].to_numpy(dtype=float)
            other &= np.isfinite(values) & (values >= band.lower) & (values <= band.upper)
        checks = {
            "Zumbach > 0.03": frame.zumbach_rank.to_numpy(dtype=float) > 0.03,
            "Zumbach > 0.05": frame.zumbach_rank.to_numpy(dtype=float) > 0.05,
            "joint with Zumbach > 0.03": other
            & (frame.zumbach_rank.to_numpy(dtype=float) > 0.03),
            "joint with Zumbach > 0.05": other
            & (frame.zumbach_rank.to_numpy(dtype=float) > 0.05),
        }
        rows.extend(
            {"seed": int(seed), "variant": variant, "check": name, "pass_rate": mask.mean()}
            for name, mask in checks.items()
        )
    per_seed = pd.DataFrame(rows)
    return per_seed.groupby(["variant", "check"], as_index=False).pass_rate.agg(
        mean="mean", between_seed_sd="std", seed_min="min", seed_max="max"
    )


def _study_comparison(
    training_path: Path,
    validation: pd.DataFrame,
    criteria: AcceptanceCriteria,
) -> pd.DataFrame:
    """Compare targeted studies under one unchanged proxy definition."""
    training = pd.read_csv(training_path)
    rows = []
    for study, source in (("Stage 2", training), ("Stage 3", validation)):
        frame = source.copy()
        frame["proxy_accepted"] = _criteria_mask(frame, criteria)
        per_seed = _per_seed_metrics(frame)
        for (variant, metric), group in per_seed.groupby(["variant", "metric"]):
            values = group.value.to_numpy(dtype=float)
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


def _selection_record(record: Mapping[str, object]) -> dict[str, object]:
    """Freeze the training evidence that defined Stage 3 before it ran."""
    return {
        "source_latent_seed": record.get("selected_seed"),
        "selected_spike_rates": list(geometric_rates(10.0, 500.0, 4)),
        "selected_spike_weight_power": -0.5,
        "training_anchor": {
            "cluster_loading": 0.24,
            "feedback_curvature": 0.073,
            "feedback_shift": 3.0,
            "spike_curvature": 0.22,
            "five_seed_pass_fraction": 0.8,
        },
        "echo_loading_fixed": 0.0,
        "primary_criteria_changed_after_training": False,
    }


def _report_supplement(
    result: FactorialExperimentResult,
    stability: pd.DataFrame,
    top: pd.DataFrame,
    strong: pd.DataFrame,
    comparison: pd.DataFrame,
) -> str:
    """Render the confirmation-specific evidence and decision boundary."""
    return "\n".join(
        (
            "## What Stage 3 changed",
            "",
            "Stage 2 selected a slower spike grid (10–500 inverse years), spike projection "
            "power -0.5, and a high-shift/low-feedback-curvature neighborhood. Stage 3 "
            "freezes those choices, uses a new Sobol design and disjoint cache seeds, and "
            "raises the Q-normalisation ensemble to 24 paths. The echo loading is fixed at "
            "zero, so M3 duplicates M2 and M4 duplicates M2+O by construction.",
            "",
            "## Same-gate Stage-2/Stage-3 comparison",
            "",
            _comparison_table(comparison),
            "",
            "## Stability of the same configuration across seeds",
            "",
            _stability_table(stability),
            "",
            "## Strong-Zumbach sensitivity",
            "",
            _strong_table(strong),
            "",
            "## Best untouched configurations",
            "",
            _top_table(top),
            "",
            "## Confirmation reading",
            "",
            _decision_text(result, stability, top),
            "",
            "The gates remain synthetic capacity proxies. They do not replace empirical "
            "estimation, option-surface calibration, conditional VIX pricing validation, "
            "or an out-of-time market-data test.",
        )
    ).rstrip() + "\n"


def _comparison_table(comparison: pd.DataFrame) -> str:
    """Show M2 and M2+O evolution on the key facts."""
    metrics = (
        "stage2_proxy_acceptance",
        "rough_hurst",
        "leverage_rank",
        "zumbach_rank",
        "excess_kurtosis",
    )
    rows = []
    for variant in ("M2", "M2+O"):
        row = [variant]
        for metric in metrics:
            subset = comparison[
                (comparison.variant == variant) & (comparison.metric == metric)
            ].set_index("study")
            row.append(
                f"{_compact(subset.loc['Stage 2', 'mean'])} → "
                f"{_compact(subset.loc['Stage 3', 'mean'])}"
            )
        rows.append(row)
    return _markdown_table(["model", *metrics], rows)


def _stability_table(stability: pd.DataFrame) -> str:
    """Show neighborhood volume that survives most or all cache seeds."""
    rows = []
    for variant in ("M2", "M2+O"):
        frame = stability[stability.variant == variant]
        rows.append(
            [
                variant,
                str(len(frame)),
                _percent((frame.seed_pass_fraction >= 0.60).mean()),
                _percent((frame.seed_pass_fraction >= 0.80).mean()),
                _percent((frame.seed_pass_fraction >= 1.00).mean()),
                _percent(frame.seed_pass_fraction.max()),
            ]
        )
    return _markdown_table(
        ["model", "configs", "pass ≥3/5", "pass ≥4/5", "pass 5/5", "best rate"],
        rows,
    )


def _strong_table(strong: pd.DataFrame) -> str:
    """Show both marginal and joint stronger-Zumbach rates."""
    checks = list(strong.check.drop_duplicates())
    rows = []
    for variant in ("M2", "M2+O"):
        subset = strong[strong.variant == variant].set_index("check")
        rows.append(
            [variant, *(_percent(subset.loc[check, "mean"]) for check in checks)]
        )
    return _markdown_table(["model", *checks], rows)


def _top_table(top: pd.DataFrame) -> str:
    """Expose the best stable M2 and M2+O configurations without hiding misses."""
    rows = []
    for variant in ("M2", "M2+O"):
        item = top[top.variant == variant].iloc[0]
        rows.append(
            [
                variant,
                str(int(item.candidate_id)),
                str(int(item.structural_scenario_id)),
                _percent(item.seed_pass_fraction),
                _compact(item.rough_hurst),
                _compact(item.leverage_rank),
                _compact(item.zumbach_rank),
                _compact(item.excess_kurtosis),
            ]
        )
    return _markdown_table(
        ["model", "candidate", "scenario", "seed pass", "H", "leverage", "Zumbach", "kurtosis"],
        rows,
    )


def _decision_text(
    result: FactorialExperimentResult,
    stability: pd.DataFrame,
    top: pd.DataFrame,
) -> str:
    """State whether the local region, optional O block, and best cell survived."""
    metrics = result.metric_summary.set_index(["variant", "metric"])
    effects = result.effect_summary.set_index(["effect", "metric"])
    lines = []
    for variant in ("M2", "M2+O"):
        acceptance = metrics.loc[(variant, "acceptance_rate"), "mean"]
        stable = stability[stability.variant == variant].passes_80pct.mean()
        best = top[top.variant == variant].seed_pass_fraction.max()
        lines.append(
            f"- {variant}: row acceptance `{_percent(acceptance)}`, neighborhood ≥4/5-seed "
            f"stability `{_percent(stable)}`, best configuration `{_percent(best)}`."
        )
    orthogonal = effects.loc[("orthogonal_at_echo0", "acceptance_rate")]
    lines.append(
        f"- The paired standardized-O acceptance effect is `{_compact(orthogonal['mean'])}` "
        f"with seed interval `[{_compact(orthogonal['ci_lower'])}, "
        f"{_compact(orthogonal['ci_upper'])}]`."
    )
    lines.append(
        "- Recommend the simpler M2 unless M2+O shows a stable joint-rate gain rather than "
        "only a marginal roughness increase. M3 remains rejected for this stage."
    )
    return "\n".join(lines)


if __name__ == "__main__":
    main()
