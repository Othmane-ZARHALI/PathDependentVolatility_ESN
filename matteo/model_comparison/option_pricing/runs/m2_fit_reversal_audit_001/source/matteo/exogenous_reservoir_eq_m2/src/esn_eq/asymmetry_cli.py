"""Test a downside-oriented shift of the fast M2 spike parabola."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .config import MarketEnvironment, ReadoutParameters, SimulationConfig
from .diagnostics import StylisedFactDiagnostics
from .experiment import AcceptanceCriteria
from .model import ExogenousReservoirVolatilityModel, SimulationResult
from .sampling import QmcParameterSampler, StructuralParameters
from .stage2_cli import _atomic_csv, _atomic_json, _atomic_text, _load_json, _markdown_table
from .stage3_cli import _stage3_architecture
from .states import ReservoirFactory


NEW_METRICS = (
    "tail_negative_day0_response",
    "tail_positive_day0_response",
    "tail_day0_asymmetry_gap",
    "tail_days1_5_asymmetry_gap",
)


def build_parser() -> argparse.ArgumentParser:
    """Configure a focused paired screen around stable Stage-3 M2 candidates."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("experiments/asymmetry_001"))
    parser.add_argument("--provenance", type=Path, default=Path("provenance"))
    parser.add_argument("--stage3", type=Path, default=Path("experiments/stage3_001"))
    parser.add_argument("--configurations", type=int, default=16)
    parser.add_argument(
        "--spike-shifts",
        type=float,
        nargs="+",
        default=[0.0, 0.10, 0.20, 0.30, 0.40, 0.50],
    )
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=[510007, 540011, 570013, 600017, 630019],
    )
    parser.add_argument("--years", type=float, default=6.0)
    parser.add_argument("--burn-years", type=float, default=1.0)
    parser.add_argument("--paths", type=int, default=24)
    parser.add_argument("--steps-per-observation", type=int, default=2)
    return parser


def main() -> None:
    """Run the paired fresh-seed screen and write an auditable result ledger."""
    args = build_parser().parse_args()
    shifts = _validate_shifts(args.spike_shifts)
    selected = _select_configurations(args.stage3, args.configurations)
    architecture = _stage3_architecture(
        _load_json(args.provenance / "latent_calibration.json")
    )
    simulation = SimulationConfig(
        years=args.years,
        burn_years=args.burn_years,
        paths=args.paths,
        steps_per_observation=args.steps_per_observation,
    )
    criteria = AcceptanceCriteria.stage2_proxy()
    runs = run_screen(
        architecture,
        simulation,
        MarketEnvironment(),
        criteria,
        selected,
        shifts,
        tuple(args.seeds),
    )
    summary = summarize_by_shift(runs, criteria)
    paired = paired_shift_effects(runs)
    manifest = {
        "schema_version": 1,
        "model_package_version": "0.4.0",
        "study_name": "Fast-factor downside-asymmetry screen",
        "selection_source": str(args.stage3 / "top_configurations.csv"),
        "selection_rule": "top Stage-3 M2 configuration-scenario rows",
        "configuration_count": len(selected),
        "spike_shifts": list(shifts),
        "seeds": list(args.seeds),
        "simulation": asdict(simulation),
        "acceptance_criteria": asdict(criteria),
        "asymmetry_metrics_are_empirical_targets": False,
        "notes": (
            "The old joint gates are unchanged. Tail-event asymmetry is reported as a "
            "mechanism diagnostic until matched SPX confidence bands are supplied."
        ),
    }
    args.output.mkdir(parents=True, exist_ok=True)
    _atomic_csv(args.output / "selected_configurations.csv", selected)
    _atomic_csv(args.output / "runs.csv", runs)
    _atomic_csv(args.output / "shift_summary.csv", summary)
    _atomic_csv(args.output / "paired_shift_effects.csv", paired)
    _atomic_json(args.output / "manifest.json", manifest)
    _atomic_text(args.output / "RESULTS.md", build_report(summary, paired, manifest))
    print(f"Completed {len(runs)} paired evaluations across {len(args.seeds)} fresh seeds.")
    print(args.output / "RESULTS.md")


def _validate_shifts(values: Iterable[float]) -> tuple[float, ...]:
    """Require a unique non-negative grid containing the exact Stage-3 control."""
    shifts = tuple(float(value) for value in values)
    if not shifts or any(not np.isfinite(value) or value < 0.0 for value in shifts):
        raise ValueError("spike shifts must be finite and non-negative")
    if len(set(shifts)) != len(shifts) or 0.0 not in shifts:
        raise ValueError("spike shifts must be unique and contain the zero control")
    return tuple(sorted(shifts))


def _select_configurations(stage3: Path, count: int) -> pd.DataFrame:
    """Select stable M2 rows while retaining their exact structural scenarios."""
    if count < 1:
        raise ValueError("configurations must be positive")
    source = pd.read_csv(stage3 / "top_configurations.csv")
    selected = source[source.variant == "M2"].head(count).copy()
    if len(selected) < count:
        raise ValueError("the Stage-3 ledger contains too few M2 configurations")
    selected.insert(0, "shortlist_id", np.arange(len(selected), dtype=int))
    return selected.reset_index(drop=True)


def run_screen(
    architecture,
    simulation: SimulationConfig,
    market: MarketEnvironment,
    criteria: AcceptanceCriteria,
    selected: pd.DataFrame,
    shifts: tuple[float, ...],
    seeds: tuple[int, ...],
) -> pd.DataFrame:
    """Reuse each seed's primitives and each scenario's features across all shifts."""
    if len(seeds) < 2 or len(set(seeds)) != len(seeds):
        raise ValueError("at least two distinct fresh seeds are required")
    diagnostics = StylisedFactDiagnostics()
    base_factory = ReservoirFactory(architecture, simulation)
    rows: list[dict[str, object]] = []
    for seed in seeds:
        randomness = base_factory.draw(seed)
        prepared_by_scenario = {}
        for item in selected.itertuples(index=False):
            scenario_id = int(item.structural_scenario_id)
            if scenario_id not in prepared_by_scenario:
                scenario = _scenario_from_row(item)
                scenario_architecture = scenario.architecture(architecture)
                reservoir = ReservoirFactory(
                    scenario_architecture, simulation
                ).from_randomness(randomness)
                model = ExogenousReservoirVolatilityModel(scenario_architecture)
                prepared = model.prepare(reservoir, scenario.risk_premium())
                prepared_by_scenario[scenario_id] = (scenario, model, prepared)
            scenario, model, prepared = prepared_by_scenario[scenario_id]
            base = _parameters_from_row(item)
            for shift in shifts:
                parameters = replace(base, spike_shift=shift)
                result = model.simulate(prepared, parameters, market, measure="P")
                summary = diagnostics.evaluate(result)
                accepted, distance = criteria.evaluate(summary.mean)
                row = {
                    "seed": seed,
                    "shortlist_id": int(item.shortlist_id),
                    "stage3_candidate_id": int(item.candidate_id),
                    "structural_scenario_id": scenario_id,
                    "spike_shift": shift,
                    "accepted": accepted,
                    "distance": distance,
                }
                row.update(scenario.to_record())
                row.update(QmcParameterSampler.to_record(parameters))
                row.update(summary.mean)
                row.update(tail_event_asymmetry(result))
                row.update(asdict(result.metadata))
                rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["seed", "shortlist_id", "spike_shift"]
    ).reset_index(drop=True)


def _scenario_from_row(row) -> StructuralParameters:
    """Restore the exact Stage-3 outer parameters for one shortlisted row."""
    return StructuralParameters(
        feedback_asset_correlation=float(row.feedback_asset_correlation),
        spike_asset_correlation=float(row.spike_asset_correlation),
        equity_risk_price=float(row.equity_risk_price),
        feedback_idiosyncratic_risk_price=float(
            row.feedback_idiosyncratic_risk_price
        ),
        spike_idiosyncratic_risk_price=float(row.spike_idiosyncratic_risk_price),
    )


def _parameters_from_row(row) -> ReadoutParameters:
    """Restore one Stage-3 M2 core and explicitly disable optional components."""
    return ReadoutParameters(
        cluster_loading=float(row.cluster_loading),
        feedback_curvature=float(row.feedback_curvature),
        feedback_shift=float(row.feedback_shift),
        spike_curvature=float(row.spike_curvature),
        spike_shift=0.0,
        orthogonal_curvature=0.0,
        echo_loading=0.0,
        cap_level=float(row.cap_level),
        cap_sharpness=float(row.cap_sharpness),
    )


def tail_event_asymmetry(
    result: SimulationResult,
    tail_probability: float = 0.025,
    baseline_days: int = 10,
    forward_days: int = 5,
) -> dict[str, float]:
    """Compare volatility around equally frequent negative and positive tail returns."""
    negative: list[np.ndarray] = []
    positive: list[np.ndarray] = []
    volatility = np.sqrt(result.realized_variance)
    for returns, path_volatility in zip(result.log_returns, volatility):
        low, high = np.quantile(returns, [tail_probability, 1.0 - tail_probability])
        valid = np.arange(baseline_days, returns.size - forward_days)
        _collect_event_responses(
            valid[returns[valid] <= low], path_volatility, baseline_days, forward_days, negative
        )
        _collect_event_responses(
            valid[returns[valid] >= high], path_volatility, baseline_days, forward_days, positive
        )
    if not negative or not positive:
        return {name: float("nan") for name in NEW_METRICS}
    negative_median = np.median(np.asarray(negative), axis=0)
    positive_median = np.median(np.asarray(positive), axis=0)
    return {
        "tail_negative_day0_response": float(negative_median[0]),
        "tail_positive_day0_response": float(positive_median[0]),
        "tail_day0_asymmetry_gap": float(negative_median[0] - positive_median[0]),
        "tail_days1_5_asymmetry_gap": float(
            np.mean(negative_median[1:] - positive_median[1:])
        ),
    }


def _collect_event_responses(
    events: np.ndarray,
    volatility: np.ndarray,
    baseline_days: int,
    forward_days: int,
    target: list[np.ndarray],
) -> None:
    """Append normalized day-zero-through-forward-day event responses."""
    for event in events:
        baseline = float(np.mean(volatility[event - baseline_days : event]))
        target.append(volatility[event : event + forward_days + 1] / max(baseline, 1e-15) - 1.0)


def summarize_by_shift(
    runs: pd.DataFrame, criteria: AcceptanceCriteria | None = None
) -> pd.DataFrame:
    """Summarize seed uncertainty, robust configurations, and individual gates."""
    criteria = criteria or AcceptanceCriteria.stage2_proxy()
    metrics = (
        "accepted",
        "distance",
        "rough_hurst",
        "leverage_rank",
        "zumbach_rank",
        "excess_kurtosis",
        "max_abs_return_acf",
        "mean_log_variance_acf",
        "taylor_gap",
        "p_cap_exceedance_fraction",
        *NEW_METRICS,
    )
    per_seed = runs.groupby(["seed", "spike_shift"], as_index=False)[list(metrics)].mean()
    rows = []
    for shift, frame in per_seed.groupby("spike_shift", sort=True):
        row: dict[str, float | int] = {
            "spike_shift": float(shift),
            "seed_count": len(frame),
        }
        for metric in metrics:
            values = frame[metric].to_numpy(dtype=float)
            row[f"{metric}_mean"] = float(np.mean(values))
            row[f"{metric}_between_seed_sd"] = float(np.std(values, ddof=1))
        shift_runs = runs[runs.spike_shift == shift]
        configuration_passes = shift_runs.groupby("shortlist_id").accepted.agg(
            ["sum", "count"]
        )
        row["robust_configuration_count"] = int(
            (configuration_passes["sum"] == configuration_passes["count"]).sum()
        )
        row["configuration_count"] = int(len(configuration_passes))
        for name, band in criteria.bands.items():
            values = shift_runs[name].to_numpy(dtype=float)
            row[f"pass_{name}_mean"] = float(
                np.mean(np.isfinite(values) & (values >= band.lower) & (values <= band.upper))
            )
        rows.append(row)
    return pd.DataFrame(rows)


def paired_shift_effects(runs: pd.DataFrame) -> pd.DataFrame:
    """Compute within-seed/configuration effects relative to the exact zero control."""
    keys = ["seed", "shortlist_id"]
    metrics = (
        "accepted",
        "distance",
        "rough_hurst",
        "leverage_rank",
        "zumbach_rank",
        "excess_kurtosis",
        "max_abs_return_acf",
        "mean_log_variance_acf",
        "taylor_gap",
        "p_cap_exceedance_fraction",
        *NEW_METRICS,
    )
    rows = []
    control = runs[runs.spike_shift == 0.0].set_index(keys)
    for shift in sorted(value for value in runs.spike_shift.unique() if value > 0.0):
        active = runs[runs.spike_shift == shift].set_index(keys)
        if not active.index.equals(control.index):
            raise ValueError("shift cells are not exactly paired with the zero control")
        for metric in metrics:
            difference = active[metric].astype(float) - control[metric].astype(float)
            per_seed = difference.groupby(level="seed").mean().to_numpy(dtype=float)
            rows.append(
                {
                    "spike_shift": float(shift),
                    "metric": "acceptance_rate" if metric == "accepted" else metric,
                    "mean_paired_effect": float(np.mean(per_seed)),
                    "between_seed_sd": float(np.std(per_seed, ddof=1)),
                    "seed_min": float(np.min(per_seed)),
                    "seed_max": float(np.max(per_seed)),
                }
            )
    return pd.DataFrame(rows)


def build_report(
    summary: pd.DataFrame,
    paired: pd.DataFrame,
    manifest: dict[str, object],
) -> str:
    """Create a compact interpretation without inventing an empirical asymmetry gate."""
    columns = (
        "accepted_mean",
        "rough_hurst_mean",
        "leverage_rank_mean",
        "zumbach_rank_mean",
        "excess_kurtosis_mean",
        "tail_negative_day0_response_mean",
        "tail_positive_day0_response_mean",
        "tail_day0_asymmetry_gap_mean",
        "tail_days1_5_asymmetry_gap_mean",
    )
    table_rows = [
        [f"{row.spike_shift:.2f}", *(f"{getattr(row, name):.5g}" for name in columns)]
        for row in summary.itertuples(index=False)
    ]
    best_acceptance = summary.iloc[summary.accepted_mean.argmax()]
    nonzero = summary[summary.spike_shift > 0.0]
    best_nonzero = nonzero.iloc[nonzero.accepted_mean.argmax()]
    best_gap = summary.iloc[summary.tail_day0_asymmetry_gap_mean.argmax()]
    gate_columns = [name for name in summary.columns if name.startswith("pass_")]
    failed_gates = [
        name.removeprefix("pass_").removesuffix("_mean")
        for name in gate_columns
        if float(best_nonzero[name]) < 1.0
    ]
    failure_reading = (
        ", ".join(f"`{name}`" for name in failed_gates)
        if failed_gates
        else "none"
    )
    positive_effects = paired[
        (paired.metric == "tail_positive_day0_response")
        & (paired.mean_paired_effect < 0.0)
    ]
    positive_reading = (
        f"Positive-tail day-zero response fell for {len(positive_effects)} of "
        f"{summary.spike_shift.gt(0).sum()} tested non-zero shifts."
    )
    headers = [
        "spike shift",
        "joint pass",
        "H",
        "leverage",
        "Zumbach",
        "kurtosis",
        "negative day 0",
        "positive day 0",
        "day-0 gap",
        "days 1-5 gap",
    ]
    return "\n".join(
        (
            "# Fast-factor downside-asymmetry screen",
            "",
            "## Interpretation boundary",
            "",
            "This is a paired fresh-seed screen of shortlisted Stage-3 M2 configurations. "
            "The original seven stylised-fact gates are unchanged. The tail-event metrics "
            "measure direction, but are not empirical SPX acceptance targets.",
            "",
            "## Results by shift",
            "",
            _markdown_table(headers, table_rows),
            "",
            "Event responses are relative changes in volatility versus the preceding ten-day "
            "mean. Joint pass is the fraction satisfying all original Stage-3 proxy gates.",
            "",
            "## Screen reading",
            "",
            f"- Highest joint pass rate: spike shift `{best_acceptance.spike_shift:.2f}` "
            f"with `{best_acceptance.accepted_mean:.1%}`.",
            f"- Best non-zero shift on the old gates: `{best_nonzero.spike_shift:.2f}` "
            f"with `{best_nonzero.accepted_mean:.1%}` joint acceptance and "
            f"`{int(best_nonzero.robust_configuration_count)}` of "
            f"`{int(best_nonzero.configuration_count)}` configurations passing on every seed.",
            f"- At that non-zero shift, the only individual gate with failures was: "
            f"{failure_reading}.",
            f"- Largest day-zero downside-minus-upside gap: spike shift "
            f"`{best_gap.spike_shift:.2f}` with `{best_gap.tail_day0_asymmetry_gap_mean:.1%}`.",
            f"- {positive_reading}",
            "- Select a shift only after considering joint acceptance, leverage, Zumbach, "
            "return tails, cap use, and both asymmetry horizons together.",
            "",
            "## Next evidential step",
            "",
            "Estimate the same event-response statistics and uncertainty from matched SPX and "
            "intraday realised-variance data. Then run an untouched joint search around the "
            "best non-zero shift rather than treating this shortlist screen as final calibration.",
            "",
            f"Package version: `{manifest['model_package_version']}`.",
            "",
        )
    )


if __name__ == "__main__":
    main()
