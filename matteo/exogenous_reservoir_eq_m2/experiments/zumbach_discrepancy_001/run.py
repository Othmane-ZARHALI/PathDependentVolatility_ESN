"""Decompose the Pearson/rank Strong-Zumbach discrepancy on frozen paths."""

from __future__ import annotations

import importlib.util
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, norm, rankdata, t

from esn_eq import (
    ArchitectureConfig,
    ExogenousReservoirVolatilityModel,
    MarketEnvironment,
    ReadoutParameters,
    ReservoirFactory,
    RiskPremium,
)
OUTPUT = Path(__file__).resolve().parent
STAGE4 = OUTPUT.parent / "stage4_fine_grid_001"
WINDOWS = (5, 10, 20)
VALIDATION_SEEDS = (1230001, 1270009, 1290017, 1310021, 1370003)
MODELS = ("primary", "anchor")


def load_stage4_module():
    """Load the frozen Stage-4 simulation specification without duplicating it."""
    path = STAGE4 / "run.py"
    spec = importlib.util.spec_from_file_location("stage4_fine_grid", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the Stage-4 experiment module.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def configurations() -> dict[str, tuple[ArchitectureConfig, ReadoutParameters]]:
    """Reconstruct the frozen primary and preceding anchor exactly."""
    stage4 = load_stage4_module()
    parameters = stage4.readout_design()
    scenarios = stage4.structural_design()
    selection = json.loads((STAGE4 / "selection_record.json").read_text())
    primary_pair = (
        int(selection["primary_candidate_id"]),
        int(selection["primary_scenario_id"]),
    )
    pairs = {"primary": primary_pair, "anchor": (0, 0)}
    result = {}
    for name, (candidate_id, scenario_id) in pairs.items():
        result[name] = (
            stage4.architecture_for(stage4.base_architecture(), scenarios[scenario_id]),
            parameters[candidate_id],
        )
    return result


def paired_variables(
    returns: np.ndarray, variance: np.ndarray, window: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return the forward and backward variable pairs in the Zumbach statistic."""
    cumulative_r = np.concatenate(([0.0], np.cumsum(returns)))
    cumulative_v = np.concatenate(([0.0], np.cumsum(variance)))
    count = returns.size - 2 * window
    center = window + np.arange(count)
    past_r = cumulative_r[center] - cumulative_r[center - window]
    future_r = cumulative_r[center + window] - cumulative_r[center]
    past_v = (cumulative_v[center] - cumulative_v[center - window]) / window
    future_v = (cumulative_v[center + window] - cumulative_v[center]) / window
    return past_r**2, future_v, past_v, future_r**2


def correlation(x: np.ndarray, y: np.ndarray) -> float:
    """Return the ordinary product-moment correlation."""
    return float(np.corrcoef(x, y)[0, 1])


def normal_scores(values: np.ndarray) -> np.ndarray:
    """Map mid-ranks to standard-normal scores."""
    ranks = rankdata(values)
    probabilities = (ranks - 0.5) / values.size
    return norm.ppf(probabilities)


def winsor(values: np.ndarray, quantile: float) -> np.ndarray:
    """Upper-winsorise a nonnegative Zumbach variable."""
    return np.minimum(values, np.quantile(values, quantile))


def estimate_pair(
    x: np.ndarray, y: np.ndarray, method: str
) -> tuple[float, int]:
    """Estimate dependence under one declared bulk/tail treatment."""
    if method == "pearson":
        return correlation(x, y), int(x.size)
    if method == "spearman":
        return correlation(rankdata(x), rankdata(y)), int(x.size)
    if method == "kendall":
        return float(kendalltau(x, y).statistic), int(x.size)
    if method == "gaussian_rank":
        return correlation(normal_scores(x), normal_scores(y)), int(x.size)
    operation, level_text = method.split("_")
    quantile = float(level_text) / 1000.0
    if operation == "winsor":
        return correlation(winsor(x, quantile), winsor(y, quantile)), int(x.size)
    if operation == "trim":
        keep = (x <= np.quantile(x, quantile)) & (y <= np.quantile(y, quantile))
        return correlation(x[keep], y[keep]), int(np.sum(keep))
    raise ValueError(f"Unknown estimator: {method}")


METHODS = (
    "pearson",
    "spearman",
    "gaussian_rank",
    "kendall",
    "winsor_995",
    "winsor_990",
    "winsor_975",
    "winsor_950",
    "trim_995",
    "trim_990",
    "trim_975",
    "trim_950",
)


def tail_concentration(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """Measure how much positive covariance is associated with upper-tail observations."""
    contribution = (x - np.mean(x)) * (y - np.mean(y))
    positive_total = float(np.sum(np.maximum(contribution, 0.0)))
    result = {}
    for level in (0.95, 0.975, 0.99, 0.995):
        tail = (x > np.quantile(x, level)) | (y > np.quantile(y, level))
        result[f"positive_covariance_share_tail_{int(level * 1000)}"] = float(
            np.sum(np.maximum(contribution[tail], 0.0)) / max(positive_total, 1e-15)
        )
    return result


def decile_profile(x: np.ndarray, y: np.ndarray) -> list[dict[str, float]]:
    """Describe the conditional response across ranks of the conditioning variable."""
    bins = np.minimum((rankdata(x, method="average") * 10 / (x.size + 1)).astype(int), 9)
    standardized_y = (y - np.mean(y)) / max(float(np.std(y)), 1e-15)
    return [
        {
            "decile": decile + 1,
            "mean_standardized_response": float(np.mean(standardized_y[bins == decile])),
            "count": int(np.sum(bins == decile)),
        }
        for decile in range(10)
    ]


def analyse_path(
    model_name: str,
    seed: int,
    path_id: int,
    log_returns: np.ndarray,
    variance: np.ndarray,
    observations_per_year: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Return estimator, tail-concentration, and decile records for one path."""
    estimator_rows: list[dict[str, object]] = []
    tail_rows: list[dict[str, object]] = []
    decile_rows: list[dict[str, object]] = []
    return_series = {
        "log": log_returns,
        "ito": log_returns + 0.5 * variance / observations_per_year,
    }
    for return_type, returns in return_series.items():
        for window in WINDOWS:
            past_r2, future_v, past_v, future_r2 = paired_variables(returns, variance, window)
            pairs = {"forward": (past_r2, future_v), "backward": (past_v, future_r2)}
            estimates: dict[tuple[str, str], tuple[float, int]] = {}
            for direction, (x, y) in pairs.items():
                for method in METHODS:
                    estimates[(direction, method)] = estimate_pair(x, y, method)
                tail_row: dict[str, object] = {
                    "model": model_name,
                    "seed": seed,
                    "path": path_id,
                    "return_type": return_type,
                    "window": window,
                    "direction": direction,
                }
                tail_row.update(tail_concentration(x, y))
                tail_rows.append(tail_row)
                for record in decile_profile(x, y):
                    decile_rows.append(
                        {
                            "model": model_name,
                            "seed": seed,
                            "path": path_id,
                            "return_type": return_type,
                            "window": window,
                            "direction": direction,
                            **record,
                        }
                    )
            for method in METHODS:
                forward, forward_n = estimates[("forward", method)]
                backward, backward_n = estimates[("backward", method)]
                estimator_rows.append(
                    {
                        "model": model_name,
                        "seed": seed,
                        "path": path_id,
                        "return_type": return_type,
                        "window": window,
                        "method": method,
                        "forward": forward,
                        "backward": backward,
                        "zumbach": forward - backward,
                        "forward_count": forward_n,
                        "backward_count": backward_n,
                    }
                )
    return estimator_rows, tail_rows, decile_rows


def seed_summary(rows: pd.DataFrame) -> pd.DataFrame:
    """Average paths and windows exactly as the project-level diagnostic does."""
    return (
        rows.groupby(["model", "seed", "return_type", "method"], as_index=False)[
            ["forward", "backward", "zumbach"]
        ]
        .mean()
        .sort_values(["model", "return_type", "method", "seed"])
    )


def overall_summary(seeds: pd.DataFrame) -> pd.DataFrame:
    """Report across-seed means and Student intervals."""
    rows = []
    for keys, group in seeds.groupby(["model", "return_type", "method"], sort=False):
        record = dict(zip(("model", "return_type", "method"), keys, strict=True))
        for metric in ("forward", "backward", "zumbach"):
            values = group[metric].to_numpy(dtype=float)
            mean = float(np.mean(values))
            half_width = float(t.ppf(0.975, values.size - 1) * np.std(values, ddof=1) / np.sqrt(values.size))
            record[f"mean_{metric}"] = mean
            record[f"ci_low_{metric}"] = mean - half_width
            record[f"ci_high_{metric}"] = mean + half_width
        rows.append(record)
    return pd.DataFrame(rows)


def main() -> None:
    """Rebuild matched validation paths and persist every discrepancy diagnostic."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    stage4 = load_stage4_module()
    configurations_by_name = configurations()
    estimator_rows: list[dict[str, object]] = []
    tail_rows: list[dict[str, object]] = []
    decile_rows: list[dict[str, object]] = []
    for model_name in MODELS:
        architecture, parameters = configurations_by_name[model_name]
        model = ExogenousReservoirVolatilityModel(architecture)
        for seed in VALIDATION_SEEDS:
            cache = ReservoirFactory(architecture, stage4.simulation("validation")).build(seed)
            result = model.simulate(model.prepare(cache, RiskPremium()), parameters, MarketEnvironment())
            for path_id, (returns, variance) in enumerate(
                zip(result.log_returns, result.realized_variance, strict=True)
            ):
                estimates, tails, deciles = analyse_path(
                    model_name,
                    seed,
                    path_id,
                    returns,
                    variance,
                    result.observations_per_year,
                )
                estimator_rows.extend(estimates)
                tail_rows.extend(tails)
                decile_rows.extend(deciles)
    estimates = pd.DataFrame(estimator_rows)
    tails = pd.DataFrame(tail_rows)
    deciles = pd.DataFrame(decile_rows)
    seeds = seed_summary(estimates)
    estimates.to_csv(OUTPUT / "path_window_estimators.csv", index=False)
    seeds.to_csv(OUTPUT / "seed_summary.csv", index=False)
    overall_summary(seeds).to_csv(OUTPUT / "estimator_summary.csv", index=False)
    tails.to_csv(OUTPUT / "tail_concentration.csv", index=False)
    (
        tails.groupby(["model", "return_type", "direction"], as_index=False)
        .mean(numeric_only=True)
        .drop(columns=["seed", "path", "window"])
        .to_csv(OUTPUT / "tail_concentration_summary.csv", index=False)
    )
    deciles.to_csv(OUTPUT / "decile_profiles.csv", index=False)
    (
        deciles.groupby(["model", "return_type", "direction", "decile"], as_index=False)[
            "mean_standardized_response"
        ]
        .mean()
        .to_csv(OUTPUT / "decile_profile_summary.csv", index=False)
    )
    manifest = {
        "study": "Pearson/rank Strong-Zumbach discrepancy decomposition",
        "source_experiment": str(STAGE4),
        "models": MODELS,
        "validation_seeds": VALIDATION_SEEDS,
        "simulation": asdict(stage4.simulation("validation")),
        "windows_days": WINDOWS,
        "return_types": ["log", "ito"],
        "methods": METHODS,
        "selection_boundary": "No parameter was fitted or selected in this diagnostic study.",
    }
    (OUTPUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
