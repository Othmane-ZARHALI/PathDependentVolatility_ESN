"""Reproducible paired Pearson-Zumbach mechanism experiment."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

from esn_eq import (
    ArchitectureConfig,
    ExogenousReservoirVolatilityModel,
    MarketEnvironment,
    ReadoutParameters,
    ReservoirFactory,
    RiskPremium,
    SimulationConfig,
)
from esn_eq.config import BankConfig
from esn_eq.diagnostics import zumbach_statistic


OUTPUT = Path(__file__).resolve().parent
SEEDS = (710003, 730003, 750019, 770011, 790021, 810013, 830003, 850009, 870007, 890023)
WINDOWS = (5, 10, 20)
SIMULATION = SimulationConfig(
    years=10.0,
    burn_years=2.0,
    paths=32,
    observations_per_year=252,
    steps_per_observation=2,
)


def selected_architecture() -> ArchitectureConfig:
    """Return the documented selected M2 architecture with no optional blocks."""
    manifest = json.loads((OUTPUT.parent / "stage3_001" / "manifest.json").read_text())
    base = ArchitectureConfig.from_mapping(manifest["architecture"])
    return replace(
        base,
        feedback=replace(base.feedback, asset_correlation=0.996259961884662),
        spike=replace(base.spike, asset_correlation=0.9790548338377681),
        orthogonal=BankConfig(()),
        echo=replace(base.echo, enabled=False),
    )


def selected_readout() -> ReadoutParameters:
    """Use the top all-seed Stage-3 M2 candidate and the selected nonzero spike shift."""
    return ReadoutParameters(
        cluster_loading=0.2333221411332488,
        feedback_curvature=0.0694007444800809,
        feedback_shift=2.9358963377773764,
        spike_curvature=0.24959932155907155,
        spike_shift=0.1,
        orthogonal_curvature=0.0,
        echo_loading=0.0,
        cap_level=6.073741978500038,
        cap_sharpness=7.956517500802875,
    )


def path_mean_zumbach(returns: np.ndarray, variance: np.ndarray) -> float:
    """Average the Pearson statistic over paths; each path averages declared horizons."""
    values = [
        zumbach_statistic(returns[index], variance[index], windows=WINDOWS)
        for index in range(returns.shape[0])
    ]
    return float(np.mean(values))


def garch_positive_control(seed: int) -> float:
    """Generate matched GARCH returns and conditional variance as a directional control."""
    rng = np.random.default_rng(seed)
    kept = int((SIMULATION.years - SIMULATION.burn_years) * SIMULATION.observations_per_year)
    burn = int(SIMULATION.burn_years * SIMULATION.observations_per_year)
    total = kept + burn
    path_values: list[float] = []
    alpha, beta = 0.08, 0.90
    omega = 1.0 - alpha - beta
    for _ in range(SIMULATION.paths):
        eps = rng.standard_normal(total)
        variance = np.empty(total)
        returns = np.empty(total)
        variance[0] = 1.0
        returns[0] = eps[0]
        for index in range(1, total):
            variance[index] = omega + alpha * returns[index - 1] ** 2 + beta * variance[index - 1]
            returns[index] = np.sqrt(variance[index]) * eps[index]
        path_values.append(zumbach_statistic(returns[burn:], variance[burn:], windows=WINDOWS))
    return float(np.mean(path_values))


def main() -> None:
    """Run paired model ablations and save raw and summary results."""
    architecture = selected_architecture()
    null_architecture = replace(
        architecture,
        feedback=replace(architecture.feedback, asset_correlation=0.0),
        spike=replace(architecture.spike, asset_correlation=0.0),
    )
    parameters = selected_readout()
    no_quadratic = replace(parameters, feedback_curvature=0.0)
    no_fast_quadratic = replace(parameters, spike_curvature=0.0)
    no_correlated_quadratics = replace(
        parameters, feedback_curvature=0.0, spike_curvature=0.0
    )
    market = MarketEnvironment()
    premium = RiskPremium()
    factory = ReservoirFactory(architecture, SIMULATION)
    null_factory = ReservoirFactory(null_architecture, SIMULATION)
    model = ExogenousReservoirVolatilityModel(architecture)
    null_model = ExogenousReservoirVolatilityModel(null_architecture)
    rows: list[dict[str, object]] = []

    for seed in SEEDS:
        randomness = factory.draw(seed)
        prepared = model.prepare(factory.from_randomness(randomness), premium)
        null_prepared = null_model.prepare(null_factory.from_randomness(randomness), premium)
        configurations = {
            "Recommended model": model.simulate(prepared, parameters, market),
            "Feedback quadratic removed": model.simulate(prepared, no_quadratic, market),
            "Fast quadratic removed": model.simulate(prepared, no_fast_quadratic, market),
            "Both correlated quadratics removed": model.simulate(
                prepared, no_correlated_quadratics, market
            ),
            "Return-state correlations removed": null_model.simulate(
                null_prepared, parameters, market
            ),
        }
        for name, result in configurations.items():
            rows.append(
                {
                    "seed": seed,
                    "configuration": name,
                    "pearson_zumbach": path_mean_zumbach(
                        result.log_returns, result.realized_variance
                    ),
                }
            )
        rows.append(
            {
                "seed": seed,
                "configuration": "GARCH(1,1) positive control",
                "pearson_zumbach": garch_positive_control(seed),
            }
        )

    with (OUTPUT / "per_seed.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("seed", "configuration", "pearson_zumbach"))
        writer.writeheader()
        writer.writerows(rows)

    summaries: list[dict[str, object]] = []
    for name in dict.fromkeys(str(row["configuration"]) for row in rows):
        values = np.asarray(
            [float(row["pearson_zumbach"]) for row in rows if row["configuration"] == name]
        )
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=1))
        half_width = float(student_t.ppf(0.975, values.size - 1) * sd / np.sqrt(values.size))
        summaries.append(
            {
                "configuration": name,
                "seed_count": values.size,
                "mean": mean,
                "between_seed_sd": sd,
                "ci95_lower": mean - half_width,
                "ci95_upper": mean + half_width,
            }
        )
    with (OUTPUT / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    manifest = {
        "study": "Paired Pearson-Zumbach mechanism ablation",
        "seeds": SEEDS,
        "windows_in_daily_observations": WINDOWS,
        "simulation": asdict(SIMULATION),
        "path_aggregation": "Pearson Z averaged over windows within path, paths within seed",
        "reported_uncertainty": "Between-seed sample SD and two-sided 95% Student-t CI",
        "variance_measure": "Model-integrated average instantaneous variance per day",
        "common_random_numbers": True,
        "readout": asdict(parameters),
        "architecture_source": "stage3_001 manifest; scenario 0; optional echo and orthogonal blocks off",
        "positive_control": {"model": "GARCH(1,1)", "alpha": 0.08, "beta": 0.90},
    }
    (OUTPUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
