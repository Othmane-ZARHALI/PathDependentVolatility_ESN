"""Nested time-step convergence audit for the selected M2 configuration."""

from __future__ import annotations

import json
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pandas as pd

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
from esn_eq.diagnostics import StylisedFactDiagnostics
from esn_eq.experiment import AcceptanceCriteria
from esn_eq.states import ReservoirCache


OUTPUT = Path(__file__).resolve().parent
STAGE3 = OUTPUT.parent / "stage3_001" / "manifest.json"
SEEDS = (910003, 930011, 950009, 970019, 990001)
LEVELS = (2, 4, 8, 16, 32, 64)
REFERENCE_LEVEL = max(LEVELS)
YEARS = 6.0
BURN_YEARS = 1.0
PATHS = 12
OBSERVATIONS_PER_YEAR = 252


def selected_architecture() -> ArchitectureConfig:
    """Return the frozen Stage-3 M2 architecture and selected correlations."""
    manifest = json.loads(STAGE3.read_text())
    base = ArchitectureConfig.from_mapping(manifest["architecture"])
    selected = manifest["structural_scenarios"][0]
    return replace(
        base,
        feedback=replace(
            base.feedback, asset_correlation=float(selected["feedback_asset_correlation"])
        ),
        spike=replace(base.spike, asset_correlation=float(selected["spike_asset_correlation"])),
        orthogonal=BankConfig(()),
        echo=replace(base.echo, enabled=False),
    )


def selected_readout() -> ReadoutParameters:
    """Use the all-seed Stage-3 configuration and documented fast-factor shift."""
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


def simulation(level: int) -> SimulationConfig:
    """Return one internal resolution with the same daily output grid."""
    return SimulationConfig(
        years=YEARS,
        burn_years=BURN_YEARS,
        paths=PATHS,
        observations_per_year=OBSERVATIONS_PER_YEAR,
        steps_per_observation=level,
    )


def coarsen_cache(fine: ReservoirCache, level: int) -> ReservoirCache:
    """Subsample exact OU endpoints and sum the matched equity Brownian increments."""
    if REFERENCE_LEVEL % level:
        raise ValueError("Every level must divide the reference level.")
    factor = REFERENCE_LEVEL // level
    target = simulation(level)
    states = fine.states[::factor]
    equity = fine.equity_increments.reshape(target.steps, factor, PATHS).sum(axis=1)
    if states.shape[0] != target.steps + 1 or equity.shape != (target.steps, PATHS):
        raise RuntimeError("Nested cache dimensions are inconsistent.")
    return ReservoirCache(
        states=states,
        equity_increments=equity,
        layout=fine.layout,
        architecture=fine.architecture,
        simulation=target,
        randomness=fine.randomness,
    )


def flattened_correlation(left: np.ndarray, right: np.ndarray) -> float:
    """Return a finite pooled correlation for matched daily arrays."""
    x, y = left.ravel(), right.ravel()
    return float(np.corrcoef(x, y)[0, 1])


def comparison_metrics(result, reference) -> dict[str, float]:
    """Compare matched daily paths with the finest-grid reference."""
    return_delta = result.log_returns - reference.log_returns
    variance_delta = result.realized_variance - reference.realized_variance
    reference_variance_mean = float(np.mean(reference.realized_variance))
    reference_return_scale = float(np.std(reference.log_returns))
    return {
        "return_path_correlation": flattened_correlation(result.log_returns, reference.log_returns),
        "variance_path_correlation": flattened_correlation(
            result.realized_variance, reference.realized_variance
        ),
        "return_rmse_in_reference_sd": float(
            np.sqrt(np.mean(return_delta**2)) / max(reference_return_scale, 1e-15)
        ),
        "variance_relative_bias": float(
            np.mean(variance_delta) / max(reference_variance_mean, 1e-15)
        ),
        "variance_relative_rmse": float(
            np.sqrt(np.mean(variance_delta**2)) / max(reference_variance_mean, 1e-15)
        ),
        "maximum_volatility": float(np.max(np.sqrt(result.realized_variance))),
        "q999_volatility": float(np.quantile(np.sqrt(result.realized_variance), 0.999)),
    }


def main() -> None:
    """Run the nested grid audit and save seed-level and aggregated evidence."""
    architecture = selected_architecture()
    parameters = selected_readout()
    market = MarketEnvironment()
    premium = RiskPremium()
    model = ExogenousReservoirVolatilityModel(architecture)
    diagnostics = StylisedFactDiagnostics()
    criteria = AcceptanceCriteria.stage2_proxy()
    rows: list[dict[str, object]] = []

    for seed in SEEDS:
        fine_factory = ReservoirFactory(architecture, simulation(REFERENCE_LEVEL))
        fine_cache = fine_factory.build(seed)
        results = {}
        summaries = {}
        for level in LEVELS:
            cache = coarsen_cache(fine_cache, level)
            result = model.simulate(model.prepare(cache, premium), parameters, market)
            results[level] = result
            summaries[level] = diagnostics.evaluate(result)
        reference = results[REFERENCE_LEVEL]
        for level in LEVELS:
            result = results[level]
            summary = summaries[level]
            accepted, distance = criteria.evaluate(summary.mean)
            row: dict[str, object] = {
                "seed": seed,
                "steps_per_day": level,
                "dt_years": simulation(level).dt,
                "maximum_lambda_dt": float(max(fine_cache.layout.rates) * simulation(level).dt),
                "accepted": accepted,
                "distance": distance,
                "p_cap_exceedance_fraction": result.metadata.p_cap_exceedance_fraction,
            }
            row.update(summary.mean)
            row.update(comparison_metrics(result, reference))
            rows.append(row)

    per_seed = pd.DataFrame(rows)
    per_seed.to_csv(OUTPUT / "per_seed.csv", index=False)
    numeric = [
        column
        for column in per_seed.columns
        if column not in {"seed", "steps_per_day", "accepted"}
        and np.issubdtype(per_seed[column].dtype, np.number)
    ]
    grouped = per_seed.groupby("steps_per_day", sort=True)
    summary_rows = []
    for level, frame in grouped:
        record: dict[str, object] = {
            "steps_per_day": int(level),
            "seed_count": int(frame.shape[0]),
            "acceptance_rate": float(frame["accepted"].mean()),
        }
        for column in numeric:
            record[f"mean_{column}"] = float(frame[column].mean())
            record[f"sd_{column}"] = float(frame[column].std(ddof=1))
        summary_rows.append(record)
    pd.DataFrame(summary_rows).to_csv(OUTPUT / "summary.csv", index=False)

    rate_audit = {}
    for name, bank in architecture.banks.items():
        if not bank.rates:
            continue
        rate_audit[name] = {
            "minimum_rate_per_year": min(bank.rates),
            "maximum_rate_per_year": max(bank.rates),
            "fastest_half_life_trading_days": OBSERVATIONS_PER_YEAR
            * np.log(2.0)
            / max(bank.rates),
            "maximum_lambda_dt_by_level": {
                str(level): max(bank.rates) * simulation(level).dt for level in LEVELS
            },
        }
    manifest = {
        "study": "Nested time-step convergence audit for selected M2",
        "seeds": SEEDS,
        "levels_steps_per_trading_day": LEVELS,
        "reference_level": REFERENCE_LEVEL,
        "simulation": asdict(simulation(REFERENCE_LEVEL)),
        "architecture": asdict(architecture),
        "readout": asdict(parameters),
        "market": asdict(market),
        "risk_premium": asdict(premium),
        "rate_audit": rate_audit,
        "nesting": (
            "Exact OU endpoints are simulated once at 64 steps/day; coarser states are exact "
            "subsamples and coarser equity increments are sums of the same fine increments."
        ),
        "selection_boundary": "No parameter is recalibrated at any resolution.",
    }
    (OUTPUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
