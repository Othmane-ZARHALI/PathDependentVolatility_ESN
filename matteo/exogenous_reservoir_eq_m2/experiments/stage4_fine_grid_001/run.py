"""Targeted fine-grid recalibration of the fixed Stage-3 M2 architecture."""

from __future__ import annotations

import json
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import qmc

from esn_eq import (
    AcceptanceCriteria,
    ArchitectureConfig,
    ExogenousReservoirVolatilityModel,
    MarketEnvironment,
    ParameterSpace,
    QmcParameterSampler,
    ReadoutParameters,
    ReservoirFactory,
    RiskPremium,
    SimulationConfig,
    StylisedFactDiagnostics,
)
from esn_eq.config import BankConfig
from esn_eq.diagnostics import ito_return_zumbach_statistic
from esn_eq.sampling import ParameterBound


OUTPUT = Path(__file__).resolve().parent
STAGE3_MANIFEST = OUTPUT.parent / "stage3_001" / "manifest.json"
SEARCH_SEEDS = (1100003, 1130011, 1170007, 1190021)
VALIDATION_SEEDS = (1230001, 1270009, 1290017, 1310021, 1370003)
STEPS_PER_DAY = 64
SEARCH_CANDIDATES = 64
STRUCTURAL_SCENARIOS = 8
SHORTLIST_SIZE = 8
DESIGN_SEED = 20260916
STRUCTURAL_SEED = 20260917

# These scales only break ties between candidates that meet the unchanged gates.
# A score of zero is a boundary; positive values are inside every gate.
MARGIN_SCALES = {
    "max_abs_return_acf": 0.025,
    "mean_log_variance_acf": 0.05,
    "rough_hurst": 0.03,
    "leverage_rank": 0.005,
    "zumbach_rank": 0.01,
    "excess_kurtosis": 2.0,
    "taylor_gap": 0.03,
}


def base_architecture() -> ArchitectureConfig:
    """Load the fitted OU rates/weights and disable rejected optional blocks."""
    manifest = json.loads(STAGE3_MANIFEST.read_text())
    base = ArchitectureConfig.from_mapping(manifest["architecture"])
    return replace(
        base,
        orthogonal=BankConfig(()),
        echo=replace(base.echo, enabled=False),
    )


def anchor_readout() -> ReadoutParameters:
    """Return the candidate used by the preceding time-step audit."""
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


def readout_design() -> tuple[ReadoutParameters, ...]:
    """Create a local Sobol design and retain the old candidate as an anchor."""
    space = ParameterSpace(
        bounds={
            "cluster_loading": ParameterBound(0.18, 0.32),
            "feedback_curvature": ParameterBound(0.058, 0.100),
            "feedback_shift": ParameterBound(2.55, 3.55),
            "spike_curvature": ParameterBound(0.14, 0.36),
            "spike_shift": ParameterBound(0.04, 0.55),
            "cap_level": ParameterBound(5.5, 7.0),
            "cap_sharpness": ParameterBound(6.0, 10.0),
        },
        fixed={"orthogonal_curvature": 0.0, "echo_loading": 0.0},
    )
    sampled = QmcParameterSampler(space).sample(SEARCH_CANDIDATES, DESIGN_SEED, "sobol")
    return (anchor_readout(), *sampled)


def structural_design() -> tuple[tuple[float, float], ...]:
    """Search only the two P-measure correlations while keeping OU rates fixed."""
    manifest = json.loads(STAGE3_MANIFEST.read_text())
    selected = manifest["structural_scenarios"][0]
    anchor = (
        float(selected["feedback_asset_correlation"]),
        float(selected["spike_asset_correlation"]),
    )
    unit = qmc.LatinHypercube(2, seed=STRUCTURAL_SEED).random(STRUCTURAL_SCENARIOS - 1)
    lower = np.asarray([0.90, 0.85])
    upper = np.asarray([0.999, 0.999])
    points = lower + unit * (upper - lower)
    return (anchor, *(tuple(float(value) for value in row) for row in points))


def simulation(stage: str) -> SimulationConfig:
    """Use a cheaper search and a longer, larger untouched validation."""
    if stage == "search":
        return SimulationConfig(4.0, 1.0, 8, 252, STEPS_PER_DAY)
    if stage == "validation":
        return SimulationConfig(6.0, 1.0, 12, 252, STEPS_PER_DAY)
    raise ValueError(f"Unknown stage: {stage}")


def architecture_for(base: ArchitectureConfig, correlations: tuple[float, float]):
    """Apply one correlation scenario without changing any OU rate or weight."""
    feedback_correlation, spike_correlation = correlations
    return replace(
        base,
        feedback=replace(base.feedback, asset_correlation=feedback_correlation),
        spike=replace(base.spike, asset_correlation=spike_correlation),
    )


def ito_metrics(result) -> dict[str, float]:
    """Average the matched Ito-relative-return Zumbach measures across paths."""
    pearson, rank = [], []
    for log_returns, variance in zip(result.log_returns, result.realized_variance, strict=True):
        pearson.append(
            ito_return_zumbach_statistic(
                log_returns, variance, result.observations_per_year, rank_based=False
            )
        )
        rank.append(
            ito_return_zumbach_statistic(
                log_returns, variance, result.observations_per_year, rank_based=True
            )
        )
    return {
        "ito_zumbach": float(np.mean(pearson)),
        "ito_zumbach_rank": float(np.mean(rank)),
    }


def run_design(
    stage: str,
    seeds: tuple[int, ...],
    candidate_ids: tuple[int, ...],
    scenario_ids: tuple[int, ...],
    parameters: tuple[ReadoutParameters, ...],
    scenarios: tuple[tuple[float, float], ...],
) -> pd.DataFrame:
    """Evaluate a declared candidate/scenario set with common paths within each cache."""
    base = base_architecture()
    market = MarketEnvironment()
    premium = RiskPremium()
    diagnostics = StylisedFactDiagnostics()
    criteria = AcceptanceCriteria.stage2_proxy()
    rows: list[dict[str, object]] = []
    total = len(seeds) * len(scenario_ids)
    completed = 0
    for seed in seeds:
        for scenario_id in scenario_ids:
            correlations = scenarios[scenario_id]
            architecture = architecture_for(base, correlations)
            model = ExogenousReservoirVolatilityModel(architecture)
            cache = ReservoirFactory(architecture, simulation(stage)).build(seed)
            prepared = model.prepare(cache, premium)
            for candidate_id in candidate_ids:
                result = model.simulate(prepared, parameters[candidate_id], market)
                summary = diagnostics.evaluate(result)
                accepted, distance = criteria.evaluate(summary.mean)
                row: dict[str, object] = {
                    "stage": stage,
                    "seed": seed,
                    "scenario_id": scenario_id,
                    "candidate_id": candidate_id,
                    "accepted": accepted,
                    "distance": distance,
                    "feedback_asset_correlation": correlations[0],
                    "spike_asset_correlation": correlations[1],
                    "p_cap_exceedance_fraction": result.metadata.p_cap_exceedance_fraction,
                    "q_forward_mean_max_error": result.metadata.q_forward_mean_max_error,
                }
                row.update(summary.mean)
                row.update(ito_metrics(result))
                rows.append(row)
            completed += 1
            (OUTPUT / "progress.json").write_text(
                json.dumps({"stage": stage, "completed_caches": completed, "total_caches": total})
                + "\n"
            )
    return pd.DataFrame(rows)


def signed_margin(metric: str, value: float, criteria: AcceptanceCriteria) -> float:
    """Measure distance to the closest gate boundary in declared practical units."""
    band = criteria.bands[metric]
    scale = MARGIN_SCALES[metric]
    return min((value - band.lower) / scale, (band.upper - value) / scale)


def summarize(frame: pd.DataFrame) -> pd.DataFrame:
    """Summarize seed stability and predeclared worst-seed gate margin."""
    criteria = AcceptanceCriteria.stage2_proxy()
    metrics = tuple(criteria.bands)
    rows = []
    for (candidate_id, scenario_id), group in frame.groupby(
        ["candidate_id", "scenario_id"], sort=True
    ):
        record: dict[str, object] = {
            "candidate_id": int(candidate_id),
            "scenario_id": int(scenario_id),
            "seed_count": int(group.shape[0]),
            "seed_pass_fraction": float(group.accepted.mean()),
            "mean_distance": float(group.distance.mean()),
            "worst_seed_distance": float(group.distance.max()),
            "worst_seed_gate_margin": min(
                signed_margin(metric, float(row[metric]), criteria)
                for _, row in group.iterrows()
                for metric in metrics
            ),
        }
        for column in (*metrics, "ito_zumbach", "ito_zumbach_rank", "p_cap_exceedance_fraction"):
            record[f"mean_{column}"] = float(group[column].mean())
            record[f"sd_{column}"] = float(group[column].std(ddof=1))
        rows.append(record)
    return pd.DataFrame(rows).sort_values(
        ["seed_pass_fraction", "worst_seed_gate_margin", "mean_distance"],
        ascending=[False, False, True],
    )


def design_frames(
    parameters: tuple[ReadoutParameters, ...], scenarios: tuple[tuple[float, float], ...]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return complete auditable design tables."""
    readout = pd.DataFrame([QmcParameterSampler.to_record(item) for item in parameters])
    readout.insert(0, "candidate_id", np.arange(len(parameters)))
    structural = pd.DataFrame(
        [
            {
                "scenario_id": index,
                "feedback_asset_correlation": values[0],
                "spike_asset_correlation": values[1],
            }
            for index, values in enumerate(scenarios)
        ]
    )
    return readout, structural


def criterion_summary(frame: pd.DataFrame) -> pd.DataFrame:
    """Expose every gate instead of hiding failures behind joint acceptance."""
    criteria = AcceptanceCriteria.stage2_proxy()
    rows = []
    for stage, stage_frame in frame.groupby("stage", sort=False):
        for metric, band in criteria.bands.items():
            passed = stage_frame[metric].between(band.lower, band.upper)
            rows.append(
                {
                    "stage": stage,
                    "criterion": metric,
                    "pass_fraction": float(passed.mean()),
                    "minimum": float(stage_frame[metric].min()),
                    "mean": float(stage_frame[metric].mean()),
                    "maximum": float(stage_frame[metric].max()),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    """Run search, freeze the primary selection, and validate on untouched seeds."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    parameters = readout_design()
    scenarios = structural_design()
    readout_frame, structural_frame = design_frames(parameters, scenarios)
    readout_frame.to_csv(OUTPUT / "readout_design.csv", index=False)
    structural_frame.to_csv(OUTPUT / "structural_design.csv", index=False)

    search = run_design(
        "search",
        SEARCH_SEEDS,
        tuple(range(len(parameters))),
        tuple(range(len(scenarios))),
        parameters,
        scenarios,
    )
    search.to_csv(OUTPUT / "search_runs.csv", index=False)
    search_summary = summarize(search)
    search_summary.to_csv(OUTPUT / "search_summary.csv", index=False)

    shortlist = search_summary.head(SHORTLIST_SIZE)
    shortlist.to_csv(OUTPUT / "frozen_shortlist.csv", index=False)
    pairs = tuple(
        (int(row.candidate_id), int(row.scenario_id)) for row in shortlist.itertuples()
    )
    validation_parts = []
    for candidate_id, scenario_id in pairs:
        validation_parts.append(
            run_design(
                "validation",
                VALIDATION_SEEDS,
                (candidate_id,),
                (scenario_id,),
                parameters,
                scenarios,
            )
        )
    validation = pd.concat(validation_parts, ignore_index=True)
    validation.to_csv(OUTPUT / "validation_runs.csv", index=False)
    validation_summary = summarize(validation)
    validation_summary.to_csv(OUTPUT / "validation_summary.csv", index=False)

    # Re-evaluate the preceding candidate on the same untouched seeds. This is a
    # benchmark only and cannot enter the frozen Stage-4 selection ordering.
    anchor_validation = run_design(
        "validation",
        VALIDATION_SEEDS,
        (0,),
        (0,),
        parameters,
        scenarios,
    )
    anchor_validation.to_csv(OUTPUT / "anchor_validation_runs.csv", index=False)
    summarize(anchor_validation).to_csv(OUTPUT / "anchor_validation_summary.csv", index=False)
    criterion_summary(pd.concat([search, validation], ignore_index=True)).to_csv(
        OUTPUT / "criterion_summary.csv", index=False
    )

    primary_candidate, primary_scenario = pairs[0]
    selection = {
        "primary_candidate_id": primary_candidate,
        "primary_scenario_id": primary_scenario,
        "selection_stage": "search only",
        "selection_order": [
            "descending seed_pass_fraction",
            "descending worst_seed_gate_margin",
            "ascending mean_distance",
        ],
        "validation_must_not_reselect_primary": True,
        "shortlist_pairs": pairs,
    }
    (OUTPUT / "selection_record.json").write_text(json.dumps(selection, indent=2) + "\n")

    criteria = AcceptanceCriteria.stage2_proxy()
    manifest = {
        "schema_version": 1,
        "study": "Stage-4 targeted M2 recalibration on a fine time grid",
        "search_seeds": SEARCH_SEEDS,
        "validation_seeds": VALIDATION_SEEDS,
        "design_seed": DESIGN_SEED,
        "structural_seed": STRUCTURAL_SEED,
        "search_candidate_count_plus_anchor": len(parameters),
        "structural_scenario_count": len(scenarios),
        "shortlist_size": SHORTLIST_SIZE,
        "search_simulation": asdict(simulation("search")),
        "validation_simulation": asdict(simulation("validation")),
        "base_architecture": asdict(base_architecture()),
        "market": asdict(MarketEnvironment()),
        "risk_premium": asdict(RiskPremium()),
        "acceptance_criteria": asdict(criteria),
        "margin_scales": MARGIN_SCALES,
        "selection_boundary": (
            "The primary candidate is selected only on search seeds. Validation results "
            "cannot be used to substitute another shortlisted candidate."
        ),
        "ito_zumbach_role": "reported robustness diagnostic, not an optimization gate",
        "anchor_validation_role": "benchmark only; excluded from Stage-4 selection",
        "next_step": "128-step nested confirmation only for a successful primary candidate",
    }
    (OUTPUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
