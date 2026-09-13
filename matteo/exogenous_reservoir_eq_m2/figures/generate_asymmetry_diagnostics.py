"""Generate paired visual diagnostics for the fast-factor downside shift.

The selected configuration (Stage-3 candidate 26, scenario 0) is one of the
six configurations that passes every original gate on all five fresh seeds at
``spike_shift=0.10``.  Control and shifted paths reuse identical reservoir
randomness.  Tail events are selected from the control returns so both models
are evaluated on exactly the same event dates.
"""

from __future__ import annotations

import copy
import json
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from esn_eq.config import (
    ArchitectureConfig,
    MarketEnvironment,
    ReadoutParameters,
    RiskPremium,
    SimulationConfig,
)
from esn_eq.model import ExogenousReservoirVolatilityModel, SimulationResult
from esn_eq.states import ReservoirFactory


CONTROL_SHIFT = 0.0
ASYMMETRIC_SHIFT = 0.10
CANDIDATE_ID = 26
SCENARIO_ID = 0
FRESH_SEEDS = (510007, 540011, 570013, 600017, 630019)
CONTROL_COLOR = "#50616f"
ASYMMETRIC_COLOR = "#c94f3d"
NEGATIVE_COLOR = "#b33c49"
POSITIVE_COLOR = "#277a63"
INK = "#18324a"
GRID = "#dfe5ea"


def _load_configuration(root: Path) -> tuple[
    ArchitectureConfig,
    SimulationConfig,
    MarketEnvironment,
    RiskPremium,
    ReadoutParameters,
]:
    """Restore one robust Stage-3 M2 configuration from its audit ledger."""
    stage3 = root / "experiments" / "stage3_001"
    manifest = json.loads((stage3 / "manifest.json").read_text(encoding="utf-8"))
    table = pd.read_csv(stage3 / "top_configurations.csv")
    selected = table[
        table.variant.eq("M2")
        & table.candidate_id.eq(CANDIDATE_ID)
        & table.structural_scenario_id.eq(SCENARIO_ID)
    ]
    if len(selected) != 1:
        raise ValueError("Expected exactly one robust Stage-3 configuration.")
    row = selected.iloc[0]
    architecture_values = copy.deepcopy(manifest["architecture"])
    architecture_values["feedback"]["asset_correlation"] = float(
        row.feedback_asset_correlation
    )
    architecture_values["spike"]["asset_correlation"] = float(
        row.spike_asset_correlation
    )
    architecture = ArchitectureConfig.from_mapping(architecture_values)
    simulation = SimulationConfig(**manifest["simulation"])
    market = MarketEnvironment(**manifest["market"])
    premium = RiskPremium(
        equity=float(row.equity_risk_price),
        feedback_idiosyncratic=float(row.feedback_idiosyncratic_risk_price),
        spike_idiosyncratic=float(row.spike_idiosyncratic_risk_price),
    )
    parameters = ReadoutParameters(
        cluster_loading=float(row.cluster_loading),
        feedback_curvature=float(row.feedback_curvature),
        feedback_shift=float(row.feedback_shift),
        spike_curvature=float(row.spike_curvature),
        spike_shift=CONTROL_SHIFT,
        orthogonal_curvature=0.0,
        echo_loading=0.0,
        cap_level=float(row.cap_level),
        cap_sharpness=float(row.cap_sharpness),
    )
    return architecture, simulation, market, premium, parameters


def _paired_simulations(root: Path) -> list[tuple[int, SimulationResult, SimulationResult]]:
    """Simulate control and asymmetric readouts on identical latent paths."""
    architecture, simulation, market, premium, parameters = _load_configuration(root)
    factory = ReservoirFactory(architecture, simulation)
    model = ExogenousReservoirVolatilityModel(architecture)
    paired = []
    for seed in FRESH_SEEDS:
        prepared = model.prepare(factory.build(seed), premium)
        control = model.simulate(prepared, parameters, market, measure="P")
        asymmetric = model.simulate(
            prepared,
            replace(parameters, spike_shift=ASYMMETRIC_SHIFT),
            market,
            measure="P",
        )
        paired.append((seed, control, asymmetric))
    return paired


def _collect_paired_events(
    paired: list[tuple[int, SimulationResult, SimulationResult]],
) -> tuple[np.ndarray, dict[str, dict[str, np.ndarray]]]:
    """Use control-return event dates for a genuinely matched event study."""
    lags = np.arange(-10, 21, dtype=int)
    values = {
        "control": {"negative": [], "positive": []},
        "asymmetric": {"negative": [], "positive": []},
    }
    for _, control, asymmetric in paired:
        control_vol = 100.0 * np.sqrt(control.realized_variance)
        asymmetric_vol = 100.0 * np.sqrt(asymmetric.realized_variance)
        for returns, vol0, vol1 in zip(control.log_returns, control_vol, asymmetric_vol):
            low, high = np.quantile(returns, [0.025, 0.975])
            valid = np.arange(10, returns.size - 20)
            for direction, events in (
                ("negative", valid[returns[valid] <= low]),
                ("positive", valid[returns[valid] >= high]),
            ):
                for event in events:
                    for name, volatility in (("control", vol0), ("asymmetric", vol1)):
                        baseline = float(np.mean(volatility[event - 10 : event]))
                        response = 100.0 * (
                            volatility[event + lags] / max(baseline, 1e-15) - 1.0
                        )
                        values[name][direction].append(response)
    arrays = {
        model: {direction: np.asarray(rows) for direction, rows in directions.items()}
        for model, directions in values.items()
    }
    return lags, arrays


def _median_band(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return event-level interquartile band and median."""
    return tuple(np.quantile(values, quantile, axis=0) for quantile in (0.25, 0.50, 0.75))


def _style_axis(axis: plt.Axes) -> None:
    """Apply the common restrained scientific style."""
    axis.grid(True, color=GRID, linewidth=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(colors="#334553")
    axis.title.set_color(INK)


def _save(figure: plt.Figure, path: Path) -> None:
    """Save a sharp publication-size PNG with deterministic metadata."""
    figure.savefig(path, dpi=180, bbox_inches="tight", facecolor="white", metadata={"Software": "esn-eq-model"})
    plt.close(figure)


def event_study_figure(
    output: Path,
    lags: np.ndarray,
    events: dict[str, dict[str, np.ndarray]],
) -> dict[str, float]:
    """Show how the shifted model changes matched tail-event responses."""
    figure, axes = plt.subplots(2, 2, figsize=(13.2, 9.2))
    summaries: dict[str, dict[str, np.ndarray]] = {}
    for column, (model, title) in enumerate(
        (("control", r"Control: $\beta_J=0$"), ("asymmetric", r"Shifted: $\beta_J=0.10$"))
    ):
        summaries[model] = {}
        for direction, color, label in (
            ("negative", NEGATIVE_COLOR, "Negative 2.5% tail"),
            ("positive", POSITIVE_COLOR, "Positive 2.5% tail"),
        ):
            lower, median, upper = _median_band(events[model][direction])
            summaries[model][direction] = median
            axes[0, column].fill_between(lags, lower, upper, color=color, alpha=0.14)
            axes[0, column].plot(lags, median, color=color, linewidth=2.2, label=label)
        axes[0, column].axvline(0, color=CONTROL_COLOR, linestyle=":", linewidth=1.2)
        axes[0, column].axhline(0, color=CONTROL_COLOR, linewidth=0.9)
        axes[0, column].set_title(title)
        axes[0, column].set_ylabel("Volatility vs 10-day baseline (%)")
        axes[0, column].set_xlim(lags[0], lags[-1])
        axes[0, column].legend(frameon=False, loc="upper right")
        _style_axis(axes[0, column])

    for model, color, label in (
        ("control", CONTROL_COLOR, r"Control $\beta_J=0$"),
        ("asymmetric", ASYMMETRIC_COLOR, r"Shifted $\beta_J=0.10$"),
    ):
        gap = summaries[model]["negative"] - summaries[model]["positive"]
        axes[1, 0].plot(lags, gap, color=color, linewidth=2.5, label=label)
    axes[1, 0].axvline(0, color=CONTROL_COLOR, linestyle=":", linewidth=1.2)
    axes[1, 0].axhline(0, color=CONTROL_COLOR, linewidth=0.9)
    axes[1, 0].set_title("Downside-minus-upside response gap")
    axes[1, 0].set_ylabel("Gap (percentage points)")
    axes[1, 0].set_xlabel("Trading days from control tail return")
    axes[1, 0].set_xlim(lags[0], lags[-1])
    axes[1, 0].legend(frameon=False)
    _style_axis(axes[1, 0])

    day0 = np.where(lags == 0)[0][0]
    days1_5 = (lags >= 1) & (lags <= 5)
    metrics = {
        "Control": (
            summaries["control"]["negative"][day0] - summaries["control"]["positive"][day0],
            float(np.mean(summaries["control"]["negative"][days1_5] - summaries["control"]["positive"][days1_5])),
        ),
        "Shifted": (
            summaries["asymmetric"]["negative"][day0] - summaries["asymmetric"]["positive"][day0],
            float(np.mean(summaries["asymmetric"]["negative"][days1_5] - summaries["asymmetric"]["positive"][days1_5])),
        ),
    }
    x = np.arange(2)
    width = 0.34
    control_bars = axes[1, 1].bar(
        x - width / 2, metrics["Control"], width, color=CONTROL_COLOR, label="Control"
    )
    shifted_bars = axes[1, 1].bar(
        x + width / 2, metrics["Shifted"], width, color=ASYMMETRIC_COLOR, label="Shifted"
    )
    axes[1, 1].bar_label(control_bars, fmt="%.1f", padding=3)
    axes[1, 1].bar_label(shifted_bars, fmt="%.1f", padding=3)
    axes[1, 1].set_xticks(x, ("Day 0", "Mean days 1-5"))
    axes[1, 1].set_ylabel("Downside-minus-upside gap (pp)")
    axes[1, 1].set_title("Shift enlarges and prolongs asymmetry")
    axes[1, 1].legend(frameon=False)
    _style_axis(axes[1, 1])
    figure.suptitle(
        "Matched tail events: a small fast-factor shift strengthens downside volatility",
        color=INK,
        fontsize=17,
        fontweight="bold",
    )
    figure.text(
        0.5,
        0.015,
        "Candidate 26, scenario 0; 5 fresh seeds x 24 paths. Events fixed from control returns; bands are event-level IQRs.",
        ha="center",
        color="#52626f",
    )
    figure.tight_layout(rect=(0, 0.04, 1, 0.95))
    _save(figure, output / "05_m2_leverage_event_study.png")
    return {
        "control_day0_gap_pp": float(metrics["Control"][0]),
        "shifted_day0_gap_pp": float(metrics["Shifted"][0]),
        "control_days1_5_gap_pp": float(metrics["Control"][1]),
        "shifted_days1_5_gap_pp": float(metrics["Shifted"][1]),
        "negative_event_count": int(events["control"]["negative"].shape[0]),
        "positive_event_count": int(events["control"]["positive"].shape[0]),
    }


def mechanism_figure(output: Path, root: Path) -> None:
    """Visualize the exact equal-shock orientation of the shifted parabola."""
    _, _, _, _, parameters = _load_configuration(root)
    factor = np.linspace(-2.5, 2.5, 501)
    figure, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))
    for shift, color, label in (
        (CONTROL_SHIFT, CONTROL_COLOR, r"Control $\beta_J=0$"),
        (ASYMMETRIC_SHIFT, ASYMMETRIC_COLOR, r"Shifted $\beta_J=0.10$"),
    ):
        relative_score = parameters.spike_curvature * ((factor - shift) ** 2 - shift**2)
        axes[0].plot(factor, relative_score, color=color, linewidth=2.5, label=label)
    axes[0].axvline(0, color=CONTROL_COLOR, linewidth=0.8)
    axes[0].axhline(0, color=CONTROL_COLOR, linewidth=0.8)
    axes[0].set_xlabel("Fast factor shock $F$")
    axes[0].set_ylabel("Fast score relative to $F=0$")
    axes[0].set_title("The shift tilts an otherwise symmetric parabola")
    axes[0].legend(frameon=False)
    _style_axis(axes[0])

    magnitudes = np.linspace(0.0, 2.5, 251)
    for shift, color, label in (
        (CONTROL_SHIFT, CONTROL_COLOR, r"Control $\beta_J=0$"),
        (ASYMMETRIC_SHIFT, ASYMMETRIC_COLOR, r"Shifted $\beta_J=0.10$"),
    ):
        log_ratio = 8.0 * parameters.spike_curvature * shift * magnitudes
        axes[1].plot(magnitudes, np.exp(log_ratio), color=color, linewidth=2.5, label=label)
    axes[1].axhline(1.0, color=CONTROL_COLOR, linewidth=0.8)
    axes[1].set_xlabel("Equal shock magnitude $f$")
    axes[1].set_ylabel(r"Variance ratio $V(F=-f)/V(F=+f)$")
    axes[1].set_title("Equal downside shocks now produce more variance")
    axes[1].legend(frameon=False)
    _style_axis(axes[1])
    figure.suptitle("How the fast-factor shift creates downside asymmetry", color=INK, fontsize=17, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    _save(figure, output / "07_m2_asymmetry_mechanism.png")


def tradeoff_figure(output: Path, root: Path) -> None:
    """Expose the empirical asymmetry-versus-Zumbach trade-off across shifts."""
    summary = pd.read_csv(root / "experiments" / "asymmetry_001" / "shift_summary.csv")
    shift = summary.spike_shift.to_numpy()
    figure, axes = plt.subplots(1, 3, figsize=(15.2, 5.1))

    axes[0].plot(shift, 100 * summary.tail_day0_asymmetry_gap_mean, "o-", color=NEGATIVE_COLOR, linewidth=2.3, label="Day 0")
    axes[0].plot(shift, 100 * summary.tail_days1_5_asymmetry_gap_mean, "s-", color=POSITIVE_COLOR, linewidth=2.3, label="Mean days 1-5")
    axes[0].axvline(ASYMMETRIC_SHIFT, color=ASYMMETRIC_COLOR, linestyle=":")
    axes[0].set_xlabel(r"Fast-factor shift $\beta_J$")
    axes[0].set_ylabel("Downside-minus-upside gap (pp)")
    axes[0].set_title("Asymmetry strengthens")
    axes[0].legend(frameon=False)
    _style_axis(axes[0])

    axes[1].plot(shift, 100 * summary.accepted_mean, "o-", color=INK, linewidth=2.3, label="Joint pass rate")
    axes[1].axvline(ASYMMETRIC_SHIFT, color=ASYMMETRIC_COLOR, linestyle=":")
    axes[1].set_xlabel(r"Fast-factor shift $\beta_J$")
    axes[1].set_ylabel("All seven gates passed (%)")
    axes[1].set_ylim(-4, 104)
    axes[1].set_title("A small shift remains viable")
    _style_axis(axes[1])
    twin = axes[1].twinx()
    twin.plot(shift, summary.robust_configuration_count, "D--", color=ASYMMETRIC_COLOR, linewidth=1.8, label="All-seed configurations")
    twin.set_ylabel("Configurations passing all 5 seeds")
    twin.set_ylim(-0.5, 16.5)
    handles, labels = axes[1].get_legend_handles_labels()
    handles2, labels2 = twin.get_legend_handles_labels()
    axes[1].legend(handles + handles2, labels + labels2, frameon=False, loc="upper right")

    zumbach = summary.zumbach_rank_mean.to_numpy()
    zumbach_sd = summary.zumbach_rank_between_seed_sd.to_numpy()
    axes[2].fill_between(shift, zumbach - zumbach_sd, zumbach + zumbach_sd, color=CONTROL_COLOR, alpha=0.16)
    axes[2].plot(shift, zumbach, "o-", color=INK, linewidth=2.3, label="Mean rank Zumbach")
    axes[2].axhline(0.015, color=NEGATIVE_COLOR, linestyle="--", linewidth=1.5, label="Existing lower gate")
    axes[2].axvline(ASYMMETRIC_SHIFT, color=ASYMMETRIC_COLOR, linestyle=":")
    axes[2].set_xlabel(r"Fast-factor shift $\beta_J$")
    axes[2].set_ylabel("Rank Zumbach statistic")
    axes[2].set_title("Zumbach is the binding trade-off")
    axes[2].legend(frameon=False)
    _style_axis(axes[2])
    figure.suptitle("Asymmetry screen: benefit, joint viability, and binding constraint", color=INK, fontsize=17, fontweight="bold")
    figure.text(0.5, 0.01, "Means across 16 shortlisted configurations and five fresh seeds; shaded Zumbach band is +/-1 between-seed SD.", ha="center", color="#52626f")
    figure.tight_layout(rect=(0, 0.04, 1, 0.94))
    _save(figure, output / "08_m2_asymmetry_tradeoff.png")


def _write_dashboard(output: Path) -> None:
    """Create a portable index combining baseline and asymmetric diagnostics."""
    figures = (
        ("Asymmetric leverage event study", "05_m2_leverage_event_study.png"),
        ("Mechanism: shifted fast parabola", "07_m2_asymmetry_mechanism.png"),
        ("Asymmetry-versus-Zumbach trade-off", "08_m2_asymmetry_tradeoff.png"),
        ("Baseline paths", "01_m2_paths_overview.png"),
        ("Baseline maximum-volatility windows", "02_m2_maximum_volatility_windows.png"),
        ("Baseline return tails", "03_m2_return_tails.png"),
        ("Baseline dependence", "04_m2_dependence.png"),
        ("Baseline drawdown and volatility", "06_m2_drawdown_volatility.png"),
    )
    sections = "\n".join(
        f'<section><h2>{title}</h2><img src="{filename}" alt="{title}"></section>'
        for title, filename in figures
    )
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>M2 baseline and asymmetric extension</title>
<style>
body{{margin:0;background:#f4f6f8;color:#263744;font-family:Arial,sans-serif}}
main{{max-width:1500px;margin:auto;padding:28px}} h1,h2{{color:{INK}}}
.note{{background:white;border-left:5px solid {ASYMMETRIC_COLOR};padding:15px 18px;line-height:1.5}}
section{{background:white;margin:22px 0;padding:16px;border-radius:8px;box-shadow:0 1px 4px rgba(24,50,74,.09)}}
section img{{display:block;width:100%;height:auto}} h2{{margin:4px 8px 14px}}
</style></head><body><main>
<h1>M2 visual diagnostics: baseline and downside-asymmetric extension</h1>
<p class="note"><b>Comparison design.</b> The first three figures compare the exact zero-shift
control with <b>&beta;<sub>J</sub> = 0.10</b> using common random numbers. Tail-event dates are
fixed from the control returns. The remaining figures document the preserved Stage-3 baseline.
Volatility is 100 &times; &radic;(daily average instantaneous variance), not formal VIX. These are
synthetic mechanism diagnostics rather than market calibration.</p>
{sections}
</main></body></html>"""
    (output / "m2_visual_diagnostics_dashboard.html").write_text(html, encoding="utf-8")


def main() -> None:
    """Regenerate the paired asymmetry figures and their audit summary."""
    root = Path(__file__).resolve().parents[1]
    output = root / "figures"
    paired = _paired_simulations(root)
    lags, events = _collect_paired_events(paired)
    event_metrics = event_study_figure(output, lags, events)
    mechanism_figure(output, root)
    tradeoff_figure(output, root)
    _write_dashboard(output)
    manifest = {
        "model": "M2 with paired fast-factor shift",
        "control_spike_shift": CONTROL_SHIFT,
        "asymmetric_spike_shift": ASYMMETRIC_SHIFT,
        "stage3_candidate_id": CANDIDATE_ID,
        "structural_scenario_id": SCENARIO_ID,
        "fresh_seeds": list(FRESH_SEEDS),
        "event_selection": "control-return 2.5% tails; identical dates in both models",
        "volatility_definition": "100 * sqrt(daily average instantaneous variance)",
        "is_formal_vix": False,
        "event_metrics": event_metrics,
        "outputs": [
            "05_m2_leverage_event_study.png",
            "07_m2_asymmetry_mechanism.png",
            "08_m2_asymmetry_tradeoff.png",
            "m2_visual_diagnostics_dashboard.html",
        ],
    }
    (output / "m2_asymmetry_visual_summary.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
