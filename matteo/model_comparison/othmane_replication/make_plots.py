"""Create offline Plotly diagnostics and a simulation-only bootstrap summary."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from compare import load_legacy


HERE = Path(__file__).resolve().parent
MODELS = ("othmane_archived_fit", "m2_stage3_anchor")
LABELS = {MODELS[0]: "Othmane saved fit", MODELS[1]: "M2 synthetic anchor"}
COLORS = {MODELS[0]: "#2563eb", MODELS[1]: "#e24a33"}
FACTS = ("Hhat", "leverage", "abs_acf", "sq_acf", "zumbach", "kurtosis")
CURVES = (
    ("leverage", list(range(1, 41)), "Leverage: corr(returnₜ, varianceₜ₊ₗ)"),
    ("abs_acf", [1, 3, 5, 10, 15, 20], "Absolute-return ACF"),
    ("sq_acf", [1, 3, 5, 10, 15, 20], "Squared-return ACF"),
    ("zumbach", [2, 10, 20, 30, 40], "Pearson Zumbach"),
)
RETURN_FACTS = ("abs_acf", "sq_acf", "kurtosis")


def bootstrap_scores(legacy, record, target, rng, reps):
    arrays = {fact: np.asarray(record["path_statistics"][fact]) for fact in FACTS}
    count = len(arrays["Hhat"])
    full = np.empty(reps)
    return_only = np.empty(reps)
    for rep in range(reps):
        indices = rng.integers(0, count, count)
        average = {fact: np.mean(values[indices], axis=0) for fact, values in arrays.items()}
        residual = legacy.joint_residual(average, target)
        full[rep] = np.sum(residual ** 2)
        return_only[rep] = np.sum(residual[41:53] ** 2) + residual[58] ** 2
    return full, return_only


def summary(data, reps):
    legacy, _ = load_legacy()
    rng = np.random.default_rng(380521)
    rows = []
    bootstrap = {}
    for asset, entry in data["assets"].items():
        target = {fact: np.asarray(entry["target"][fact]) for fact in FACTS}
        bootstrap[asset] = {}
        for model in MODELS:
            record = entry["models"][model]
            full, return_only = bootstrap_scores(legacy, record, target, rng, reps)
            bootstrap[asset][model] = {"full": full, "return_only": return_only}
            blocks = record["weighted_squared_error_by_block"]
            row = {
                "asset": asset, "model": model,
                "weighted_sse": record["weighted_squared_error_total"],
                "weighted_sse_boot_p5": float(np.percentile(full, 5)),
                "weighted_sse_boot_p95": float(np.percentile(full, 95)),
                "return_only_sse": sum(blocks[fact] for fact in RETURN_FACTS),
                "return_only_sse_boot_p5": float(np.percentile(return_only, 5)),
                "return_only_sse_boot_p95": float(np.percentile(return_only, 95)),
                "mean_pathwise_sse": record["mean_pathwise_weighted_squared_error"],
                "Hhat": record["statistics"]["Hhat"],
                "Hhat_common_lags_2_40": record["Hhat_common_lags_2_40"],
                "kurtosis": record["statistics"]["kurtosis"],
            }
            row.update({f"sse_{fact}": blocks[fact] for fact in FACTS})
            rows.append(row)
    return rows, bootstrap


def save_html(fig, path):
    fig.write_html(path, include_plotlyjs="directory", full_html=True,
                   config={"responsive": True, "displaylogo": False})


def overview(data, rows, output):
    assets = list(data["assets"])
    lookup = {(row["asset"], row["model"]): row for row in rows}
    fig = make_subplots(rows=2, cols=2, vertical_spacing=0.15,
                        subplot_titles=("Full 59-term weighted SSE",
                                        "Return-only weighted SSE: ACFs + kurtosis",
                                        "Roughness H estimate", "Excess return kurtosis"))
    for model in MODELS:
        color = COLORS[model]
        values = [lookup[(asset, model)]["weighted_sse"] for asset in assets]
        upper = [lookup[(asset, model)]["weighted_sse_boot_p95"] - value
                 for asset, value in zip(assets, values)]
        lower = [value - lookup[(asset, model)]["weighted_sse_boot_p5"]
                 for asset, value in zip(assets, values)]
        fig.add_trace(go.Bar(x=assets, y=values, name=LABELS[model], marker_color=color,
                             error_y=dict(type="data", symmetric=False, array=upper,
                                          arrayminus=lower), legendgroup=model), 1, 1)
        values = [lookup[(asset, model)]["return_only_sse"] for asset in assets]
        upper = [lookup[(asset, model)]["return_only_sse_boot_p95"] - value
                 for asset, value in zip(assets, values)]
        lower = [value - lookup[(asset, model)]["return_only_sse_boot_p5"]
                 for asset, value in zip(assets, values)]
        fig.add_trace(go.Bar(x=assets, y=values, name=LABELS[model], marker_color=color,
                             error_y=dict(type="data", symmetric=False, array=upper,
                                          arrayminus=lower), legendgroup=model,
                             showlegend=False), 1, 2)
    for fact, col in (("Hhat", 1), ("kurtosis", 2)):
        fig.add_trace(go.Scatter(x=assets,
                                 y=[data["assets"][a]["target"][fact] for a in assets],
                                 mode="lines+markers", name="Historical target",
                                 line=dict(color="#1f2937", width=2.5),
                                 legendgroup="target", showlegend=(col == 1)), 2, col)
        for model in MODELS:
            fig.add_trace(go.Scatter(x=assets,
                                     y=[lookup[(a, model)][fact] for a in assets],
                                     mode="lines+markers", name=LABELS[model],
                                     line=dict(color=COLORS[model], width=2),
                                     legendgroup=model, showlegend=False), 2, col)
    fig.update_layout(title={"text": "Nine-asset diagnostic: saved ESN fit vs unfitted M2 anchor"
                             "<br><sup>20 fresh paths × 4,000 retained days; 5–95% bars bootstrap simulation paths only."
                             " Historical targets were used to fit ESN; variance proxies differ.</sup>"},
                      template="plotly_white", height=950, width=1450,
                      barmode="group", legend=dict(orientation="h", y=1.06),
                      margin=dict(t=145, b=90))
    fig.update_xaxes(tickangle=-35, row=1, col=1)
    fig.update_xaxes(tickangle=-35, row=1, col=2)
    fig.update_xaxes(tickangle=-35, row=2, col=1)
    fig.update_xaxes(tickangle=-35, row=2, col=2)
    save_html(fig, output / "overview.html")


def block_differences(data, output):
    assets = list(data["assets"])
    z = []
    text = []
    for fact in FACTS:
        values = []
        labels = []
        for asset in assets:
            models = data["assets"][asset]["models"]
            esn = models[MODELS[0]]["weighted_squared_error_by_block"][fact]
            m2 = models[MODELS[1]]["weighted_squared_error_by_block"][fact]
            values.append(m2 - esn)
            labels.append(f"ESN {esn:.3f}<br>M2 {m2:.3f}")
        z.append(values)
        text.append(labels)
    limit = max(abs(np.min(z)), abs(np.max(z)))
    fig = go.Figure(go.Heatmap(x=assets, y=list(FACTS), z=z, customdata=text,
                               zmin=-limit, zmax=limit, zmid=0, colorscale="RdBu_r",
                               colorbar=dict(title="M2 − ESN SSE"),
                               hovertemplate="%{x} · %{y}<br>%{customdata}"
                                             "<br>Difference %{z:.3f}<extra></extra>"))
    fig.update_layout(title={"text": "Which weighted facts drive the score difference?"
                             "<br><sup>Positive means the saved ESN fit has lower error."
                             " M2 is the same synthetic anchor for every asset.</sup>"},
                      template="plotly_white", height=620, width=1250,
                      xaxis=dict(tickangle=-35), margin=dict(t=100, b=90))
    save_html(fig, output / "block_differences.html")


def curves(data, output):
    for asset, entry in data["assets"].items():
        fig = make_subplots(rows=2, cols=2,
                            subplot_titles=[item[2] for item in CURVES],
                            vertical_spacing=0.13, horizontal_spacing=0.09)
        for index, (fact, lags, _) in enumerate(CURVES):
            row, col = divmod(index, 2)
            row += 1
            col += 1
            fig.add_trace(go.Scatter(x=lags, y=entry["target"][fact],
                                     mode="lines+markers", name="Historical target",
                                     line=dict(color="#1f2937", width=2.8),
                                     marker=dict(size=5), legendgroup="target",
                                     showlegend=(index == 0)), row, col)
            for model in MODELS:
                result = entry["models"][model]
                color = COLORS[model]
                lower = result["path_band_p5"][fact]
                upper = result["path_band_p95"][fact]
                fig.add_trace(go.Scatter(x=lags, y=lower, mode="lines",
                                         line=dict(width=0), hoverinfo="skip",
                                         showlegend=False, legendgroup=model), row, col)
                fig.add_trace(go.Scatter(x=lags, y=upper, mode="lines",
                                         line=dict(width=0), fill="tonexty",
                                         fillcolor=("rgba(37,99,235,0.13)" if model == MODELS[0]
                                                    else "rgba(226,74,51,0.13)"),
                                         hoverinfo="skip", showlegend=False,
                                         legendgroup=model), row, col)
                fig.add_trace(go.Scatter(x=lags, y=result["statistics"][fact],
                                         mode="lines+markers", name=LABELS[model],
                                         line=dict(color=color, width=2),
                                         marker=dict(size=4), legendgroup=model,
                                         showlegend=(index == 0)), row, col)
            fig.update_xaxes(title_text="Lag (trading days)", row=row, col=col)
        target = entry["target"]
        esn = entry["models"][MODELS[0]]["statistics"]
        m2 = entry["models"][MODELS[1]]["statistics"]
        scalar_note = (f"H: target {target['Hhat']:.3f}, ESN {esn['Hhat']:.3f}, "
                       f"M2 {m2['Hhat']:.3f} · Excess kurtosis: target "
                       f"{target['kurtosis']:.2f}, ESN {esn['kurtosis']:.2f}, "
                       f"M2 {m2['kurtosis']:.2f}")
        fig.update_layout(title={"text": f"{asset}: historical facts and fresh model paths"
                                 f"<br><sup>{scalar_note}</sup>"},
                          template="plotly_white", height=900, width=1350,
                          hovermode="closest", legend=dict(orientation="h", y=1.08),
                          margin=dict(t=145, b=75))
        save_html(fig, output / f"curves_{asset.lower()}.html")


def late_return_screen(data, screen, output, reps):
    """Show the early-selected readout on later historical return targets."""
    legacy, _ = load_legacy()
    rng = np.random.default_rng(480522)
    assets = list(data["assets"])
    series = {"othmane_archived_fit": [], "m2_stage3_anchor": [],
              "m2_early_selected": []}
    comparison = {}
    rows = []
    for asset in assets:
        target = {fact: np.asarray(screen["targets"][asset]["late"][fact])
                  for fact in RETURN_FACTS}
        selected_id = str(screen["comparison"][asset]["selected_candidate_id"])
        path_data = {
            "othmane_archived_fit": data["assets"][asset]["models"][MODELS[0]]["path_statistics"],
            "m2_stage3_anchor": data["assets"][asset]["models"][MODELS[1]]["path_statistics"],
            "m2_early_selected": screen["evaluation_path_statistics"][selected_id],
        }
        boot = {}
        for model, all_stats in path_data.items():
            arrays = {fact: np.asarray(all_stats[fact]) for fact in RETURN_FACTS}
            average = {fact: np.mean(values, axis=0) for fact, values in arrays.items()}
            tol = legacy.TOLERANCES

            def score(stats):
                return float(
                    np.sum(((stats["abs_acf"] - target["abs_acf"]) /
                            (tol["abs_acf"] * np.sqrt(6))) ** 2)
                    + np.sum(((stats["sq_acf"] - target["sq_acf"]) /
                              (tol["sq_acf"] * np.sqrt(6))) ** 2)
                    + ((stats["kurtosis"] - target["kurtosis"]) /
                       tol["kurtosis"]) ** 2
                )

            value = score(average)
            draws = np.empty(reps)
            count = len(arrays["kurtosis"])
            for rep in range(reps):
                ix = rng.integers(0, count, count)
                draws[rep] = score({fact: np.mean(values[ix], axis=0)
                                    for fact, values in arrays.items()})
            boot[model] = draws
            entry = {"asset": asset, "model": model, "late_return_sse": value,
                     "bootstrap_p5": float(np.percentile(draws, 5)),
                     "bootstrap_p95": float(np.percentile(draws, 95)),
                     "selected_candidate_id": selected_id}
            series[model].append(entry)
            rows.append(entry)
        comparison[asset] = {
            "bootstrap_probability_selected_below_anchor": float(np.mean(
                boot["m2_early_selected"] < boot["m2_stage3_anchor"])),
            "bootstrap_probability_selected_below_archived_esn": float(np.mean(
                boot["m2_early_selected"] < boot["othmane_archived_fit"])),
        }

    with (output / "late_return_scores.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output / "late_return_bootstrap.json").write_text(
        json.dumps({"replicates": reps, "seed": 480522,
                    "method": "resample simulated paths independently; later historical target fixed",
                    "probabilities": comparison}, indent=2) + "\n")

    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("Before 2015: candidate screening score",
                                        "After 2014: fresh-seed return-only SSE"),
                        horizontal_spacing=0.10)
    anchor_early = [screen["comparison"][a]["early_screen_score_anchor"] for a in assets]
    selected_early = [screen["comparison"][a]["early_screen_score_selected"] for a in assets]
    fig.add_trace(go.Bar(x=assets, y=anchor_early, name="M2 synthetic anchor",
                         marker_color=COLORS[MODELS[1]], legendgroup="anchor"), 1, 1)
    fig.add_trace(go.Bar(x=assets, y=selected_early, name="M2 early selected",
                         marker_color="#16a085", legendgroup="selected"), 1, 1)
    for model, label, color in ((MODELS[0], "ESN full-sample fit", COLORS[MODELS[0]]),
                                (MODELS[1], "M2 synthetic anchor", COLORS[MODELS[1]]),
                                ("m2_early_selected", "M2 early selected", "#16a085")):
        entries = series[model]
        fig.add_trace(go.Bar(x=assets, y=[e["late_return_sse"] for e in entries],
                             name=label, marker_color=color, legendgroup=model,
                             showlegend=(model == MODELS[0]),
                             error_y=dict(type="data", symmetric=False,
                                          array=[e["bootstrap_p95"] - e["late_return_sse"]
                                                 for e in entries],
                                          arrayminus=[e["late_return_sse"] - e["bootstrap_p5"]
                                                      for e in entries])), 1, 2)
    fig.update_layout(title={"text": "M2 readout selection on early historical returns"
                             "<br><sup>32 existing readouts; pre-2015 selection, later historical targets, independent simulation seeds."
                             " ESN parameters used the full historical sample. Bars show simulation bootstrap 5–95%.</sup>"},
                      barmode="group", template="plotly_white", height=620, width=1550,
                      legend=dict(orientation="h", y=1.10), margin=dict(t=140, b=100))
    fig.update_xaxes(tickangle=-35)
    save_html(fig, output / "readout_screen.html")


def lag_sensitivity(lag_data, output):
    assets = list(lag_data["assets"])
    fig = make_subplots(rows=2, cols=2,
                        subplot_titles=("H on lags 1–40", "H on lags 2–40",
                                        "Full SSE with common lags 1–40",
                                        "Full SSE with common lags 2–40"),
                        vertical_spacing=0.15, horizontal_spacing=0.09)
    for col, band in ((1, "1_40"), (2, "2_40")):
        fig.add_trace(go.Scatter(x=assets,
                                 y=[lag_data["assets"][a][f"historical_H_lags_{band}"]
                                    for a in assets], mode="lines+markers",
                                 name="Historical GK", line=dict(color="#1f2937", width=2.5),
                                 legendgroup="target", showlegend=(col == 1)), 1, col)
        for model in MODELS:
            fig.add_trace(go.Scatter(
                x=assets,
                y=[lag_data["assets"][a]["models"][model][f"H_lags_{band}"]
                   for a in assets], mode="lines+markers", name=LABELS[model],
                line=dict(color=COLORS[model], width=2), legendgroup=model,
                showlegend=(col == 1)), 1, col)
            fig.add_trace(go.Bar(
                x=assets,
                y=[lag_data["assets"][a]["models"][model][f"full_sse_common_{band}"]
                   for a in assets], name=LABELS[model], marker_color=COLORS[model],
                legendgroup=model, showlegend=False), 2, col)
    fig.update_layout(title={"text": "Roughness lag-band sensitivity"
                             "<br><sup>Both historical GK and model variance use the displayed lag band."
                             " Variance observations still differ across data and models.</sup>"},
                      barmode="group", template="plotly_white", width=1450,
                      height=930, margin=dict(t=140, b=95),
                      legend=dict(orientation="h", y=1.07))
    fig.update_xaxes(tickangle=-35)
    save_html(fig, output / "roughness_lag.html")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=HERE / "nine_assets_fresh_001.json")
    parser.add_argument("--output-dir", type=Path, default=HERE / "plots")
    parser.add_argument("--bootstrap-reps", type=int, default=400)
    parser.add_argument("--screen", type=Path, default=HERE / "m2_readout_screen_001.json")
    parser.add_argument("--lag-audit", type=Path,
                        default=HERE / "roughness_lag_sensitivity_001.json")
    args = parser.parse_args()
    if args.bootstrap_reps < 100:
        parser.error("Use at least 100 bootstrap replicates")
    data = json.loads(args.input.read_text())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows, bootstrap = summary(data, args.bootstrap_reps)
    csv_path = args.output_dir / "nine_asset_scores.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    comparisons = {}
    for asset, scores in bootstrap.items():
        esn = scores[MODELS[0]]["full"]
        m2 = scores[MODELS[1]]["full"]
        er = scores[MODELS[0]]["return_only"]
        mr = scores[MODELS[1]]["return_only"]
        comparisons[asset] = {
            "simulation_bootstrap_probability_m2_lower_full_sse": float(np.mean(m2 < esn)),
            "simulation_bootstrap_probability_m2_lower_return_only_sse": float(np.mean(mr < er)),
        }
    (args.output_dir / "bootstrap_comparisons.json").write_text(
        json.dumps({"replicates": args.bootstrap_reps, "seed": 380521,
                    "method": "resample model paths independently, keep historical target fixed",
                    "comparisons": comparisons}, indent=2) + "\n")
    overview(data, rows, args.output_dir)
    block_differences(data, args.output_dir)
    curves(data, args.output_dir)
    if args.screen.exists():
        late_return_screen(data, json.loads(args.screen.read_text()),
                           args.output_dir, args.bootstrap_reps)
    if args.lag_audit.exists():
        lag_sensitivity(json.loads(args.lag_audit.read_text()), args.output_dir)
    print(args.output_dir)


if __name__ == "__main__":
    main()
