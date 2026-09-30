"""Plot rolling physical-statistic estimates from physical_stability.py.

Requires plotly in addition to the preprocessing dependencies. The HTML output
embeds Plotly so it can be opened without a network connection.
"""

from __future__ import annotations

import argparse
from html import escape
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .market_data import ASSET_NAMES
from .physical_stability import OUT_DIR


PANELS = (
    ("hurst", "H from log GK volatility"),
    ("clustering_1_20", "Log-volatility ACF, mean lags 1-20"),
    ("leverage_1_5", "Leverage, mean lags 1-5"),
    ("zumbach_10", "10-day Zumbach proxy"),
    ("kurtosis", "Daily excess kurtosis"),
    ("taylor_gap_1_20", "Taylor gap, mean lags 1-20"),
    ("return_acf_1", "Lag-1 return ACF"),
    ("open_equal_previous_close_share", "Open = previous close, share of days"),
)


def plot_asset(input_dir: Path, output: Path, asset: str, interval_style: str = "bars",
               additional_inputs: tuple[Path, ...] = ()):
    estimate_frames = []
    check_frames = []
    seen_windows = set()
    for root in (input_dir, *additional_inputs):
        asset_dir = root / asset
        estimate_frame = pd.read_csv(asset_dir / "rolling_estimates.csv")
        check_frame = pd.read_csv(asset_dir / "rolling_source_checks.csv")
        window_years = set(estimate_frame.window_years.unique())
        if seen_windows & window_years:
            raise ValueError(f"Duplicate rolling window lengths in {asset_dir}")
        if window_years != set(check_frame.window_years.unique()):
            raise ValueError(f"Estimate/source-check window lengths differ in {asset_dir}")
        seen_windows.update(window_years)
        estimate_frames.append(estimate_frame)
        check_frames.append(check_frame)
    estimates = pd.concat(estimate_frames, ignore_index=True)
    checks = pd.concat(check_frames, ignore_index=True)
    windows = sorted(estimates.window_years.unique())
    if not windows:
        raise ValueError("No rolling windows found")
    palette = ("#1766a3", "#db7415", "#669933", "#8c55aa")
    colors = {window: palette[i % len(palette)] for i, window in enumerate(windows)}
    fig = make_subplots(rows=4, cols=2, subplot_titles=[label for _, label in PANELS],
                        vertical_spacing=0.085, horizontal_spacing=0.12)
    for panel, (metric, _) in enumerate(PANELS):
        row, col = divmod(panel, 2)
        has_interval = metric != "open_equal_previous_close_share"
        table = estimates[estimates.metric == metric] if has_interval else checks
        value = "estimate" if has_interval else metric
        for window in windows:
            data = table[table.window_years == window].sort_values("window_end_year")
            error = None
            if has_interval:
                if interval_style == "band":
                    rgb = tuple(int(colors[window][i:i + 2], 16) for i in (1, 3, 5))
                    fig.add_trace(go.Scatter(
                        x=[*data.window_end_year, *data.window_end_year.iloc[::-1]],
                        y=[*data.ci_high, *data.ci_low.iloc[::-1]],
                        mode="lines", fill="toself",
                        fillcolor=f"rgba({rgb[0]},{rgb[1]},{rgb[2]},0.14)",
                        line=dict(color="rgba(0,0,0,0)"),
                        name="Pointwise 95% interval", legendgroup=str(window),
                        showlegend=False, hoverinfo="skip",
                    ), row=row + 1, col=col + 1)
                else:
                    error = dict(type="data", symmetric=False,
                                 array=(data.ci_high - data.estimate).to_numpy(),
                                 arrayminus=(data.estimate - data.ci_low).to_numpy(),
                                 color=colors[window], thickness=1, width=2)
            hover = ("%{customdata[0]}-%{customdata[1]}<br>"
                     "Estimate: %{y:.3f}<br>Returns: %{customdata[2]}")
            fields = ["window_start_year", "window_end_year", "n_returns"]
            if has_interval:
                fields.extend(("ci_low", "ci_high"))
                hover += ("<br>95% interval: [%{customdata[3]:.3f}, "
                          "%{customdata[4]:.3f}]")
            fig.add_trace(go.Scatter(
                x=data.window_end_year, y=data[value], mode="lines+markers",
                name=f"{window}-year window", legendgroup=str(window),
                showlegend=panel == 0, line=dict(color=colors[window], width=2),
                marker=dict(size=5), error_y=error,
                customdata=data[fields],
                hovertemplate=hover + "<extra></extra>",
            ), row=row + 1, col=col + 1)
        fig.update_xaxes(title_text="Window end year", row=row + 1, col=col + 1)
    fig.update_layout(
        template="plotly_white", height=1450, width=1300,
        margin=dict(t=115, b=60, l=60, r=30),
        legend=dict(orientation="h", x=0.5, xanchor="center",
                    y=1.04, yanchor="bottom", font=dict(size=16)),
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    chart = fig.to_html(include_plotlyjs=True, full_html=False)
    title = escape(f"{asset}: rolling physical statistics from daily OHLC data")
    interval_description = ("Shaded areas" if interval_style == "band" else "Bars")
    output.write_text(f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  body {{ margin: 0; background: white; color: #26384d; font-family: Arial, sans-serif; }}
  .report-header {{ max-width: 1240px; margin: 28px auto 0; padding: 0 24px; }}
  .report-header h1 {{ margin: 0 0 12px; font-size: 27px; line-height: 1.2; }}
  .report-header p {{ margin: 0; font-size: 16px; line-height: 1.45; }}
  .chart {{ width: 1300px; max-width: 100%; margin: 0 auto; }}
</style>
</head>
<body>
<header class="report-header">
  <h1>{title}</h1>
  <p>{interval_description} are pointwise 95% block-bootstrap intervals; overlapping windows are dependent.
     Check the recorded-open panel for source changes.</p>
</header>
<main class="chart">{chart}</main>
</body>
</html>
""", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset", choices=ASSET_NAMES, default="SP500")
    parser.add_argument("--input", type=Path, default=OUT_DIR)
    parser.add_argument("--additional-input", type=Path, action="append", default=[],
                        help="Another experiment root with distinct rolling window lengths")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--interval-style", choices=("bars", "band"), default="bars")
    args = parser.parse_args()
    output = args.output or args.input / args.asset / "rolling_estimates.html"
    plot_asset(args.input, output, args.asset, args.interval_style,
               tuple(args.additional_input))
    print(output)


if __name__ == "__main__":
    main()
