"""Generate reproducible Plotly diagnostics for the selected Stage-3 M2 model.

The figures are deliberately qualitative.  The plotted volatility series is
100 * sqrt(daily average instantaneous variance), not a model-implied VIX.
"""

from __future__ import annotations

import copy
import json
import warnings
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from numpy.typing import NDArray
from plotly.subplots import make_subplots
from scipy.stats import kurtosis, norm, skew

from esn_eq.config import (
    ArchitectureConfig,
    MarketEnvironment,
    ReadoutParameters,
    RiskPremium,
    SimulationConfig,
)
from esn_eq.model import ExogenousReservoirVolatilityModel, SimulationResult
from esn_eq.states import ReservoirFactory

FloatArray = NDArray[np.float64]

ASSET_COLOR = "#17324D"
VOL_COLOR = "#D55E00"
NEGATIVE_COLOR = "#C43B45"
POSITIVE_COLOR = "#2A7F62"
NEUTRAL_COLOR = "#68727D"
LIGHT_GRID = "#E5E9ED"


@dataclass(frozen=True, slots=True)
class SelectedM2Configuration:
    """Fully specified Stage-3 M2 selection and its simulation protocol."""

    architecture: ArchitectureConfig
    simulation: SimulationConfig
    parameters: ReadoutParameters
    premium: RiskPremium
    market: MarketEnvironment
    seeds: tuple[int, ...]
    candidate_id: int
    scenario_id: int


@dataclass(frozen=True, slots=True)
class SeedSimulation:
    """One independent cache seed and its full collection of physical paths."""

    seed: int
    result: SimulationResult


@dataclass(frozen=True, slots=True)
class EventStudy:
    """Robust volatility response around positive and negative return extremes."""

    lags: NDArray[np.int64]
    negative: FloatArray
    positive: FloatArray


class Stage3M2Loader:
    """Recover the untouched Stage-3 winner from the experiment audit files."""

    def __init__(self, experiment_dir: Path) -> None:
        """Store the directory containing the manifest and ranked configurations."""
        self.experiment_dir = experiment_dir

    def load(self, candidate_id: int = 18, scenario_id: int = 0) -> SelectedM2Configuration:
        """Load M2 and explicitly zero both optional readout channels."""
        manifest = self._read_manifest()
        selected = self._read_selected_row(candidate_id, scenario_id)
        architecture = self._architecture(manifest, selected)
        return SelectedM2Configuration(
            architecture=architecture,
            simulation=SimulationConfig(**manifest["simulation"]),
            parameters=self._parameters(selected),
            premium=self._premium(selected),
            market=MarketEnvironment(**manifest["market"]),
            seeds=tuple(int(seed) for seed in manifest["seeds"]),
            candidate_id=candidate_id,
            scenario_id=scenario_id,
        )

    def _read_manifest(self) -> dict[str, Any]:
        """Read the complete Stage-3 audit manifest."""
        path = self.experiment_dir / "manifest.json"
        return json.loads(path.read_text(encoding="utf-8"))

    def _read_selected_row(self, candidate_id: int, scenario_id: int) -> pd.Series:
        """Select the predeclared M2 candidate/scenario without visual screening."""
        table = pd.read_csv(self.experiment_dir / "top_configurations.csv")
        mask = (
            table["variant"].eq("M2")
            & table["candidate_id"].eq(candidate_id)
            & table["structural_scenario_id"].eq(scenario_id)
        )
        rows = table.loc[mask]
        if len(rows) != 1:
            raise ValueError("Expected exactly one selected M2 configuration row.")
        return rows.iloc[0]

    @staticmethod
    def _architecture(manifest: dict[str, Any], selected: pd.Series) -> ArchitectureConfig:
        """Apply the selected structural correlations to the frozen architecture."""
        values = copy.deepcopy(manifest["architecture"])
        values["feedback"]["asset_correlation"] = float(
            selected["feedback_asset_correlation"]
        )
        values["spike"]["asset_correlation"] = float(selected["spike_asset_correlation"])
        return ArchitectureConfig.from_mapping(values)

    @staticmethod
    def _parameters(selected: pd.Series) -> ReadoutParameters:
        """Build the selected readout while enforcing the M2 cell definition."""
        names = {field.name for field in fields(ReadoutParameters)}
        values = {name: float(selected[name]) for name in names if name in selected.index}
        values.setdefault("spike_shift", 0.0)
        values["orthogonal_curvature"] = 0.0
        values["echo_loading"] = 0.0
        return ReadoutParameters.from_mapping(values)

    @staticmethod
    def _premium(selected: pd.Series) -> RiskPremium:
        """Recover the deterministic primitive-driver market prices of risk."""
        return RiskPremium(
            equity=float(selected["equity_risk_price"]),
            feedback_idiosyncratic=float(
                selected["feedback_idiosyncratic_risk_price"]
            ),
            spike_idiosyncratic=float(selected["spike_idiosyncratic_risk_price"]),
        )


class M2SimulationSuite:
    """Simulate the selected configuration on every declared independent seed."""

    def __init__(self, configuration: SelectedM2Configuration) -> None:
        """Construct one reusable model and reservoir factory."""
        self.configuration = configuration
        self.factory = ReservoirFactory(configuration.architecture, configuration.simulation)
        self.model = ExogenousReservoirVolatilityModel(configuration.architecture)

    def run(self) -> tuple[SeedSimulation, ...]:
        """Generate physical-measure paths without selecting on their appearance."""
        simulations: list[SeedSimulation] = []
        for seed in self.configuration.seeds:
            cache = self.factory.build(seed)
            prepared = self.model.prepare(cache, self.configuration.premium)
            result = self.model.simulate(
                prepared,
                self.configuration.parameters,
                self.configuration.market,
                measure="P",
            )
            simulations.append(SeedSimulation(seed, result))
        return tuple(simulations)


class PlotlyVisualDiagnostics:
    """Build a compact visual audit of path shape, tails, memory, and asymmetry."""

    def __init__(
        self,
        simulations: Sequence[SeedSimulation],
        output_dir: Path,
        display_path: int = 0,
    ) -> None:
        """Use a fixed path index for displays and all paths for pooled diagnostics."""
        if not simulations:
            raise ValueError("At least one seed simulation is required.")
        self.simulations = tuple(simulations)
        self.output_dir = output_dir
        self.display_path = display_path
        self.observations_per_year = simulations[0].result.observations_per_year
        self.figures: list[tuple[str, go.Figure]] = []
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate(self) -> dict[str, Path]:
        """Create figures, dashboard, plotted-path data, and numerical summaries."""
        self.figures = [
            ("Five-seed paths", self.paths_overview()),
            ("Maximum-volatility windows", self.spike_windows()),
            ("Return tails", self.return_tails()),
            ("Dependence", self.dependence()),
            ("Leverage event study", self.leverage_event_study()),
            ("Drawdown and volatility", self.drawdown_volatility()),
        ]
        outputs = self._export_figures()
        outputs["dashboard"] = self._write_dashboard()
        outputs["display_data"] = self._write_display_data()
        outputs["summary"] = self._write_summary()
        return outputs

    def paths_overview(self) -> go.Figure:
        """Show one preselected asset/volatility path for every independent seed."""
        titles = [f"Seed {item.seed} · fixed path {self.display_path}" for item in self.simulations]
        specs = [[{"secondary_y": True}] for _ in self.simulations]
        figure = make_subplots(
            rows=len(self.simulations), cols=1, specs=specs, subplot_titles=titles,
            vertical_spacing=0.055,
        )
        for row, item in enumerate(self.simulations, start=1):
            years, asset, volatility = self._display_series(item)
            figure.add_trace(self._line(years, asset, "Asset index", ASSET_COLOR, row == 1),
                             row=row, col=1, secondary_y=False)
            figure.add_trace(self._line(years, volatility, "Volatility proxy", VOL_COLOR, row == 1),
                             row=row, col=1, secondary_y=True)
            figure.update_yaxes(title_text="Asset (log)", type="log", row=row, col=1,
                                secondary_y=False)
            figure.update_yaxes(title_text="Vol. (%)", rangemode="tozero", row=row, col=1,
                                secondary_y=True)
        figure.update_xaxes(title_text="Simulated years after burn-in", row=len(titles), col=1)
        return self._style(
            figure,
            "M2 asset and daily volatility proxy across five independent seeds",
            "One fixed path index per seed; asset starts at 100. Volatility is 100√RV, not VIX.",
            height=1900,
        )

    def spike_windows(self) -> go.Figure:
        """Align each displayed path on its objectively selected maximum-volatility day."""
        titles = [f"Seed {item.seed} · maximum eligible volatility event" for item in self.simulations]
        specs = [[{"secondary_y": True}] for _ in self.simulations]
        figure = make_subplots(
            rows=len(self.simulations), cols=1, specs=specs, subplot_titles=titles,
            vertical_spacing=0.055,
        )
        for row, item in enumerate(self.simulations, start=1):
            event_days, asset, volatility = self._spike_series(item)
            figure.add_trace(self._line(event_days, asset, "Reindexed asset", ASSET_COLOR, row == 1),
                             row=row, col=1, secondary_y=False)
            figure.add_trace(self._line(event_days, volatility, "Volatility proxy", VOL_COLOR, row == 1),
                             row=row, col=1, secondary_y=True)
            figure.add_vline(x=0, line_dash="dot", line_color=NEUTRAL_COLOR, row=row, col=1)
            figure.update_yaxes(title_text="Asset", row=row, col=1, secondary_y=False)
            figure.update_yaxes(title_text="Vol. (%)", rangemode="tozero", row=row, col=1,
                                secondary_y=True)
        figure.update_xaxes(title_text="Trading days from volatility maximum", row=len(titles), col=1)
        return self._style(
            figure,
            "Maximum-volatility episodes: rise, asset move, and decay",
            "Maximum with a complete 126-day pre/252-day post window on each fixed path; no visual selection.",
            height=1900,
        )

    def return_tails(self) -> go.Figure:
        """Contrast pooled standardized returns with a Gaussian density and QQ line."""
        standardized = self._standardized_returns()
        bound = float(np.clip(np.ceil(np.quantile(np.abs(standardized), 0.9995)), 6, 16))
        edges = np.linspace(-bound, bound, 161)
        density, _ = np.histogram(standardized, bins=edges, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        probabilities = np.linspace(0.0005, 0.9995, 501)
        theoretical = norm.ppf(probabilities)
        empirical = np.quantile(standardized, probabilities)
        figure = make_subplots(rows=1, cols=2, subplot_titles=("Log-density", "Normal QQ"))
        plotted_density = np.where(density > 0.0, density, np.nan)
        figure.add_trace(go.Scatter(x=centers, y=plotted_density, mode="lines", name="M2 returns",
                                    line={"color": ASSET_COLOR, "width": 2}), row=1, col=1)
        figure.add_trace(go.Scatter(x=centers, y=norm.pdf(centers), mode="lines", name="Gaussian",
                                    line={"color": NEUTRAL_COLOR, "dash": "dash"}), row=1, col=1)
        figure.add_trace(go.Scatter(x=theoretical, y=empirical, mode="markers", name="Empirical quantiles",
                                    marker={"color": VOL_COLOR, "size": 5}), row=1, col=2)
        line_bound = max(abs(theoretical[0]), abs(empirical[0]), abs(empirical[-1]))
        figure.add_trace(go.Scatter(x=[-line_bound, line_bound], y=[-line_bound, line_bound],
                                    mode="lines", name="Gaussian reference",
                                    line={"color": NEUTRAL_COLOR, "dash": "dash"}), row=1, col=2)
        figure.update_yaxes(type="log", title_text="Density (log scale)", row=1, col=1)
        figure.update_xaxes(title_text="Return / path standard deviation", row=1, col=1)
        figure.update_xaxes(title_text="Gaussian quantile", row=1, col=2)
        figure.update_yaxes(title_text="M2 quantile", row=1, col=2)
        return self._style(
            figure,
            "Heavy tails in pooled daily M2 returns",
            f"{standardized.size:,} observations; standardized separately within each model path.",
            height=650,
        )

    def dependence(self) -> go.Figure:
        """Display return efficiency and volatility/absolute-return persistence."""
        returns, variance = self._pooled_matrices()
        return_acf = self._acf_matrix(returns, 30)
        absolute_acf = self._acf_matrix(np.abs(returns), 60)
        squared_acf = self._acf_matrix(returns**2, 60)
        log_variance_acf = self._acf_matrix(np.log(np.maximum(variance, 1e-15)), 60)
        figure = make_subplots(rows=1, cols=2, subplot_titles=("Return autocorrelation", "Volatility persistence"))
        self._add_acf_band(figure, return_acf, "Returns", ASSET_COLOR, row=1, col=1)
        self._add_acf_band(figure, absolute_acf, "Absolute returns", POSITIVE_COLOR, row=1, col=2)
        self._add_acf_band(figure, squared_acf, "Squared returns", VOL_COLOR, row=1, col=2)
        self._add_acf_band(figure, log_variance_acf, "Log variance", ASSET_COLOR, row=1, col=2)
        figure.add_hline(y=0, line_color=NEUTRAL_COLOR, line_width=1, row=1, col=1)
        figure.add_hline(y=0, line_color=NEUTRAL_COLOR, line_width=1, row=1, col=2)
        figure.update_xaxes(title_text="Lag (trading days)")
        figure.update_yaxes(title_text="Autocorrelation")
        return self._style(
            figure,
            "Weak return memory with persistent volatility proxies",
            f"Lines are means over {returns.shape[0]} paths; shaded bands span the 10th–90th path percentiles.",
            height=650,
        )

    def leverage_event_study(self) -> go.Figure:
        """Compare volatility around equally extreme negative and positive returns."""
        study = self._event_study()
        negative = self._quantile_curves(study.negative)
        positive = self._quantile_curves(study.positive)
        gap = negative[1] - positive[1]
        figure = make_subplots(rows=2, cols=1, shared_xaxes=True,
                               subplot_titles=("Tail-event volatility response", "Negative-minus-positive median"),
                               row_heights=[0.72, 0.28], vertical_spacing=0.12)
        self._add_event_band(figure, study.lags, negative, "Negative 2.5% tail", NEGATIVE_COLOR, 1)
        self._add_event_band(figure, study.lags, positive, "Positive 2.5% tail", POSITIVE_COLOR, 1)
        figure.add_trace(go.Scatter(x=study.lags, y=gap, mode="lines", name="Asymmetry gap",
                                    line={"color": ASSET_COLOR, "width": 2}), row=2, col=1)
        figure.add_hline(y=0, line_color=NEUTRAL_COLOR, line_width=1, row=2, col=1)
        figure.add_vline(x=0, line_dash="dot", line_color=NEUTRAL_COLOR, row=1, col=1)
        figure.add_vline(x=0, line_dash="dot", line_color=NEUTRAL_COLOR, row=2, col=1)
        figure.update_yaxes(title_text="Vol. vs pre-event baseline (%)", row=1, col=1)
        figure.update_yaxes(title_text="Gap (pp)", row=2, col=1)
        figure.update_xaxes(title_text="Trading days from tail return", row=2, col=1)
        return self._style(
            figure,
            "Descriptive leverage event study",
            f"Negative events: {study.negative.shape[0]:,}; positive events: {study.positive.shape[0]:,}. "
            "Bands are event-level IQRs and may contain overlapping windows.",
            height=850,
        )

    def drawdown_volatility(self) -> go.Figure:
        """Summarise how the volatility proxy changes across drawdown regimes."""
        drawdowns, volatility = self._drawdown_observations()
        edges = np.array([0.0, 0.02, 0.05, 0.10, 0.20, 0.30, np.inf])
        labels = ["0–2%", "2–5%", "5–10%", "10–20%", "20–30%", ">30%"]
        median, lower, upper, counts = [], [], [], []
        depths = -drawdowns
        for left, right in zip(edges[:-1], edges[1:]):
            sample = volatility[(depths >= left) & (depths < right)]
            median.append(float(np.median(sample)))
            lower.append(float(np.quantile(sample, 0.10)))
            upper.append(float(np.quantile(sample, 0.90)))
            counts.append(int(sample.size))
        figure = go.Figure()
        figure.add_trace(go.Scatter(x=labels, y=upper, mode="lines", line={"width": 0},
                                    showlegend=False, hoverinfo="skip"))
        figure.add_trace(go.Scatter(x=labels, y=lower, mode="lines", fill="tonexty",
                                    fillcolor="rgba(213,94,0,0.18)", line={"width": 0},
                                    name="10th–90th percentile"))
        figure.add_trace(go.Scatter(x=labels, y=median, mode="lines+markers", name="Median volatility",
                                    line={"color": VOL_COLOR, "width": 3},
                                    marker={"size": 8}, customdata=np.asarray(counts)[:, None],
                                    hovertemplate="Drawdown %{x}<br>Median vol %{y:.1f}%<br>n=%{customdata[0]:,}<extra></extra>"))
        figure.update_xaxes(title_text="Current drawdown from running peak")
        figure.update_yaxes(title_text="Daily annualized volatility proxy (%)", rangemode="tozero")
        return self._style(
            figure,
            "Volatility across asset drawdown regimes",
            f"Pooled over all {self._pooled_matrices()[0].shape[0]} paths; bands show within-regime "
            "dispersion, not confidence intervals.",
            height=650,
        )

    def _display_series(self, item: SeedSimulation) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Return years, a 100-based asset index, and annualized volatility percent."""
        returns = item.result.log_returns[self.display_path]
        variance = item.result.realized_variance[self.display_path]
        years = np.arange(returns.size, dtype=float) / self.observations_per_year
        asset = self._asset_index(returns)
        volatility = 100.0 * np.sqrt(variance)
        return years, asset, volatility

    def _spike_series(self, item: SeedSimulation) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Extract a fixed-length window around the displayed path's volatility maximum."""
        _, asset, volatility = self._display_series(item)
        first_eligible = 126
        last_eligible = volatility.size - 252
        peak = first_eligible + int(np.argmax(volatility[first_eligible:last_eligible]))
        start = peak - 126
        stop = peak + 253
        days = np.arange(start, stop, dtype=float) - peak
        reindexed = 100.0 * asset[start:stop] / asset[start]
        return days, reindexed, volatility[start:stop]

    def _pooled_matrices(self) -> tuple[FloatArray, FloatArray]:
        """Stack all paths while retaining path boundaries for diagnostics."""
        returns = np.vstack([item.result.log_returns for item in self.simulations])
        variance = np.vstack([item.result.realized_variance for item in self.simulations])
        return returns, variance

    def _standardized_returns(self) -> FloatArray:
        """Standardize within paths before pooling to avoid scale-mixture artefacts."""
        returns, _ = self._pooled_matrices()
        centered = returns - np.mean(returns, axis=1, keepdims=True)
        scales = np.std(centered, axis=1, ddof=1, keepdims=True)
        return (centered / np.maximum(scales, 1e-15)).ravel()

    def _event_study(self) -> EventStudy:
        """Collect normalized volatility windows around within-path tail events."""
        returns, variance = self._pooled_matrices()
        volatility = 100.0 * np.sqrt(variance)
        lags = np.arange(-10, 41, dtype=int)
        negative: list[FloatArray] = []
        positive: list[FloatArray] = []
        for path_returns, path_volatility in zip(returns, volatility):
            low, high = np.quantile(path_returns, [0.025, 0.975])
            valid = np.arange(10, path_returns.size - 40)
            self._collect_events(valid[path_returns[valid] <= low], path_volatility, lags, negative)
            self._collect_events(valid[path_returns[valid] >= high], path_volatility, lags, positive)
        return EventStudy(lags, np.vstack(negative), np.vstack(positive))

    @staticmethod
    def _collect_events(
        events: Iterable[int], volatility: FloatArray, lags: NDArray[np.int64], target: list[FloatArray]
    ) -> None:
        """Normalize each event window by its preceding ten-day mean."""
        for event in events:
            baseline = float(np.mean(volatility[event - 10 : event]))
            response = 100.0 * (volatility[event + lags] / max(baseline, 1e-15) - 1.0)
            target.append(response)

    def _drawdown_observations(self) -> tuple[FloatArray, FloatArray]:
        """Pool pointwise drawdown and volatility observations across all paths."""
        returns, variance = self._pooled_matrices()
        drawdown_rows: list[FloatArray] = []
        for path_returns in returns:
            asset = self._asset_index(path_returns)
            drawdown_rows.append(asset / np.maximum.accumulate(asset) - 1.0)
        return np.concatenate(drawdown_rows), 100.0 * np.sqrt(variance).ravel()

    @staticmethod
    def _asset_index(returns: FloatArray) -> FloatArray:
        """Exponentiate cumulative log returns and rebase the first plotted point to 100."""
        log_asset = np.cumsum(returns)
        log_asset -= log_asset[0]
        return 100.0 * np.exp(log_asset)

    @staticmethod
    def _acf_matrix(values: FloatArray, maximum_lag: int) -> FloatArray:
        """Compute pathwise sample autocorrelations for positive lags."""
        rows = np.empty((values.shape[0], maximum_lag), dtype=float)
        for index, series in enumerate(values):
            centered = series - np.mean(series)
            denominator = float(np.dot(centered, centered))
            for lag in range(1, maximum_lag + 1):
                rows[index, lag - 1] = np.dot(centered[:-lag], centered[lag:]) / denominator
        return rows

    @staticmethod
    def _quantile_curves(values: FloatArray) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Return the 25th, 50th, and 75th event-response percentiles."""
        lower, median, upper = np.quantile(values, [0.25, 0.50, 0.75], axis=0)
        return lower, median, upper

    @staticmethod
    def _line(x: FloatArray, y: FloatArray, name: str, color: str, legend: bool) -> go.Scatter:
        """Construct a consistent line trace for path panels."""
        return go.Scatter(x=x, y=y, mode="lines", name=name, showlegend=legend,
                          line={"color": color, "width": 1.5})

    @staticmethod
    def _add_acf_band(
        figure: go.Figure, values: FloatArray, name: str, color: str, row: int, col: int
    ) -> None:
        """Add a path-quantile ribbon and mean ACF line."""
        x = np.arange(1, values.shape[1] + 1)
        lower, upper = np.quantile(values, [0.10, 0.90], axis=0)
        rgba = PlotlyVisualDiagnostics._rgba(color, 0.12)
        figure.add_trace(go.Scatter(x=x, y=upper, mode="lines", line={"width": 0},
                                    showlegend=False, hoverinfo="skip"), row=row, col=col)
        figure.add_trace(go.Scatter(x=x, y=lower, mode="lines", line={"width": 0},
                                    fill="tonexty", fillcolor=rgba, showlegend=False,
                                    hoverinfo="skip"), row=row, col=col)
        figure.add_trace(go.Scatter(x=x, y=np.mean(values, axis=0), mode="lines", name=name,
                                    line={"color": color, "width": 2}), row=row, col=col)

    @staticmethod
    def _add_event_band(
        figure: go.Figure,
        lags: NDArray[np.int64],
        curves: tuple[FloatArray, FloatArray, FloatArray],
        name: str,
        color: str,
        row: int,
    ) -> None:
        """Add an interquartile response ribbon and median event line."""
        lower, median, upper = curves
        figure.add_trace(go.Scatter(x=lags, y=upper, mode="lines", line={"width": 0},
                                    showlegend=False, hoverinfo="skip"), row=row, col=1)
        figure.add_trace(go.Scatter(x=lags, y=lower, mode="lines", line={"width": 0},
                                    fill="tonexty", fillcolor=PlotlyVisualDiagnostics._rgba(color, 0.16),
                                    showlegend=False, hoverinfo="skip"), row=row, col=1)
        figure.add_trace(go.Scatter(x=lags, y=median, mode="lines", name=name,
                                    line={"color": color, "width": 2.5}), row=row, col=1)

    @staticmethod
    def _rgba(hex_color: str, alpha: float) -> str:
        """Convert a hexadecimal color to an RGBA string."""
        red, green, blue = (int(hex_color[index : index + 2], 16) for index in (1, 3, 5))
        return f"rgba({red},{green},{blue},{alpha})"

    @staticmethod
    def _style(figure: go.Figure, title: str, subtitle: str, height: int) -> go.Figure:
        """Apply one publication-oriented Plotly theme."""
        top_margin = 155
        plot_height = max(300, height - top_margin - 70)
        legend_y = 1.0 + 32.0 / plot_height
        figure.update_layout(
            template="plotly_white",
            title={"text": f"{title}<br><sup>{subtitle}</sup>", "x": 0.02, "xanchor": "left"},
            height=height,
            width=1500,
            margin={"l": 80, "r": 90, "t": top_margin, "b": 70},
            hovermode="x unified",
            legend={"orientation": "h", "yanchor": "bottom", "y": legend_y, "x": 0.5,
                    "xanchor": "center"},
            font={"family": "Arial, sans-serif", "size": 14, "color": "#25313B"},
        )
        figure.update_xaxes(showgrid=True, gridcolor=LIGHT_GRID, zeroline=False)
        figure.update_yaxes(showgrid=True, gridcolor=LIGHT_GRID, zeroline=False)
        return figure

    def _export_figures(self) -> dict[str, Path]:
        """Write every Plotly figure as a high-resolution PNG."""
        names = (
            "01_m2_paths_overview.png",
            "02_m2_maximum_volatility_windows.png",
            "03_m2_return_tails.png",
            "04_m2_dependence.png",
            "05_m2_leverage_event_study.png",
            "06_m2_drawdown_volatility.png",
        )
        outputs: dict[str, Path] = {}
        for (label, figure), filename in zip(self.figures, names):
            path = self.output_dir / filename
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                figure.write_image(path, width=figure.layout.width, height=figure.layout.height, scale=1.25)
            outputs[label] = path
        return outputs

    def _write_dashboard(self) -> Path:
        """Write a self-contained interactive HTML dashboard with all figures."""
        cards = self._dashboard_cards()
        sections: list[str] = []
        for index, (label, figure) in enumerate(self.figures):
            fragment = pio.to_html(
                figure,
                full_html=False,
                include_plotlyjs=True if index == 0 else False,
                config={"responsive": True, "displaylogo": False},
            )
            sections.append(f"<section><h2>{label}</h2>{fragment}</section>")
        html = self._dashboard_shell(cards, "\n".join(sections))
        path = self.output_dir / "m2_visual_diagnostics_dashboard.html"
        path.write_text(html, encoding="utf-8")
        return path

    def _dashboard_cards(self) -> str:
        """Create compact headline cards from the pooled numerical summary."""
        summary = self._summary_payload()["pooled"]
        values = (
            ("Paths", f"{sum(item.result.log_returns.shape[0] for item in self.simulations)}"),
            ("Daily observations", f"{summary['observations']:,}"),
            ("Excess kurtosis", f"{summary['excess_kurtosis']:.2f}"),
            ("Median volatility", f"{summary['median_volatility_pct']:.1f}%"),
            ("99% volatility", f"{summary['q99_volatility_pct']:.1f}%"),
        )
        return "".join(f"<div class='card'><b>{value}</b><span>{label}</span></div>" for label, value in values)

    @staticmethod
    def _dashboard_shell(cards: str, sections: str) -> str:
        """Wrap Plotly fragments in a portable, styled HTML document."""
        return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>M2 visual diagnostics</title>
<style>
body{{margin:0;background:#f4f6f8;color:#25313b;font-family:Arial,sans-serif}}
main{{max-width:1560px;margin:auto;padding:28px}} h1{{color:{ASSET_COLOR};margin-bottom:8px}}
.note{{background:#fff4e8;border-left:5px solid {VOL_COLOR};padding:14px 18px;line-height:1.5}}
.cards{{display:flex;gap:12px;flex-wrap:wrap;margin:20px 0}} .card{{background:white;padding:14px 20px;
border:1px solid #dfe4e8;border-radius:8px;min-width:150px}} .card b{{display:block;font-size:24px;color:{ASSET_COLOR}}}
.card span{{color:#5b6670}} section{{background:white;margin:22px 0;padding:16px;border-radius:8px;
box-shadow:0 1px 4px rgba(23,50,77,.08)}} h2{{margin:5px 12px;color:{ASSET_COLOR}}}
</style></head><body><main><h1>Stage-3 M2 visual diagnostics</h1>
<p class="note"><b>Interpretation boundary.</b> These are synthetic physical-measure paths from the selected
M2 capacity experiment. The orange series is 100 × √(daily average instantaneous variance): it is VIX-like
only in appearance and is <b>not</b> the 30-day risk-neutral conditional expectation defining VIX. No market
data were used, so resemblance is qualitative rather than validation.</p>
<div class="cards">{cards}</div>{sections}</main></body></html>"""

    def _write_display_data(self) -> Path:
        """Export the exact five displayed paths for independent reproduction."""
        frames: list[pd.DataFrame] = []
        for item in self.simulations:
            years, asset, volatility = self._display_series(item)
            returns = item.result.log_returns[self.display_path]
            frames.append(pd.DataFrame({
                "seed": item.seed,
                "path_index": self.display_path,
                "observation": np.arange(returns.size),
                "simulated_year": years,
                "asset_index": asset,
                "log_return_pct": 100.0 * returns,
                "annualized_volatility_proxy_pct": volatility,
            }))
        path = self.output_dir / "m2_display_paths.csv"
        pd.concat(frames, ignore_index=True).to_csv(path, index=False)
        return path

    def _write_summary(self) -> Path:
        """Export pooled and displayed-path statistics used in the assessment."""
        path = self.output_dir / "m2_visual_summary.json"
        path.write_text(json.dumps(self._summary_payload(), indent=2), encoding="utf-8")
        return path

    def _summary_payload(self) -> dict[str, Any]:
        """Compute audit statistics without treating them as empirical estimates."""
        returns, variance = self._pooled_matrices()
        standardized = self._standardized_returns()
        volatility = 100.0 * np.sqrt(variance)
        return {
            "method": {
                "measure": "P",
                "seeds": [item.seed for item in self.simulations],
                "paths_per_seed": int(returns.shape[0] / len(self.simulations)),
                "display_path_index": self.display_path,
                "volatility_definition": "100 * sqrt(daily average instantaneous variance)",
                "is_formal_vix": False,
            },
            "pooled": {
                "observations": int(returns.size),
                "mean_daily_log_return_pct": 100.0 * float(np.mean(returns)),
                "daily_return_skewness": float(skew(standardized, bias=False)),
                "excess_kurtosis": float(kurtosis(standardized, fisher=True, bias=False)),
                "absolute_return_q99_in_path_sd": float(np.quantile(np.abs(standardized), 0.99)),
                "mean_volatility_pct": float(np.mean(volatility)),
                "median_volatility_pct": float(np.median(volatility)),
                "root_mean_variance_pct": 100.0 * float(np.sqrt(np.mean(variance))),
                "q95_volatility_pct": float(np.quantile(volatility, 0.95)),
                "q99_volatility_pct": float(np.quantile(volatility, 0.99)),
                "maximum_volatility_pct": float(np.max(volatility)),
                "fraction_volatility_above_40pct": float(np.mean(volatility > 40.0)),
                "fraction_volatility_above_60pct": float(np.mean(volatility > 60.0)),
            },
            "displayed_paths": [self._display_path_summary(item) for item in self.simulations],
            "event_study": self._event_summary(),
            "dependence": self._dependence_summary(),
            "drawdown_regimes": self._drawdown_regime_summary(),
        }

    def _display_path_summary(self, item: SeedSimulation) -> dict[str, Any]:
        """Summarise the fixed path shown for one seed."""
        _, asset, volatility = self._display_series(item)
        drawdown = asset / np.maximum.accumulate(asset) - 1.0
        years = (asset.size - 1) / self.observations_per_year
        annualized_growth = (asset[-1] / asset[0]) ** (1.0 / years) - 1.0
        return {
            "seed": item.seed,
            "path_index": self.display_path,
            "terminal_asset_index": float(asset[-1]),
            "annualized_compound_growth_pct": 100.0 * float(annualized_growth),
            "maximum_drawdown_pct": 100.0 * float(np.min(drawdown)),
            "median_volatility_pct": float(np.median(volatility)),
            "q99_volatility_pct": float(np.quantile(volatility, 0.99)),
            "maximum_volatility_pct": float(np.max(volatility)),
        }

    def _event_summary(self) -> dict[str, Any]:
        """Report the robust asymmetry visible in the event-study figure."""
        study = self._event_study()
        negative = np.median(study.negative, axis=0)
        positive = np.median(study.positive, axis=0)
        future = (study.lags >= 1) & (study.lags <= 5)
        first_six = (study.lags >= 0) & (study.lags <= 5)
        return {
            "negative_event_count": int(study.negative.shape[0]),
            "positive_event_count": int(study.positive.shape[0]),
            "day0_negative_response_pct": float(negative[study.lags == 0][0]),
            "day0_positive_response_pct": float(positive[study.lags == 0][0]),
            "days1_to_5_negative_minus_positive_median_pp": float(
                np.mean(negative[future] - positive[future])
            ),
            "days0_to_5_negative_median_response_pct": negative[first_six].tolist(),
            "days0_to_5_positive_median_response_pct": positive[first_six].tolist(),
        }

    def _drawdown_regime_summary(self) -> list[dict[str, Any]]:
        """Return the numerical values behind the drawdown-regime plot."""
        drawdowns, volatility = self._drawdown_observations()
        edges = np.array([0.0, 0.02, 0.05, 0.10, 0.20, 0.30, np.inf])
        labels = ["0-2%", "2-5%", "5-10%", "10-20%", "20-30%", ">30%"]
        depths = -drawdowns
        rows: list[dict[str, Any]] = []
        for label, left, right in zip(labels, edges[:-1], edges[1:]):
            sample = volatility[(depths >= left) & (depths < right)]
            rows.append({
                "drawdown": label,
                "observations": int(sample.size),
                "median_volatility_pct": float(np.median(sample)),
                "q10_volatility_pct": float(np.quantile(sample, 0.10)),
                "q90_volatility_pct": float(np.quantile(sample, 0.90)),
            })
        return rows

    def _dependence_summary(self) -> dict[str, Any]:
        """Report mean pathwise ACF values corresponding to the visual curves."""
        returns, variance = self._pooled_matrices()
        return_acf = self._acf_matrix(returns, 30)
        log_variance_acf = self._acf_matrix(np.log(np.maximum(variance, 1e-15)), 60)
        absolute_acf = self._acf_matrix(np.abs(returns), 60)
        squared_acf = self._acf_matrix(returns**2, 60)
        return {
            "mean_max_absolute_return_acf_lags_1_30": float(
                np.mean(np.max(np.abs(return_acf), axis=1))
            ),
            "mean_log_variance_acf_lag_1": float(np.mean(log_variance_acf[:, 0])),
            "mean_log_variance_acf_lag_20": float(np.mean(log_variance_acf[:, 19])),
            "mean_absolute_return_acf_lags_1_20": float(np.mean(absolute_acf[:, :20])),
            "mean_squared_return_acf_lags_1_20": float(np.mean(squared_acf[:, :20])),
        }


def _bundle_manifest(configuration: SelectedM2Configuration, outputs: dict[str, Path]) -> dict[str, Any]:
    """Create a compact provenance record for the visual suite."""
    return {
        "model": "M2",
        "stage3_candidate_id": configuration.candidate_id,
        "stage3_structural_scenario_id": configuration.scenario_id,
        "seeds": list(configuration.seeds),
        "display_path_index": 0,
        "measure": "P",
        "volatility_is_formal_vix": False,
        "outputs": {name: path.name for name, path in outputs.items()},
    }


def main() -> None:
    """Load, simulate, visualize, and write an auditable output manifest."""
    model_root = Path(__file__).resolve().parents[1]
    experiment_dir = model_root / "experiments" / "stage3_001"
    output_dir = model_root / "visual_diagnostics"
    configuration = Stage3M2Loader(experiment_dir).load(candidate_id=18, scenario_id=0)
    simulations = M2SimulationSuite(configuration).run()
    diagnostics = PlotlyVisualDiagnostics(simulations, output_dir, display_path=0)
    outputs = diagnostics.generate()
    manifest = _bundle_manifest(configuration, outputs)
    path = output_dir / "visual_diagnostics_manifest.json"
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"manifest": str(path), **{key: str(value) for key, value in outputs.items()}}, indent=2))


if __name__ == "__main__":
    main()
