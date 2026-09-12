"""Matched-observation diagnostics for the targeted physical facts."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import kurtosis, rankdata

from .config import FloatArray
from .model import SimulationResult


@dataclass(frozen=True, slots=True)
class DiagnosticConfig:
    """Lag conventions expressed in units of the output observation frequency."""

    maximum_dependence_lag: int = 20
    roughness_maximum_lag: int = 40
    zumbach_windows: tuple[int, ...] = (5, 10, 20)

    def __post_init__(self) -> None:
        """Require strictly positive diagnostic horizons."""
        values = (self.maximum_dependence_lag, self.roughness_maximum_lag, *self.zumbach_windows)
        if any(value < 1 for value in values):
            raise ValueError("Diagnostic lags and windows must be positive.")


def _correlation(left: FloatArray, right: FloatArray) -> float:
    """Return NaN when a sample is too short or effectively constant."""
    if left.size < 3 or np.std(left) < 1e-14 or np.std(right) < 1e-14:
        return float("nan")
    return float(np.corrcoef(left, right)[0, 1])


def _finite_mean(values: list[float]) -> float:
    """Average only well-defined diagnostics."""
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    return float(np.mean(finite)) if finite.size else float("nan")


def rough_hurst(variance: FloatArray, maximum_lag: int = 40) -> float:
    """Estimate finite-band H from the log-variance structure function."""
    log_variance = np.log(np.maximum(variance, 1e-15))
    largest = min(maximum_lag, max(2, log_variance.size // 4))
    lags = np.unique(np.round(np.geomspace(1, largest, 12)).astype(int))
    x_values, y_values = [], []
    for lag in lags:
        increment = log_variance[lag:] - log_variance[:-lag]
        structure = float(np.mean(increment**2))
        if structure > 0.0:
            x_values.append(np.log(float(lag)))
            y_values.append(np.log(structure))
    if len(x_values) < 3:
        return float("nan")
    return 0.5 * float(np.polyfit(x_values, y_values, 1)[0])


def hill_tail_summary(returns: FloatArray) -> tuple[float, float]:
    """Average several Hill estimates and report threshold instability."""
    ordered = np.sort(np.abs(returns))[::-1]
    estimates: list[float] = []
    for fraction in (0.02, 0.05, 0.10):
        order = min(max(20, int(fraction * ordered.size)), ordered.size - 2)
        if order < 2 or ordered[order] <= 0.0:
            continue
        denominator = float(np.mean(np.log(ordered[:order] / ordered[order])))
        if denominator > 0.0:
            estimates.append(1.0 / denominator)
    return _finite_mean(estimates), float(np.std(estimates)) if estimates else float("nan")


def zumbach_statistic(
    returns: FloatArray,
    variance: FloatArray,
    rank_based: bool = False,
    windows: tuple[int, ...] = (5, 10, 20),
) -> float:
    """Measure causal past-return/future-RV asymmetry on matched windows."""
    values: list[float] = []
    cumulative_r = np.concatenate(([0.0], np.cumsum(returns)))
    cumulative_v = np.concatenate(([0.0], np.cumsum(variance)))
    for window in windows:
        count = returns.size - 2 * window
        if count < 30:
            continue
        center = window + np.arange(count)
        past_r = cumulative_r[center] - cumulative_r[center - window]
        future_r = cumulative_r[center + window] - cumulative_r[center]
        past_v = (cumulative_v[center] - cumulative_v[center - window]) / window
        future_v = (cumulative_v[center + window] - cumulative_v[center]) / window
        transform = rankdata if rank_based else np.asarray
        forward = _correlation(transform(past_r**2), transform(future_v))
        backward = _correlation(transform(past_v), transform(future_r**2))
        values.append(forward - backward)
    return _finite_mean(values)


def lagged_effects(
    returns: FloatArray,
    variance: FloatArray,
    maximum_lag: int = 20,
) -> dict[str, float]:
    """Compute return ACF, leverage, clustering, and Taylor diagnostics."""
    lags = range(1, min(maximum_lag, returns.size - 2) + 1)
    return_acf, leverage, leverage_rank = [], [], []
    vol_acf, taylor = [], []
    for lag in lags:
        return_acf.append(_correlation(returns[:-lag], returns[lag:]))
        leverage.append(_correlation(returns[:-lag], variance[lag:]))
        leverage_rank.append(_correlation(rankdata(returns[:-lag]), rankdata(variance[lag:])))
        vol_acf.append(_correlation(np.log(variance[:-lag]), np.log(variance[lag:])))
        absolute = _correlation(np.abs(returns[:-lag]), np.abs(returns[lag:]))
        squared = _correlation(returns[:-lag] ** 2, returns[lag:] ** 2)
        taylor.append(absolute - squared)
    return {
        "max_abs_return_acf": float(np.nanmax(np.abs(return_acf))),
        "leverage": _finite_mean(leverage),
        "leverage_rank": _finite_mean(leverage_rank),
        "mean_log_variance_acf": _finite_mean(vol_acf),
        "taylor_gap": _finite_mean(taylor),
        "taylor_fraction": float(np.mean(np.asarray(taylor) > 0.0)),
    }


@dataclass(frozen=True, slots=True)
class DiagnosticSummary:
    """Across-path diagnostic means, uncertainty, and raw path estimates."""

    mean: dict[str, float]
    standard_error: dict[str, float]
    per_path: tuple[dict[str, float], ...]


class StylisedFactDiagnostics:
    """Evaluate every requested P-measure mechanism with matched data."""

    def __init__(self, config: DiagnosticConfig | None = None) -> None:
        """Make all frequency-dependent lag conventions explicit."""
        self.config = config or DiagnosticConfig()

    def evaluate(self, result: SimulationResult) -> DiagnosticSummary:
        """Compute pathwise estimates before aggregating simulation uncertainty."""
        per_path = tuple(
            self._single_path(result.log_returns[index], result.realized_variance[index])
            for index in range(result.log_returns.shape[0])
        )
        names = tuple(per_path[0])
        means = {name: _finite_mean([row[name] for row in per_path]) for name in names}
        errors = {
            name: self._standard_error([row[name] for row in per_path]) for name in names
        }
        return DiagnosticSummary(means, errors, per_path)

    def _single_path(self, returns: FloatArray, variance: FloatArray) -> dict[str, float]:
        """Compute one internally consistent set of stylised-fact estimates."""
        tail_index, tail_instability = hill_tail_summary(returns)
        scale = float(np.std(returns))
        result = {
            "mean_volatility": float(np.mean(np.sqrt(variance))),
            "q995_volatility": float(np.quantile(np.sqrt(variance), 0.995)),
            "excess_kurtosis": float(kurtosis(returns, fisher=True, bias=False)),
            "absolute_return_tail_ratio": float(
                np.quantile(np.abs(returns), 0.99) / max(scale, 1e-15)
            ),
            "hill_tail_index": tail_index,
            "hill_threshold_instability": tail_instability,
            "rough_hurst": rough_hurst(variance, self.config.roughness_maximum_lag),
            "zumbach": zumbach_statistic(
                returns, variance, windows=self.config.zumbach_windows
            ),
            "zumbach_rank": zumbach_statistic(
                returns,
                variance,
                rank_based=True,
                windows=self.config.zumbach_windows,
            ),
        }
        result.update(lagged_effects(returns, variance, self.config.maximum_dependence_lag))
        return result

    @staticmethod
    def _standard_error(values: list[float]) -> float:
        """Return across-path Monte Carlo standard error."""
        array = np.asarray(values, dtype=float)
        array = array[np.isfinite(array)]
        if array.size < 2:
            return float("nan")
        return float(np.std(array, ddof=1) / np.sqrt(array.size))


def aggregate_result(result: SimulationResult, factor: int) -> SimulationResult:
    """Aggregate matched returns and variance to a coarser frequency."""
    if factor < 1 or result.observations_per_year % factor:
        raise ValueError("factor must divide observations_per_year.")
    complete = (result.log_returns.shape[1] // factor) * factor
    if complete < factor:
        raise ValueError("The result is too short for the requested aggregation.")
    path_count = result.log_returns.shape[0]
    returns = result.log_returns[:, :complete].reshape(path_count, -1, factor).sum(axis=2)
    variance = result.realized_variance[:, :complete].reshape(path_count, -1, factor).mean(axis=2)
    return SimulationResult(
        returns,
        variance,
        result.observations_per_year // factor,
        result.measure,
        result.metadata,
    )
