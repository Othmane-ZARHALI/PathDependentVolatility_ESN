"""Deterministic calibration of latent roughness and clustering features."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import differential_evolution

from .config import ArchitectureConfig, FloatArray, geometric_rates


@dataclass(frozen=True, slots=True)
class LatentTargets:
    """Targets that can be assessed before choosing nonlinear readout parameters."""

    rough_hurst: float = 0.10
    rough_lags_days: tuple[float, ...] = (1 / 26, 1 / 13, 1 / 6.5, 1 / 3.25, 1, 2, 5, 10)
    cluster_lags_days: tuple[float, ...] = (1, 5, 20, 60, 120, 252)
    cluster_acf: tuple[float, ...] = (0.98, 0.92, 0.78, 0.58, 0.42, 0.25)
    trading_days: int = 252

    def __post_init__(self) -> None:
        """Validate target curves and the roughness domain."""
        if not 0.0 < self.rough_hurst < 0.5:
            raise ValueError("rough_hurst must lie in (0, 0.5).")
        if len(self.cluster_lags_days) != len(self.cluster_acf):
            raise ValueError("Clustering lags and ACF targets must have equal length.")
        if any(value <= 0.0 for value in (*self.rough_lags_days, *self.cluster_lags_days)):
            raise ValueError("All calibration lags must be positive.")


@dataclass(frozen=True, slots=True)
class ArchitectureSearchBounds:
    """Search box for rate endpoints and fixed projection powers."""

    feedback_low: tuple[float, float] = (0.05, 5.0)
    feedback_high: tuple[float, float] = (100.0, 5000.0)
    feedback_power: tuple[float, float] = (-0.5, 1.0)
    cluster_low: tuple[float, float] = (0.005, 0.5)
    cluster_high: tuple[float, float] = (2.0, 100.0)
    cluster_power: tuple[float, float] = (-2.0, 1.0)


@dataclass(frozen=True, slots=True)
class LatentCalibrationResult:
    """Fitted architecture and transparent latent-curve residuals."""

    architecture: ArchitectureConfig
    objective: float
    rough_variogram: FloatArray
    cluster_acf: FloatArray
    optimizer_message: str


@dataclass(frozen=True, slots=True)
class LatentCalibrationEnsemble:
    """Replicated stochastic-optimizer fits and their best objective."""

    seeds: tuple[int, ...]
    results: tuple[LatentCalibrationResult, ...]

    def __post_init__(self) -> None:
        """Require one labelled fit for every distinct optimizer seed."""
        if not self.seeds or len(self.seeds) != len(self.results):
            raise ValueError("Seeds and calibration results must be non-empty and aligned.")
        if len(set(self.seeds)) != len(self.seeds):
            raise ValueError("Optimizer seeds must be distinct.")

    @property
    def best_index(self) -> int:
        """Return the index of the smallest replicated objective."""
        return int(np.argmin([item.objective for item in self.results]))

    @property
    def best(self) -> LatentCalibrationResult:
        """Return the selected fit without discarding the alternatives."""
        return self.results[self.best_index]

    @property
    def best_seed(self) -> int:
        """Return the optimizer seed attached to the selected fit."""
        return self.seeds[self.best_index]


def projection_covariance(
    rates: FloatArray,
    weights: FloatArray,
    lags: FloatArray,
    shared_driver: bool,
) -> FloatArray:
    """Return scalar lag covariance for a stationary OU projection."""
    driver = np.ones((rates.size, rates.size)) if shared_driver else np.eye(rates.size)
    stationary = driver * (
        2.0 * np.sqrt(np.outer(rates, rates)) / (rates[:, None] + rates[None, :])
    )
    covariance = []
    for lag in lags:
        lagged = stationary * np.exp(-rates[None, :] * lag)
        covariance.append(float(weights @ lagged @ weights))
    return np.asarray(covariance)


def power_weights(rates: FloatArray, power: float, shared_driver: bool) -> FloatArray:
    """Normalise rate-power weights to variance one."""
    raw = rates**power
    variance = projection_covariance(rates, raw, np.asarray([0.0]), shared_driver)[0]
    return raw / np.sqrt(variance)


class LatentArchitectureCalibrator:
    """Fit only the latent feature shapes that precede the nonlinear readout."""

    def __init__(self, bounds: ArchitectureSearchBounds | None = None) -> None:
        """Store explicit broad bounds to avoid hidden optimisation choices."""
        self.bounds = bounds or ArchitectureSearchBounds()

    def fit(
        self,
        base: ArchitectureConfig,
        targets: LatentTargets,
        seed: int = 2026,
        max_iterations: int = 60,
    ) -> LatentCalibrationResult:
        """Calibrate fixed dimensions, rate endpoints, and projection powers."""
        bounds = self._optimizer_bounds()
        result = differential_evolution(
            lambda point: self._objective(point, base, targets),
            bounds,
            seed=seed,
            maxiter=max_iterations,
            polish=True,
            updating="immediate",
        )
        architecture = self._decode(result.x, base)
        rough, cluster = self.curves(architecture, targets)
        return LatentCalibrationResult(
            architecture, float(result.fun), rough, cluster, str(result.message)
        )

    def fit_replicated(
        self,
        base: ArchitectureConfig,
        targets: LatentTargets,
        seeds: tuple[int, ...],
        max_iterations: int = 60,
    ) -> LatentCalibrationEnsemble:
        """Repeat the stochastic optimiser and select only after retaining all fits."""
        if not seeds or len(set(seeds)) != len(seeds):
            raise ValueError("Provide one or more distinct optimizer seeds.")
        results = tuple(
            self.fit(base, targets, seed=seed, max_iterations=max_iterations)
            for seed in seeds
        )
        return LatentCalibrationEnsemble(seeds, results)

    def curves(
        self,
        architecture: ArchitectureConfig,
        targets: LatentTargets,
    ) -> tuple[FloatArray, FloatArray]:
        """Return normalised rough variogram and clustering ACF curves."""
        rough_lags = np.asarray(targets.rough_lags_days) / targets.trading_days
        feedback_rates = np.asarray(architecture.feedback.rates)
        feedback_weights = power_weights(
            feedback_rates, architecture.feedback_weight_power, True
        )
        covariance = projection_covariance(feedback_rates, feedback_weights, rough_lags, True)
        rough = np.maximum(2.0 * (1.0 - covariance), 1e-15)
        cluster_lags = np.asarray(targets.cluster_lags_days) / targets.trading_days
        cluster_rates = np.asarray(architecture.cluster.rates)
        cluster_weights = power_weights(
            cluster_rates, architecture.cluster_weight_power, False
        )
        cluster = projection_covariance(cluster_rates, cluster_weights, cluster_lags, False)
        return rough, cluster

    def _objective(
        self,
        point: FloatArray,
        base: ArchitectureConfig,
        targets: LatentTargets,
    ) -> float:
        """Compare curve shapes without pretending to fit full-model facts."""
        architecture = self._decode(point, base)
        rough, cluster = self.curves(architecture, targets)
        rough_lags = np.asarray(targets.rough_lags_days, dtype=float)
        rough_target = rough_lags ** (2.0 * targets.rough_hurst)
        log_error = np.log(rough) - np.log(rough_target)
        log_error -= np.mean(log_error)
        cluster_error = cluster - np.asarray(targets.cluster_acf)
        return float(np.mean(log_error**2) + 2.0 * np.mean(cluster_error**2))

    def _decode(self, point: FloatArray, base: ArchitectureConfig) -> ArchitectureConfig:
        """Decode log-rate endpoints and powers into a valid architecture."""
        feedback_low, feedback_high = np.exp(point[0]), np.exp(point[1])
        cluster_low, cluster_high = np.exp(point[3]), np.exp(point[4])
        feedback = replace(
            base.feedback,
            rates=geometric_rates(feedback_low, feedback_high, len(base.feedback.rates)),
        )
        cluster = replace(
            base.cluster,
            rates=geometric_rates(cluster_low, cluster_high, len(base.cluster.rates)),
        )
        return replace(
            base,
            feedback=feedback,
            cluster=cluster,
            feedback_weight_power=float(point[2]),
            cluster_weight_power=float(point[5]),
        )

    def _optimizer_bounds(self) -> list[tuple[float, float]]:
        """Transform only positive rate endpoints to log coordinates."""
        return [
            tuple(np.log(self.bounds.feedback_low)),
            tuple(np.log(self.bounds.feedback_high)),
            self.bounds.feedback_power,
            tuple(np.log(self.bounds.cluster_low)),
            tuple(np.log(self.bounds.cluster_high)),
            self.bounds.cluster_power,
        ]
