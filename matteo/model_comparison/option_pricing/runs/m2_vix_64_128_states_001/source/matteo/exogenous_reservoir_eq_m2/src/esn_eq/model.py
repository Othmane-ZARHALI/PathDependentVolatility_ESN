"""P/Q simulation on prepared exogenous feature paths."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .config import (
    ArchitectureConfig,
    FloatArray,
    MarketEnvironment,
    ReadoutParameters,
    RiskPremium,
)
from .features import FeatureEngine, FeaturePaths
from .readout import TemperedExponentialQuadratic
from .states import ReservoirCache, q_state_paths

Measure = Literal["P", "Q"]


@dataclass(frozen=True, slots=True)
class PreparedPathCache:
    """Expensive state features shared by every cheap readout evaluation."""

    reservoir: ReservoirCache
    physical_features: FeaturePaths
    pricing_features: FeaturePaths
    risk_premium: RiskPremium


@dataclass(frozen=True, slots=True)
class SimulationMetadata:
    """Numerical checks needed to audit normalisation and tempering."""

    q_forward_mean_max_error: float
    p_cap_exceedance_fraction: float
    q_cap_exceedance_fraction: float
    normalizer_min: float
    normalizer_max: float


@dataclass(frozen=True, slots=True)
class SimulationResult:
    """Matched output returns and average instantaneous variance."""

    log_returns: FloatArray
    realized_variance: FloatArray
    observations_per_year: int
    measure: Measure
    metadata: SimulationMetadata


class ExogenousReservoirVolatilityModel:
    """Current tempered exponential-quadratic descendant of ESN-A2."""

    def __init__(self, architecture: ArchitectureConfig) -> None:
        """Fix the architecture while leaving all readout parameters free."""
        self.architecture = architecture
        self.readout = TemperedExponentialQuadratic()

    def prepare(self, cache: ReservoirCache, premium: RiskPremium) -> PreparedPathCache:
        """Compute P/Q features once; no readout coefficient enters this method."""
        self._validate_cache(cache)
        engine = FeatureEngine(self.architecture, cache.layout)
        physical = engine.build(cache.states[:-1])
        pricing = engine.build(q_state_paths(cache, premium)[:-1])
        return PreparedPathCache(cache, physical, pricing, premium)

    def simulate(
        self,
        prepared: PreparedPathCache,
        parameters: ReadoutParameters,
        market: MarketEnvironment,
        measure: Measure = "P",
    ) -> SimulationResult:
        """Recompute only score, variance, returns, and aggregation."""
        if measure not in ("P", "Q"):
            raise ValueError("measure must be 'P' or 'Q'.")
        q_score = self.readout.score(prepared.pricing_features, parameters)
        q_multiplier = self.readout.multiplier(prepared.pricing_features, parameters)
        normalizer = self.readout.normalizer(q_multiplier)
        p_score = self.readout.score(prepared.physical_features, parameters)
        features = prepared.physical_features if measure == "P" else prepared.pricing_features
        multiplier = self.readout.multiplier(features, parameters)
        curve = market.curve(prepared.reservoir.simulation.steps)
        variance = self._variance(multiplier, normalizer, curve, market.variance_floor)
        q_variance = self._variance(q_multiplier, normalizer, curve, market.variance_floor)
        returns, realized = self._observations(prepared, variance, market, measure)
        metadata = self._metadata(p_score, q_score, q_variance, curve, normalizer, parameters)
        kept = slice(prepared.reservoir.simulation.burn_observations, None)
        return SimulationResult(
            returns[:, kept],
            realized[:, kept],
            prepared.reservoir.simulation.observations_per_year,
            measure,
            metadata,
        )

    @staticmethod
    def _variance(
        multiplier: FloatArray,
        normalizer: FloatArray,
        curve: FloatArray,
        floor: float,
    ) -> FloatArray:
        """Apply the forward-variance-normalised positive map."""
        return floor + (curve[:, None] - floor) * multiplier / normalizer[:, None]

    @staticmethod
    def _observations(
        prepared: PreparedPathCache,
        variance: FloatArray,
        market: MarketEnvironment,
        measure: Measure,
    ) -> tuple[FloatArray, FloatArray]:
        """Integrate adapted fine-step returns and aggregate matched variance."""
        cache = prepared.reservoir
        premium_drift = prepared.risk_premium.equity * np.sqrt(variance) if measure == "P" else 0.0
        drift = market.risk_free_rate - market.dividend_yield + premium_drift - 0.5 * variance
        fine_returns = drift * cache.simulation.dt + np.sqrt(variance) * cache.equity_increments
        shape = (
            cache.simulation.observations,
            cache.simulation.steps_per_observation,
            cache.simulation.paths,
        )
        returns = fine_returns.reshape(shape).sum(axis=1).T
        realized = variance.reshape(shape).mean(axis=1).T
        return returns, realized

    @staticmethod
    def _metadata(
        p_score: FloatArray,
        q_score: FloatArray,
        q_variance: FloatArray,
        curve: FloatArray,
        normalizer: FloatArray,
        parameters: ReadoutParameters,
    ) -> SimulationMetadata:
        """Summarise checks without retaining all fine-step variance paths."""
        mean_error = float(np.max(np.abs(np.mean(q_variance, axis=1) - curve)))
        return SimulationMetadata(
            q_forward_mean_max_error=mean_error,
            p_cap_exceedance_fraction=float(np.mean(p_score > parameters.cap_level)),
            q_cap_exceedance_fraction=float(np.mean(q_score > parameters.cap_level)),
            normalizer_min=float(np.min(normalizer)),
            normalizer_max=float(np.max(normalizer)),
        )

    def _validate_cache(self, cache: ReservoirCache) -> None:
        """Prevent feature reuse across a different state architecture."""
        if cache.architecture != self.architecture:
            raise ValueError("Reservoir cache was built for another architecture.")

