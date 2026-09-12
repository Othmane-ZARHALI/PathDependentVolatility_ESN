"""Space-filling sampling of cheap readout parameters."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Literal, Mapping

import numpy as np
import pandas as pd
from scipy.stats import qmc

from .config import ArchitectureConfig, FloatArray, ReadoutParameters, RiskPremium


@dataclass(frozen=True, slots=True)
class ParameterBound:
    """One continuous box bound with linear or logarithmic scaling."""

    lower: float
    upper: float
    scale: Literal["linear", "log"] = "linear"

    def __post_init__(self) -> None:
        """Reject empty bounds and invalid logarithmic domains."""
        if not np.isfinite(self.lower) or not np.isfinite(self.upper) or self.upper <= self.lower:
            raise ValueError("Require finite lower < upper.")
        if self.scale == "log" and self.lower <= 0.0:
            raise ValueError("Logarithmic bounds must be strictly positive.")

    def transform(self, unit_values: FloatArray) -> FloatArray:
        """Map unit-cube coordinates into the requested physical scale."""
        if self.scale == "linear":
            return self.lower + unit_values * (self.upper - self.lower)
        low, high = np.log(self.lower), np.log(self.upper)
        return np.exp(low + unit_values * (high - low))


@dataclass(frozen=True, slots=True)
class ParameterSpace:
    """Named readout bounds, deliberately excluding state-changing choices."""

    bounds: Mapping[str, ParameterBound]
    fixed: Mapping[str, float]

    @classmethod
    def plausible_m2(cls) -> "ParameterSpace":
        """Return broad research bounds, not market-calibrated priors."""
        return cls(
            bounds={
                "cluster_loading": ParameterBound(0.05, 0.80),
                "feedback_curvature": ParameterBound(0.015, 0.24, "log"),
                "feedback_shift": ParameterBound(0.05, 1.00),
                "spike_curvature": ParameterBound(0.002, 0.25, "log"),
                "cap_level": ParameterBound(2.0, 5.5),
                "cap_sharpness": ParameterBound(2.0, 12.0),
            },
            fixed={"orthogonal_curvature": 0.0, "echo_loading": 0.0},
        )

    @classmethod
    def plausible_m3(cls) -> "ParameterSpace":
        """Activate one fixed echo direction for an explicit M2/M3 ablation."""
        base = cls.plausible_m2()
        bounds = dict(base.bounds)
        bounds["echo_loading"] = ParameterBound(-0.50, 0.50)
        return cls(bounds, {"orthogonal_curvature": 0.0})

    @classmethod
    def plausible_m4(cls) -> "ParameterSpace":
        """Add the orthogonal dispersion loading after creating that bank."""
        base = cls.plausible_m3()
        bounds = dict(base.bounds)
        bounds["orthogonal_curvature"] = ParameterBound(0.001, 0.10, "log")
        return cls(bounds, {})

    @classmethod
    def targeted_stage2(cls) -> "ParameterSpace":
        """Focus on the pilot-supported M2 region while retesting both optional channels."""
        return cls(
            bounds={
                "cluster_loading": ParameterBound(0.10, 0.55),
                "feedback_curvature": ParameterBound(0.10, 0.50, "log"),
                "feedback_shift": ParameterBound(0.25, 0.75),
                "spike_curvature": ParameterBound(0.04, 0.50, "log"),
                "orthogonal_curvature": ParameterBound(0.02, 0.50, "log"),
                "echo_loading": ParameterBound(-0.35, 0.35),
                "cap_level": ParameterBound(3.0, 6.5),
                "cap_sharpness": ParameterBound(5.0, 12.0),
            },
            fixed={},
        )

    @classmethod
    def stage3_local(cls) -> "ParameterSpace":
        """Validate the slower-spike/high-shift ridge selected in Stage 2."""
        return cls(
            bounds={
                "cluster_loading": ParameterBound(0.20, 0.30),
                "feedback_curvature": ParameterBound(0.060, 0.085),
                "feedback_shift": ParameterBound(2.70, 3.30),
                "spike_curvature": ParameterBound(0.18, 0.26),
                "orthogonal_curvature": ParameterBound(0.02, 0.25, "log"),
                "cap_level": ParameterBound(5.0, 6.5),
                "cap_sharpness": ParameterBound(6.0, 9.0),
            },
            fixed={"echo_loading": 0.0},
        )

    @classmethod
    def asymmetric_stage3_local(cls) -> "ParameterSpace":
        """Add a downside-oriented fast-factor shift to the frozen Stage-3 box."""
        base = cls.stage3_local()
        bounds = dict(base.bounds)
        bounds["spike_shift"] = ParameterBound(0.05, 0.50)
        return cls(bounds, dict(base.fixed))


class QmcParameterSampler:
    """Use Latin hypercube or Sobol designs instead of naive IID draws."""

    def __init__(self, space: ParameterSpace) -> None:
        """Keep the exploration domain immutable across replications."""
        self.space = space

    def sample(
        self,
        count: int,
        seed: int,
        method: Literal["latin_hypercube", "sobol"] = "latin_hypercube",
    ) -> tuple[ReadoutParameters, ...]:
        """Return space-filling parameter objects in deterministic order."""
        if count < 1:
            raise ValueError("count must be positive.")
        names = tuple(self.space.bounds)
        unit = self._unit_design(count, len(names), seed, method)
        columns = {
            name: self.space.bounds[name].transform(unit[:, index])
            for index, name in enumerate(names)
        }
        parameters = []
        for row in range(count):
            values = {name: float(columns[name][row]) for name in names}
            values.update(self.space.fixed)
            parameters.append(ReadoutParameters.from_mapping(values))
        return tuple(parameters)

    def frame(self, parameters: tuple[ReadoutParameters, ...]) -> pd.DataFrame:
        """Convert sampled objects into an auditable tabular design."""
        return pd.DataFrame([self.to_record(item) for item in parameters])

    @staticmethod
    def to_record(parameters: ReadoutParameters) -> dict[str, float]:
        """Expose every cheap parameter without private dataclass machinery."""
        return {
            "cluster_loading": parameters.cluster_loading,
            "feedback_curvature": parameters.feedback_curvature,
            "feedback_shift": parameters.feedback_shift,
            "spike_curvature": parameters.spike_curvature,
            "spike_shift": parameters.spike_shift,
            "orthogonal_curvature": parameters.orthogonal_curvature,
            "echo_loading": parameters.echo_loading,
            "cap_level": parameters.cap_level,
            "cap_sharpness": parameters.cap_sharpness,
        }

    @staticmethod
    def _unit_design(
        count: int,
        dimension: int,
        seed: int,
        method: Literal["latin_hypercube", "sobol"],
    ) -> FloatArray:
        """Generate a scrambled low-discrepancy design on the unit cube."""
        if method == "latin_hypercube":
            return qmc.LatinHypercube(dimension, seed=seed).random(count)
        if method == "sobol":
            exponent = int(np.ceil(np.log2(count)))
            return qmc.Sobol(dimension, scramble=True, seed=seed).random_base2(exponent)[:count]
        raise ValueError("Unknown QMC method.")


@dataclass(frozen=True, slots=True)
class StructuralParameters:
    """Outer-loop choices that alter states or the deterministic Q shift."""

    feedback_asset_correlation: float = 0.90
    spike_asset_correlation: float = 0.65
    equity_risk_price: float = 0.35
    feedback_idiosyncratic_risk_price: float = 0.0
    spike_idiosyncratic_risk_price: float = 0.0

    def architecture(self, base: ArchitectureConfig) -> ArchitectureConfig:
        """Change Brownian loadings without changing state dimensions."""
        return replace(
            base,
            feedback=replace(
                base.feedback, asset_correlation=self.feedback_asset_correlation
            ),
            spike=replace(base.spike, asset_correlation=self.spike_asset_correlation),
        )

    def risk_premium(self) -> RiskPremium:
        """Build the deterministic primitive-driver P-to-Q kernel."""
        return RiskPremium(
            equity=self.equity_risk_price,
            feedback_idiosyncratic=self.feedback_idiosyncratic_risk_price,
            spike_idiosyncratic=self.spike_idiosyncratic_risk_price,
        )

    def to_record(self) -> dict[str, float]:
        """Expose outer parameters in result tables."""
        return {
            "feedback_asset_correlation": self.feedback_asset_correlation,
            "spike_asset_correlation": self.spike_asset_correlation,
            "equity_risk_price": self.equity_risk_price,
            "feedback_idiosyncratic_risk_price": self.feedback_idiosyncratic_risk_price,
            "spike_idiosyncratic_risk_price": self.spike_idiosyncratic_risk_price,
        }


class QmcStructuralSampler:
    """Generate an outer design for state-affecting and measure-change choices."""

    def __init__(self, bounds: Mapping[str, ParameterBound] | None = None) -> None:
        """Use broad sign-aware ranges that remain economically interpretable."""
        self.bounds = bounds or {
            "feedback_asset_correlation": ParameterBound(0.50, 0.99),
            "spike_asset_correlation": ParameterBound(0.10, 0.90),
            "equity_risk_price": ParameterBound(0.10, 0.60),
            "feedback_idiosyncratic_risk_price": ParameterBound(-0.50, 0.50),
            "spike_idiosyncratic_risk_price": ParameterBound(-0.50, 0.50),
        }

    def sample(self, count: int, seed: int) -> tuple[StructuralParameters, ...]:
        """Return a reproducible Latin-hypercube outer design."""
        if count < 1:
            raise ValueError("count must be positive.")
        names = tuple(self.bounds)
        unit = qmc.LatinHypercube(len(names), seed=seed).random(count)
        columns = {
            name: self.bounds[name].transform(unit[:, index])
            for index, name in enumerate(names)
        }
        return tuple(
            StructuralParameters(**{name: float(columns[name][row]) for name in names})
            for row in range(count)
        )
