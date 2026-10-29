"""Typed model, architecture, and simulation configuration."""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Literal, Mapping, cast

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def geometric_rates(low: float, high: float, count: int) -> tuple[float, ...]:
    """Create a positive geometric rate grid in inverse years."""
    if count < 0 or low <= 0.0 or high < low:
        raise ValueError("Require count >= 0 and 0 < low <= high.")
    if count == 0:
        return ()
    return tuple(float(value) for value in np.geomspace(low, high, count))


@dataclass(frozen=True, slots=True)
class BankConfig:
    """One exogenous OU bank and its Brownian-driver convention."""

    rates: tuple[float, ...]
    asset_correlation: float = 0.0
    shared_driver: bool = False

    def __post_init__(self) -> None:
        """Reject invalid decay rates and Brownian correlations."""
        if any(rate <= 0.0 or not np.isfinite(rate) for rate in self.rates):
            raise ValueError("Every OU decay rate must be finite and positive.")
        if abs(self.asset_correlation) > 1.0:
            raise ValueError("asset_correlation must lie in [-1, 1].")


@dataclass(frozen=True, slots=True)
class EchoStateConfig:
    """Fixed optional nonlinear second layer; disabled in recommended M2."""

    enabled: bool = False
    units: int = 12
    leak: float = 0.35
    recurrent_norm: float = 0.75
    input_scale: float = 0.20
    bias_scale: float = 0.05
    seed: int = 20260901
    input_banks: tuple[str, ...] = ("feedback", "spike", "cluster")

    def __post_init__(self) -> None:
        """Enforce the sufficient contractive echo-state condition."""
        if self.units < 1:
            raise ValueError("units must be positive.")
        if not 0.0 < self.leak <= 1.0:
            raise ValueError("leak must lie in (0, 1].")
        if not 0.0 <= self.recurrent_norm < 1.0:
            raise ValueError("recurrent_norm must lie in [0, 1).")
        if self.input_scale < 0.0 or self.bias_scale < 0.0:
            raise ValueError("Echo-state scales must be non-negative.")
        allowed = {"feedback", "spike", "cluster", "orthogonal"}
        if not self.input_banks or len(set(self.input_banks)) != len(self.input_banks):
            raise ValueError("input_banks must be non-empty and contain no duplicates.")
        if any(name not in allowed for name in self.input_banks):
            raise ValueError("input_banks contains an unknown OU bank.")


def _feedback_bank() -> BankConfig:
    """Return the recommended shared-driver rough lift."""
    return BankConfig(geometric_rates(0.5, 1000.0, 12), 0.90, True)


def _spike_bank() -> BankConfig:
    """Return a fast, partially equity-correlated spike bank."""
    return BankConfig(geometric_rates(20.0, 2000.0, 4), 0.65, True)


def _cluster_bank() -> BankConfig:
    """Return independent slow modes for exogenous clustering."""
    return BankConfig(geometric_rates(0.05, 25.0, 8), 0.0, False)


@dataclass(frozen=True, slots=True)
class ArchitectureConfig:
    """State-changing hyperparameters and fixed feature projections."""

    feedback: BankConfig = field(default_factory=_feedback_bank)
    spike: BankConfig = field(default_factory=_spike_bank)
    cluster: BankConfig = field(default_factory=_cluster_bank)
    orthogonal: BankConfig = field(default_factory=lambda: BankConfig(()))
    feedback_weight_power: float = 0.40
    spike_weight_power: float = 0.0
    cluster_weight_power: float = -0.50
    echo: EchoStateConfig = field(default_factory=EchoStateConfig)
    orthogonal_energy_mode: Literal["raw_mean", "standardized_centered"] = "raw_mean"

    def __post_init__(self) -> None:
        """Require the two load-bearing feedback and clustering banks."""
        if not self.feedback.rates or not self.cluster.rates:
            raise ValueError("Feedback and clustering banks cannot be empty.")
        powers = (
            self.feedback_weight_power,
            self.spike_weight_power,
            self.cluster_weight_power,
        )
        if any(not np.isfinite(value) for value in powers):
            raise ValueError("Projection powers must be finite.")
        if self.orthogonal_energy_mode not in {"raw_mean", "standardized_centered"}:
            raise ValueError("Unknown orthogonal_energy_mode.")

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> "ArchitectureConfig":
        """Rebuild a nested architecture from an auditable JSON record."""
        banks = {
            name: cls._bank_from_mapping(cast(Mapping[str, object], values[name]))
            for name in ("feedback", "spike", "cluster", "orthogonal")
        }
        echo_values = cast(Mapping[str, object], values["echo"])
        echo = EchoStateConfig(
            enabled=bool(echo_values["enabled"]),
            units=int(cast(int, echo_values["units"])),
            leak=float(cast(float, echo_values["leak"])),
            recurrent_norm=float(cast(float, echo_values["recurrent_norm"])),
            input_scale=float(cast(float, echo_values["input_scale"])),
            bias_scale=float(cast(float, echo_values["bias_scale"])),
            seed=int(cast(int, echo_values["seed"])),
            input_banks=tuple(cast(tuple[str, ...], echo_values["input_banks"])),
        )
        return cls(
            **banks,
            feedback_weight_power=float(cast(float, values["feedback_weight_power"])),
            spike_weight_power=float(cast(float, values["spike_weight_power"])),
            cluster_weight_power=float(cast(float, values["cluster_weight_power"])),
            echo=echo,
            orthogonal_energy_mode=cast(
                Literal["raw_mean", "standardized_centered"],
                values.get("orthogonal_energy_mode", "raw_mean"),
            ),
        )

    @staticmethod
    def _bank_from_mapping(values: Mapping[str, object]) -> BankConfig:
        """Convert one serialised OU-bank record back to typed configuration."""
        return BankConfig(
            rates=tuple(float(value) for value in cast(tuple[float, ...], values["rates"])),
            asset_correlation=float(cast(float, values["asset_correlation"])),
            shared_driver=bool(values["shared_driver"]),
        )

    @property
    def factor_count(self) -> int:
        """Return the total number of OU coordinates."""
        return sum(len(bank.rates) for bank in self.banks.values())

    @property
    def banks(self) -> dict[str, BankConfig]:
        """Expose banks in the canonical state-vector order."""
        return {
            "feedback": self.feedback,
            "spike": self.spike,
            "cluster": self.cluster,
            "orthogonal": self.orthogonal,
        }


@dataclass(frozen=True, slots=True)
class SimulationConfig:
    """Horizon and output frequency; OU rates remain in inverse years."""

    years: float = 8.0
    burn_years: float = 2.0
    paths: int = 8
    observations_per_year: int = 252
    steps_per_observation: int = 4

    def __post_init__(self) -> None:
        """Validate dimensions before allocating the reusable path cache."""
        if self.years <= self.burn_years or self.burn_years < 0.0:
            raise ValueError("Require years > burn_years >= 0.")
        if self.paths < 2:
            raise ValueError("At least two paths are needed for Q normalisation.")
        if self.observations_per_year < 1 or self.steps_per_observation < 1:
            raise ValueError("Discrete frequencies must be positive.")

    @property
    def dt(self) -> float:
        """Return the simulation step in years."""
        return 1.0 / (self.observations_per_year * self.steps_per_observation)

    @property
    def observations(self) -> int:
        """Return the rounded number of output observations."""
        return int(round(self.years * self.observations_per_year))

    @property
    def burn_observations(self) -> int:
        """Return the number of discarded output observations."""
        return int(round(self.burn_years * self.observations_per_year))

    @property
    def steps(self) -> int:
        """Return the total number of fine simulation steps."""
        return self.observations * self.steps_per_observation


@dataclass(frozen=True, slots=True)
class ReadoutParameters:
    """Cheap parameters evaluated on fixed P/Q feature paths."""

    cluster_loading: float = 0.35
    feedback_curvature: float = 0.12
    feedback_shift: float = 0.35
    spike_curvature: float = 0.05
    spike_shift: float = 0.0
    orthogonal_curvature: float = 0.0
    echo_loading: float = 0.0
    cap_level: float = 3.5
    cap_sharpness: float = 6.0

    def __post_init__(self) -> None:
        """Keep convex channels and smooth tempering in their valid domains."""
        values = tuple(getattr(self, item.name) for item in fields(self))
        if any(not np.isfinite(value) for value in values):
            raise ValueError("All readout parameters must be finite.")
        curvatures = (
            self.feedback_curvature,
            self.spike_curvature,
            self.orthogonal_curvature,
        )
        if any(value < 0.0 for value in curvatures):
            raise ValueError("Quadratic curvatures must be non-negative.")
        if self.cap_sharpness <= 0.0:
            raise ValueError("cap_sharpness must be positive.")

    @classmethod
    def from_mapping(cls, values: Mapping[str, float]) -> "ReadoutParameters":
        """Build a parameter object from a sampled name-to-value mapping."""
        return cls(**{name: float(value) for name, value in values.items()})


@dataclass(frozen=True, slots=True)
class RiskPremium:
    """Constant primitive-driver prices for the deterministic P-to-Q shift."""

    equity: float = 0.35
    feedback_idiosyncratic: float = 0.0
    spike_idiosyncratic: float = 0.0
    cluster_idiosyncratic: tuple[float, ...] = ()
    orthogonal_idiosyncratic: tuple[float, ...] = ()

    def __post_init__(self) -> None:
        """Reject non-finite market prices of risk."""
        values = (
            self.equity,
            self.feedback_idiosyncratic,
            self.spike_idiosyncratic,
            *self.cluster_idiosyncratic,
            *self.orthogonal_idiosyncratic,
        )
        if any(not np.isfinite(value) for value in values):
            raise ValueError("Risk-premium entries must be finite.")


@dataclass(frozen=True, slots=True)
class MarketEnvironment:
    """Daily inputs kept separate from structural model parameters."""

    forward_variance: float | tuple[float, ...] = 0.20**2
    variance_floor: float = 0.01**2
    risk_free_rate: float = 0.0
    dividend_yield: float = 0.0

    def curve(self, steps: int) -> FloatArray:
        """Return one strictly-above-floor forward variance per fine step."""
        if np.isscalar(self.forward_variance):
            values = np.full(steps, float(self.forward_variance), dtype=float)
        else:
            values = np.asarray(self.forward_variance, dtype=float)
            if values.size != steps:
                raise ValueError("A forward-variance tuple must have one value per step.")
        if np.any(~np.isfinite(values)) or np.any(values <= self.variance_floor):
            raise ValueError("Forward variance must be finite and above the floor.")
        return values
