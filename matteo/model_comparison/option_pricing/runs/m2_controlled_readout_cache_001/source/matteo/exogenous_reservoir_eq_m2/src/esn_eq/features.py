"""Fixed linear projections and the optional contractive echo-state layer."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import ArchitectureConfig, EchoStateConfig, FloatArray
from .states import ModeLayout, stationary_covariance


@dataclass(frozen=True, slots=True)
class FeaturePaths:
    """Readout-ready features with time-by-path orientation."""

    feedback: FloatArray
    spike: FloatArray
    cluster: FloatArray
    orthogonal_energy: FloatArray
    echo: FloatArray


def normalized_power_weights(
    layout: ModeLayout,
    bank_slice: slice,
    power: float,
) -> FloatArray:
    """Normalise a rate-power projection to stationary variance one."""
    rates = layout.rates[bank_slice]
    if rates.size == 0:
        return np.empty(0, dtype=float)
    raw = rates**power
    covariance = stationary_covariance(layout)[bank_slice, bank_slice]
    variance = float(raw @ covariance @ raw)
    if variance <= 0.0:
        raise ValueError("A projection has zero stationary variance.")
    return raw / np.sqrt(variance)


class FixedEchoState:
    """Random but fixed bounded reservoir driven only by exogenous OU states."""

    def __init__(self, config: EchoStateConfig, factor_count: int) -> None:
        """Draw matrices once and scale the recurrent operator by spectral norm."""
        self.config = config
        self.factor_count = factor_count
        generator = np.random.default_rng(config.seed)
        raw = generator.standard_normal((config.units, config.units))
        norm = float(np.linalg.svd(raw, compute_uv=False)[0])
        self.recurrent = raw * (config.recurrent_norm / max(norm, 1e-15))
        input_count = 2 * factor_count
        self.input = generator.standard_normal((config.units, input_count))
        self.input *= config.input_scale / np.sqrt(max(input_count, 1))
        self.bias = generator.standard_normal(config.units) * config.bias_scale
        direction = generator.standard_normal(config.units)
        self.direction = direction / np.linalg.norm(direction)

    @property
    def contraction_bound(self) -> float:
        """Return the documented sufficient Lipschitz bound."""
        return (1.0 - self.config.leak) + self.config.leak * self.config.recurrent_norm

    def project(self, states: FloatArray) -> FloatArray:
        """Return a fixed scalar echo feature for every adapted state."""
        if not self.config.enabled:
            return np.zeros(states.shape[:2], dtype=float)
        output = np.empty((*states.shape[:2], self.config.units), dtype=float)
        output[0] = 0.0
        for step in range(states.shape[0] - 1):
            features = np.concatenate((states[step], states[step] ** 2 - 1.0), axis=1)
            drive = output[step] @ self.recurrent.T + features @ self.input.T + self.bias
            candidate = np.tanh(drive)
            output[step + 1] = (
                (1.0 - self.config.leak) * output[step] + self.config.leak * candidate
            )
        return output @ self.direction


class FeatureEngine:
    """Convert raw OU states into the identified score directions."""

    def __init__(self, architecture: ArchitectureConfig, layout: ModeLayout) -> None:
        """Precompute projections because their profiles are hyperparameters."""
        self.architecture = architecture
        self.layout = layout
        self.weights = {
            "feedback": normalized_power_weights(
                layout, layout.slices["feedback"], architecture.feedback_weight_power
            ),
            "spike": normalized_power_weights(
                layout, layout.slices["spike"], architecture.spike_weight_power
            ),
            "cluster": normalized_power_weights(
                layout, layout.slices["cluster"], architecture.cluster_weight_power
            ),
        }
        self.echo_slices = tuple(
            layout.slices[name] for name in architecture.echo.input_banks
        )
        echo_factor_count = sum(item.stop - item.start for item in self.echo_slices)
        if echo_factor_count < 1:
            raise ValueError("The selected echo input banks contain no OU factors.")
        self.echo = FixedEchoState(architecture.echo, echo_factor_count)

    def build(self, states: FloatArray) -> FeaturePaths:
        """Build feature paths without any readout parameter dependence."""
        return FeaturePaths(
            feedback=self._projection(states, "feedback"),
            spike=self._projection(states, "spike"),
            cluster=self._projection(states, "cluster"),
            orthogonal_energy=self._orthogonal_energy(states),
            echo=self.echo.project(self._echo_inputs(states)),
        )

    def _projection(self, states: FloatArray, bank: str) -> FloatArray:
        """Apply one stationary-normalised bank projection."""
        bank_slice = self.layout.slices[bank]
        weights = self.weights[bank]
        if weights.size == 0:
            return np.zeros(states.shape[:2], dtype=float)
        return states[..., bank_slice] @ weights

    def _orthogonal_energy(self, states: FloatArray) -> FloatArray:
        """Return raw or stationary-standardised independent quadratic energy."""
        block = states[..., self.layout.slices["orthogonal"]]
        if block.shape[-1] == 0:
            return np.zeros(states.shape[:2], dtype=float)
        raw_mean = np.mean(block**2, axis=-1)
        if self.architecture.orthogonal_energy_mode == "raw_mean":
            return raw_mean
        stationary_scale = np.sqrt(2.0 / block.shape[-1])
        return (raw_mean - 1.0) / stationary_scale

    def _echo_inputs(self, states: FloatArray) -> FloatArray:
        """Select fixed banks so an orthogonal ablation need not alter M3."""
        blocks = tuple(states[..., item] for item in self.echo_slices if item.stop > item.start)
        if len(blocks) == 1:
            return blocks[0]
        return np.concatenate(blocks, axis=-1)
