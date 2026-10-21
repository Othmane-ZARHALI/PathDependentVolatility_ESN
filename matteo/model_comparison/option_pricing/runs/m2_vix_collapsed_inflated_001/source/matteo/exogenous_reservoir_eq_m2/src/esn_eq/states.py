"""Exact joint OU/equity transitions and common-random-number caches."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import ArchitectureConfig, FloatArray, RiskPremium, SimulationConfig


@dataclass(frozen=True, slots=True)
class ModeLayout:
    """Flattened OU metadata and Brownian loadings."""

    rates: FloatArray
    driver_loadings: FloatArray
    primitive_names: tuple[str, ...]
    slices: dict[str, slice]


@dataclass(frozen=True, slots=True)
class PrimitiveRandomness:
    """Parameter-free standard normals that can serve many architectures."""

    transition_normals: FloatArray
    initial_normals: FloatArray


@dataclass(frozen=True, slots=True)
class ReservoirCache:
    """Adapted OU states and matched equity Brownian increments."""

    states: FloatArray
    equity_increments: FloatArray
    layout: ModeLayout
    architecture: ArchitectureConfig
    simulation: SimulationConfig
    randomness: PrimitiveRandomness


def build_mode_layout(architecture: ArchitectureConfig) -> ModeLayout:
    """Represent each bank driver as loadings on independent primitives."""
    names = ["equity"]
    bank_primitives: dict[str, list[int]] = {}
    for bank_name, bank in architecture.banks.items():
        count = 1 if bank.shared_driver and bank.rates else len(bank.rates)
        bank_primitives[bank_name] = list(range(len(names), len(names) + count))
        names.extend(f"{bank_name}:{index}" for index in range(count))
    rates, rows, slices = _layout_rows(architecture, names, bank_primitives)
    return ModeLayout(np.asarray(rates), np.asarray(rows), tuple(names), slices)


def _layout_rows(
    architecture: ArchitectureConfig,
    names: list[str],
    bank_primitives: dict[str, list[int]],
) -> tuple[list[float], list[FloatArray], dict[str, slice]]:
    """Build factor rows while preserving canonical bank slices."""
    rates: list[float] = []
    rows: list[FloatArray] = []
    slices: dict[str, slice] = {}
    for bank_name, bank in architecture.banks.items():
        start = len(rates)
        for index, rate in enumerate(bank.rates):
            row = np.zeros(len(names), dtype=float)
            row[0] = bank.asset_correlation
            primitive = bank_primitives[bank_name][0 if bank.shared_driver else index]
            row[primitive] = np.sqrt(max(0.0, 1.0 - bank.asset_correlation**2))
            rates.append(rate)
            rows.append(row)
        slices[bank_name] = slice(start, len(rates))
    return rates, rows, slices


def driver_correlation(layout: ModeLayout) -> FloatArray:
    """Return instantaneous correlations among all factor drivers."""
    return layout.driver_loadings @ layout.driver_loadings.T


def stationary_covariance(layout: ModeLayout) -> FloatArray:
    """Return the exact stationary covariance of unit-variance OU modes."""
    rates = layout.rates
    kernel = 2.0 * np.sqrt(np.outer(rates, rates)) / (rates[:, None] + rates[None, :])
    return driver_correlation(layout) * kernel


def joint_transition_covariance(layout: ModeLayout, dt: float) -> FloatArray:
    """Couple one equity Brownian increment to exact OU innovations."""
    rates = layout.rates
    sums = rates[:, None] + rates[None, :]
    mode = driver_correlation(layout) * (
        2.0 * np.sqrt(np.outer(rates, rates)) * (-np.expm1(-sums * dt)) / sums
    )
    asset = layout.driver_loadings[:, 0] * (
        np.sqrt(2.0 * rates) * (-np.expm1(-rates * dt)) / rates
    )
    covariance = np.empty((rates.size + 1, rates.size + 1), dtype=float)
    covariance[0, 0] = dt
    covariance[0, 1:] = asset
    covariance[1:, 0] = asset
    covariance[1:, 1:] = mode
    return covariance


def symmetric_square_root(covariance: FloatArray) -> FloatArray:
    """Compute a stable square root of a positive-semidefinite matrix."""
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    tolerance = 1e-11 * max(1.0, float(np.max(eigenvalues)))
    if float(np.min(eigenvalues)) < -tolerance:
        raise ValueError("Innovation covariance is not positive semidefinite.")
    return (eigenvectors * np.sqrt(np.maximum(eigenvalues, 0.0))) @ eigenvectors.T


class ReservoirFactory:
    """Build immutable exact-transition paths from reusable standard normals."""

    def __init__(self, architecture: ArchitectureConfig, simulation: SimulationConfig) -> None:
        """Fix state-changing choices without drawing any randomness."""
        self.architecture = architecture
        self.simulation = simulation
        self.layout = build_mode_layout(architecture)

    def draw(self, seed: int) -> PrimitiveRandomness:
        """Draw primitive normals once for common-random-number experiments."""
        generator = np.random.default_rng(seed)
        joint_shape = (self.simulation.steps, self.simulation.paths, self.layout.rates.size + 1)
        initial_shape = (self.simulation.paths, self.layout.rates.size)
        return PrimitiveRandomness(
            generator.standard_normal(joint_shape),
            generator.standard_normal(initial_shape),
        )

    def build(self, seed: int) -> ReservoirCache:
        """Generate an exact stationary P cache from a seed."""
        return self.from_randomness(self.draw(seed))

    def from_randomness(self, randomness: PrimitiveRandomness) -> ReservoirCache:
        """Rebuild states deterministically when structural choices change."""
        self._validate_randomness(randomness)
        stationary_root = symmetric_square_root(stationary_covariance(self.layout))
        equity, innovations = self._transition_innovations(randomness.transition_normals)
        initial = randomness.initial_normals @ stationary_root.T
        states = self._propagate(initial, innovations)
        return ReservoirCache(
            states, equity, self.layout, self.architecture, self.simulation, randomness
        )

    def _transition_innovations(self, normals: FloatArray) -> tuple[FloatArray, FloatArray]:
        """Condition on the first primitive so equity shocks remain bitwise fixed."""
        covariance = joint_transition_covariance(self.layout, self.simulation.dt)
        equity = np.sqrt(self.simulation.dt) * normals[..., 0]
        asset_covariance = covariance[1:, 0]
        conditional = covariance[1:, 1:] - (
            np.outer(asset_covariance, asset_covariance) / self.simulation.dt
        )
        conditional_root = symmetric_square_root(conditional)
        innovations = normals[..., 1:] @ conditional_root.T
        innovations += equity[..., None] * asset_covariance / self.simulation.dt
        return equity, innovations

    def _propagate(self, initial: FloatArray, innovations: FloatArray) -> FloatArray:
        """Store start-of-step states so the volatility integrand is adapted."""
        decay = np.exp(-self.layout.rates * self.simulation.dt)
        states = np.empty((self.simulation.steps + 1, *initial.shape), dtype=float)
        states[0] = initial
        for step in range(self.simulation.steps):
            states[step + 1] = decay * states[step] + innovations[step]
        return states

    def _validate_randomness(self, randomness: PrimitiveRandomness) -> None:
        """Prevent accidental reuse across incompatible dimensions."""
        expected_joint = (self.simulation.steps, self.simulation.paths, self.layout.rates.size + 1)
        expected_initial = (self.simulation.paths, self.layout.rates.size)
        if randomness.transition_normals.shape != expected_joint:
            raise ValueError("Transition normals have an incompatible shape.")
        if randomness.initial_normals.shape != expected_initial:
            raise ValueError("Initial normals have an incompatible shape.")


def primitive_risk_prices(layout: ModeLayout, premium: RiskPremium) -> FloatArray:
    """Resolve named deterministic prices on independent Brownian primitives."""
    values = np.zeros(len(layout.primitive_names), dtype=float)
    values[0] = premium.equity
    cluster = list(premium.cluster_idiosyncratic)
    orthogonal = list(premium.orthogonal_idiosyncratic)
    _validate_optional_prices(layout, "cluster", cluster)
    _validate_optional_prices(layout, "orthogonal", orthogonal)
    for index, name in enumerate(layout.primitive_names[1:], start=1):
        bank, position = name.split(":")
        if bank == "feedback":
            values[index] = premium.feedback_idiosyncratic
        elif bank == "spike":
            values[index] = premium.spike_idiosyncratic
        elif bank == "cluster" and cluster:
            values[index] = cluster[int(position)]
        elif bank == "orthogonal" and orthogonal:
            values[index] = orthogonal[int(position)]
    return values


def _validate_optional_prices(layout: ModeLayout, bank: str, values: list[float]) -> None:
    """Require either no prices or one per primitive of an independent bank."""
    required = sum(name.startswith(f"{bank}:") for name in layout.primitive_names)
    if values and len(values) != required:
        raise ValueError(f"{bank} risk prices need {required} entries.")


def q_state_paths(cache: ReservoirCache, premium: RiskPremium) -> FloatArray:
    """Apply the exact deterministic Q mean shift from the common initial state."""
    primitive = primitive_risk_prices(cache.layout, premium)
    effective = cache.layout.driver_loadings @ primitive
    rates = cache.layout.rates
    times = np.arange(cache.simulation.steps + 1, dtype=float) * cache.simulation.dt
    response = -np.sqrt(2.0 * rates) * effective / rates
    shift = (-np.expm1(-times[:, None] * rates[None, :])) * response[None, :]
    return cache.states + shift[:, None, :]
