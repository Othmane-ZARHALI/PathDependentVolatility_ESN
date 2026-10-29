"""Unit tests for exact transitions and reusable randomness."""

from __future__ import annotations

import unittest
from dataclasses import replace

import numpy as np

from esn_eq.config import ArchitectureConfig, BankConfig, RiskPremium, SimulationConfig
from esn_eq.states import (
    ReservoirFactory,
    build_mode_layout,
    joint_transition_covariance,
    q_state_paths,
    stationary_covariance,
)


def small_architecture() -> ArchitectureConfig:
    """Return a low-dimensional architecture for deterministic unit tests."""
    return ArchitectureConfig(
        feedback=BankConfig((1.0, 8.0), 0.75, True),
        spike=BankConfig((20.0,), 0.40, True),
        cluster=BankConfig((0.2, 2.0), 0.0, False),
        orthogonal=BankConfig((3.0,), 0.0, False),
    )


class ExactTransitionTests(unittest.TestCase):
    """Check analytical covariances and the P/Q deterministic shift."""

    def test_covariances_are_symmetric_positive_semidefinite(self) -> None:
        """The exact Gaussian construction must define valid covariances."""
        layout = build_mode_layout(small_architecture())
        covariances = (
            stationary_covariance(layout),
            joint_transition_covariance(layout, 0.01),
        )
        for covariance in covariances:
            np.testing.assert_allclose(covariance, covariance.T, atol=1e-14)
            self.assertGreaterEqual(float(np.min(np.linalg.eigvalsh(covariance))), -1e-11)

    def test_stationary_marginal_variances_are_one(self) -> None:
        """Every standardised OU coordinate has unit stationary variance."""
        layout = build_mode_layout(small_architecture())
        np.testing.assert_allclose(np.diag(stationary_covariance(layout)), 1.0)

    def test_q_shift_starts_at_common_state_and_is_deterministic(self) -> None:
        """Changing deterministic risk premia must not draw new noise."""
        simulation = SimulationConfig(
            years=0.2,
            burn_years=0.0,
            paths=3,
            observations_per_year=20,
            steps_per_observation=2,
        )
        cache = ReservoirFactory(small_architecture(), simulation).build(8)
        shifted = q_state_paths(cache, RiskPremium(equity=0.4, feedback_idiosyncratic=0.2))
        np.testing.assert_array_equal(shifted[0], cache.states[0])
        self.assertFalse(np.array_equal(shifted[-1], cache.states[-1]))
        difference = shifted - cache.states
        np.testing.assert_allclose(difference[:, 0], difference[:, 1])

    def test_equity_shocks_are_fixed_when_correlation_changes(self) -> None:
        """Outer structural sweeps must retain the same primitive equity path."""
        simulation = SimulationConfig(
            years=0.2,
            burn_years=0.0,
            paths=3,
            observations_per_year=20,
            steps_per_observation=2,
        )
        first_architecture = small_architecture()
        first_factory = ReservoirFactory(first_architecture, simulation)
        randomness = first_factory.draw(13)
        first = first_factory.from_randomness(randomness)
        changed = replace(
            first_architecture,
            feedback=replace(first_architecture.feedback, asset_correlation=0.20),
        )
        second = ReservoirFactory(changed, simulation).from_randomness(randomness)
        np.testing.assert_array_equal(first.equity_increments, second.equity_increments)
        self.assertFalse(np.array_equal(first.states, second.states))


if __name__ == "__main__":
    unittest.main()
