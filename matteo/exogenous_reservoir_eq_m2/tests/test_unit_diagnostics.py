"""Unit tests for physical-measure stylised-fact estimators."""

from __future__ import annotations

import unittest

import numpy as np

from esn_eq.diagnostics import aggregate_result, lagged_effects, rough_hurst, zumbach_statistic
from esn_eq.model import SimulationMetadata, SimulationResult


class DiagnosticTests(unittest.TestCase):
    """Use synthetic positive controls to check statistic orientation."""

    def test_brownian_log_variance_has_half_scaling(self) -> None:
        """A random walk positive control should have H near one half."""
        generator = np.random.default_rng(22)
        log_variance = np.cumsum(0.02 * generator.standard_normal(6000))
        estimate = rough_hurst(np.exp(log_variance), maximum_lag=80)
        self.assertGreater(estimate, 0.40)
        self.assertLess(estimate, 0.60)

    def test_causal_variance_filter_has_positive_zumbach(self) -> None:
        """Past squared shocks feeding future variance should break time reversal."""
        generator = np.random.default_rng(31)
        returns = generator.standard_normal(12000)
        variance = np.empty_like(returns)
        variance[0] = 1.0
        for index in range(1, returns.size):
            variance[index] = 0.92 * variance[index - 1] + 0.08 * returns[index - 1] ** 2
        self.assertGreater(zumbach_statistic(returns, variance, rank_based=True), 0.01)

    def test_negative_shocks_create_negative_rank_leverage(self) -> None:
        """A leverage positive control must have the documented negative sign."""
        generator = np.random.default_rng(9)
        returns = generator.standard_normal(8000)
        variance = np.ones_like(returns)
        for index in range(1, returns.size):
            shock = max(-returns[index - 1], 0.0)
            variance[index] = 0.8 * variance[index - 1] + 0.2 * (1.0 + shock)
        effects = lagged_effects(returns, variance)
        self.assertLess(effects["leverage_rank"], -0.01)

    def test_frequency_aggregation_keeps_returns_and_variance_matched(self) -> None:
        """Coarser returns sum while average instantaneous variance averages."""
        metadata = SimulationMetadata(0.0, 0.0, 0.0, 1.0, 1.0)
        result = SimulationResult(
            np.asarray([[1.0, 2.0, 3.0, 4.0]]),
            np.asarray([[2.0, 4.0, 6.0, 8.0]]),
            252,
            "P",
            metadata,
        )
        aggregated = aggregate_result(result, 2)
        np.testing.assert_array_equal(aggregated.log_returns, [[3.0, 7.0]])
        np.testing.assert_array_equal(aggregated.realized_variance, [[3.0, 7.0]])
        self.assertEqual(aggregated.observations_per_year, 126)


if __name__ == "__main__":
    unittest.main()
