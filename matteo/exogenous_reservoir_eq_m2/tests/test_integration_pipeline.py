"""Integration tests for preparation, simulation, and hierarchical exploration."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from esn_eq import (
    AcceptanceCriteria,
    ArchitectureConfig,
    ExogenousReservoirVolatilityModel,
    HierarchicalExplorer,
    MarketEnvironment,
    ParameterSpace,
    QmcParameterSampler,
    ReadoutParameters,
    ReservoirFactory,
    RiskPremium,
    SimulationConfig,
    StructuralParameters,
)


def research_simulation() -> SimulationConfig:
    """Return a small but statistically usable integration horizon."""
    return SimulationConfig(
        years=1.6,
        burn_years=0.2,
        paths=4,
        observations_per_year=100,
        steps_per_observation=2,
    )


class PipelineIntegrationTests(unittest.TestCase):
    """Exercise the complete common-random-number research workflow."""

    def test_readout_changes_do_not_mutate_prepared_features(self) -> None:
        """Path reuse means bitwise-stable features, not identical final returns."""
        architecture = ArchitectureConfig()
        cache = ReservoirFactory(architecture, research_simulation()).build(101)
        model = ExogenousReservoirVolatilityModel(architecture)
        prepared = model.prepare(cache, RiskPremium())
        saved_feedback = prepared.physical_features.feedback.copy()
        first = model.simulate(prepared, ReadoutParameters(), MarketEnvironment())
        second = model.simulate(
            prepared,
            ReadoutParameters(feedback_curvature=0.22, cluster_loading=0.60),
            MarketEnvironment(),
        )
        np.testing.assert_array_equal(saved_feedback, prepared.physical_features.feedback)
        self.assertFalse(np.array_equal(first.log_returns, second.log_returns))
        self.assertLess(first.metadata.q_forward_mean_max_error, 1e-12)
        self.assertLess(second.metadata.q_forward_mean_max_error, 1e-12)

    def test_hierarchical_experiment_saves_every_candidate(self) -> None:
        """Two outer by three inner candidates must produce six retained rows."""
        architecture = ArchitectureConfig()
        readouts = QmcParameterSampler(ParameterSpace.plausible_m2()).sample(3, 202)
        structures = (
            StructuralParameters(),
            StructuralParameters(feedback_asset_correlation=0.70, equity_risk_price=0.20),
        )
        explorer = HierarchicalExplorer(
            architecture,
            research_simulation(),
            MarketEnvironment(),
            AcceptanceCriteria.illustrative(),
        )
        result = explorer.run(structures, readouts, seed=203)
        self.assertEqual(len(result.candidates), 6)
        self.assertIn("feedback_asset_correlation", result.candidates)
        self.assertTrue(result.candidates["distance"].is_monotonic_increasing)
        with tempfile.TemporaryDirectory() as directory:
            paths = result.save(directory)
            self.assertTrue(all(Path(path).is_file() for path in paths))


if __name__ == "__main__":
    unittest.main()

