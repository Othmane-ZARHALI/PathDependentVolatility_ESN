"""Unit tests for space-filling designs and latent calibration utilities."""

from __future__ import annotations

import unittest

import numpy as np

from esn_eq.calibration import LatentArchitectureCalibrator, LatentTargets
from esn_eq.config import ArchitectureConfig
from esn_eq.sampling import (
    ParameterSpace,
    QmcParameterSampler,
    QmcStructuralSampler,
)


class SamplingAndCalibrationTests(unittest.TestCase):
    """Check reproducibility, parameter bounds, and latent curve construction."""

    def test_latin_hypercube_is_reproducible_and_within_bounds(self) -> None:
        """A fixed seed must produce exactly the same candidate design."""
        space = ParameterSpace.plausible_m2()
        sampler = QmcParameterSampler(space)
        first = sampler.frame(sampler.sample(16, 17))
        second = sampler.frame(sampler.sample(16, 17))
        np.testing.assert_array_equal(first.to_numpy(), second.to_numpy())
        for name, bound in space.bounds.items():
            self.assertTrue(first[name].between(bound.lower, bound.upper).all())

    def test_targeted_stage2_design_respects_every_declared_bound(self) -> None:
        """The validation sampler must not silently fall back to the broad pilot box."""
        space = ParameterSpace.targeted_stage2()
        sampler = QmcParameterSampler(space)
        frame = sampler.frame(sampler.sample(32, 18, method="sobol"))
        for name, bound in space.bounds.items():
            self.assertTrue(frame[name].between(bound.lower, bound.upper).all())

    def test_stage3_design_fixes_the_rejected_echo_loading(self) -> None:
        """The confirmation design must not spend samples on the null M3 effect."""
        space = ParameterSpace.stage3_local()
        frame = QmcParameterSampler(space).frame(
            QmcParameterSampler(space).sample(16, 20, method="sobol")
        )
        self.assertTrue((frame.echo_loading == 0.0).all())
        for name, bound in space.bounds.items():
            self.assertTrue(frame[name].between(bound.lower, bound.upper).all())

    def test_asymmetric_stage3_design_samples_only_positive_spike_shifts(self) -> None:
        """The asymmetric extension keeps the zero-shift model as a separate control."""
        space = ParameterSpace.asymmetric_stage3_local()
        frame = QmcParameterSampler(space).frame(
            QmcParameterSampler(space).sample(16, 208, method="sobol")
        )
        self.assertTrue(frame.spike_shift.between(0.05, 0.50).all())
        self.assertTrue((frame.echo_loading == 0.0).all())

    def test_architecture_json_round_trip_preserves_nested_configuration(self) -> None:
        """A frozen pilot architecture must be reusable without recalibration."""
        from dataclasses import asdict

        architecture = ArchitectureConfig(orthogonal_energy_mode="standardized_centered")
        rebuilt = ArchitectureConfig.from_mapping(asdict(architecture))
        self.assertEqual(rebuilt, architecture)

    def test_structural_sampler_respects_correlation_domain(self) -> None:
        """Outer designs must remain valid Brownian loading scenarios."""
        scenarios = QmcStructuralSampler().sample(12, 19)
        for scenario in scenarios:
            self.assertLess(abs(scenario.feedback_asset_correlation), 1.0)
            self.assertLess(abs(scenario.spike_asset_correlation), 1.0)

    def test_latent_curves_are_finite_and_normalised(self) -> None:
        """The architecture-only diagnostics must be deterministic and finite."""
        calibrator = LatentArchitectureCalibrator()
        rough, cluster = calibrator.curves(ArchitectureConfig(), LatentTargets())
        self.assertTrue(np.all(np.isfinite(rough)))
        self.assertTrue(np.all(np.isfinite(cluster)))
        self.assertTrue(np.all(rough > 0.0))
        self.assertTrue(np.all((cluster > 0.0) & (cluster < 1.0)))

    def test_replicated_calibration_retains_every_optimizer_seed(self) -> None:
        """Stochastic calibration must not report only a favourable single run."""
        ensemble = LatentArchitectureCalibrator().fit_replicated(
            ArchitectureConfig(),
            LatentTargets(),
            seeds=(23, 29),
            max_iterations=1,
        )
        self.assertEqual(ensemble.seeds, (23, 29))
        self.assertEqual(len(ensemble.results), 2)
        self.assertEqual(
            ensemble.best.objective,
            min(item.objective for item in ensemble.results),
        )


if __name__ == "__main__":
    unittest.main()
