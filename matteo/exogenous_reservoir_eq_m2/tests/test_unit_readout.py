"""Unit tests for fixed features and the tempered readout."""

from __future__ import annotations

import unittest
from dataclasses import replace

import numpy as np

from esn_eq.config import EchoStateConfig, ReadoutParameters
from esn_eq.features import FeaturePaths, FixedEchoState
from esn_eq.readout import TemperedExponentialQuadratic


def feature_paths(feedback: np.ndarray) -> FeaturePaths:
    """Create otherwise-zero features for isolated feedback tests."""
    zeros = np.zeros_like(feedback, dtype=float)
    return FeaturePaths(feedback, zeros, zeros, zeros, zeros)


class TemperedReadoutTests(unittest.TestCase):
    """Check monotonicity, leverage orientation, and Q normalisation."""

    def test_tempering_is_monotone_and_upper_bounded(self) -> None:
        """The smooth cap must stay stable on extreme scores."""
        parameters = ReadoutParameters(cap_level=3.0, cap_sharpness=5.0)
        scores = np.linspace(-100.0, 100.0, 1001)
        tempered = TemperedExponentialQuadratic.tempered_score(scores, parameters)
        self.assertTrue(np.all(np.diff(tempered) >= -1e-14))
        self.assertLessEqual(float(np.max(tempered)), parameters.cap_level)
        self.assertTrue(np.all(np.isfinite(tempered)))

    def test_shifted_parabola_has_negative_local_orientation(self) -> None:
        """For q,beta>0, a positive feedback shock near zero lowers eta."""
        parameters = ReadoutParameters(
            cluster_loading=0.0,
            feedback_curvature=0.2,
            feedback_shift=0.5,
            spike_curvature=0.0,
        )
        feedback = np.asarray([[0.0, 0.1]])
        score = TemperedExponentialQuadratic.score(feature_paths(feedback), parameters)
        self.assertLess(score[0, 1], score[0, 0])

    def test_spike_shift_breaks_equal_magnitude_symmetry_toward_downside(self) -> None:
        """A positive shift must make a negative fast shock more volatile."""
        parameters = ReadoutParameters(
            cluster_loading=0.0,
            feedback_curvature=0.0,
            spike_curvature=0.2,
            spike_shift=0.4,
        )
        spike = np.asarray([[-2.0, 2.0]])
        features = replace(feature_paths(np.zeros_like(spike)), spike=spike)
        score = TemperedExponentialQuadratic.score(features, parameters)
        self.assertGreater(score[0, 0], score[0, 1])
        self.assertAlmostEqual(float(score[0, 0] - score[0, 1]), 0.64)

    def test_zero_spike_shift_preserves_the_original_readout(self) -> None:
        """The new parameter must provide an exact nested Stage-3 control."""
        spike = np.linspace(-3.0, 3.0, 13)[None, :]
        features = replace(feature_paths(np.zeros_like(spike)), spike=spike)
        parameters = ReadoutParameters(spike_curvature=0.23, spike_shift=0.0)
        score = TemperedExponentialQuadratic.score(features, parameters)
        expected = (
            parameters.feedback_curvature
            * (features.feedback - parameters.feedback_shift) ** 2
            + parameters.spike_curvature * spike**2
        )
        np.testing.assert_array_equal(score, expected)

    def test_q_normalizer_reproduces_cross_sectional_forward_mean(self) -> None:
        """Dividing by M_Q must produce a unit cross-sectional multiplier."""
        multipliers = np.asarray([[1.0, 2.0, 3.0], [0.5, 1.5, 4.0]])
        normalizer = TemperedExponentialQuadratic.normalizer(multipliers)
        np.testing.assert_allclose(np.mean(multipliers / normalizer[:, None], axis=1), 1.0)

    def test_echo_layer_is_contractive_and_bounded(self) -> None:
        """The optional layer cannot violate its sufficient echo-state bound."""
        config = EchoStateConfig(enabled=True, units=6, recurrent_norm=0.7, leak=0.4)
        echo = FixedEchoState(config, factor_count=3)
        states = np.random.default_rng(4).standard_normal((40, 5, 3))
        projection = echo.project(states)
        self.assertLess(echo.contraction_bound, 1.0)
        self.assertLessEqual(float(np.max(np.abs(projection))), np.sqrt(config.units) + 1e-12)


if __name__ == "__main__":
    unittest.main()
