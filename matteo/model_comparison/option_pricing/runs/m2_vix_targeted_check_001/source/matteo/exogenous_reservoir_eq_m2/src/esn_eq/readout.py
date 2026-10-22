"""Tempered exponential-quadratic variance readout."""

from __future__ import annotations

import numpy as np

from .config import FloatArray, ReadoutParameters
from .features import FeaturePaths


class TemperedExponentialQuadratic:
    """Implement the recommended M2 score, cap, and Q normalisation."""

    @staticmethod
    def score(features: FeaturePaths, parameters: ReadoutParameters) -> FloatArray:
        """Evaluate the identified linear and quadratic score channels."""
        return (
            parameters.cluster_loading * features.cluster
            + parameters.feedback_curvature * (features.feedback - parameters.feedback_shift) ** 2
            + parameters.spike_curvature * (features.spike - parameters.spike_shift) ** 2
            + parameters.orthogonal_curvature * features.orthogonal_energy
            + parameters.echo_loading * features.echo
        )

    @staticmethod
    def tempered_score(score: FloatArray, parameters: ReadoutParameters) -> FloatArray:
        """Apply a numerically stable smooth upper cap h(L, kappa)."""
        scaled_gap = parameters.cap_sharpness * (parameters.cap_level - score)
        softplus = np.logaddexp(0.0, scaled_gap)
        return parameters.cap_level - softplus / parameters.cap_sharpness

    def multiplier(self, features: FeaturePaths, parameters: ReadoutParameters) -> FloatArray:
        """Return G(eta)=exp(2h(eta)) without unsafe overflow."""
        tempered = self.tempered_score(self.score(features, parameters), parameters)
        return np.exp(np.clip(2.0 * tempered, -700.0, 700.0))

    @staticmethod
    def normalizer(q_multiplier: FloatArray) -> FloatArray:
        """Estimate M_Q(t) cross-sectionally from matched Q paths."""
        normalizer = np.mean(q_multiplier, axis=1)
        if np.any(~np.isfinite(normalizer)) or np.any(normalizer <= 0.0):
            raise FloatingPointError("The Q normaliser is not finite and positive.")
        return normalizer
