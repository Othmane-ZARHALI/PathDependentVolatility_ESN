import unittest

import numpy as np

from esn_eq.config import ArchitectureConfig, BankConfig, RiskPremium
from esn_eq.outer_sampling import GaussianOuterSampler, pca_root, replication_statistics
from esn_eq.pricing import DatedM2Pricer


def fixture():
    architecture = ArchitectureConfig(feedback=BankConfig((1., 8.), .7, True),
        spike=BankConfig((3., 25.), .4, True), cluster=BankConfig((.2, 2.), 0., False))
    return DatedM2Pricer(architecture, RiskPremium(equity=.2, feedback_idiosyncratic=-.1))


class OuterSamplingTests(unittest.TestCase):
    def test_factor_reconstructs_covariance_without_rank_truncation(self):
        dated = fixture()
        sample = GaussianOuterSampler(dated, .2, np.arange(6)/10)
        np.testing.assert_allclose(sample.root @ sample.root.T, sample.covariance, atol=2e-15)
        self.assertEqual(sample.root.shape, (6, 6))
        self.assertTrue(np.all(np.diff(sample.coordinate_scores) <= 0))
        mean, cov = dated.gaussian_projection_law(.2, np.arange(6)/10)
        np.testing.assert_allclose(dated.projection @ sample.mean, mean, atol=1e-15)
        np.testing.assert_allclose(dated.projection @ sample.covariance @ dated.projection.T, cov, atol=1e-15)

    def test_prefixes_reproducibility_and_distinct_scrambles(self):
        sample = GaussianOuterSampler(fixture(), .2, np.zeros(6))
        for method in ("iid", "sobol"):
            small = sample.draw(5, 12, method)
            np.testing.assert_array_equal(small, sample.draw(6, 12, method)[:32])
            self.assertFalse(np.array_equal(small, sample.draw(5, 13, method)))

    def test_zero_expiry_returns_initial_state(self):
        initial = np.arange(6)/10
        sample = GaussianOuterSampler(fixture(), 0., initial)
        np.testing.assert_array_equal(sample.draw(3, 12), np.broadcast_to(initial, (8, 6)))

    def test_gaussian_mean_and_covariance_from_independent_sobol_scrambles(self):
        sample = GaussianOuterSampler(fixture(), .2, np.zeros(6))
        draws = [sample.draw(10, seed) for seed in range(8)]
        means = np.array([x.mean(axis=0) for x in draws])
        second = np.array([x.T @ x/len(x) for x in draws])
        # A scale-aware numerical correspondence check, not a confidence certificate.
        np.testing.assert_allclose(means.mean(axis=0), sample.mean, atol=.002)
        np.testing.assert_allclose(second.mean(axis=0), sample.covariance+np.outer(sample.mean, sample.mean), atol=.005)

    def test_standard_error_uses_complete_replications(self):
        mean, se = replication_statistics([[1., 0.], [2., 0.], [3., 0.]])
        np.testing.assert_allclose(mean, [2., 0.])
        np.testing.assert_allclose(se, [1/np.sqrt(3), 0.])
        with self.assertRaises(ValueError):
            replication_statistics([[1.]])

    def test_rejects_invalid_laws_and_methods(self):
        with self.assertRaises(ValueError):
            pca_root([[1., 0.], [0., -1.]], np.eye(2))
        sample = GaussianOuterSampler(fixture(), .2, np.zeros(6))
        for method, power in (("unknown", 3), ("sobol", 0), ("iid", 1.5)):
            with self.assertRaises(ValueError):
                sample.draw(power, 1, method)


if __name__ == "__main__":
    unittest.main()
