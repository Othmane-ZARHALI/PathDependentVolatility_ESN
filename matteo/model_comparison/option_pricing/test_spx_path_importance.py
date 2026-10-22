"""Target-zero weighted controls and independent-pilot estimator correspondence."""
import unittest
import numpy as np
from .spx_path_importance import weighted_controls
from .spx_hedging_study import beta_fit
from esn_eq.importance_sampling import DefensiveGaussianMixture


class WeightedControlTests(unittest.TestCase):
    def test_target_zero_controls_require_full_likelihood_and_no_mass_rescale(self):
        rng = np.random.default_rng(94)
        mixture = DefensiveGaussianMixture([[0.], [2.], [-2.]], [.5,.25,.25])
        def sample():
            comp = rng.choice(3, 40000, p=mixture.probabilities)
            z = rng.standard_normal(40000)+mixture.shifts[comp,0]
            w = mixture.weights(z[:,None])
            # E_f[z²+z]=1, target-zero controls z and z²-1.
            return weighted_controls(z*z+z, np.c_[z,z*z-1.], w)
        y, c = sample()
        beta = beta_fit(y, c)
        check_y, check_c = sample()
        estimate = check_y-check_c@beta
        # Algebraically exact quadratic oracle, evaluated on independent draws.
        np.testing.assert_allclose(estimate, 1., atol=2e-13)
        np.testing.assert_allclose(beta, [1.,1.,1.], atol=2e-13)
        self.assertGreater(np.max(np.abs(check_c.mean(axis=0))), 1e-4)
        # Values are unnormalized w*control and w-1, including finite sample mass.
        raw, control = weighted_controls(np.array([2.,3.]),np.array([[4.],[5.]]),np.array([.2,1.5]))
        np.testing.assert_array_equal(raw, [.4,4.5])
        np.testing.assert_array_equal(control, [[.8,-.8],[7.5,.5]])
