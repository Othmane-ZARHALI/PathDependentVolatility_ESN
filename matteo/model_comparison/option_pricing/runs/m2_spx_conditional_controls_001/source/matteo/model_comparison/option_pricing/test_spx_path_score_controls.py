"""Analytic Gaussian quadratic oracle for independently frozen ridge controls."""
import unittest
import numpy as np
from .spx_path_score_controls import dictionary, ridge_fit
from .spx_volatility_sampling import proposal


class PathScoreControlTests(unittest.TestCase):
    def test_independent_quadratic_oracle_and_known_zero_dictionary(self):
        mixture=proposal(1);rng=np.random.default_rng(451)
        def draw(n):
            component=rng.choice(len(mixture.shifts),n,p=mixture.probabilities)
            z=rng.standard_normal((n,4))+mixture.shifts[component]
            w=mixture.weights(z)
            return z,w,dictionary(np.zeros((n,4)),w,z)
        z,w,c=draw(4096)
        y=w*(1.+.3*z[:,0]+.7*z[:,2]**2)
        beta=ridge_fit(y,c)
        z,w,c=draw(60000)
        zero_means=c.mean(axis=0)
        se=c.std(axis=0,ddof=1)/np.sqrt(len(c))
        self.assertTrue(np.all(np.abs(zero_means)<6*se+1e-12))
        raw=w*(1.+.3*z[:,0]+.7*z[:,2]**2)
        adjusted=raw-c@beta
        self.assertLess(abs(adjusted.mean()-1.7),6*adjusted.std(ddof=1)/np.sqrt(len(adjusted)))
        self.assertLess(adjusted.var()/raw.var(),1e-4)
        np.testing.assert_array_equal(c[:,:4],0.)
        self.assertEqual(c.shape[1],13)
