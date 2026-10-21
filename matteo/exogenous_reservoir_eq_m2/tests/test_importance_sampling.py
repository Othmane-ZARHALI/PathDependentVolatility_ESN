import unittest
import numpy as np
from scipy.special import ndtr
from esn_eq.importance_sampling import DefensiveGaussianMixture


class ImportanceSamplingTests(unittest.TestCase):
    def test_density_ratio_and_bound(self):
        mixture = DefensiveGaussianMixture([[0.],[2.],[-2.]],[.5,.25,.25])
        z = np.array([[-100.],[-2.],[0.],[3.],[100.]])
        expected = 1/(.5+.25*np.exp(2*z[:,0]-2)+.25*np.exp(-2*z[:,0]-2))
        np.testing.assert_allclose(mixture.weights(z),expected,rtol=1e-13)
        self.assertTrue(np.all(mixture.weights(z) <= 2.))

    def test_gaussian_moments_and_rare_call_oracle(self):
        mixture = DefensiveGaussianMixture.symmetric_leading(3)
        rows = []
        for seed in range(8):
            z,w,_ = mixture.draw(14,seed)
            rows.append(np.r_[w.mean(),(w[:,None]*z).mean(axis=0),
                (w[:,None]*z**2).mean(axis=0),np.mean(w*np.maximum(z[:,0]-3,0))])
        estimate = np.mean(rows,axis=0)
        oracle = np.exp(-4.5)/np.sqrt(2*np.pi)-3*ndtr(-3.)
        np.testing.assert_allclose(estimate[:7],[1,0,0,0,1,1,1],atol=.002)
        self.assertAlmostEqual(estimate[7],oracle,delta=2e-5)

    def test_prefix_and_rejection(self):
        mixture = DefensiveGaussianMixture.symmetric_leading(4)
        small = mixture.draw(5,13)
        large = mixture.draw(6,13)
        for a,b in zip(small,large):
            np.testing.assert_array_equal(a,b[:32])
        for shifts,probabilities in (([[1.]],[1.]),([[0.]],[.9]),([[0.],[1.]],[1.,0.])):
            with self.assertRaises(ValueError):
                DefensiveGaussianMixture(shifts,probabilities)
