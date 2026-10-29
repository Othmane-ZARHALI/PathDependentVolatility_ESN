import unittest
import numpy as np
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.stratified_sampling import StratifiedGaussianMixture


class StratifiedSamplingTests(unittest.TestCase):
    def test_allocation_correction_recovers_gaussian_moments(self):
        proposal=DefensiveGaussianMixture([[0.],[2.],[-2.]],[.5,.25,.25])
        sampler=StratifiedGaussianMixture(proposal,[1,1,14],4)
        rows=[]
        for seed in range(8):
            d=sampler.draw(14,seed);z=d.normals[:,0];w=d.weights
            rows.append([w.mean(),(w*z).mean(),(w*z*z).mean()])
        np.testing.assert_allclose(np.mean(rows,axis=0),[1,0,1],atol=.002)
        self.assertGreater(np.max(d.allocation_factors),1.)
        # Deliberately unequal allocation makes an uncorrected plain mean wrong.
        self.assertGreater(abs(np.mean(d.likelihood_weights*z)),.1)

    def test_prefix_preserves_each_component_and_its_weights(self):
        proposal=DefensiveGaussianMixture.symmetric_leading(2)
        sampler=StratifiedGaussianMixture(proposal,[8,2,2,2,2],4)
        small=sampler.draw(4,17);fine=sampler.draw(5,17)
        mask=sampler.prefix_mask(fine,4)
        for a,b in ((small.normals,fine.normals),(small.weights,fine.weights),
                    (small.components,fine.components),(small.component_points,fine.component_points)):
            np.testing.assert_array_equal(a,b[mask])
        self.assertFalse(np.array_equal(np.flatnonzero(mask),np.arange(16)))

    def test_identical_components_and_invalid_allocations(self):
        proposal=DefensiveGaussianMixture([[0.],[0.],[0.]],[.5,.25,.25])
        sampler=StratifiedGaussianMixture(proposal,[1,7,8],4)
        d=sampler.draw(4,12)
        self.assertAlmostEqual(d.weights.mean(),1.,places=14)
        for counts in ([1,0,15],[1.,7.,8.],[1,2,3]):
            with self.assertRaises(ValueError):
                StratifiedGaussianMixture(proposal,counts,4)
        with self.assertRaises(ValueError):
            sampler.prefix_mask(d,5)
