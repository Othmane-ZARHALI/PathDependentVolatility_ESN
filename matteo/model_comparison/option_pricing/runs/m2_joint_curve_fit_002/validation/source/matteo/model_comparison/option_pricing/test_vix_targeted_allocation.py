import unittest
import numpy as np
import pandas as pd
from scipy.special import ndtr
from .vix_targeted_allocation import prefix_mask
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.stratified_sampling import StratifiedGaussianMixture


class TargetedAllocationTests(unittest.TestCase):
    def proposal(self):
        angles=np.arange(16)*2*np.pi/16
        shifts=np.r_[np.zeros((1,2)),3*np.column_stack((np.cos(angles),np.sin(angles))) ]
        return DefensiveGaussianMixture(shifts,[.5]+[1/24]*6+[.025]*10)

    def test_seventeen_component_prefix_and_effective_weights(self):
        p=self.proposal();s=StratifiedGaussianMixture(p,[32]+[2]*16)
        fine=s.draw(7,81);coarse=s.draw(6,81)
        frame=pd.DataFrame(dict(method='fixed_targeted',point=np.arange(128),component=fine.components,component_point=fine.component_points))
        mask=prefix_mask(frame,6)
        np.testing.assert_array_equal(fine.normals[mask],coarse.normals)
        np.testing.assert_array_equal(fine.weights[mask],coarse.weights)
        np.testing.assert_array_equal(np.bincount(fine.components[mask]),[32]+[2]*16)
        np.testing.assert_allclose(coarse.allocation_factors,np.r_[np.ones(32),np.full(12,4/3),np.full(20,.8)],rtol=1e-14)
        frame.method='random_targeted'
        np.testing.assert_array_equal(np.flatnonzero(prefix_mask(frame,6)),np.arange(64))
        self.assertFalse(np.array_equal(np.flatnonzero(mask),np.arange(64)))

    def test_target_density_oracles_under_actual_allocation_pattern(self):
        p=self.proposal();s=StratifiedGaussianMixture(p,[32]+[2]*16);rows=[]
        for seed in range(8):
            d=s.draw(14,seed);z=d.normals;w=d.weights
            rows.append([w.mean(),(w*z[:,0]).mean(),(w*z[:,0]**2).mean(),(w*np.maximum(z[:,0]-3,0)).mean()])
            self.assertLessEqual(w.max(),8/3+1e-12)
        mean=np.mean(rows,axis=0);oracle=np.exp(-4.5)/np.sqrt(2*np.pi)-3*ndtr(-3.)
        np.testing.assert_allclose(mean[:3],[1,0,1],atol=.002)
        self.assertAlmostEqual(mean[3],oracle,delta=3e-5)
