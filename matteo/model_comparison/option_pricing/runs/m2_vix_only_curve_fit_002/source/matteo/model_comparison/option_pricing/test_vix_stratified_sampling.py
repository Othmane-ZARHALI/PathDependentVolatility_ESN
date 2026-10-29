import unittest
import numpy as np
import pandas as pd
from .vix_stratified_sampling import prefix_mask
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.stratified_sampling import StratifiedGaussianMixture


class AllocationPrefixTests(unittest.TestCase):
    def test_stratified_driver_preserves_all_component_prefixes(self):
        proposal=DefensiveGaussianMixture(np.zeros((7,2)),[.5]+[1/12]*6)
        sampler=StratifiedGaussianMixture(proposal)
        draw=sampler.draw(7,123)
        frame=pd.DataFrame(dict(method='stratified',point=np.arange(128),
            component=draw.components,component_point=draw.component_points))
        mask=prefix_mask(frame,6)
        np.testing.assert_array_equal(mask,sampler.prefix_mask(draw,6))
        np.testing.assert_array_equal(np.bincount(draw.components[mask]),[32,6,6,5,5,5,5])
        self.assertFalse(np.array_equal(np.flatnonzero(mask),np.arange(64)))

    def test_random_mixture_uses_global_prefix(self):
        frame=pd.DataFrame(dict(method='random_mixture',point=np.arange(128),
            component=np.zeros(128,int),component_point=np.full(128,-1)))
        np.testing.assert_array_equal(np.flatnonzero(prefix_mask(frame,6)),np.arange(64))
