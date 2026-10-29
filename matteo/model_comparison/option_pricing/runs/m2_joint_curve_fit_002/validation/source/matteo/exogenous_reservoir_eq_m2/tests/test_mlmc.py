"""Allocation optimum, independent variance addition and analytic target oracles."""
import unittest
import numpy as np
from esn_eq.mlmc import allocate_levels,independent_level_statistics


class MLMCTests(unittest.TestCase):
    def test_allocation_respects_budget_and_continuous_optimum(self):
        # sqrt(V/C)=(8,2,1), exact continuous solution is integral in blocks.
        v=np.array([64.,16.,16.]);c=np.array([1.,4.,16.]);budget=2048.
        n=allocate_levels(v,c,budget,64,64)
        np.testing.assert_array_equal(n,[512,128,64])
        self.assertLessEqual(n@c,budget)
        best=sum(v/n)
        for a in range(64,1025,64):
            for b in range(64,513,64):
                for d in range(64,129,64):
                    if np.array([a,b,d])@c<=budget:
                        self.assertLessEqual(best,sum(v/np.array([a,b,d]))+1e-14)

    def test_minimum_zero_variance_and_invalid_budgets(self):
        n=allocate_levels([100.,.001,0.],[1.,2.,4.],1024,64,64)
        self.assertEqual(n[-1],64);self.assertEqual(n[1],64)
        self.assertLessEqual(n@np.array([1.,2.,4.]),1024)
        np.testing.assert_array_equal(allocate_levels([0.,0.],[1.,2.],1024),[64,64])
        with self.assertRaises(ValueError):allocate_levels([1.,1.],[1.,2.],100)
        with self.assertRaises(ValueError):allocate_levels([1.,-1.],[1.,2.],1000)

    def test_independent_mean_and_variance_with_unequal_counts(self):
        base=np.array([[1.,2.],[3.,4.],[5.,6.]])
        correction=np.array([[-2.,1.],[0.,3.]])
        mean,se=independent_level_statistics([base,correction])
        np.testing.assert_allclose(mean,[2.,6.])
        np.testing.assert_allclose(se,np.sqrt([4/3+1.,4/3+1.]))
        with self.assertRaises(ValueError):independent_level_statistics([np.array([1.])])
