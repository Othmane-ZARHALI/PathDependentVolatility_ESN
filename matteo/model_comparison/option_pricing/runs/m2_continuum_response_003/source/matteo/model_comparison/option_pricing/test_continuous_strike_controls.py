import unittest
import numpy as np
from .continuous_strike_controls import correction_draws,affine_correction_draws,interpolation_weights,warm_start

class ContinuousControlTests(unittest.TestCase):
 def test_exact_known_zero_mean_preserved(self):
  draws=np.array([[1.,3.,-2.],[-1.,-3.,2.]])
  got=correction_draws(np.linspace(.5,4,100),[1,2,3],draws)
  np.testing.assert_allclose(got.mean(0),0,atol=1e-15)
 def test_continuity_across_anchor(self):
  draws=np.array([[2.,-3.,7.],[1.,8.,-4.]])
  for k in [1.,2.,3.]:
   y=correction_draws([k-1e-9,k,k+1e-9],[1,2,3],draws)
   self.assertLess(np.max(np.abs(np.diff(y,axis=1))),2e-8)
 def test_partition_and_frozen_anchor_replay(self):
  k=np.array([1.,2.,3.]);np.testing.assert_allclose(interpolation_weights(k,k),np.eye(3))
  w=interpolation_weights([.5,1.5,2.5,4.],k);np.testing.assert_allclose(w.sum(1),1);self.assertTrue((w>=0).all())
 def test_affine_known_mean_and_convexity(self):
  draws=np.array([[1.,3.,-2.],[-1.,-3.,2.]])
  strikes=np.linspace(.5,4,100);c=affine_correction_draws(strikes,[1,2,3],draws)
  np.testing.assert_allclose(c.mean(0),0,atol=1e-15)
  np.testing.assert_allclose(np.diff(c,n=2,axis=1),0,atol=3e-15)
  convex=np.maximum(2-strikes,0)[None,:]-c
  self.assertGreaterEqual(np.diff(convex,n=2,axis=1).min(),-4e-15)
 def test_start_only_changes_beta_j(self):
  x=np.array([.2,0.,.06]);got=warm_start(x,[.1,.01,.01],[1,2,.2],['cluster_loading','spike_shift','xi_1']);np.testing.assert_array_equal(got,[.2,.01,.06])
 def test_unrelated_bound_violation_is_not_silently_projected(self):
  with self.assertRaises(ValueError):warm_start([.05,0],[.1,.01],[1,1],['cluster_loading','spike_shift'])
if __name__=='__main__':unittest.main()
