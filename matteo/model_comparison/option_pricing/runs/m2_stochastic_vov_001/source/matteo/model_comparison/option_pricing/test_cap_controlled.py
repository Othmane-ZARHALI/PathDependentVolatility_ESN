"""Meaningful constraints for the clean cap experiment; standard-library runner."""
import unittest
import json
import numpy as np
from .full_readout_fit import HERE
from .cap_controlled_core import window_basis,adjusted,curve_jac

class WindowTests(unittest.TestCase):
    def setUp(self):
        d=json.loads((HERE/'runs/m2_maturity_curve_fit_001/run.json').read_text());self.p=np.array(d['selected']);self.B=window_basis(d['ends_days'])
    def test_known_windows_and_overlap(self):
        np.testing.assert_array_equal(np.flatnonzero(self.B.sum(1)),[1,2,3,4,7,8,10])
        np.testing.assert_array_equal(self.B[3],[.5,.5,0,0]);np.testing.assert_array_equal(self.B[10],[0,0,0,1])
        np.testing.assert_array_equal(adjusted(self.p,np.zeros(4),self.B),self.p)
    def test_all_corners_preserve_frozen_parameters_and_floor(self):
        for bits in range(16):
            s=np.array([.15 if bits&(1<<j) else -.15 for j in range(4)]);x=adjusted(self.p,s,self.B);off=self.B.sum(1)==0
            np.testing.assert_array_equal(x[7:-1][off],self.p[7:-1][off]);np.testing.assert_array_equal(x[:7],self.p[:7]);self.assertEqual(x[-1],self.p[-1]);self.assertGreater(x[4],0)
            self.assertTrue(np.all(x[7:-1]>=x[-1]));ratio=(x[7:-1]-x[-1])/(self.p[7:-1]-x[-1]);self.assertTrue(np.all(ratio>=.85-1e-14));self.assertTrue(np.all(ratio<=1.15+1e-14))
    def test_exact_curve_sensitivities(self):
        D=curve_jac(self.p,self.B)
        for j in range(4):
            s=np.zeros(4);s[j]=1e-6;a=adjusted(self.p,s,self.B);s[j]=-1e-6;b=adjusted(self.p,s,self.B);np.testing.assert_allclose((a-b)/2e-6,D[:,j],atol=5e-12,rtol=1e-8)
        with self.assertRaises(ValueError):adjusted(self.p,np.full(4,.16),self.B)

if __name__=='__main__':unittest.main()
