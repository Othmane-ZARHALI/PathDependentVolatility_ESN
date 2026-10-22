import unittest
import numpy as np
from esn_eq.gaussian_controls import gaussian_control_rows,fit_replication_controls
from esn_eq.importance_sampling import DefensiveGaussianMixture


class GaussianControlTests(unittest.TestCase):
    def test_exact_polynomial_oracle_with_importance_weights(self):
        U=np.linalg.qr(np.random.default_rng(12).normal(size=(6,6)))[0][:,:2]
        p=DefensiveGaussianMixture(np.r_[np.zeros((1,6)),(2*U[:,0])[None,:]],[.5,.5])
        z,w,_=p.draw(8,14);y=z@U;C=gaussian_control_rows(z,w,U)
        h=2+y[:,0]+.3*(y[:,0]**2-1)+.2*y[:,1]
        beta=np.array([2,1,.2,.3,0.])
        np.testing.assert_allclose(w*h-C@beta,2,rtol=1e-13,atol=1e-13)
        # No clipping: arbitrary fixed coefficients may produce negative samples.
        self.assertTrue(np.any(w*h-C@(beta*100)<0))

    def test_gaussian_means_under_shifted_sampling(self):
        U=np.linalg.qr(np.random.default_rng(5).normal(size=(5,5)))[0][:,:3]
        shifts=np.r_[np.zeros((1,5)),(2*U[:,0])[None,:],(-2*U[:,0])[None,:]]
        p=DefensiveGaussianMixture(shifts,[.5,.25,.25]);means=[]
        for seed in range(8):
            z,w,_=p.draw(14,seed);means.append(gaussian_control_rows(z,w,U).mean(axis=0))
        np.testing.assert_allclose(np.mean(means,axis=0),0,atol=.002)

    def test_replication_ridge_and_constant_controls(self):
        x=np.arange(8.)[:,None];C=np.column_stack((x,np.ones(8)))
        Y=4+3*x;beta=fit_replication_controls(Y,C,ridge=1.)
        np.testing.assert_allclose(beta[:,0],[1.5,0],rtol=1e-14,atol=1e-14)
        self.assertAlmostEqual(float(fit_replication_controls(Y*7,C)[0,0]),10.5)
        for bad in (0.,-1.,np.nan):
            with self.assertRaises(ValueError):fit_replication_controls(Y,C,bad)
        with self.assertRaises(ValueError):gaussian_control_rows(np.zeros((2,3)),np.ones(2),np.ones((3,2)))
