import unittest
import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
from .vix_targeted_pilots import blend,fit_candidate,second_moment_contributions
from esn_eq.importance_sampling import DefensiveGaussianMixture


class TargetedPilotTests(unittest.TestCase):
    def test_defensive_blend_has_exact_mixture_density(self):
        a=DefensiveGaussianMixture([[0.,0.],[2.,0.]],[.5,.5])
        b=DefensiveGaussianMixture([[0.,0.],[0.,-3.]],[.5,.5])
        c=blend(a,b);z=np.array([[0.,0.],[2.,-3.],[-1.,4.]])
        np.testing.assert_allclose(1/c.weights(z),.5/a.weights(z)+.5/b.weights(z),rtol=1e-14)
        self.assertEqual(c.probabilities[0],.5)
        self.assertTrue(np.all(c.weights(z)<=2))

    def test_second_moment_proxy_corrects_both_source_and_candidate(self):
        source=DefensiveGaussianMixture([[0.],[1.5]],[.5,.5])
        target=DefensiveGaussianMixture([[0.]],[1.])
        x,q=hermgauss(96);z=np.r_[np.sqrt(2)*x,np.sqrt(2)*x+1.5][:,None]
        mass=np.r_[q,q]/(2*np.sqrt(np.pi));h=np.exp(.3*z)
        value=second_moment_contributions(z,source.weights(z),h,target)[:,0]
        self.assertAlmostEqual(float(mass@value),np.exp(.18),places=12)
        # The quantity is a second moment, not variance (which subtracts E[h]^2).
        self.assertGreater(float(mass@value),1.)

    def test_score_lift_retains_orthogonal_state_and_target_scale_invariance(self):
        rng=np.random.default_rng(91);z=rng.normal(size=(1024,6));w=np.ones(len(z));v=20+3*z[:,0]**2
        q=pd.DataFrame(dict(option_type=['put','call','call'],strike=[21.,21.,30.],price_scale_U=[.1,.1,.1]))
        base=DefensiveGaussianMixture(np.zeros((1,6)),[1.]);U=np.eye(6)[:,:4]
        a,diag=fit_candidate(z,w,v,q,U,base,'score4_balanced',4.)
        self.assertEqual(len(diag),5)
        self.assertTrue(np.all(a.shifts[:,4:]==0))
        self.assertTrue(np.all(np.linalg.norm(a.shifts,axis=1)<=4.+1e-12))
        self.assertEqual(a.draw(4,3)[0].shape,(16,6))
        # Scaling quote-family targets cannot alter their weighted law.
        q.price_scale_U*=100
        b,_=fit_candidate(z,w,v,q,U,base,'score4_balanced',4.)
        np.testing.assert_allclose(a.shifts,b.shifts,atol=1e-12)
