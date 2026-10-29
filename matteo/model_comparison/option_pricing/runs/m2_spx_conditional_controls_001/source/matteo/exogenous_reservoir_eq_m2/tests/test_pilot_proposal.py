import unittest
import numpy as np
from esn_eq.pilot_proposal import fit_payoff_mixture
from esn_eq.importance_sampling import DefensiveGaussianMixture


class PilotProposalTests(unittest.TestCase):
    def test_symmetric_tail_detects_direction_despite_zero_mean(self):
        rng=np.random.default_rng(71)
        z=rng.standard_normal((20000,4)); z=np.r_[z,-z]
        h=np.maximum(z[:,0]**2-4,0)[:,None]
        proposal,records=fit_payoff_mixture(z,np.ones(len(z)),h,coordinates=2)
        self.assertGreater(abs(proposal.shifts[1,0]),1.)
        np.testing.assert_allclose(proposal.shifts[1],-proposal.shifts[2],atol=1e-12)
        self.assertTrue(np.all(proposal.shifts[:,2:]==0))
        self.assertGreater(records[0]['positive_pilot_count'],100)

    def test_payoff_scaling_invariance_and_shift_cap(self):
        z=np.array([[0.,0.],[10.,0.],[11.,1.],[12.,0.]])
        w=np.array([1.,.5,.5,.25]); h=np.array([[0.],[1.],[2.],[3.]])
        a,_=fit_payoff_mixture(z,w,h,max_shift=2.)
        b,_=fit_payoff_mixture(z,w,h*123,max_shift=2.)
        np.testing.assert_allclose(a.shifts,b.shifts,atol=1e-14)
        self.assertTrue(np.all(np.linalg.norm(a.shifts,axis=1)<=2.+1e-14))
        self.assertEqual(a.probabilities[0],.5)
        # Learned shifts are still corrected with the ordinary exact density ratio.
        zz,weights,_=a.draw(12,12)
        self.assertTrue(np.all(weights<=2.+1e-14))

    def test_zero_mass_and_invalid_weights_rejected(self):
        for w,h in (([1.,1.],[[0.],[0.]]),([1.,-1.],[[1.],[2.]])):
            with self.assertRaises(ValueError):
                fit_payoff_mixture(np.zeros((2,3)),w,h)
