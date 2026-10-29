import unittest
import numpy as np
from scipy.stats import qmc
from .vix_composed_moment import discrete_gaussian,composed_integrands
from . import test_vix_collapsed_moment as fixtures
from esn_eq.config import ReadoutParameters
from esn_eq.pricing import ForwardVarianceCurve


class ComposedMomentTests(unittest.TestCase):
    def test_discrete_rule_recovers_gaussian_polynomial_moments(self):
        u=qmc.Sobol(3,scramble=True,bits=52,rng=31).random_base2(18)
        z,mass=discrete_gaussian(u,8)
        self.assertAlmostEqual(mass,1.,places=14)
        np.testing.assert_allclose(z.mean(axis=0),0.,atol=1e-4)
        np.testing.assert_allclose((z*z).mean(axis=0),1.,atol=3e-4)
        np.testing.assert_allclose((z**4).mean(axis=0),3.,atol=.003)

    def test_linear_exponential_readout_matches_known_gaussian_expectation(self):
        dated,initial=fixtures.CollapsedMomentTests().fixture()
        parameters=ReadoutParameters(cluster_loading=.15,feedback_curvature=0.,spike_curvature=0.,cap_level=30.)
        curve=ForwardVarianceCurve((1.,),(.04,))
        times=np.array([.21,.24,.27]);weights=np.full(3,(30/365)/3)
        denominators=[]
        for t in times:
            mean,cov=dated.gaussian_projection_law(t,initial)
            denominators.append(np.exp(.3*mean[0]+.5*.3**2*cov[0,0]))
        u=qmc.Sobol(6,scramble=True,bits=52,rng=14).random_base2(15)
        values=composed_integrands(dated,parameters,initial,curve,.2,times,weights,np.array([denominators]),u,orders=(8,))
        np.testing.assert_allclose(values.mean(axis=0),400.,atol=.08)

    def test_inflated_discrete_proposal_recovers_exact_rule_moments(self):
        rows=[]
        for seed in range(4):
            u=qmc.Sobol(4,scramble=True,bits=52,rng=seed+30).random_base2(16)
            z,mass,log_ratio=discrete_gaussian(u[:,1:],8,2.5,u[:,0]>.5)
            weight=np.exp(-np.logaddexp(-np.log(2.),-np.log(2.)+log_ratio))
            rows.append(np.r_[weight.mean(),(weight[:,None]*z).mean(axis=0),(weight[:,None]*z*z).mean(axis=0)])
            self.assertLessEqual(weight.max(),2.)
        np.testing.assert_allclose(np.mean(rows,axis=0),[1,0,0,0,1,1,1],atol=.002)

    def test_inflated_joint_proposal_preserves_linear_gaussian_oracle(self):
        dated,initial=fixtures.CollapsedMomentTests().fixture()
        parameters=ReadoutParameters(cluster_loading=.15,feedback_curvature=0.,spike_curvature=0.,cap_level=30.)
        curve=ForwardVarianceCurve((1.,),(.04,));times=np.array([.21,.24,.27]);weights=np.full(3,(30/365)/3)
        denominators=[]
        for t in times:
            mean,cov=dated.gaussian_projection_law(t,initial)
            denominators.append(np.exp(.3*mean[0]+.5*.3**2*cov[0,0]))
        results=[]
        for seed in range(4):
            u=qmc.Sobol(7,scramble=True,bits=52,rng=seed+70).random_base2(16)
            results.append(composed_integrands(dated,parameters,initial,curve,.2,times,weights,np.array([denominators]),u,orders=(8,),proposal_scale=2.5).mean(axis=0))
        np.testing.assert_allclose(np.mean(results,axis=0),400.,atol=.5)
