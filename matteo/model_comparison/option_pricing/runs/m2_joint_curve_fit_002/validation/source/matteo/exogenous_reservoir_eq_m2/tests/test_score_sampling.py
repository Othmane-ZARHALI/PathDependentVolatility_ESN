import unittest
import numpy as np
from scipy.special import ndtr
from esn_eq.score_sampling import LowRankInflatedGaussian,score_subspace
from esn_eq.config import ArchitectureConfig,BankConfig,RiskPremium,ReadoutParameters
from esn_eq.pricing import DatedM2Pricer,ForwardVarianceCurve
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.vix_pricing import ConditionalVIXPricer


class ScoreSamplingTests(unittest.TestCase):
    def test_density_ratio_and_unaffected_complement(self):
        U=np.array([[1.],[0.],[0.]])
        p=LowRankInflatedGaussian(U,2.)
        z=np.array([[0.,0.,0.],[3.,0.,0.],[3.,-40.,70.]])
        expected=1/(.5+.25*np.exp(.375*z[:,0]**2))
        np.testing.assert_allclose(p.weights(z),expected,rtol=1e-14)
        self.assertAlmostEqual(p.weights(z)[1],p.weights(z)[2])
        np.testing.assert_allclose(p.rotation@p.rotation.T,np.eye(3),atol=1e-14)

    def test_gaussian_oracles_in_oblique_subspace_and_prefix(self):
        Q=np.linalg.qr(np.random.default_rng(19).normal(size=(6,6)))[0]
        p=LowRankInflatedGaussian(Q[:,:2]);rows=[]
        for seed in range(8):
            d=p.draw(14,seed);z=d.normals;w=d.weights
            h=np.maximum(z[:,0]-3.,0.)
            rows.append(np.r_[w.mean(),(w[:,None]*z).mean(axis=0),(w[:,None]*z*z).mean(axis=0),(w*h).mean()])
            self.assertTrue(np.all(w<=2.+1e-14))
        oracle=np.exp(-4.5)/np.sqrt(2*np.pi)-3*ndtr(-3.)
        np.testing.assert_allclose(np.mean(rows,axis=0)[:13],np.r_[1,np.zeros(6),np.ones(6)],atol=.003)
        self.assertAlmostEqual(np.mean(rows,axis=0)[13],oracle,delta=5e-5)
        small=p.draw(6,47);large=p.draw(7,47)
        np.testing.assert_array_equal(small.normals,large.normals[:64]);np.testing.assert_array_equal(small.weights,large.weights[:64])

    def test_score_geometry_is_deterministic_and_keeps_full_state(self):
        arch=ArchitectureConfig(feedback=BankConfig((1.,8.),.7,True),spike=BankConfig((3.,25.),.4,True),cluster=BankConfig((.2,2.),0.,False))
        d=DatedM2Pricer(arch,RiskPremium(equity=.2));initial=np.zeros(6);params=ReadoutParameters()
        e=ConditionalVIXPricer(d,params,ForwardVarianceCurve((1.,),(.04,)),initial,normalizer_order=6)
        w=e.gaussian_window(.2,time_order=4);sampler=GaussianOuterSampler(d,.2,initial)
        U,values=score_subspace(d,params,sampler,w,rank=3);U2,_=score_subspace(d,params,sampler,w,rank=3)
        np.testing.assert_array_equal(U,U2);np.testing.assert_allclose(U.T@U,np.eye(3),atol=1e-14)
        self.assertTrue(np.all(np.diff(values)<=0));self.assertEqual(LowRankInflatedGaussian(U).draw(4,12).normals.shape,(16,6))
        with self.assertRaises(ValueError):LowRankInflatedGaussian(np.ones((3,2)))
