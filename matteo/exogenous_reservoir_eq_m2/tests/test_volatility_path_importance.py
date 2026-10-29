"""Full joint Gaussian and backward-correspondence oracles for OU path shifts."""
import unittest
import numpy as np
from scipy.special import ndtr
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.path_importance import prepare_importance_paths
from esn_eq.volatility_path_importance import volatility_path_basis, prepare_volatility_paths, GaussianPathTilt
from esn_eq.pricing import ForwardVarianceCurve, coarsen_paths
from esn_eq.states import joint_transition_covariance
from test_pricing import small_pricer, constant_readout


class VolatilityPathTests(unittest.TestCase):
    def assert_mc(self, draws, target):
        delta = np.abs(draws.mean(axis=0)-target)
        se = draws.std(axis=0, ddof=1)/np.sqrt(len(draws))
        self.assertTrue(np.all(delta < 6*se+2e-12), (delta,se))

    def test_joint_score_and_state_oracles_with_asymmetric_volatility_shift(self):
        p = small_pricer()
        times = np.array([0., .03, .08, .2])
        basis = volatility_path_basis(p, times, .08)
        shifts = np.zeros((4,4)); shifts[1,0]=2.; shifts[2,2]=-2.; shifts[3,3]=2.
        mixture = DefensiveGaussianMixture(shifts, [.5,.2,.2,.1])
        initial = np.linspace(-.1,.2,len(p.layout.rates))
        x = prepare_volatility_paths(p,times,initial,120000,527,basis,mixture,checkpoints=[.08,.2])
        w,z = x.weights,x.score
        self.assertLessEqual(w.max(),2.)
        self.assert_mc(w,1.)
        self.assert_mc(w[:,None]*z,np.zeros(4))
        self.assert_mc(w[:,None,None]*z[:,:,None]*z[:,None,:],np.eye(4))
        for t in (.08,.2):
            mean = np.exp(-p.layout.rates*t)*initial+p.q_drift*(-np.expm1(-p.layout.rates*t))/p.layout.rates
            centered = x.cache.checkpoint_states[t]-mean
            self.assert_mc(w[:,None]*centered,np.zeros(len(initial)))
            cov = joint_transition_covariance(p.layout,t)[1:,1:]
            self.assert_mc(w[:,None]*centered**2,np.diag(cov))
            self.assert_mc(w[:,None,None]*centered[:,:,None]*centered[:,None,:],cov)
        curve = ForwardVarianceCurve((.2,),(.04,))
        stock = p.spx(x.cache,constant_readout(),curve,[.08,.2],order=4).terminal_ratios
        for t in (.08,.2):
            self.assert_mc(w*stock[t],1.)
            k=1.1; sigma=np.sqrt(.04*t); d1=(-np.log(k)+.5*sigma**2)/sigma
            exact=ndtr(d1)-k*ndtr(d1-sigma)
            self.assert_mc(w*np.maximum(stock[t]-k,0),exact)
        coarse = coarsen_paths(x.cache,[0.,.08,.2])
        y = p.spx(coarse,constant_readout(),curve,[.08,.2],order=4).terminal_ratios
        np.testing.assert_allclose(y[.2],stock[.2],atol=1e-14)

    def test_equity_only_sampler_correspondence_and_single_use_guard(self):
        p=small_pricer(); times=np.array([0.,.03,.1]); initial=np.zeros(len(p.layout.rates))
        basis=volatility_path_basis(p,times,.03)
        shifts=np.zeros((3,4)); shifts[1,0]=2.; shifts[2,0]=-2.
        mixture=DefensiveGaussianMixture(shifts,[.5,.25,.25])
        old=prepare_importance_paths(p,times,initial,64,29)
        new=prepare_volatility_paths(p,times,initial,64,29,basis,mixture)
        np.testing.assert_array_equal(old.components,new.components)
        np.testing.assert_allclose(old.cache.projections,new.cache.projections,atol=1e-14)
        np.testing.assert_allclose(old.cache.equity_increments,new.cache.equity_increments,atol=1e-14)
        np.testing.assert_allclose(old.weights,new.weights,atol=1e-14)
        tilt=GaussianPathTilt(basis,mixture,64,29)
        with self.assertRaises(ValueError):
            tilt.transform(1,np.zeros((64,len(initial)+1)))
        with self.assertRaises(ValueError):
            GaussianPathTilt(basis*2,mixture,64,29)
        with self.assertRaises(ValueError):
            volatility_path_basis(p,times,.02)
        with self.assertRaises(ValueError):
            p.prepare(times,initial,64,29,equity_drift=np.zeros(64),gaussian_tilt=tilt)
