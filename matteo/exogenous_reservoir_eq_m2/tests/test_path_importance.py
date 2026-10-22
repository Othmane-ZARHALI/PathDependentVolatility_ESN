"""Independent Gaussian and Black oracles for full-path likelihood weighting."""
import unittest
import numpy as np
from scipy.special import ndtr
from esn_eq.path_importance import prepare_importance_paths
from esn_eq.pricing import coarsen_paths, ForwardVarianceCurve
from test_pricing import small_pricer, constant_readout


class PathImportanceTests(unittest.TestCase):
    @staticmethod
    def assert_mc(test, draws, exact):
        error = np.abs(draws.mean(axis=0)-exact)
        se = draws.std(axis=0, ddof=1)/np.sqrt(len(draws))
        test.assertTrue(np.all(error < 6*se+1e-12), (error, se))

    def test_weighted_joint_gaussian_moments_and_black_tails(self):
        p = small_pricer()
        initial = np.linspace(-.2, .3, len(p.layout.rates))
        t = .3
        x = prepare_importance_paths(p, [0., .1, t], initial, 160000, 412,
                                     checkpoints=[t])
        w, z = x.weights, x.score
        self.assertLessEqual(w.max(), 2.)
        self.assert_mc(self, np.c_[w, w*z, w*z*z], [1., 0., 1.])
        # Joint OU/equity oracle, independently derived from exact transitions.
        state = x.cache.checkpoint_states[t]
        decay = np.exp(-p.layout.rates*t)
        mean = decay*initial+p.q_drift*(-np.expm1(-p.layout.rates*t))/p.layout.rates
        centered = state-mean
        self.assert_mc(self, w[:, None]*centered, np.zeros(len(initial)))
        cross = p.layout.driver_loadings[:, 0]*np.sqrt(2*p.layout.rates)*(-np.expm1(-p.layout.rates*t))/p.layout.rates/np.sqrt(t)
        self.assert_mc(self, w[:, None]*centered*z[:, None], cross)
        # Constant-volatility pricer, including rare calls and puts.
        curve = ForwardVarianceCurve((t,), (.04,))
        stock = p.spx(x.cache, constant_readout(), curve, [t], order=4).terminal_ratios[t]
        self.assert_mc(self, w*stock, 1.)
        for k, kind in ((1., 'call'), (1.5, 'call'), (.65, 'put')):
            sigma = np.sqrt(.04*t)
            d1 = (-np.log(k)+.5*sigma*sigma)/sigma
            d2 = d1-sigma
            exact = ndtr(d1)-k*ndtr(d2) if kind == 'call' else k*ndtr(-d2)-ndtr(-d1)
            payoff = np.maximum((1 if kind == 'call' else -1)*(stock-k), 0.)
            self.assertGreater(np.count_nonzero(payoff), 30)
            self.assert_mc(self, w*payoff, exact)
        # Full-path weight survives coarsening; score and constant-vol stock agree.
        coarse = coarsen_paths(x.cache, [0., t])
        np.testing.assert_allclose(coarse.equity_increments.sum(axis=0)/np.sqrt(t), z, atol=1e-14)
        np.testing.assert_allclose(p.spx(coarse, constant_readout(), curve, [t], order=4).terminal_ratios[t], stock, atol=1e-14)
        # Earlier-expiry expectation also valid with full-horizon weight.
        early = np.exp(-.5*.04*.1+np.sqrt(.04)*x.cache.equity_increments[0])
        self.assert_mc(self, w*early, 1.)

    def test_shift_coupling_and_zero_shift_reproduce_original(self):
        p = small_pricer()
        initial = np.zeros(len(p.layout.rates))
        plain = p.prepare([0., .1], initial, 8, 7, checkpoints=[.1])
        zero = p.prepare([0., .1], initial, 8, 7, checkpoints=[.1], equity_drift=np.zeros(8))
        np.testing.assert_array_equal(plain.equity_increments, zero.equity_increments)
        np.testing.assert_array_equal(plain.checkpoint_states[.1], zero.checkpoint_states[.1])
        shifted = p.prepare([0., .1], initial, 8, 7, checkpoints=[.1], equity_drift=np.full(8, 2.))
        cross = p.layout.driver_loadings[:, 0]*np.sqrt(2*p.layout.rates)*(-np.expm1(-p.layout.rates*.1))/p.layout.rates
        np.testing.assert_allclose(shifted.checkpoint_states[.1]-plain.checkpoint_states[.1], np.broadcast_to(2*cross, (8, len(cross))), atol=1e-14)
        with self.assertRaises(ValueError):
            p.prepare([0., .1], initial, 8, 7, equity_drift=[1.])
        with self.assertRaises(ValueError):
            prepare_importance_paths(p, [0., .1], initial, 8, 7, magnitude=0)
