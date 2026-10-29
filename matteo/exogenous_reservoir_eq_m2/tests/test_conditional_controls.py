"""Independent Gaussian integration, Q controls and retained-law oracles."""
from dataclasses import replace
import unittest

import numpy as np
from scipy.special import ndtr, roots_hermitenorm

from esn_eq.conditional_controls import stock_delta_mean, retained_path_basis, simulate_conditional_controls
from esn_eq.conditional_spx import lognormal_payoffs
from esn_eq.pricing import DatedM2Pricer, ForwardVarianceCurve
from esn_eq.importance_sampling import DefensiveGaussianMixture
from test_pricing import small_pricer, constant_readout


class ConditionalControlTests(unittest.TestCase):
    def assert_mc(self, x, target=0.):
        self.assertTrue(np.all(abs(x.mean(axis=0)-target) <
            6*x.std(axis=0, ddof=1)/np.sqrt(len(x))+2e-11))

    def test_stock_delta_against_independent_gaussian_quadrature(self):
        z, w = roots_hermitenorm(128)
        w /= np.sqrt(2*np.pi)
        for kind in ['put', 'call']:
            for m, v, k, remaining in [(.1, .03, 1.1, .08), (-.3, .2, .8, .15), (.01, 0., 1., .02)]:
                l = m+np.sqrt(v)*z
                d1 = (l-np.log(k)+.5*remaining)/np.sqrt(remaining)
                actual = np.sum(w*np.exp(l)*(ndtr(d1)-(kind == 'put')))
                self.assertAlmostEqual(float(stock_delta_mean(m, v, k, remaining, kind)), actual, places=12)

    def test_fixed_retained_path_hedge_by_two_dimensional_quadrature(self):
        # At step two the past eliminated equity shock randomizes the delta.
        # Directly integrate both eliminated shocks without using the formula.
        z, w = roots_hermitenorm(100)
        w /= np.sqrt(2*np.pi)
        m, v, d, q, k, remaining = -.03, .025, .012, .008, 1.04, .04
        l = m+np.sqrt(v)*z[:, None]
        following = l+d+np.sqrt(q)*z[None, :]
        for kind in ['put', 'call']:
            delta = ndtr((l-np.log(k)+.5*remaining)/np.sqrt(remaining))-(kind == 'put')
            actual = np.sum(w[:, None]*w[None, :]*delta*(np.exp(following)-np.exp(l)))
            expected = np.expm1(d+.5*q)*stock_delta_mean(m, v, k, remaining, kind)
            self.assertAlmostEqual(actual, float(expected), places=13)
            self.assertGreater(abs(expected), 1e-4)  # E[gain|U] is not zero.

    def test_independent_constant_volatility_black_and_zero_hedges(self):
        p = small_pricer()
        arch = replace(p.architecture,
            feedback=replace(p.architecture.feedback, asset_correlation=0.),
            spike=replace(p.architecture.spike, asset_correlation=0.))
        p = DatedM2Pricer(arch, p.premium)
        times = np.array([0., .02, .05, .1])
        result = simulate_conditional_controls(p, times, np.zeros(len(p.layout.rates)),
            constant_readout(), ForwardVarianceCurve((.05, .1), (.04, .09)),
            [.05, .1], [.98, 1.04], ['put', 'call'], 128, 310, order=4)
        total = np.array([.002, .0065])
        expected = lognormal_payoffs(-.5*total, total, [.98, 1.04], ['put', 'call'])
        np.testing.assert_allclose(result.payoff, np.broadcast_to(expected, result.payoff.shape), atol=1e-14)
        np.testing.assert_allclose(result.controls, 0., atol=1e-14)
        np.testing.assert_allclose(result.stock, 1., atol=1e-14)

    def test_correlated_nonlinear_raw_conditioning_and_forward_controls(self):
        p = small_pricer()
        times = np.array([0., .02, .05, .1])
        params = replace(constant_readout(), cluster_loading=.15, feedback_curvature=.03, spike_curvature=.02)
        result = simulate_conditional_controls(p, times, np.zeros(len(p.layout.rates)),
            params, ForwardVarianceCurve((.1,), (.04,)), [.05, .1], [.97, 1.04],
            ['put', 'call'], 50000, 311, order=12, include_raw=True)
        self.assert_mc(result.raw_payoff-result.payoff)
        self.assert_mc(result.raw_controls-result.controls)
        self.assert_mc(result.controls)
        self.assert_mc(result.stock, 1.)
        self.assert_mc(result.reference_stock, 1.)
        # Fixed coefficients isolate conditioning from independently trained fits.
        beta = np.array([[.1, .4, .3, .2], [.2, .4, .2, .3]])
        raw = result.raw_payoff-np.einsum('njk,jk->nj', result.raw_controls, beta)
        conditional = result.payoff-np.einsum('njk,jk->nj', result.controls, beta)
        self.assert_mc(raw-conditional)
        self.assertTrue(np.all(conditional.var(axis=0) < 1.03*raw.var(axis=0)))

    def test_retained_importance_law_mass_black_and_controls(self):
        p = small_pricer()
        times = np.array([0., .02, .05, .1])
        basis = retained_path_basis(p, times, .02)
        flat = basis.reshape(-1, 4)
        np.testing.assert_allclose(flat.T@flat, np.eye(4), atol=1e-14)
        mixture = DefensiveGaussianMixture([[0., 0., 0., 0.], [1.5, 0., 0., 0.], [0., -1.5, 0., 0.]], [.5, .25, .25])
        r = simulate_conditional_controls(p, times, np.zeros(len(p.layout.rates)),
            constant_readout(), ForwardVarianceCurve((.1,), (.04,)), [.1, .1], [.98, 1.04],
            ['put', 'call'], 60000, 312, order=4, basis=basis, mixture=mixture, include_raw=True)
        self.assert_mc(r.weights, 1.)
        self.assert_mc(r.weights[:, None]*r.scores)
        self.assert_mc(r.weights[:, None]*r.scores**2, 1.)
        self.assertLessEqual(r.weights.max(), 2.)
        expected = lognormal_payoffs(-.002, .004, [.98, 1.04], ['put', 'call'])
        self.assert_mc(r.weights[:, None]*r.payoff, expected)
        self.assert_mc(r.weights[:, None, None]*r.controls)
        self.assert_mc(r.weights[:, None]*r.stock, 1.)
        self.assert_mc(r.weights[:, None]*(r.raw_payoff-r.payoff))
        self.assert_mc(r.weights[:, None, None]*(r.raw_controls-r.controls))

    def test_perfect_equity_driver_recovers_raw_controls(self):
        p = small_pricer()
        arch = replace(p.architecture, feedback=replace(p.architecture.feedback, asset_correlation=1.))
        p = DatedM2Pricer(arch, p.premium)
        r = simulate_conditional_controls(p, [0., .02, .1], np.zeros(len(p.layout.rates)),
            constant_readout(), ForwardVarianceCurve((.1,), (.04,)), [.1, .1], [.98, 1.04],
            ['put', 'call'], 64, 313, order=4, include_raw=True)
        np.testing.assert_allclose(r.raw_payoff, r.payoff, atol=1e-14)
        np.testing.assert_allclose(r.raw_controls, r.controls, atol=1e-14)

    def test_invalid_proposals_and_delta_inputs(self):
        with self.assertRaises(ValueError):
            stock_delta_mean(0., -.1, 1., .1, 'call')
        with self.assertRaises(ValueError):
            retained_path_basis(small_pricer(), [0., .1], .02)
