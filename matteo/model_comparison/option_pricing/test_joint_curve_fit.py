"""Exact repricing bridge, analytical derivatives and unchanged loss weights."""
import unittest
import numpy as np
import pandas as pd
from dataclasses import replace
from .joint_curve_core import CurvePathBatch, path_payoffs, donor_indices, ExploratoryJointLoss
from .joint_curve_fit import VIXData
from .m2_spx import ROOT
from esn_eq.market_targets import MarketTargets
from esn_eq.calibration_loss import JointPriceObjective
from esn_eq.pricing import ForwardVarianceCurve, DatedM2Pricer
from esn_eq.config import ArchitectureConfig, BankConfig, ReadoutParameters, RiskPremium


class JointCurveFitTests(unittest.TestCase):
    def test_exponential_euler_bridge_and_derivative(self):
        p = DatedM2Pricer(ArchitectureConfig(feedback=BankConfig((.5, 10.), .7, True),
            spike=BankConfig((3., 25.), .4, True), cluster=BankConfig((.2, 2.), 0., False)), RiskPremium(equity=.35))
        times = np.array([0., .02, .05, .08, .1]); params = ReadoutParameters(cluster_loading=.1, feedback_curvature=.03, spike_curvature=.02)
        initial = np.zeros(len(p.layout.rates)); cache = p.prepare(times, initial, 128, 400)
        ends = (.02, .08, .1); levels = np.array([.04, .06, .05]); curve = ForwardVarianceCurve(ends, tuple(levels))
        a = p.multiplier(cache.projections, params)/p.normalizer(times[:-1], initial, params, 12)[:, None]
        batch = CurvePathBatch(a, cache.equity_increments, np.diff(times), np.minimum(np.searchsorted(ends, times[:-1], side='right'), 2))
        stock, derivative = batch.terminal(levels, chunk=2)
        direct = p.spx(cache, params, curve, [.1], order=12).terminal_ratios[.1]
        np.testing.assert_allclose(stock, direct, atol=1e-14)
        for j in range(3):
            d = np.eye(3)[j]*1e-6
            plus = batch.terminal(levels+d)[0]; minus = batch.terminal(levels-d)[0]
            np.testing.assert_allclose(derivative[:, j], (plus-minus)/2e-6, rtol=1e-7, atol=2e-9)

    def test_signed_extension_matches_production_on_nonnegative_prices(self):
        t = MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
        quotes = pd.concat([t.select('SPX').groupby('expiry').head(3), t.select('VIX')], ignore_index=True)
        levels = np.array([.04, .06, .05]); prices = (quotes.price_bid_U+quotes.price_ask_U).to_numpy()/2+.03
        f = t.vix_forwards.sort_values('expiry').reset_index(drop=True)
        forwards = f.forward_vix_points.to_numpy()+.1
        expected = JointPriceObjective(quotes, f).evaluate(prices, forwards, curve_penalty=.1*np.sum(np.diff(levels)**2))
        actual = ExploratoryJointLoss(quotes, f).evaluate(prices, forwards, levels)
        self.assertAlmostEqual(expected['total'], actual['total'], places=11)
        for key in expected['components']:
            self.assertAlmostEqual(expected['components'][key], actual['components'][key], places=12)
        prices[0] = -1.
        with self.assertRaises(ValueError):
            JointPriceObjective(quotes, f).evaluate(prices, forwards)
        actual = ExploratoryJointLoss(quotes, f).evaluate(prices, forwards, levels)
        self.assertEqual(actual['negative_count'], 1)
        self.assertGreater(actual['total'], expected['total'])

    def test_loss_and_vix_derivatives(self):
        q = pd.DataFrame(dict(option_type=['put', 'call'], strike=[18., 24.]))
        rng = np.random.default_rng(401); a = rng.uniform(.01, .04, (64, 3)); control = rng.normal(size=(64, 2))
        data = VIXData(q, a, .00001, .08, np.ones(64), control, np.zeros((2, 3)), 8)
        levels = np.array([.04, .06, .05]); mean, jac, _, _ = data.price(levels)
        for j in range(3):
            d = np.eye(3)[j]*1e-6
            np.testing.assert_allclose(jac[:, j], (data.price(levels+d)[0]-data.price(levels-d)[0])/2e-6, rtol=1e-6, atol=1e-7)
        t = MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
        quotes = pd.concat([t.select('SPX').groupby('expiry').head(1), t.select('VIX').groupby('expiry').head(1)], ignore_index=True)
        f = t.vix_forwards.sort_values('expiry').reset_index(drop=True)
        loss = ExploratoryJointLoss(quotes, f)
        prices = quotes.price_bid_U.to_numpy()*.5-1.
        forward = f.forward_vix_points.to_numpy()+.2
        pj = rng.normal(size=(len(quotes), 3)); fj = rng.normal(size=(len(f), 3))
        actual = loss.evaluate(prices, forward, levels, pj, fj)['jacobian']
        for j in range(3):
            delta = 1e-6
            plus = loss.evaluate(prices+delta*pj[:, j], forward+delta*fj[:, j], levels+delta*np.eye(3)[j])['residuals']
            minus = loss.evaluate(prices-delta*pj[:, j], forward-delta*fj[:, j], levels-delta*np.eye(3)[j])['residuals']
            np.testing.assert_allclose(actual[:, j], (plus-minus)/(2*delta), rtol=1e-6, atol=1e-7)

    def test_market_only_hedge_assignment(self):
        q = pd.DataFrame(dict(option_type=['put', 'call'], K_over_F=[.95, 1.1]))
        r = pd.DataFrame(dict(option_type=['put', 'put', 'call', 'call'], K_over_F=[.8, .98, 1.02, 1.2]))
        np.testing.assert_array_equal(donor_indices(q, r), [1, 2])
