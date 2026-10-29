"""Dated-pricing correspondence and analytic-oracle tests."""
from dataclasses import replace
import unittest

import numpy as np
from scipy.stats import norm

from esn_eq.config import ArchitectureConfig, BankConfig, ReadoutParameters, RiskPremium
from esn_eq.pricing import DatedM2Pricer, ForwardVarianceCurve, coarsen_paths, dated_grid, payoff_statistics
from esn_eq.states import joint_transition_covariance


def small_pricer():
    return DatedM2Pricer(ArchitectureConfig(
        feedback=BankConfig((0.5, 10.), 0.7, True),
        spike=BankConfig((3., 25.), 0.4, True),
        cluster=BankConfig((0.2, 2.), 0., False)), RiskPremium(equity=0.35))


def constant_readout():
    return ReadoutParameters(cluster_loading=0, feedback_curvature=0, spike_curvature=0)


class DatedPricingTests(unittest.TestCase):
    def test_event_grid_preserves_exact_expiry_and_curve_jumps(self):
        maturity = 29 / 365
        times = dated_grid([maturity, 0.03], 1 / 504)
        self.assertIn(maturity, times)
        self.assertIn(0.03, times)
        self.assertLessEqual(np.max(np.diff(times)), 1 / 504 + 1e-15)
        self.assertAlmostEqual(np.diff(times).sum(), maturity, places=15)

    def test_curve_integral_and_right_continuity(self):
        curve = ForwardVarianceCurve((0.1, 0.3), (0.04, 0.09))
        self.assertEqual(float(curve(0.1)), 0.09)
        self.assertAlmostEqual(curve.integral(0.05, 0.2), 0.05 * 0.04 + 0.1 * 0.09)
        with self.assertRaises(ValueError):
            curve(0.31)

    def test_projection_law_matches_full_exact_transition_and_risk_shift(self):
        p = small_pricer()
        initial = np.linspace(-0.3, 0.4, len(p.layout.rates))
        t = 0.17
        mean, covariance = p.gaussian_projection_law(t, initial)
        exact = joint_transition_covariance(p.layout, t)[1:, 1:]
        np.testing.assert_allclose(covariance, p.projection @ exact @ p.projection.T, atol=1e-14)
        # Integrate constant Q drift analytically, using a nonzero full initial state.
        full_mean = np.exp(-p.layout.rates * t) * initial + p.q_drift * (1 - np.exp(-p.layout.rates * t)) / p.layout.rates
        np.testing.assert_allclose(mean, p.projection @ full_mean, atol=1e-14)
        self.assertLess(p.q_drift[0], 0)

    def test_projected_initial_state_recovers_requested_coordinates(self):
        p = small_pricer()
        requested = np.array([0.3, -0.4, 0.2])
        full, penalty = p.initial_state(requested)
        np.testing.assert_allclose(p.projection @ full, requested, atol=1e-14)
        self.assertGreater(penalty, 0)
        np.testing.assert_array_equal(p.initial_state(np.zeros(3))[0], 0)

    def test_gaussian_normalizer_against_independent_monte_carlo(self):
        p = small_pricer()
        initial, _ = p.initial_state([0.1, -0.2, 0.1])
        parameters = ReadoutParameters(cluster_loading=0.15, feedback_curvature=0.04,
                                      feedback_shift=0.5, spike_curvature=0.03)
        mean, cov = p.gaussian_projection_law(0.4, initial)
        value = float(p.normalizer([0.4], initial, parameters, order=20)[0])
        draws = np.random.default_rng(49).multivariate_normal(mean, cov, size=50000)
        multipliers = p.multiplier(draws, parameters)
        se = multipliers.std(ddof=1) / np.sqrt(len(multipliers))
        self.assertLess(abs(value - multipliers.mean()), 5 * se)

    def test_zero_time_normalizer_is_the_deterministic_multiplier(self):
        p = small_pricer()
        initial, _ = p.initial_state([0.2, -0.5, 0.1])
        parameters = ReadoutParameters()
        np.testing.assert_allclose(p.normalizer([0], initial, parameters, 12),
                                   p.multiplier(p.projection @ initial, parameters))

    def test_constant_volatility_black_prices_and_stock_forward(self):
        p = small_pricer()
        maturity = 0.5
        times = dated_grid([maturity, 0.2], 0.1)
        cache = p.prepare(times, np.zeros(len(p.layout.rates)), 32768, seed=104)
        curve = ForwardVarianceCurve((0.2, maturity), (0.04, 0.09))
        result = p.spx(cache, constant_readout(), curve, [maturity], order=8)
        y = result.terminal_ratios[maturity]
        prices, errors = payoff_statistics(y, [0.9, 1.0, 1.1], ["put", "call", "call"])
        variance = curve.integral(0, maturity)
        for strike, kind, price, error in zip([0.9, 1., 1.1], ["put", "call", "call"], prices, errors):
            d1 = (-np.log(strike) + 0.5 * variance) / np.sqrt(variance)
            d2 = d1 - np.sqrt(variance)
            oracle = (norm.cdf(d1) - strike * norm.cdf(d2) if kind == "call"
                      else strike * norm.cdf(-d2) - norm.cdf(-d1))
            self.assertLess(abs(price - oracle), 5 * error)
        self.assertLess(abs(y.mean() - 1), 5 * y.std(ddof=1) / np.sqrt(len(y)))
        np.testing.assert_allclose(result.variance_mean, curve(times[:-1]), rtol=0, atol=1e-12)
        # Empirical parity is exact on shared paths, even when sample mean differs from 1.
        both, _ = payoff_statistics(y, [1., 1.], ["call", "put"])
        self.assertAlmostEqual(both[0] - both[1], y.mean() - 1, places=14)

    def test_exact_state_distribution_and_checkpoint(self):
        p = small_pricer()
        initial = np.zeros(len(p.layout.rates))
        cache = p.prepare([0, 0.4], initial, 16384, seed=50, checkpoints=[0.4])
        draws = cache.checkpoint_states[0.4] @ p.projection.T
        mean, cov = p.gaussian_projection_law(0.4, initial)
        np.testing.assert_allclose(draws.mean(axis=0), mean, atol=0.03)
        np.testing.assert_allclose(np.cov(draws.T), cov, atol=0.03)

    def test_readout_reuse_and_nonconstant_martingale(self):
        p = small_pricer()
        cache = p.prepare(dated_grid([0.1], 0.01), np.zeros(len(p.layout.rates)), 8192, seed=12)
        original = cache.projections.copy()
        parameters = ReadoutParameters(cluster_loading=0.2, feedback_curvature=0.05,
                                      feedback_shift=0.8, spike_curvature=0.03)
        result = p.spx(cache, parameters, ForwardVarianceCurve((0.1,), (0.04,)), [0.1], order=12)
        y = result.terminal_ratios[0.1]
        self.assertLess(abs(y.mean() - 1), 5 * y.std(ddof=1) / np.sqrt(len(y)))
        np.testing.assert_array_equal(cache.projections, original)
        other = p.spx(cache, replace(parameters, feedback_shift=1.2),
                      ForwardVarianceCurve((0.1,), (0.04,)), [0.1], order=12)
        self.assertFalse(np.array_equal(y, other.terminal_ratios[0.1]))

    def test_reject_incompatible_cache_and_missing_curve_event(self):
        p = small_pricer()
        cache = p.prepare([0., 0.1], np.zeros(len(p.layout.rates)), 2, seed=12)
        with self.assertRaises(ValueError):
            p.spx(cache, constant_readout(), ForwardVarianceCurve((0.05, 0.1), (0.04, 0.05)), [0.1])
        other = DatedM2Pricer(p.architecture, RiskPremium(equity=0.1))
        with self.assertRaises(ValueError):
            other.spx(cache, constant_readout(), ForwardVarianceCurve((0.1,), (0.04,)), [0.1])
        with self.assertRaises(ValueError):
            p.spx(cache, ReadoutParameters(echo_loading=0.2), ForwardVarianceCurve((0.1,), (0.04,)), [0.1])

    def test_coupled_coarsening_preserves_constant_volatility_stock_paths(self):
        p = small_pricer()
        fine = p.prepare([0., 0.01, 0.02, 0.03, 0.04], np.zeros(len(p.layout.rates)), 64, seed=9)
        coarse = coarsen_paths(fine, [0., 0.02, 0.04])
        np.testing.assert_allclose(coarse.equity_increments.sum(axis=0), fine.equity_increments.sum(axis=0), atol=1e-16)
        curve = ForwardVarianceCurve((0.04,), (0.04,))
        yf = p.spx(fine, constant_readout(), curve, [0.04], order=4).terminal_ratios[0.04]
        yc = p.spx(coarse, constant_readout(), curve, [0.04], order=4).terminal_ratios[0.04]
        np.testing.assert_allclose(yf, yc, rtol=0, atol=3e-16)
        with self.assertRaises(ValueError):
            coarsen_paths(fine, [0., 0.015, 0.04])


if __name__ == "__main__":
    unittest.main()
