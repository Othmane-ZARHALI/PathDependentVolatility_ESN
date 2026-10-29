"""Conditional VIX checks against constant-volatility and nested-MC oracles."""
import unittest
from dataclasses import replace

import numpy as np

from esn_eq.config import ArchitectureConfig, BankConfig, ReadoutParameters, RiskPremium
from esn_eq.pricing import DatedM2Pricer, ForwardVarianceCurve, payoff_statistics
from esn_eq.vix_pricing import ConditionalVIXPricer


def fixture(parameters=None):
    dated = DatedM2Pricer(ArchitectureConfig(feedback=BankConfig((1., 8.), .7, True),
        spike=BankConfig((3., 25.), .4, True), cluster=BankConfig((.2, 2.), 0., False)), RiskPremium())
    parameters = parameters or ReadoutParameters(cluster_loading=.15,
        feedback_curvature=.04, feedback_shift=.5, spike_curvature=.03)
    curve = ForwardVarianceCurve((.21, 1.), (.04, .06))
    vix = ConditionalVIXPricer(dated, parameters, curve, np.zeros(6), normalizer_order=16)
    return dated, vix


class VIXPricingTests(unittest.TestCase):
    def test_deterministic_variance_and_payoffs_are_exact_for_both_engines(self):
        dated, vix = fixture(ReadoutParameters(cluster_loading=0., feedback_curvature=0., spike_curvature=0.))
        rule = vix.window(.2, max_step=.03)
        states = np.random.default_rng(1).normal(size=(4, 6))
        expected = 100*np.sqrt(vix.curve.integral(.2, .2+30/365)/(30/365))
        for engine in ("gaussian_quadrature", "nested_mc"):
            result = vix.evaluate(states, rule, engine=engine, order=8, inner_paths=8, seed=2)
            np.testing.assert_allclose(result.vix, expected, rtol=0, atol=1e-12)
            price, se = payoff_statistics(result.vix, [expected-1, expected+1], ["call", "put"])
            np.testing.assert_allclose(price, [1., 1.], atol=1e-12)
            np.testing.assert_allclose(se, 0., atol=1e-12)

    def test_state_changes_values_without_changing_date_zero_normalizer(self):
        dated, vix = fixture()
        rule = vix.window(.2, max_step=.04)
        saved = rule.normalizer.copy()
        low, _ = dated.initial_state([-1., 0., 0.])
        high, _ = dated.initial_state([1., 0., 0.])
        result = vix.evaluate(np.stack([low, high]), rule, order=12)
        self.assertGreater(result.vix[1], result.vix[0]+1.)
        np.testing.assert_array_equal(rule.normalizer, saved)

    def test_date_zero_tower_endpoint_recovers_curve_integral(self):
        dated, vix = fixture()
        rule = vix.window(0., max_step=.025)
        result = vix.evaluate(np.zeros((2, 6)), rule, order=16)
        np.testing.assert_allclose(result.vix**2, vix.moment_target(rule), rtol=1e-13)

    def test_nested_full_paths_agree_with_conditional_quadrature(self):
        dated, vix = fixture()
        rule = vix.window(.2, max_step=.025)
        states = np.array([dated.initial_state(u)[0] for u in ([0., 0., 0.], [.5, -.3, .7], [-.5, .8, -.4])])
        quadrature = vix.evaluate(states, rule, order=20)
        nested = vix.evaluate(states, rule, engine="nested_mc", inner_paths=8192, seed=71)
        self.assertTrue(np.all(np.abs(nested.integrated_variance-quadrature.integrated_variance)
                               < 6*nested.integrated_variance_se))

    def test_inner_streams_are_invariant_to_outer_batching(self):
        dated, vix = fixture()
        rule = vix.window(.2, max_step=.04)
        states = np.random.default_rng(9).normal(size=(3, 6))
        a = vix.evaluate(states, rule, engine="nested_mc", inner_paths=64, seed=22, batch_size=1)
        b = vix.evaluate(states, rule, engine="nested_mc", inner_paths=64, seed=22, batch_size=3)
        np.testing.assert_allclose(a.integrated_variance, b.integrated_variance, rtol=1e-14)

    def test_outer_second_moment_matches_curve_with_independent_sampling_error(self):
        dated, vix = fixture()
        rule = vix.window(.2, max_step=.03)
        cache = dated.prepare([0., .2], np.zeros(6), 4096, seed=94, checkpoints=[.2])
        result = vix.evaluate(cache.checkpoint_states[.2], rule, order=12)
        squared = result.vix**2
        se = squared.std(ddof=1)/np.sqrt(len(squared))
        self.assertLess(abs(squared.mean()-vix.moment_target(rule)), 5*se)
        self.assertLessEqual(result.vix.mean()**2, squared.mean()+1e-12)

    def test_window_from_another_date_zero_model_is_rejected(self):
        _, first = fixture()
        _, other = fixture()
        with self.assertRaises(ValueError):
            other.evaluate(np.zeros((2, 6)), first.window(.2))

    def test_stale_normalizer_after_parameter_change_is_rejected(self):
        _, vix = fixture()
        rule = vix.window(.2)
        vix.parameters = replace(vix.parameters, feedback_curvature=.1)
        with self.assertRaises(ValueError):
            vix.evaluate(np.zeros((2, 6)), rule)

    def test_gaussian_time_rule_preserves_curve_jumps_and_constant_variance(self):
        _, vix = fixture(ReadoutParameters(cluster_loading=0., feedback_curvature=0., spike_curvature=0.))
        rule = vix.gaussian_window(.2, time_order=6)
        self.assertAlmostEqual(rule.weights.sum(), rule.delta)
        self.assertTrue(np.all(rule.weights > 0))
        self.assertAlmostEqual(rule.weights[rule.times < .21].sum(), .01)
        result = vix.evaluate(np.zeros((1, 6)), rule, order=8)
        self.assertAlmostEqual(result.vix[0]**2, vix.moment_target(rule), places=10)
        with self.assertRaises(ValueError):
            vix.gaussian_window(.2, time_order=1)

    def test_gaussian_time_rule_matches_fine_midpoint_on_nonconstant_states(self):
        dated, vix = fixture()
        states = np.array([dated.initial_state(u)[0] for u in ([0., 0., 0.], [1., -2., -2.])])
        reference = vix.evaluate(states, vix.window(.2, max_step=1/10000), order=16)
        gaussian = vix.evaluate(states, vix.gaussian_window(.2, time_order=8), order=16)
        np.testing.assert_allclose(gaussian.vix, reference.vix, rtol=0, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
