"""Independent-pilot controls preserve known expectations and reduce sampling error."""
import unittest

import numpy as np
from scipy.special import ndtr

from esn_eq.config import ArchitectureConfig, BankConfig, RiskPremium, ReadoutParameters
from esn_eq.pricing import DatedM2Pricer, ForwardVarianceCurve, coarsen_paths
from esn_eq.variance_reduction import (fit_payoff_controls, controlled_payoff_statistics,
                                      deterministic_stock_controls)


class StockControlTests(unittest.TestCase):
    def test_deterministic_control_matches_constant_variance_pricer_and_coupled_grid(self):
        dated = DatedM2Pricer(ArchitectureConfig(feedback=BankConfig((1.,)),
            spike=BankConfig((5.,)), cluster=BankConfig((.2,))), RiskPremium())
        fine = dated.prepare([0., .1, .2, .3, .4], np.zeros(3), 64, seed=3)
        coarse = coarsen_paths(fine, [0., .2, .4])
        curve = ForwardVarianceCurve((.2, .4), (.04, .06))
        first = deterministic_stock_controls(fine, curve, [.4])[.4]
        second = deterministic_stock_controls(coarse, curve, [.4])[.4]
        params = ReadoutParameters(cluster_loading=0., feedback_curvature=0., spike_curvature=0.)
        model = dated.spx(fine, params, curve, [.4], order=4).terminal_ratios[.4]
        # Coupled summation changes grouping; allow four machine epsilons.
        np.testing.assert_allclose(first, second, atol=4*np.finfo(float).eps, rtol=0)
        np.testing.assert_allclose(first, model, atol=4*np.finfo(float).eps, rtol=0)

    def test_constant_volatility_control_recovers_black_with_negligible_variance(self):
        pilot = np.exp(-.5*.04+.2*np.random.default_rng(1).normal(size=2048))
        evaluation = np.exp(-.5*.04+.2*np.random.default_rng(2).normal(size=4096))
        strikes, kinds = [.9, 1., 1.1], ["put", "call", "call"]
        beta = fit_payoff_controls(pilot, pilot, strikes, kinds, .04)
        prices, se = controlled_payoff_statistics(evaluation, evaluation, strikes, kinds, .04, beta)
        for k, kind, price, error in zip(strikes, kinds, prices, se):
            d1 = (-np.log(k)+.02)/.2
            d2 = d1-.2
            known = ndtr(d1)-k*ndtr(d2) if kind == "call" else k*ndtr(-d2)-ndtr(-d1)
            self.assertAlmostEqual(price, known, places=14)
            self.assertLess(error, 1e-15)

    def test_independent_pilot_preserves_expectation_for_different_model_variance(self):
        zp = np.random.default_rng(10).normal(size=4096)
        ze = np.random.default_rng(11).normal(size=16384)
        model_p, reference_p = np.exp(-.045+.3*zp), np.exp(-.02+.2*zp)
        model_e, reference_e = np.exp(-.045+.3*ze), np.exp(-.02+.2*ze)
        beta = fit_payoff_controls(model_p, reference_p, [1.], ["call"], .04)
        prices, se = controlled_payoff_statistics(model_e, reference_e, [1.], ["call"], .04, beta)
        oracle = 2*ndtr(.15)-1
        self.assertLess(abs(prices[0]-oracle), 5*se[0])
        plain_se = np.maximum(model_e-1, 0).std(ddof=1)/np.sqrt(len(ze))
        self.assertLess(se[0], plain_se/3)


if __name__ == "__main__":
    unittest.main()
