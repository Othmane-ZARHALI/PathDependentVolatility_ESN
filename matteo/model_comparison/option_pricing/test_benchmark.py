"""Focused correspondence checks for the legacy pricing benchmark."""
import unittest

import numpy as np

from .benchmark import fit_metrics, legacy_modules, scalar_batch_check


class BenchmarkTests(unittest.TestCase):
    def test_scalar_and_batched_pathwise_correspondence(self):
        self.assertTrue(scalar_batch_check()["passed_at_1e-12"])

    def test_conditional_restart_uses_all_state_coordinates_without_mutation(self):
        m = legacy_modules()
        c, v = m["calibrate_smile"], m["vectorized_sim"]
        mat = c.build_matrices_raw(0.386, 0.672, 1e-5, 4.8809, 0.1745)
        _, _, state = v.simulate_batch(c.OPT_ARCH, mat, 0.22, 20, 3, seed=72)
        original = tuple(x.copy() for x in state)
        x, var, final = v.simulate_batch(c.OPT_ARCH, mat, 0.22, 2, 3,
                                       *state, seed=83, b0_delta=-0.4)
        rng = np.random.default_rng(83)
        p = v.prepare(c.OPT_ARCH, mat, b0_delta=-0.4)
        expected, dx0, var0 = v.step(original, p, 0.22, rng)
        expected, dx1, var1 = v.step(expected, p, 0.22, rng)
        np.testing.assert_allclose(x, np.stack([dx0, dx1]), rtol=0, atol=1e-14)
        np.testing.assert_allclose(var, np.stack([var0, var1]), rtol=0, atol=1e-14)
        for got, want, before, after in zip(final, expected, original, state):
            np.testing.assert_array_equal(got, want)
            np.testing.assert_array_equal(before, after)

    def test_missing_iv_is_visible_even_if_surviving_fit_is_perfect(self):
        result = fit_metrics(np.array([0.2, np.nan]), np.array([0.2, 0.9]),
                             np.array([0.25, 0.75]))
        self.assertEqual(result["finite_only_weighted_rmse"], 0)
        self.assertEqual(result["failed_target_weight"], 0.75)
        self.assertEqual(result["n_failed_model_iv"], 1)
        self.assertFalse(result["complete_iv_comparison"])

    def test_all_failed_iv_has_no_reported_rmse(self):
        result = fit_metrics(np.array([np.nan]), np.array([0.2]), np.array([1.0]))
        self.assertIsNone(result["finite_only_weighted_rmse"])
        self.assertIsNone(result["finite_only_rmse"])

    def test_constant_volatility_black_roundtrip(self):
        m = legacy_modules()
        for strike in (0.8, 1.0, 1.2):
            for is_call in (False, True):
                price = m["market_vix"].black76_undisc(1, strike, 0.5, 0.3, is_call)
                iv = m["model_pricing"].implied_vol_from_price(1, strike, 0.5, price, is_call)
                self.assertAlmostEqual(iv, 0.3, places=6)

    def test_constant_volatility_variance_replication_refines_toward_oracle(self):
        v = legacy_modules()["market_vix"]
        errors = []
        for atm, wing in ((0.02, 0.06), (0.002, 0.006)):
            strikes = v.build_moneyness_grid(atm_step=atm, wing_step=wing)
            prices = [v.black76_undisc(1, k, 30 / 365, 0.2, k >= 1) for k in strikes]
            estimate = v.cboe_variance_from_grid(1, 30 / 365, strikes, prices)["sigma2"]
            errors.append(abs(estimate - 0.04))
        self.assertLess(errors[1], errors[0])
        self.assertLess(errors[1], 1e-4)


if __name__ == "__main__":
    unittest.main()
