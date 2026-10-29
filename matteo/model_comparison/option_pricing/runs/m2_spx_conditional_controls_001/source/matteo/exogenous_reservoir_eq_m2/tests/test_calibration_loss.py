import unittest

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.special import ndtr

from esn_eq.calibration_loss import JointPriceObjective, LossConfig


def targets():
    rows = []
    for instrument, expiry, strikes in [("SPX", 1, [90., 100.]), ("SPX", 2, [100.]),
                                        ("VIX", 3, [20., 21., 22.])]:
        for strike in strikes:
            rows.append(dict(instrument=instrument, expiry=expiry, strike=strike,
                forward=100. if instrument == "SPX" else 20., T_years=.2,
                option_type="call", split="train", price_bid_U=1., price_ask_U=2., price_scale_U=.5))
    return pd.DataFrame(rows), pd.DataFrame([dict(expiry=3, forward_vix_points=20.)])


class CalibrationLossTests(unittest.TestCase):
    def test_family_and_expiry_weights_and_residual_norm(self):
        q, f = targets()
        result = JointPriceObjective(q, f).evaluate([3., 3., 4., 2.5, 2.5, 2.5], [21.], state_penalty=.1)
        # SPX: average expiry means (4+16)/2; VIX: 1; forward: 4.
        self.assertEqual(result["components"]["SPX"], 10.)
        self.assertEqual(result["components"]["VIX"], 1.)
        self.assertAlmostEqual(result["total"], 15.1)
        self.assertAlmostEqual(result["residuals"] @ result["residuals"], result["total"])

    def test_quote_count_does_not_change_family_or_expiry_weight(self):
        q, f = targets()
        objective = JointPriceObjective(q, f)
        base = objective.evaluate(np.full(len(q), 3.), [20.])["total"]
        extra = q.iloc[[0]].copy(); extra["strike"] = 95.
        expanded = pd.concat([q, extra], ignore_index=True)
        self.assertEqual(JointPriceObjective(expanded, f).evaluate(np.full(len(expanded), 3.), [20.])["total"], base)

    def test_frozen_targets_and_explicit_heldout_boundary(self):
        q, f = targets()
        objective = JointPriceObjective(q, f)
        q.loc[0, "price_bid_U"] = 999.
        self.assertEqual(objective.evaluate(np.full(len(q), 1.5), [20.])["total"], 0.)
        q, f = targets(); q.loc[0, "split"] = "final_test"
        with self.assertRaises(ValueError):
            JointPriceObjective(q, f)

    def test_failed_prices_and_missing_predictions_cannot_disappear(self):
        q, f = targets(); objective = JointPriceObjective(q, f)
        for bad in (np.nan, np.inf, -1.):
            prices = np.full(len(q), 1.5); prices[0] = bad
            with self.assertRaises(ValueError):
                objective.evaluate(prices, [20.])
        with self.assertRaises(ValueError):
            objective.evaluate(np.ones(len(q)-1), [20.])
        with self.assertRaises(ValueError):
            objective.evaluate(np.ones(len(q)), [])

    def test_configured_sensitivity_and_regularisation_are_separate(self):
        q, f = targets()
        cfg = LossConfig(forward_weight=0., midpoint_weight=.05)
        result = JointPriceObjective(q, f, cfg).evaluate(np.ones(len(q)), [21.], curve_penalty=.2)
        self.assertAlmostEqual(result["components"]["midpoint"], .1)
        self.assertAlmostEqual(result["total"], .3)
        self.assertEqual(result["forwards"].error_points.iloc[0], 1.)
        for kwargs in (dict(forward_scale=0.), dict(vix_weight=-1.), dict(spx_weight=np.nan)):
            with self.assertRaises(ValueError):
                LossConfig(**kwargs)

    def test_synthetic_constant_variance_fit_recovers_shared_spx_vix_parameter(self):
        rows = []
        for instrument, expiry in (("SPX", 1), ("SPX", 2), ("VIX", 3)):
            for strike in ([90., 100., 110.] if instrument == "SPX" else [18., 22., 26.]):
                rows.append(dict(instrument=instrument, expiry=expiry, strike=strike,
                    forward=100. if instrument == "SPX" else 22., T_years=.1*expiry,
                    option_type="call", split="train", price_scale_U=.05))
        q = pd.DataFrame(rows)

        def predictions(sigma):
            std = sigma*np.sqrt(q.T_years.to_numpy())
            strike = q.strike.to_numpy()
            d1 = (np.log(100/strike)+.5*std**2)/std
            black = 100*ndtr(d1)-strike*ndtr(d1-std)
            return np.where(q.instrument == "SPX", black, np.maximum(100*sigma-strike, 0.))

        truth = predictions(.22)
        q["price_bid_U"], q["price_ask_U"] = np.maximum(truth-.01, 0.), truth+.01
        f = pd.DataFrame([dict(expiry=3, forward_vix_points=22.)])
        objective = JointPriceObjective(q, f)
        for initial in (.1, .4):
            fit = least_squares(lambda x: objective.evaluate(predictions(x[0]), [100*x[0]])["residuals"],
                                [initial], bounds=(.05, .6), xtol=1e-12, ftol=1e-12, gtol=1e-12)
            self.assertTrue(fit.success)
            self.assertAlmostEqual(fit.x[0], .22, places=8)
            self.assertLess(objective.evaluate(predictions(fit.x[0]), [100*fit.x[0]])["total"], 1e-12)


if __name__ == "__main__":
    unittest.main()
