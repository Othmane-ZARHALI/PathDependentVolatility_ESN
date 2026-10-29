"""Common price-target validation, split boundaries and failed-price handling."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from esn_eq.market_targets import MarketTargets, interval_residuals


def snapshot(directory):
    rows = []
    for strike, split in [(90., "train"), (100., "validation"), (110., "final_test")]:
        rows.append(dict(instrument="SPX", expiry=123, expiry_days=30,
                         T_years=30/365, strike=strike, forward=100., K_over_F=strike/100,
                         option_type="call", split=split, price_bid_U=1., price_ask_U=2., price_scale_U=0.1))
    pd.DataFrame(rows).to_csv(directory / "quotes_calibration.csv", index=False)
    pd.DataFrame([dict(expiry=124, T_years=0.2, forward_vix_points=20.)]).to_csv(
        directory / "vix_forwards.csv", index=False)
    pd.DataFrame([dict(curve_id="reference", start_years=0., end_years=1., xi_annual_variance=0.04)]).to_csv(
        directory / "forward_variance.csv", index=False)
    (directory / "audit.json").write_text(json.dumps(dict(reference_date="2026-02-19",
        variance_audit=dict(overall_compatibility_status="unresolved_input_conflict_requires_nuisance_curve_profile"))))
    (directory / "manifest.json").write_text(json.dumps(dict(reference_date="2026-02-19",
        price_unit="undiscounted_forward_value")))


class CommonTargetTests(unittest.TestCase):
    def test_interval_distance_preserves_zero_only_inside_quote_interval(self):
        residual = interval_residuals([0.5, 1.5, 3.], [1., 1., 1.], [2., 2., 2.], [0.5, 0.5, 0.5])
        np.testing.assert_array_equal(residual, [-1., 0., 2.])

    def test_nan_price_and_invalid_scale_are_never_silently_ignored(self):
        for price, scale in ((np.nan, 0.1), (np.inf, 0.1), (1., 0.), (-1., 0.1)):
            with self.assertRaises(ValueError):
                interval_residuals([price], [1.], [2.], [scale])

    def test_default_selection_preserves_training_boundary_and_audit_flags(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            snapshot(directory)
            targets = MarketTargets.load(directory)
            selected = targets.select("SPX")
            self.assertEqual(selected.strike.tolist(), [90.])
            self.assertEqual(len(targets.artifact_sha256), 5)
            self.assertEqual(targets.audit["variance_audit"]["overall_compatibility_status"],
                             "unresolved_input_conflict_requires_nuisance_curve_profile")
            self.assertAlmostEqual(targets.curve("reference").integral(0., 1.), 0.04)
            selected.loc[:, "strike"] = 1.
            self.assertEqual(targets.select("SPX").strike.tolist(), [90.])

    def test_loader_rejects_moneyness_mismatch_and_duplicate_contract(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            snapshot(directory)
            file = directory / "quotes_calibration.csv"
            frame = pd.read_csv(file)
            frame.loc[0, "K_over_F"] = 0.5
            frame.to_csv(file, index=False)
            with self.assertRaises(ValueError):
                MarketTargets.load(directory)
            snapshot(directory)
            frame = pd.read_csv(file)
            pd.concat([frame, frame.iloc[:1]]).to_csv(file, index=False)
            with self.assertRaises(ValueError):
                MarketTargets.load(directory)


if __name__ == "__main__":
    unittest.main()
