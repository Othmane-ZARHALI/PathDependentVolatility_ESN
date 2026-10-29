"""Tests for the downside-oriented fast-factor extension."""

from __future__ import annotations

import unittest

import pandas as pd

from esn_eq.asymmetry_cli import _validate_shifts, summarize_by_shift


class AsymmetryScreenTests(unittest.TestCase):
    """Protect the nested control and shift-grid interpretation."""

    def test_shift_grid_is_sorted_and_contains_exact_control(self) -> None:
        """The original M2 must be present as the zero-shift paired control."""
        self.assertEqual(_validate_shifts((0.3, 0.0, 0.1)), (0.0, 0.1, 0.3))

    def test_shift_grid_rejects_negative_duplicate_or_missing_control(self) -> None:
        """Ambiguous grids would invalidate the paired comparison."""
        for values in ((0.0, -0.1), (0.0, 0.1, 0.1), (0.1, 0.2)):
            with self.assertRaises(ValueError):
                _validate_shifts(values)

    def test_summary_counts_configurations_that_pass_on_every_seed(self) -> None:
        """Robustness means joint acceptance for all seeds, not a high pooled rate."""
        rows = []
        for seed in (1, 2):
            for shortlist_id in (0, 1):
                row = {
                    "seed": seed,
                    "shortlist_id": shortlist_id,
                    "spike_shift": 0.1,
                    "accepted": shortlist_id == 0 or seed == 1,
                    "distance": 0.0,
                    "rough_hurst": 0.1,
                    "leverage_rank": -0.02,
                    "zumbach_rank": 0.03,
                    "excess_kurtosis": 5.0,
                    "max_abs_return_acf": 0.02,
                    "mean_log_variance_acf": 0.5,
                    "taylor_gap": 0.02,
                    "p_cap_exceedance_fraction": 0.0,
                }
                row.update({name: 0.1 for name in (
                    "tail_negative_day0_response",
                    "tail_positive_day0_response",
                    "tail_day0_asymmetry_gap",
                    "tail_days1_5_asymmetry_gap",
                )})
                rows.append(row)
        summary = summarize_by_shift(pd.DataFrame(rows)).iloc[0]
        self.assertEqual(summary.robust_configuration_count, 1)
        self.assertEqual(summary.configuration_count, 2)
        self.assertEqual(summary.pass_zumbach_rank_mean, 1.0)


if __name__ == "__main__":
    unittest.main()
