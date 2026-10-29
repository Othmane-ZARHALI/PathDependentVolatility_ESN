import unittest

from .m2_spx import ROOT
from .spx_hedging_study import choose_quotes
from esn_eq.market_targets import MarketTargets


class RepresentativeContractsTests(unittest.TestCase):
    def test_available_upper_wing_does_not_duplicate_a_contract(self):
        quotes = choose_quotes(MarketTargets.load(ROOT/"data/processed/options/2026-02-19"))
        self.assertEqual(len(quotes),12)
        self.assertFalse(quotes.duplicated(["expiry","strike","option_type"]).any())
        self.assertEqual(set(quotes.expiry_days),{29,57})
        self.assertTrue(quotes.split.eq("train").all())
        self.assertTrue((quotes.groupby(["expiry","option_type"]).size() == 3).all())


if __name__ == "__main__":
    unittest.main()
