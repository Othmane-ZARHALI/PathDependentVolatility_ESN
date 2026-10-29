import json
import unittest

import numpy as np

from .m2_spx import ROOT
from .vix_near_fit import active_targets,prefix_cache
from esn_eq.market_targets import MarketTargets
from esn_eq.vix_curve_profile import VIXCurveExposures


class ScopeTests(unittest.TestCase):
    def test_full_scope_restored_and_historical_near_scope_explicit(self):
        targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
        scope = json.loads((ROOT/"matteo/model_comparison/option_pricing/active_scope.json").read_text())
        self.assertEqual(scope["vix_fit_expiry_days"],[27,55,118,153])
        self.assertEqual(scope["vix_deferred_expiry_days"],[])
        self.assertEqual(len(targets.select("VIX",expiry_days=scope["vix_fit_expiry_days"])),175)
        restored_forwards = targets.vix_forwards[targets.vix_forwards.expiry_days.isin(scope["vix_fit_expiry_days"])]
        self.assertEqual(len(restored_forwards),4)
        self.assertEqual(set(restored_forwards.expiry_days),{27,55,118,153})
        with self.assertRaises(ValueError):
            active_targets(targets,scope)
        historical = json.loads((ROOT/"matteo/model_comparison/option_pricing/near_scope_historical.json").read_text())
        quotes,forwards = active_targets(targets,historical)
        self.assertEqual(len(quotes),87)
        self.assertEqual(set(quotes.expiry_days),{27,55})
        self.assertEqual(set(forwards.expiry_days),{27,55})
        self.assertEqual(len(forwards),2)
        self.assertEqual(set(targets.select("VIX").expiry_days),{27,55,118,153})
        self.assertTrue(quotes.split.eq("train").all())

    def test_prefix_preserves_each_independent_scramble(self):
        coefficients = np.arange(8*3).reshape(8,3)
        cache = VIXCurveExposures(coefficients,0.,(.25,.75,1.),.1,32,None,())
        prefix = prefix_cache(cache,2,4,2)
        np.testing.assert_array_equal(prefix.coefficients,coefficients[[0,1,4,5]])
        self.assertFalse(prefix.coefficients.flags.writeable)


if __name__ == "__main__":
    unittest.main()
