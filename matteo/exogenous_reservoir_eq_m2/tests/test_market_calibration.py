from pathlib import Path
import unittest

import numpy as np

from esn_eq.config import ArchitectureConfig, BankConfig, RiskPremium
from esn_eq.market_calibration import NuisanceCurveFamily, restricted_initial_state, risk_price_scenarios
from esn_eq.market_targets import MarketTargets
from esn_eq.pricing import DatedM2Pricer


def fixture():
    return DatedM2Pricer(ArchitectureConfig(feedback=BankConfig((1.,8.),.7,True),
        spike=BankConfig((3.,25.),.4,True),cluster=BankConfig((.2,2.),0.,False)),RiskPremium())


class MarketCalibrationInputTests(unittest.TestCase):
    def test_declared_snapshot_family_bounds_penalty_and_exact_integral(self):
        root = Path(__file__).resolve().parents[3]
        targets = MarketTargets.load(root/"data/processed/options/2026-02-19")
        family = NuisanceCurveFamily.from_targets(targets)
        self.assertEqual(family.ends[:2],(.25,.75))
        curve = family.build([.04,.05,.06])
        self.assertAlmostEqual(curve.integral(0.,1.),.04*.25+.05*.5+.06*.25)
        self.assertAlmostEqual(family.penalty([.04,.05,.06]),.1*(.01**2+.01**2))
        for levels in ([.04,.05], [0.,.05,.06], [.04,.05,.251], [.04,np.nan,.06]):
            with self.assertRaises(ValueError):
                family.build(levels)
        # Constructing a candidate must not clear the preprocessing conflict.
        self.assertIn("unresolved",targets.audit["variance_audit"]["overall_compatibility_status"])

    def test_restricted_state_order_and_penalty_are_explicit(self):
        dated = fixture()
        projections = [.4,-.7,1.2]
        result = restricted_initial_state(dated,projections,weight=.1)
        np.testing.assert_allclose(dated.projection @ result.full_state,projections,atol=1e-14)
        covariance = dated.projection @ dated.stationary @ dated.projection.T
        norm = np.array(projections) @ np.linalg.solve(covariance,projections)
        self.assertAlmostEqual(result.prior_squared_norm,norm)
        self.assertAlmostEqual(result.penalty,.1*norm/3)
        self.assertEqual(restricted_initial_state(dated,[0.,0.,0.]).penalty,0.)
        self.assertFalse(result.full_state.flags.writeable)

    def test_zero_unspanned_scenario_preserves_equity_induced_drift(self):
        scenarios = risk_price_scenarios(.35)
        self.assertEqual(len(scenarios),5)
        dated = DatedM2Pricer(fixture().architecture,scenarios[0])
        expected = -np.sqrt(2*dated.layout.rates)*dated.layout.driver_loadings[:,0]*.35
        np.testing.assert_allclose(dated.q_drift,expected,atol=1e-15)
        self.assertTrue(np.any(dated.q_drift != 0))
        self.assertTrue(all(s.equity == .35 for s in scenarios))
        self.assertEqual({(s.feedback_idiosyncratic,s.spike_idiosyncratic) for s in scenarios},
                         {(0.,0.),(.25,0.),(-.25,0.),(0.,.25),(0.,-.25)})


if __name__ == "__main__":
    unittest.main()
