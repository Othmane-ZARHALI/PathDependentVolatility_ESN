import unittest

import numpy as np

from esn_eq.config import ArchitectureConfig,BankConfig,ReadoutParameters,RiskPremium
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.pricing import DatedM2Pricer,ForwardVarianceCurve,coarsen_paths


def fixture():
    dated = DatedM2Pricer(ArchitectureConfig(feedback=BankConfig((1.,8.),.7,True),
        spike=BankConfig((3.,25.),.4,True),cluster=BankConfig((.2,2.),0.,False)),RiskPremium())
    return dated,ForwardVarianceCurve((1.,),(.04,))


class HedgingControlTests(unittest.TestCase):
    def test_constant_variance_path_and_dynamic_control_refinement(self):
        dated,curve = fixture()
        parameters = ReadoutParameters(cluster_loading=0.,feedback_curvature=0.,spike_curvature=0.)
        fine = np.linspace(0.,.25,129)
        cache = dated.prepare(fine,np.zeros(6),4096,91)
        fine_result = spx_hedging_controls(dated,cache,parameters,curve,[.25],[1.],["call"],order=8)
        oracle = dated.spx(cache,parameters,curve,[.25],order=8)
        np.testing.assert_allclose(fine_result.terminal_ratios[.25],oracle.terminal_ratios[.25],rtol=2e-15)
        coarse = spx_hedging_controls(dated,coarsen_paths(cache,fine[::8]),parameters,curve,[.25],[1.],["call"],order=8)
        payoff = np.maximum(fine_result.terminal_ratios[.25]-1.,0.)
        adjusted = payoff-fine_result.hedges[:,0,0]
        self.assertGreater(payoff.var()/adjusted.var(),20.)
        self.assertLess(adjusted.var(),(payoff-coarse.hedges[:,0,0]).var())
        np.testing.assert_allclose(fine_result.hedges[:,:,0],fine_result.hedges[:,:,1],atol=2e-15)

    def test_nonconstant_variance_control_means_are_consistent_with_zero(self):
        dated,curve = fixture()
        cache = dated.prepare(np.linspace(0.,.2,33),np.zeros(6),8192,72)
        parameters = ReadoutParameters(cluster_loading=.2,feedback_curvature=.04,spike_curvature=.03)
        result = spx_hedging_controls(dated,cache,parameters,curve,[.2,.2],[.9,1.1],["put","call"],order=12)
        se = result.hedges.std(axis=0,ddof=1)/np.sqrt(8192)
        self.assertTrue(np.all(np.abs(result.hedges.mean(axis=0)) < 6*se))
        oracle = dated.spx(cache,parameters,curve,[.2],order=12)
        np.testing.assert_allclose(result.terminal_ratios[.2],oracle.terminal_ratios[.2],rtol=2e-15)

    def test_invalid_contracts_fail(self):
        dated,curve = fixture()
        cache = dated.prepare([0.,.1],np.zeros(6),4,1)
        with self.assertRaises(ValueError):
            spx_hedging_controls(dated,cache,ReadoutParameters(),curve,[.2],[1.],["call"])


if __name__ == "__main__":
    unittest.main()
