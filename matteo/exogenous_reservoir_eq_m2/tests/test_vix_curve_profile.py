from dataclasses import replace
import unittest

import numpy as np

from esn_eq.config import ArchitectureConfig, BankConfig, ReadoutParameters, RiskPremium
from esn_eq.pricing import DatedM2Pricer, ForwardVarianceCurve
from esn_eq.vix_pricing import ConditionalVIXPricer
from esn_eq.vix_curve_profile import prepare_curve_exposures


def fixture():
    dated = DatedM2Pricer(ArchitectureConfig(feedback=BankConfig((1.,8.),.7,True),
        spike=BankConfig((3.,25.),.4,True),cluster=BankConfig((.2,2.),0.,False)),RiskPremium())
    parameters = ReadoutParameters(cluster_loading=.15,feedback_curvature=.04,spike_curvature=.03)
    curve = ForwardVarianceCurve((.21,.3,1.),(.04,.06,.05))
    pricer = ConditionalVIXPricer(dated,parameters,curve,np.zeros(6),normalizer_order=16)
    states = np.array([dated.initial_state(u)[0] for u in ([0.,0.,0.],[1.,-2.,-2.])])
    return pricer,states


class VIXCurveProfileTests(unittest.TestCase):
    def test_cached_values_match_direct_repricing_for_different_levels(self):
        pricer,states = fixture()
        window = pricer.gaussian_window(.2,time_order=6)
        cache = prepare_curve_exposures(pricer,states,window,order=16)
        for levels in ((.04,.06,.05),(.025,.09,.07)):
            other = ConditionalVIXPricer(pricer.dated,pricer.parameters,
                ForwardVarianceCurve(pricer.curve.ends,levels),pricer.initial_state,normalizer_order=16)
            direct = other.evaluate(states,other.gaussian_window(.2,time_order=6),order=16)
            np.testing.assert_allclose(cache.evaluate(levels).vix,direct.vix,rtol=2e-15)
        self.assertTrue(np.all(cache.coefficients[:,2] == 0))
        self.assertFalse(cache.coefficients.flags.writeable)

    def test_integrated_variance_is_affine_with_nonzero_floor(self):
        pricer,states = fixture()
        cache = prepare_curve_exposures(pricer,states,pricer.gaussian_window(.2),order=12)
        a,b = np.array([.02,.05,.06]),np.array([.08,.07,.03])
        lhs = cache.evaluate(.3*a+.7*b).integrated_variance
        rhs = .3*cache.evaluate(a).integrated_variance+.7*cache.evaluate(b).integrated_variance
        np.testing.assert_allclose(lhs,rhs,rtol=1e-15)

    def test_constant_variance_exposures_equal_segment_lengths(self):
        pricer,states = fixture()
        pricer.parameters = ReadoutParameters(cluster_loading=0.,feedback_curvature=0.,spike_curvature=0.)
        window = pricer.gaussian_window(.2)
        cache = prepare_curve_exposures(pricer,states,window,order=8)
        np.testing.assert_allclose(cache.coefficients,np.tile([.01,window.delta-.01,0.],(2,1)),atol=1e-16)
        levels = [.04,.06,.07]
        expected = 1e4*ForwardVarianceCurve(pricer.curve.ends,tuple(levels)).integral(.2,.2+window.delta)/window.delta
        np.testing.assert_allclose(cache.evaluate(levels).vix**2,expected,atol=1e-11)

    def test_invalid_levels_and_stale_parameters_are_rejected(self):
        pricer,states = fixture()
        cache = prepare_curve_exposures(pricer,states,pricer.gaussian_window(.2),order=8)
        for levels in ([.04,.06],[.04,np.nan,.06],[.04,1e-4,.06]):
            with self.assertRaises(ValueError):
                cache.evaluate(levels)
        pricer.parameters = replace(pricer.parameters,feedback_curvature=.1)
        with self.assertRaises(ValueError):
            cache.evaluate([.04,.05,.06])


if __name__ == "__main__":
    unittest.main()
