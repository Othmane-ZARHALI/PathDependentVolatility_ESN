"""Validate cached normalization and coupled batching before high-resolution runs."""
import contextlib
import io
import unittest

import numpy as np

from .fine_grid import CachedNormalizerPricer, bundle
from esn_eq.config import ArchitectureConfig, BankConfig, RiskPremium, ReadoutParameters
from esn_eq.pricing import DatedM2Pricer, ForwardVarianceCurve


def fixture():
    architecture = ArchitectureConfig(feedback=BankConfig((1.,)),
        spike=BankConfig((4.,)),cluster=BankConfig((.2,)))
    return architecture,RiskPremium()


class FineGridTests(unittest.TestCase):
    def test_cached_normalizer_matches_direct_value_and_changes_with_readout(self):
        architecture,premium = fixture()
        direct = DatedM2Pricer(architecture,premium)
        cached = CachedNormalizerPricer(architecture,premium)
        parameters = ReadoutParameters()
        expected = direct.normalizer([0.,.1],np.zeros(3),parameters,8)
        actual = cached.normalizer([0.,.1],np.zeros(3),parameters,8)
        np.testing.assert_array_equal(actual,expected)
        self.assertIs(actual,cached.normalizer([0.,.1],np.zeros(3),parameters,8))
        other = cached.normalizer([0.,.1],np.zeros(3),ReadoutParameters(feedback_shift=2.),8)
        self.assertFalse(np.array_equal(actual,other))

    def test_batched_constant_volatility_has_coupled_paths_and_distinct_streams(self):
        architecture,premium = fixture()
        pricer = CachedNormalizerPricer(architecture,premium)
        parameters = ReadoutParameters(cluster_loading=0.,feedback_curvature=0.,spike_curvature=0.)
        curve = ForwardVarianceCurve((.2,),(.04,))
        with contextlib.redirect_stdout(io.StringIO()):
            outputs,seeds = bundle(pricer,np.array([0.,.05,.1,.15,.2]),np.array([0.,.1,.2]),
                np.zeros(3),parameters,curve,[.2],8,4,61,1,4,True)
        self.assertEqual(len(set(seeds)),2)
        np.testing.assert_allclose(outputs['fine'][.2],outputs['coarse'][.2],rtol=0,atol=1e-15)
        np.testing.assert_allclose(outputs['fine'][.2],outputs['reference'][.2],rtol=0,atol=1e-15)
        self.assertFalse(np.array_equal(outputs['fine'][.2][:4],outputs['fine'][.2][4:]))


if __name__ == '__main__':
    unittest.main()
