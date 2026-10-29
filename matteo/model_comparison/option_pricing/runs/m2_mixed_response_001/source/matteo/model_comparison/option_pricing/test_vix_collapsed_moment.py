import unittest
import numpy as np
from .vix_collapsed_moment import tower_errors,collapsed_integrands
from esn_eq.config import ArchitectureConfig,BankConfig,RiskPremium,ReadoutParameters
from esn_eq.pricing import DatedM2Pricer,ForwardVarianceCurve
from esn_eq.vix_pricing import ConditionalVIXPricer


class CollapsedMomentTests(unittest.TestCase):
    def fixture(self):
        dated=DatedM2Pricer(ArchitectureConfig(feedback=BankConfig((1.,8.),.7,True),spike=BankConfig((3.,25.),.4,True),cluster=BankConfig((.2,2.),0.,False)),RiskPremium(equity=.2,feedback_idiosyncratic=-.1))
        return dated,np.arange(6)/10

    def test_gaussian_tower_at_multiple_lags_and_nonzero_initial_state(self):
        dated,initial=self.fixture()
        for lag in (0.,.001,.05,.4):
            self.assertLess(max(tower_errors(dated,initial,.2,lag)),2e-14)

    def test_constant_readout_recovers_piecewise_variance_integral_per_draw(self):
        dated,initial=self.fixture()
        curve=ForwardVarianceCurve((.22,1.),(.04,.09))
        parameters=ReadoutParameters(cluster_loading=0.,feedback_curvature=0.,spike_curvature=0.)
        engine=ConditionalVIXPricer(dated,parameters,curve,initial,normalizer_order=8)
        window=engine.gaussian_window(.2,time_order=4)
        values=collapsed_integrands(engine,window,np.random.default_rng(42).normal(size=(32,3)),[window.normalizer])
        truth=10000*curve.integral(.2,.2+window.delta)/window.delta
        np.testing.assert_allclose(values,truth,rtol=2e-14)

    def test_node_cache_reproduces_direct_conditional_pricer_with_floor(self):
        from .vix_matched_moment_refinement import values_from_nodes
        dated,initial=self.fixture();curve=ForwardVarianceCurve((.22,1.),(.04,.09))
        engine=ConditionalVIXPricer(dated,ReadoutParameters(),curve,initial,normalizer_order=8)
        window=engine.gaussian_window(.2,time_order=4)
        states=np.random.default_rng(17).normal(size=(8,len(initial)))
        expected=[]
        for t in window.times:
            mean,cov=dated.gaussian_projection_law(float(t-.2),states)
            expected.append(dated.gaussian_expectation(mean,cov,engine.parameters,6))
        cached=values_from_nodes(engine,window,np.column_stack(expected),window.normalizer)
        direct=engine.evaluate(states,window,order=6).vix
        np.testing.assert_allclose(cached,direct,rtol=2e-14)
