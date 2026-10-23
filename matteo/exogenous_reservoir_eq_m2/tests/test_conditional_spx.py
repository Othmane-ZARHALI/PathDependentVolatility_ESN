"""Gaussian split correspondence, nonnegative Black oracles and conditioning."""
from dataclasses import replace
import unittest
import numpy as np
from scipy.special import ndtr
from esn_eq.conditional_spx import independent_equity_fraction,partial_transition,lognormal_payoffs,simulate_conditional_spx
from esn_eq.pricing import DatedM2Pricer,ForwardVarianceCurve
from esn_eq.states import joint_transition_covariance
from test_pricing import small_pricer,constant_readout


class ConditionalSPXTests(unittest.TestCase):
    def assert_mc(self,x,target):
        self.assertTrue(np.all(abs(x.mean(axis=0)-target)<6*x.std(axis=0,ddof=1)/np.sqrt(len(x))+1e-11))

    def test_primitive_fraction_and_covariance(self):
        p=small_pricer();expected=1/(1+.7**2/(1-.7**2)+.4**2/(1-.4**2))
        self.assertAlmostEqual(independent_equity_fraction(p.layout),expected,places=14)
        for steps in [2,64,512]:
            dt=1/(252*steps);root,independent,audit=partial_transition(p,dt)
            q=joint_transition_covariance(p.layout,dt);reconstructed=root@root.T;reconstructed[0,0]+=independent
            np.testing.assert_allclose(reconstructed,q,atol=1e-14,rtol=1e-10)
            self.assertLess(audit['max_standardized_covariance_error'],2e-9)

    def test_black_nonnegative_parity_extreme_and_zero_variance(self):
        for total_vol in [.001,.2,1.]:
            k=np.array([.1,.9,1.,1.1,10.]);v=total_vol**2;m=-.5*v
            call=lognormal_payoffs(m,v,k,'call');put=lognormal_payoffs(m,v,k,'put')
            d1=(-np.log(k)+.5*v)/total_vol;d2=d1-total_vol
            np.testing.assert_allclose(call,ndtr(d1)-k*ndtr(d2),atol=1e-14)
            np.testing.assert_allclose(call-put,1-k,atol=2e-14)
            self.assertTrue(np.all(call>=0)&np.all(put>=0))
        self.assertEqual(float(lognormal_payoffs(0.,0.,1.2,'call')),0.)
        self.assertAlmostEqual(float(lognormal_payoffs(0.,0.,1.2,'put')),.2)
        self.assertTrue(np.isfinite(lognormal_payoffs(-1000.,.00001,[.1,1.,10.],'call')).all())

    def test_independent_constant_volatility_prices_exactly(self):
        p=small_pricer();arch=replace(p.architecture,feedback=replace(p.architecture.feedback,asset_correlation=0),spike=replace(p.architecture.spike,asset_correlation=0))
        p=DatedM2Pricer(arch,p.premium);times=np.array([0.,.05,.1]);k=np.array([.9,1.05]);kinds=np.array(['put','call'])
        result=simulate_conditional_spx(p,times,np.zeros(len(p.layout.rates)),constant_readout(),ForwardVarianceCurve((.1,),(.04,)),[.1,.1],k,kinds,2048,91,order=4)
        expected=lognormal_payoffs(-.5*.004,.004,k,kinds)
        np.testing.assert_allclose(result.conditional_payoffs,np.broadcast_to(expected,result.conditional_payoffs.shape),atol=1e-14)
        np.testing.assert_allclose(result.conditional_stock,1.,atol=1e-14)

    def test_nonconstant_correlated_conditional_mean_and_legacy_stock_rule(self):
        p=small_pricer();times=np.array([0.,.02,.05,.1]);initial=np.zeros(len(p.layout.rates))
        parameters=replace(constant_readout(),cluster_loading=.15,feedback_curvature=.03,spike_curvature=.02)
        curve=ForwardVarianceCurve((.1,),(.04,));k=np.array([.95,1.05]);kinds=np.array(['put','call'])
        result=simulate_conditional_spx(p,times,initial,parameters,curve,[.1,.1],k,kinds,60000,92,order=12,keep_cache=True)
        self.assert_mc(result.raw_payoffs-result.conditional_payoffs,np.zeros(2))
        self.assert_mc(result.raw_stock,np.ones(2));self.assert_mc(result.conditional_stock,np.ones(2))
        self.assertTrue(np.all(result.conditional_payoffs.var(axis=0)<1.02*result.raw_payoffs.var(axis=0)))
        old=p.spx(result.cache,parameters,curve,[.1],order=12)
        np.testing.assert_allclose(result.raw_stock[:,0],old.terminal_ratios[.1],rtol=1e-14,atol=1e-14)

    def test_perfect_equity_driver_has_no_independent_noise(self):
        p=small_pricer();arch=replace(p.architecture,feedback=replace(p.architecture.feedback,asset_correlation=1.))
        p=DatedM2Pricer(arch,p.premium)
        self.assertLess(independent_equity_fraction(p.layout),1e-25)
