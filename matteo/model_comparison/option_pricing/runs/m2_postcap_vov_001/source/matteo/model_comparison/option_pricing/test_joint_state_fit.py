"""Representation checks needed before using the state-repricing bridge."""
import numpy as np
import unittest
from functools import lru_cache

from .m2_spx import frozen_stage3
from .joint_state_core import AdaptiveDatedNormalizer, StateDatedPricer, projected_initial_shift
from .joint_curve_core import CurvePathBatch
from esn_eq.pricing import DatedM2Pricer, ForwardVarianceCurve


@lru_cache(maxsize=1)
def fixture():
    architecture,parameters,premium,_=frozen_stage3()
    p=DatedM2Pricer(architecture,premium)
    u=np.array([.3,-.4,.2]);initial,_=p.initial_state(u)
    return p,parameters,initial


class StateRepricingTests(unittest.TestCase):
    def test_exact_initial_translation_keeps_same_equity_and_gaussian_noise(self):
        fixture = globals()["fixture"]()
        p,parameters,initial=fixture
        times=np.r_[0.,np.geomspace(1e-6,.02,25)]
        a=p.prepare(times,np.zeros_like(initial),12,731)
        b=p.prepare(times,initial,12,731)
        np.testing.assert_array_equal(a.equity_increments,b.equity_increments)
        np.testing.assert_allclose(b.projections,a.projections+projected_initial_shift(p,times[:-1],initial)[:,None,:],rtol=2e-13,atol=2e-13)


    def test_rebuilt_ratio_matches_original_candidate_stock_pricer(self):
        fixture = globals()["fixture"]()
        p,parameters,initial=fixture
        times=np.unique(np.r_[0.,np.geomspace(1e-6,.02,25),.01,.015])
        curve=ForwardVarianceCurve((.01,.015,.03),(.05,.08,.04))
        zero=p.prepare(times,np.zeros_like(initial),12,183)
        shifted=p.prepare(times,initial,12,183)
        normalizer=p.normalizer(times[:-1],initial,parameters,32)
        projections=zero.projections+projected_initial_shift(p,times[:-1],initial)[:,None,:]
        ratio=p.multiplier(projections,parameters)/normalizer[:,None]
        batch=CurvePathBatch(ratio,zero.equity_increments,np.diff(times),np.minimum(np.searchsorted(curve.ends,times[:-1],side='right'),2))
        actual,_=batch.terminal(curve.levels)
        reference=p.spx(shifted,parameters,curve,[times[-1]],order=32).terminal_ratios[times[-1]]
        np.testing.assert_allclose(actual,reference,rtol=3e-13,atol=3e-13)


    def test_candidate_normalizer_and_owner_guards(self):
        fixture = globals()["fixture"]()
        p,parameters,initial=fixture
        g=AdaptiveDatedNormalizer(p,parameters,initial,.05,order=32)
        times=np.geomspace(1e-7,.05,173)
        direct=DatedM2Pricer.normalizer(p,times,initial,parameters,32)
        np.testing.assert_allclose(g(times),direct,rtol=2e-7,atol=0)
        shifted=StateDatedPricer(p.architecture,p.premium);shifted.state_normalizer=g
        with self.assertRaises(ValueError):shifted.normalizer(times,np.zeros_like(initial),parameters,32)
        with self.assertRaises(ValueError):g([.051])


    def test_whitened_coordinates_reproduce_declared_prior_penalty(self):
        fixture = globals()["fixture"]()
        p,parameters,initial=fixture
        L=np.linalg.cholesky(p.projection@p.stationary@p.projection.T)
        for z in [np.zeros(3),np.array([1.,-.3,.7]),np.array([0.,0.,-1.])]:
            x,norm=p.initial_state(L@z)
            np.testing.assert_allclose(p.projection@x,L@z,rtol=1e-13,atol=1e-13)
            np.testing.assert_allclose(norm,z@z,rtol=1e-13,atol=1e-13)

    def test_importance_weights_are_unchanged_by_initial_state(self):
        from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths
        from .spx_volatility_sampling import proposal
        p,parameters,initial=globals()['fixture']()
        times=np.unique(np.r_[0.,np.linspace(0,.04,20),7/365])
        basis=volatility_path_basis(p,times,7/365)
        a=prepare_volatility_paths(p,times,np.zeros_like(initial),16,592,basis,proposal(1))
        b=prepare_volatility_paths(p,times,initial,16,592,basis,proposal(1))
        np.testing.assert_array_equal(a.weights,b.weights)
        np.testing.assert_array_equal(a.cache.equity_increments,b.cache.equity_increments)
        np.testing.assert_allclose(b.cache.projections,a.cache.projections+projected_initial_shift(p,times[:-1],initial)[:,None,:],rtol=1e-12,atol=1e-12)

    def test_native_forward_recovers_common_call_put_iv(self):
        from .plot_joint_curve_fit import black_price,implied_vol
        forward,wrong,k,t,vol=23.,21.,22.,.3,.8
        call=black_price(vol,forward,k,t,'call');put=black_price(vol,forward,k,t,'put')
        self.assertAlmostEqual(implied_vol(call,forward,k,t,'call'),80.,places=9)
        self.assertAlmostEqual(implied_vol(put,forward,k,t,'put'),80.,places=9)
        self.assertGreater(abs(implied_vol(call,wrong,k,t,'call')-implied_vol(put,wrong,k,t,'put')),1.)

if __name__ == "__main__":
    unittest.main()
