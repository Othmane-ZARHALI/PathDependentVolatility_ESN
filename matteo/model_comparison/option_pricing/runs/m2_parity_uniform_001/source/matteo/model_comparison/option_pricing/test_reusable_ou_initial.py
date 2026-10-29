import unittest
import numpy as np
from .m2_spx import frozen_stage3
from .reusable_ou_initial import state_translation,projection_translation
from esn_eq.pricing import DatedM2Pricer

class InitialTranslationTests(unittest.TestCase):
    def test_translation_matches_fresh_simulation_with_same_increments(self):
        a,_,q,_=frozen_stage3();p=DatedM2Pricer(a,q)
        old,_=p.initial_state([1,0,0]);new,_=p.initial_state([.4,-.3,.7])
        times=np.array([0.,.0002,.017,.1,.6]);checkpoints=times[1:]
        x=p.prepare(times,old,19,7183,checkpoints);y=p.prepare(times,new,19,7183,checkpoints)
        np.testing.assert_array_equal(x.equity_increments,y.equity_increments)
        shift=projection_translation(p,times[:-1],old,new)
        np.testing.assert_allclose(x.projections+shift[:,None,:],y.projections,atol=3e-14,rtol=3e-14)
        for t in checkpoints:
            np.testing.assert_allclose(x.checkpoint_states[t]+state_translation(p,t,old,new),y.checkpoint_states[t],atol=3e-14,rtol=3e-14)
    def test_importance_weights_and_scores_are_initial_state_independent(self):
        from .spx_volatility_sampling import proposal
        from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths
        a,_,q,_=frozen_stage3();p=DatedM2Pricer(a,q)
        old,_=p.initial_state([1,0,0]);new,_=p.initial_state([.4,-.3,.7])
        times=np.array([0.,.005,7/365,.04,.08]);basis=volatility_path_basis(p,times,7/365)
        x=prepare_volatility_paths(p,times,old,23,716,basis,proposal(1))
        y=prepare_volatility_paths(p,times,new,23,716,basis,proposal(1))
        np.testing.assert_array_equal(x.weights,y.weights)
        np.testing.assert_array_equal(x.score,y.score)
        np.testing.assert_array_equal(x.components,y.components)
        np.testing.assert_array_equal(x.cache.equity_increments,y.cache.equity_increments)
        np.testing.assert_allclose(x.cache.projections+projection_translation(p,times[:-1],old,new)[:,None,:],y.cache.projections,atol=3e-14,rtol=3e-14)
    def test_q_drift_and_covariance_do_not_change(self):
        a,_,q,_=frozen_stage3();p=DatedM2Pricer(a,q);old,_=p.initial_state([1,0,0]);new,_=p.initial_state([0,.2,-.5])
        for t in [0.,.02,.8]:
            m0,c0=p.gaussian_projection_law(t,old);m1,c1=p.gaussian_projection_law(t,new)
            np.testing.assert_array_equal(c0,c1);np.testing.assert_allclose(m1-m0,projection_translation(p,t,old,new),atol=2e-14)

if __name__=='__main__':unittest.main()
