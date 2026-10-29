"""Independent target-expectation oracle for frozen full/half pilot coefficients."""
import unittest
import numpy as np
from .spx_volatility_sampling import fit_variants, adjusted, proposal


class VolatilitySamplingDriverTests(unittest.TestCase):
    def test_weighted_quadratic_exact_for_all_pilot_halves_and_new_draws(self):
        mix=proposal(1);rng=np.random.default_rng(170)
        def sample(n):
            comp=rng.choice(len(mix.shifts),n,p=mix.probabilities)
            z=rng.standard_normal((n,4))+mix.shifts[comp]
            w=mix.weights(z)
            y=np.c_[z[:,0]**2+z[:,2],2*z[:,0]**2-z[:,2]]
            c=np.broadcast_to(np.c_[z[:,0]**2-1,z[:,2],np.zeros((n,2))][:,None,:],(n,2,4)).copy()
            return y,c,w
        y,c,w=sample(1024)
        fits=fit_variants(y,c,w)
        y,c,w=sample(4096)
        x=adjusted(y,c,w,fits,np.array([3.,5.]))
        for variant in range(1,4):
            np.testing.assert_allclose(x[:,:,variant],np.broadcast_to([3.,10.],(len(x),2)),atol=1e-12)
        self.assertGreater(x[:,:,0].var(),1.)

    def test_small_end_to_end_collection_includes_pilot_and_coupled_checks(self):
        from types import SimpleNamespace
        import pandas as pd
        from .spx_volatility_sampling import collect
        from .fine_grid import CachedNormalizerPricer
        from esn_eq.config import ArchitectureConfig, BankConfig, RiskPremium, ReadoutParameters
        from esn_eq.pricing import ForwardVarianceCurve
        from esn_eq.volatility_path_importance import volatility_path_basis
        p=CachedNormalizerPricer(ArchitectureConfig(
            feedback=BankConfig((1.,2.),.3,True),spike=BankConfig((3.,4.),.2,True),cluster=BankConfig((.5,1.),0.,False)),RiskPremium())
        grids=[np.array([0.,.1]),np.array([0.,.05,.1]),np.array([0.,.025,.05,.075,.1])]
        basis=volatility_path_basis(p,grids[-1],.05)
        quotes=pd.DataFrame(dict(T_years=[.1,.1],K_over_F=[.9,1.1],option_type=['put','call']))
        parameters=ReadoutParameters(cluster_loading=0,feedback_curvature=0,spike_curvature=0)
        args=SimpleNamespace(pilot=8,selection=8,paths=8,batch_size=8,seed=182)
        initial=np.zeros(len(p.layout.rates));curve=ForwardVarianceCurve((.1,),(.04,))
        for role in (0,2):
            data,w,z,comp,monitor,seeds,costs=collect(p,parameters,curve,initial,quotes,grids,basis,proposal(1),args,1,role,0)
            self.assertEqual(set(data),{2} if role==0 else {0,1,2})
            self.assertEqual(data[2][0].shape,(8,2))
            self.assertEqual(data[2][1].shape,(8,2,4))
            self.assertEqual(monitor.shape,(8,5))
            self.assertEqual(z.shape,(8,4))
            self.assertEqual(len(seeds),1)
            self.assertEqual(len(costs),1 if role==0 else 3)
            self.assertTrue(np.isfinite(monitor).all())
