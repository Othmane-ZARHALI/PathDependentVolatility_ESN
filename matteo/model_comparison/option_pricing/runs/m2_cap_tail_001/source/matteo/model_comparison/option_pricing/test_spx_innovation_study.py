"""Full collection, centered correction mass and frozen ridge correspondence."""
import unittest
from types import SimpleNamespace
import numpy as np
import pandas as pd
from .spx_innovation_study import collect,datasets,train,adjusted
from .spx_hedging_study import refinement
from .fine_grid import CachedNormalizerPricer
from esn_eq.config import ArchitectureConfig,BankConfig,RiskPremium,ReadoutParameters
from esn_eq.pricing import ForwardVarianceCurve
from esn_eq.volatility_path_importance import volatility_path_basis


class InnovationStudyTests(unittest.TestCase):
    def test_constant_readout_new_columns_zero_and_ridge_estimates_unchanged(self):
        p=CachedNormalizerPricer(ArchitectureConfig(feedback=BankConfig((1.,2.),.3,True),
            spike=BankConfig((3.,4.),.2,True),cluster=BankConfig((.5,1.),0.,False)),RiskPremium())
        grids=[np.array([0.,.05,.1])]
        for _ in range(3):grids.append(refinement(grids[-1]))
        basis=volatility_path_basis(p,grids[-1],.05)
        quotes=pd.DataFrame(dict(T_years=[.1,.1],K_over_F=[.95,1.05],option_type=['put','call']))
        args=SimpleNamespace(pilot=16,paths=16,batch_size=8,seed=184,direct_only=False,linear_only=False)
        params=ReadoutParameters(cluster_loading=0,feedback_curvature=0,spike_curvature=0)
        data,seeds,costs=collect(p,params,ForwardVarianceCurve((.1,),(.04,)),np.zeros(len(p.layout.rates)),quotes,grids,basis,args,0,0)
        self.assertEqual(len(seeds),2);self.assertEqual({c['execution_order'] for c in costs},{'ordinary_first','augmented_first'})
        for name,(y,b,v) in datasets(data).items():
            np.testing.assert_array_equal(v,0.)
            if name.startswith('corr'):
                np.testing.assert_allclose(y,0.,atol=1e-13)
                np.testing.assert_array_equal(b[:,:,-1],np.broadcast_to(data['w'][:,None]-1.,(16,2)))
            fits=train(y,b,v);x=adjusted(y,b,v,fits,np.ones(2))
            np.testing.assert_allclose(x[:,:,2],x[:,:,3],atol=1e-12)
            self.assertEqual(x.shape,(16,2,8))

    def test_linear_direct_collection(self):
        p=CachedNormalizerPricer(ArchitectureConfig(feedback=BankConfig((1.,2.),.3,True),
            spike=BankConfig((3.,4.),.2,True),cluster=BankConfig((.5,1.),0.,False)),RiskPremium())
        grids=[np.array([0.,.05,.1])]
        for _ in range(3):grids.append(refinement(grids[-1]))
        basis=volatility_path_basis(p,grids[-1],.05)
        quotes=pd.DataFrame(dict(T_years=[.1,.1],K_over_F=[.95,1.05],option_type=['put','call']))
        args=SimpleNamespace(pilot=16,paths=16,batch_size=8,seed=185,direct_only=True,linear_only=True)
        data,_,costs=collect(p,ReadoutParameters(),ForwardVarianceCurve((.1,),(.04,)),np.zeros(len(p.layout.rates)),quotes,grids,basis,args,0,0)
        self.assertEqual(list(datasets(data)),['price512'])
        self.assertEqual(data['g3_vol'].shape,(16,2,3))
        self.assertNotIn('ordinary_grid2_seconds',costs[0])
        y,b,v=datasets(data)['price512'];fits=train(y,b,v)
        self.assertEqual(fits['augmented'].shape,(2,8))
        self.assertEqual(adjusted(y,b,v,fits,np.ones(2)).shape,(16,2,8))
