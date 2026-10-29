"""End-to-end reflected/ordinary correction and actual-timing correspondence."""
import unittest
from types import SimpleNamespace
import numpy as np
import pandas as pd
from .spx_bridge_coupling import collect,fit_variants
from .fine_grid import CachedNormalizerPricer
from esn_eq.config import ArchitectureConfig,BankConfig,RiskPremium,ReadoutParameters
from esn_eq.pricing import ForwardVarianceCurve
from esn_eq.bridge_coupling import EndpointBridge
from esn_eq.volatility_path_importance import volatility_path_basis


class BridgeDriverTests(unittest.TestCase):
    def test_constant_volatility_correction_and_frozen_fit_shapes(self):
        p=CachedNormalizerPricer(ArchitectureConfig(feedback=BankConfig((1.,2.),.3,True),
            spike=BankConfig((3.,4.),.2,True),cluster=BankConfig((.5,1.),0.,False)),RiskPremium())
        grids=[np.array([0.,.05,.1]),np.array([0.,.025,.05,.075,.1])]
        basis=volatility_path_basis(p,grids[-1],.05);kernel=EndpointBridge(p,grids[-1],basis)
        quotes=pd.DataFrame(dict(T_years=[.1,.1],K_over_F=[.95,1.05],option_type=['put','call']))
        parameters=ReadoutParameters(cluster_loading=0,feedback_curvature=0,spike_curvature=0)
        args=SimpleNamespace(pilot=16,paths=16,batch_size=8,seed=183)
        for role in (0,1):
            data,seeds,costs=collect(p,parameters,ForwardVarianceCurve((.1,),(.04,)),np.zeros(len(p.layout.rates)),quotes,grids,basis,kernel,args,1,role,0)
            np.testing.assert_allclose(data['plain_y'],0.,atol=1e-13)
            np.testing.assert_allclose(data['bridge_y'],0.,atol=1e-13)
            self.assertEqual(data['monitors'].shape,(16,7));self.assertEqual(len(seeds),2)
            self.assertTrue(all(c['bridge_seconds']>0 and c['ordinary_seconds']>0 for c in costs))
            self.assertEqual({c['execution_order'] for c in costs},{'ordinary_first','bridge_first'})
            self.assertEqual(fit_variants(data,'bridge').shape,(3,2,5))
