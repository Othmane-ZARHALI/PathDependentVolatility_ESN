"""Small full-stack MLMC correction collection and analytic Black target."""
import unittest
from types import SimpleNamespace
import numpy as np
import pandas as pd
from scipy.special import ndtr
from .spx_mlmc import collect,fit_controls,draw_estimators
from .fine_grid import CachedNormalizerPricer
from esn_eq.config import ArchitectureConfig,BankConfig,RiskPremium,ReadoutParameters
from esn_eq.pricing import ForwardVarianceCurve
from esn_eq.mlmc import independent_level_statistics
from esn_eq.volatility_path_importance import volatility_path_basis


class SPXMLMCTests(unittest.TestCase):
    def test_constant_volatility_zero_corrections_and_black_expectation(self):
        p=CachedNormalizerPricer(ArchitectureConfig(feedback=BankConfig((1.,2.),.3,True),
            spike=BankConfig((3.,4.),.2,True),cluster=BankConfig((.5,1.),0.,False)),RiskPremium())
        grids=[np.array([0.,.05,.1]),np.array([0.,.025,.05,.075,.1])]
        bases=[volatility_path_basis(p,g,.05) for g in grids]
        initial=np.zeros(len(p.layout.rates));curve=ForwardVarianceCurve((.1,),(.04,))
        parameters=ReadoutParameters(cluster_loading=0,feedback_curvature=0,spike_curvature=0)
        quotes=pd.DataFrame(dict(T_years=[.1,.1],K_over_F=[.95,1.05],option_type=['put','call']))
        args=SimpleNamespace(seed=183,batch_size=20000)
        base,_,_=collect(p,parameters,curve,initial,quotes,grids,bases,0,20000,args,2,0)
        fine,_,cost=collect(p,parameters,curve,initial,quotes,grids,bases,1,128,args,2,0)
        direct,_,_=collect(p,parameters,curve,initial,quotes,grids,bases,1,128,args,3,0,direct=True)
        np.testing.assert_allclose(fine['y'],0.,atol=1e-14)
        self.assertEqual(fine['c'].shape,(128,2,5));self.assertEqual(base['monitors'].shape,(20000,5))
        self.assertGreater(cost[0]['seconds'],cost[0]['fine_seconds'])
        self.assertTrue(np.any(direct['y']>0))
        mean,se=independent_level_statistics([base['y'],fine['y']])
        sigma=np.sqrt(.04*.1)
        for j,(k,kind) in enumerate(zip(quotes.K_over_F,quotes.option_type)):
            d1=(-np.log(k)+.5*sigma**2)/sigma;d2=d1-sigma
            exact=ndtr(d1)-k*ndtr(d2) if kind=='call' else k*ndtr(-d2)-ndtr(-d1)
            self.assertLess(abs(mean[j]-exact),6*se[j])
        beta=fit_controls(base['y'],base['c'])
        x=draw_estimators(direct['y'],direct['c'],[beta],np.ones(2))
        self.assertEqual(x.shape,(128,2,2));self.assertTrue(np.isfinite(x).all())
