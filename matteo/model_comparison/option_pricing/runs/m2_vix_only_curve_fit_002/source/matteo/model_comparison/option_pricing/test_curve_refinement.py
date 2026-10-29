import unittest
import numpy as np
import pandas as pd

from .curve_refinement_core import FlexibleCurvePathBatch,FlexibleLoss
from .joint_curve_core import CurvePathBatch,ExploratoryJointLoss
from esn_eq.market_targets import MarketTargets
from .m2_spx import ROOT


class CurveRefinementTests(unittest.TestCase):
    def test_nested_curve_stock_and_aggregated_derivative(self):
        times=np.array([0.,.05,.15,.25,.4,.5,.75,.8])
        oldseg=np.searchsorted([.25,.75,2.],times[:-1],side='right')
        newseg=np.searchsorted([.15,.25,.5,.75,2.],times[:-1],side='right')
        rng=np.random.default_rng(13);ratio=np.exp(rng.normal(0,.3,(7,16)))
        dw=rng.normal(size=(7,16))*np.sqrt(np.diff(times))[:,None]
        levels=np.array([.04,.08,.03]);embedded=levels[[0,0,1,1,2]]
        a,da=CurvePathBatch(ratio,dw,np.diff(times),oldseg).terminal(levels)
        b,db=FlexibleCurvePathBatch(ratio,dw,np.diff(times),newseg,5).terminal(embedded)
        np.testing.assert_allclose(a,b,rtol=1e-14,atol=1e-14)
        np.testing.assert_allclose(da,np.column_stack([db[:,0]+db[:,1],db[:,2]+db[:,3],db[:,4]]),rtol=1e-13,atol=1e-13)

    def test_five_level_stock_derivative(self):
        rng=np.random.default_rng(19);dt=np.full(30,.01);segment=np.repeat(np.arange(5),6)
        ratio=np.exp(rng.normal(0,.5,(30,8)));dw=rng.normal(0,.1,(30,8))
        batch=FlexibleCurvePathBatch(ratio,dw,dt,segment,5);x=np.array([.04,.06,.03,.08,.07])
        _,jac=batch.terminal(x)
        for j in range(5):
            step=np.eye(5)[j]*1e-7
            finite=(batch.terminal(x+step)[0]-batch.terminal(x-step)[0])/2e-7
            np.testing.assert_allclose(jac[:,j],finite,rtol=2e-7,atol=2e-8)

    def test_embedded_penalty_and_general_loss_jacobian(self):
        t=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
        q=pd.concat([t.select('SPX').head(3),t.select('VIX').head(3)],ignore_index=True)
        f=t.vix_forwards[t.vix_forwards.expiry==q[q.instrument=='VIX'].expiry.iloc[0]]
        old=ExploratoryJointLoss(q,f);new=FlexibleLoss(q,f)
        levels=np.array([.05,.08,.03]);embedded=levels[[0,0,1,1,2]]
        prices=q.price_ask_U.to_numpy()+2;forwards=f.forward_vix_points.to_numpy()+.2
        a=old.evaluate(prices,forwards,levels);b=new.evaluate(prices,forwards,embedded)
        self.assertAlmostEqual(a['total'],b['total'],places=12)
        rng=np.random.default_rng(3);pj=rng.normal(size=(len(q),5));fj=rng.normal(size=(len(f),5))
        jac=new.evaluate(prices,forwards,embedded,pj,fj)['jacobian']
        for j in range(5):
            h=1e-7;d=np.eye(5)[j]*h
            finite=(new.evaluate(prices+h*pj[:,j],forwards+h*fj[:,j],embedded+d)['residuals']-
                    new.evaluate(prices-h*pj[:,j],forwards-h*fj[:,j],embedded-d)['residuals'])/(2*h)
            np.testing.assert_allclose(jac[:,j],finite,rtol=1e-5,atol=1e-6)


if __name__=='__main__':unittest.main()
