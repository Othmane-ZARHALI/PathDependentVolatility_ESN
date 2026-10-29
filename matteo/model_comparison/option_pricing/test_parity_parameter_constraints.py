import unittest
from types import SimpleNamespace
import numpy as np
import pandas as pd
from .parity_parameter_constraints import constrained_local,anchor_criteria

class ConstrainedMeritTests(unittest.TestCase):
    def test_both_conventions_criteria_and_merit_derivatives(self):
        quote=pd.DataFrame(dict(price_bid_U=[0.]*4,price_ask_U=[0.]*4,price_scale_U=[1.]*4))
        vq=pd.DataFrame(dict(expiry_days=[55,55],option_type=['put','call'],strike=[10.,11.],price_bid_U=[0.,0.],price_ask_U=[0.,0.],price_scale_U=[1.,1.]))
        class Loss:
            sqrtweights=np.ones(4)
            def evaluate(self,price,forward,levels):
                residual=np.r_[price,(forward-10)/.5,np.sqrt(.1)*np.diff(levels)]
                return dict(total=float(residual@residual),residuals=residual,components=dict(SPX=float(np.sum(price[:2]**2)),VIX=float(np.sum(price[2:]**2)),forwards=float(np.sum(((forward-10)/.5)**2)),state_penalty=0.))
        def vix(p):
            prices=np.array([p[0]**2+2*p[1],np.exp(p[2])+.2*p[3]]);jac=np.zeros((2,11));jac[0,0]=2*p[0];jac[0,1]=2;jac[1,2]=np.exp(p[2]);jac[1,3]=.2
            fj=np.zeros((1,11));fj[0,4]=1
            return dict(prices=prices,jacobian=jac,forwards=np.array([10+p[4]]),forward_jacobian=fj,normalizer={})
        loss=Loss();e=SimpleNamespace(spx=[0,1],vix=vq,quotes=quote,loss=loss,otm_loss=loss,vix_prediction=vix)
        x=np.array([.1,.2,.3,.15,.2,5,2,.05,.06,.055,.001]);rng=np.random.default_rng(72)
        anchor=dict(prices=np.array([1.,2.,3.,4.]),price_jacobian=rng.normal(size=(4,11)))
        r=constrained_local(e,x,x,anchor);h=1e-6
        for j in range(11):
            plus=x.copy();minus=x.copy();plus[j]+=h;minus[j]-=h
            a=constrained_local(e,plus,x,anchor);b=constrained_local(e,minus,x,anchor)
            np.testing.assert_allclose((a['merit']-b['merit'])/(2*h),r['merit_gradient'][j],rtol=2e-6,atol=1e-8)
            np.testing.assert_allclose((a['criteria']-b['criteria'])/(2*h),r['criteria_gradients'][:,j],rtol=2e-6,atol=1e-8)
        p=vix(x)['prices'];f=vix(x)['forwards'][0];native=p.copy();native[0]-=f-10
        np.testing.assert_allclose(r['native_vix'],np.sum(native*native)+((f-10)/.5)**2)
    def test_anchor_criteria_order(self):
        a=dict(total=10,original_otm_total=12,components=dict(VIX=3,forwards=2),original_otm_components=dict(VIX=4,forwards=2))
        np.testing.assert_array_equal(anchor_criteria(a),[10,12,5,6])

if __name__=='__main__':unittest.main()
