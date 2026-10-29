import unittest
import numpy as np
import pandas as pd
from .parity_adaptive import (canonical_market, tier_weights, sample_coefficients,
    iid_statistics, scrambled_statistics, reference_calls, gaussian_controls)

class ParityAdaptiveTests(unittest.TestCase):
    def test_otm_conversion_preserves_spreads_and_parity(self):
        q=pd.DataFrame(dict(option_type=['put','call'],forward=[100.,100.],strike=[90.,110.],
                            price_bid_U=[1.,2.],price_mid_U=[1.5,2.5],price_ask_U=[2.,3.]))
        c=canonical_market(q)
        np.testing.assert_allclose(c.price_mid_U,[11.5,2.5])
        np.testing.assert_allclose(c.price_ask_U-c.price_bid_U,q.price_ask_U-q.price_bid_U)
        self.assertTrue(c.option_type.eq('call').all());self.assertTrue(q.option_type.iloc[0]=='put')
    def test_budget_continuity_and_increasing_effective_paths(self):
        d=np.linspace(0,4,1001);a=sample_coefficients(d)
        np.testing.assert_allclose(a.sum(0),1,atol=1e-13)
        self.assertTrue(np.all(np.diff(1/(a*a).sum(0))>=-1e-9))
        for z in [1.,2.,3.]:
            np.testing.assert_allclose(tier_weights(z-1e-7),tier_weights(z+1e-7),atol=1e-12)
        self.assertEqual(np.sum(a[:,0]>0),1024);self.assertEqual(np.sum(a[:,-1]>0),4096)
    def test_shared_prefix_variance(self):
        # Exact variance for .5*mean(first 4)+.5*mean(first 8):
        # 4*(3/16)^2+4*(1/16)^2 = 5/32, not independent-prefix 3/32.
        a=sample_coefficients([1.5],(4,8,16))[:,0]
        self.assertAlmostEqual(np.sum(a*a),5/32)
        rng=np.random.default_rng(718);y=rng.normal(size=(16000,16))
        self.assertAlmostEqual(np.var(y@a,ddof=1),5/32,delta=.005)
    def test_iid_se_replays_and_deterministic_parity(self):
        rng=np.random.default_rng(18);y=rng.normal(size=(16,3));d=[0.,1.5,3.]
        s=iid_statistics(y,d,(4,8,16))
        a=sample_coefficients(d,(4,8,16));np.testing.assert_allclose(s['price'],(a*y).sum(0))
        np.testing.assert_allclose(s['se'],[np.std(y[:4,0],ddof=1)/2,np.std(y[:8,1],ddof=1)*np.sqrt(5/32),np.std(y[:,2],ddof=1)/4])
        call=s['price'];forward=100.;strikes=np.array([90.,100.,110.]);put=call-forward+strikes
        np.testing.assert_allclose(call-put,forward-strikes,atol=1e-14)
    def test_scramble_se_includes_whole_prefix_blend(self):
        y=np.arange(4*16*3,dtype=float).reshape(4,16,3)
        s=scrambled_statistics(y,[0,1.5,3],(16,32,64))
        a=sample_coefficients([0,1.5,3],(4,8,16));b=np.einsum('nq,rnq->rq',a,y)
        np.testing.assert_allclose(s['price'],b.mean(0));np.testing.assert_allclose(s['se'],b.std(0,ddof=1)/2)
    def test_lognormal_oracle_is_exact_on_every_path(self):
        rng=np.random.default_rng(4);r=np.exp(-.02+.2*rng.normal(size=32));dr=r[:,None]*rng.normal(size=(32,3))
        y,j=reference_calls(r,dr,r,dr,np.linspace(.2,2,32),[.6,1.,1.4],.04,np.zeros(3))
        np.testing.assert_allclose(y,np.broadcast_to(y[0],y.shape),atol=1e-15);np.testing.assert_allclose(j,0,atol=1e-15)
    def test_reference_derivative_includes_moving_control(self):
        rng=np.random.default_rng(8);z=rng.normal(size=30);x=rng.normal(size=30);w=np.exp(.1*z);ks=np.array([.7,1.,1.3])
        def f(a):
            r=np.exp(-a+.5*np.sqrt(a)*x);dr=r*(-1+.25*x/np.sqrt(a))
            ref=np.exp(-.5*a+np.sqrt(a)*z);dref=ref*(-.5+.5*z/np.sqrt(a))
            return reference_calls(r,dr[:,None],ref,dref[:,None],w,ks,a,np.ones(1))
        a=.13;h=1e-6;y,j=f(a);fd=(f(a+h)[0]-f(a-h)[0])/(2*h)
        np.testing.assert_allclose(j[:,:,0],fd,rtol=2e-6,atol=2e-9)
    def test_gaussian_controls_known_zero_means(self):
        # Deterministic symmetric points exactly match mean and variance.
        z=np.array([[-1.,1.],[1.,-1.]])
        c=gaussian_controls(z,np.ones(2));np.testing.assert_allclose(c.mean(0),0)

if __name__=='__main__':unittest.main()
