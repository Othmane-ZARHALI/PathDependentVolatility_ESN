"""Target-law/Black oracle, off-policy mass accounting, frozen controls."""
import unittest
import numpy as np
from matteo.preprocessing.options import black76_forward
from .spx_sampling_design import call_proposal,tilted_second_moment_ratio
from .spx_call_tilt_study import adjusted,METHODS


class CallTiltTests(unittest.TestCase):
    def draw(self,mixture,n,seed):
        rng=np.random.default_rng(seed);k=rng.choice(len(mixture.probabilities),n,p=mixture.probabilities)
        z=rng.normal(size=(n,4))+mixture.shifts[k]
        return z,mixture.weights(z)

    def assert_mc(self,x,target):
        self.assertLess(abs(x.mean()-target),6*x.std(ddof=1)/np.sqrt(len(x))+1e-12)

    def test_new_tilt_target_gaussian_and_black_oracle(self):
        mixture=call_proposal();z,w=self.draw(mixture,120000,821)
        self.assertTrue(np.all((w>0)&(w<=2)))
        self.assert_mc(w,1.);self.assert_mc(w*z[:,0],0.);self.assert_mc(w*z[:,0]**2,1.)
        stock=np.exp(-.5*.2**2*.5+.2*np.sqrt(.5)*z[:,0])
        self.assert_mc(w*stock,1.)
        self.assert_mc(w*np.maximum(stock-1.3,0),black76_forward(1.,1.3,.5,.2,'call'))

    def test_off_policy_second_moment_subtracts_mass_constant(self):
        old,new=call_proposal(2.),call_proposal(3.)
        z,w=self.draw(old,180000,822);h=1+.1*z[:,0]
        b=np.zeros((len(w),1,5));b[:,0,4]=w-1;beta=np.zeros((1,5));beta[0,4]=7.
        y=(w*(h+7))[:,None]
        ratio=tilted_second_moment_ratio(y,b,beta,w,z,new)[0]
        estimate=ratio*np.mean((w*h)**2)
        z1,w1=self.draw(new,180000,823);direct=(w1*(1+.1*z1[:,0]))**2
        transformed=(w*h)**2*(new.weights(z)/w)
        se=np.sqrt(transformed.var(ddof=1)/len(w)+direct.var(ddof=1)/len(w1))
        self.assertLess(abs(estimate-direct.mean()),6*se)

    def test_constant_payoff_mass_control_preserves_known_value(self):
        _,w=self.draw(call_proposal(),128,824);b=np.zeros((128,2,5));b[:,:,-1]=w[:,None]-1
        data=dict(y=2.5*w[:,None]*np.ones((1,2)),b=b);beta=np.zeros((2,5));beta[:,-1]=2.5
        fits={name:beta for name in METHODS[1:]};x=adjusted(data,fits,np.array([1.,2.]))
        np.testing.assert_allclose(x[:,:,1:],np.broadcast_to(np.array([2.5,5.])[None,:,None],(128,2,4)),atol=1e-14)
