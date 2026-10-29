"""Numerical checks for the retained experimental reference components."""
import unittest
import numpy as np
from matteo.model_comparison.option_pricing.current.two_speed import exponential, two_speed_response
from matteo.model_comparison.option_pricing.current.rotation import stock


class ReferenceTests(unittest.TestCase):
    def test_two_speed_derivatives(self):
        rng = np.random.default_rng(42)
        t = 55 / 365
        x = np.array([.2]*7+[.4,.8,.5,.9,.7,.6,.65])
        d = {'times': t + np.array([.005,.03,.07]), 'a': {'z_normals': rng.normal(size=12), 'fast_normals': rng.normal(size=12)}, 'ratios': rng.uniform(.8,1.2,(3,12))}
        value, derivatives = two_speed_response(d,t,x)
        self.assertTrue(np.all(value > 0))
        for k, deriv in zip([7,8,9,11,12,13], derivatives):
            h=1e-6; xp=x.copy(); xm=x.copy(); xp[k]+=h; xm[k]-=h
            finite=(two_speed_response(d,t,xp)[0]-two_speed_response(d,t,xm)[0])/(2*h)
            np.testing.assert_allclose(deriv,finite,rtol=1e-5,atol=1e-8)

    def test_rotation(self):
        stats=np.array([[.04,.1,.03,-.02],[.09,-.2,.04,.05]])
        np.testing.assert_allclose(stock(stats,0,.5,.5),np.exp(-.5*stats[:,0]+stats[:,1]))
        theta = .8
        cs, cf = np.sqrt(.34), np.sqrt(.66)
        expected = np.exp(-.5*stats[:,0] + np.sqrt(1-theta**2)*stats[:,1] - theta*(cs*stats[:,2]+cf*stats[:,3]))
        np.testing.assert_allclose(stock(stats,theta,cs,cf),expected)

if __name__ == '__main__': unittest.main()
