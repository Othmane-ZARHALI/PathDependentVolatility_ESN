"""Independent integral and stock-stepping checks for the post-cap extension."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from .postcap_vov_core import PostKernel,conditional_multiplier
from .stochastic_vov_core import seed_for
from .test_stochastic_vov import python_innovations
from .full_readout_fit import HERE,arr

class PostcapTests(unittest.TestCase):
    def test_conditional_expectation_and_mean_one(self):
        x=np.linspace(-10,10,200001);density=np.exp(-x*x/2)/np.sqrt(2*np.pi)
        for rate in [2.,36.]:
            for gamma in [.25,1.]:
                T=.17;lag=.04;zT=1.3;qU=-np.expm1(-2*rate*(T+lag))
                zU=np.exp(-rate*lag)*zT+np.sqrt(-np.expm1(-2*rate*lag))*x
                truth=np.trapezoid(np.exp(gamma*zU-.5*gamma**2*qU)*density,x)
                self.assertAlmostEqual(truth,float(conditional_multiplier(zT,T,lag,gamma,rate)),places=12)
                qT=-np.expm1(-2*rate*T);mean=np.trapezoid(conditional_multiplier(np.sqrt(qT)*x,T,lag,gamma,rate)*density,x)
                self.assertAlmostEqual(mean,1.,places=12)
                np.testing.assert_allclose(conditional_multiplier(x,T,lag,0,rate),1.)

    def test_stock_step_and_zero_strength_nesting(self):
        with tempfile.TemporaryDirectory(prefix='m2-post-test-',dir='/private/tmp') as temp:
            temp=Path(temp);lib=temp/'test.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'postcap_vov_kernel.cpp'),'-o',str(lib)],check=True)
            d=json.loads((HERE/'runs/m2_maturity_curve_fit_001/run.json').read_text());p=np.array(d['selected']);k=PostKernel(lib,d['ends_years'])
            rng=np.random.default_rng(293);nt=43;n=12;folder=temp/'d29'/'batch000';folder.mkdir(parents=True)
            x=rng.normal(size=(nt,n,3));dw=rng.normal(size=(nt,n))*.01;dt=np.linspace(.0001,.0002,nt);grid=np.r_[0,np.cumsum(dt)]
            for key,value in [('features',x),('equity',dw),('grid',grid)]:np.save(folder/f'{key}.npy',value)
            nn=np.zeros((nt,8));nn[:,0]=.7;rate=8.;seed=seed_for(2,29,0);eps=python_innovations(seed,nt*n).reshape(nt,n)
            for gamma in [0.,.7]:
                k.spec=('test',gamma,rate,30,np.exp(gamma));actual=k.stock(folder,nn,p,False)[0,:,0]
                logs=np.zeros(n);z=np.zeros(n)
                for t in range(nt):
                    eta=p[0]*x[t,:,0]+p[1]*(x[t,:,1]-p[2])**2+p[3]*(x[t,:,2]-p[4])**2
                    h=p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6]
                    L=np.exp(gamma*z-.5*gamma**2*(-np.expm1(-2*rate*grid[t])))
                    v=p[-1]+(p[7]-p[-1])*np.exp(2*h-nn[t,0])*L
                    logs+=-.5*v*dt[t]+np.sqrt(v)*dw[t]
                    z=np.exp(-rate*dt[t])*z+np.sqrt(-np.expm1(-2*rate*dt[t]))*eps[t]
                np.testing.assert_allclose(actual,np.exp(logs),rtol=3e-14,atol=3e-14)
                if gamma==0:
                    old=np.zeros((n,len(p)+1));k.curve_lib.stock_curve(nt,n,len(k.ends),0,arr(x),arr(dw),arr(dt),np.zeros(nt,dtype='int32'),arr(nn),arr(p),old)
                    np.testing.assert_allclose(actual,old[:,0],rtol=3e-14)

if __name__=='__main__':unittest.main()
