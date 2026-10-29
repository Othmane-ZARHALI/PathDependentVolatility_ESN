"""Independent oracles for the new stochastic response integration and stock step."""
import ctypes as ct
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np

from .stochastic_vov_core import LayerKernel, coefficients, seed_for, SPECS
from .full_readout_fit import HERE,arr


def python_innovations(seed,n):
    mask=(1<<64)-1
    def mix(x):
        x&=mask;x=((x^(x>>30))*0xbf58476d1ce4e5b9)&mask
        x=((x^(x>>27))*0x94d049bb133111eb)&mask;return x^(x>>31)
    out=[]
    for j in range(n):
        a=mix(seed+2*j*0x9e3779b97f4a7c15);b=mix(seed+(2*j+1)*0x9e3779b97f4a7c15)
        u=((a>>11)+.5)*2**-53;v=((b>>11)+.5)*2**-53
        out.append(np.sqrt(-2*np.log(u))*np.cos(2*np.pi*v))
    return np.array(out)


class LayerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='m2-z-test-',dir='/private/tmp');cls.path=Path(cls.temp.name)
        lib=cls.path/'test.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'stochastic_vov_kernel.cpp'),'-o',str(lib)],check=True)
        d=json.loads((HERE/'runs/m2_maturity_curve_fit_001/run.json').read_text());cls.p=np.array(d['selected']);cls.k=LayerKernel(lib,d['ends_years'])
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def test_integrated_polynomial_has_known_gaussian_mean(self):
        # H(A)=3+2A+A^2 and s=0 have exact elementary values, regardless of Z.
        n=32;A=1.25+.75*np.cos(np.pi*(np.arange(n)+.5)/n)
        c=coefficients((3+2*A+A*A)[:,None]);r=self.k.response(c,[1.7],.9,0.,128)
        np.testing.assert_allclose(r,[6.],rtol=3e-14)
        # A(Z)+A(-Z)=2cosh(s tanhZ), independently integrate on a dense real grid.
        z=np.linspace(-9,9,200001);density=np.exp(-z*z/2)/np.sqrt(2*np.pi);s=np.log(2)
        truth=np.trapezoid((3+2*np.exp(s*np.tanh(z))+np.exp(2*s*np.tanh(z)))*density,z)
        actual=self.k.response(c,[0.],1.,s,192)[0];self.assertLess(abs(actual/truth-1),1e-9)

    def test_counter_normals_and_ou_variance(self):
        seed=seed_for(99);out=np.empty(100000);self.k.lib.innovations(seed,len(out),out)
        np.testing.assert_allclose(out[:128],python_innovations(seed,128),atol=5e-15,rtol=5e-15)
        self.assertLess(abs(out.mean()),.012);self.assertLess(abs(out.var()-1),.025)
        for _,_,rate,_,_ in SPECS:
            sd=np.sqrt(-np.expm1(-2*rate*27/365));sample=sd*out
            self.assertLess(abs(sample.var()/(-np.expm1(-2*rate*27/365))-1),.025)

    def test_streamed_stock_matches_python_and_nests_baseline(self):
        rng=np.random.default_rng(287);nt=41;n=12;folder=self.path/'d29'/'batch000';folder.mkdir(parents=True)
        x=rng.normal(size=(nt,n,3));dw=rng.normal(size=(nt,n))*.01;dt=np.linspace(.0001,.0002,nt);grid=np.r_[0,np.cumsum(dt)]
        np.save(folder/'features.npy',x);np.save(folder/'equity.npy',dw);np.save(folder/'grid.npy',grid)
        nn=np.zeros((nt,8));nn[:,0]=.7;seg=np.zeros(nt,dtype='int32');p=self.p;k=self.k;seed=seed_for(2,29,0)
        eps=python_innovations(seed,nt*n).reshape(nt,n)
        for continuation in [0.,.5]:
            k.configure(continuation)
            for sensitivity in [0.,np.log(2)]:
                rate=np.log(2)*365/30;k.spec=('test',sensitivity,rate,30,2)
                actual=k.stock(folder,nn,p,False)[0,:,0];logs=np.zeros(n);z=np.zeros(n)
                for t in range(nt):
                    eta=p[0]*x[t,:,0]+p[1]*(x[t,:,1]-p[2])**2+p[3]*np.exp(sensitivity*np.tanh(z))*(x[t,:,2]-p[4])**2
                    h=p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6]
                    if continuation:h+=continuation*np.log1p(np.logaddexp(0,p[6]*(eta-p[5]))/p[6])
                    v=p[-1]+(p[7]-p[-1])*np.exp(2*h-nn[t,0]);logs+=-.5*v*dt[t]+np.sqrt(v)*dw[t]
                    z=np.exp(-rate*dt[t])*z+np.sqrt(-np.expm1(-2*rate*dt[t]))*eps[t]
                np.testing.assert_allclose(actual,np.exp(logs),rtol=3e-14,atol=3e-14)
                if sensitivity==0:
                    old=np.zeros((n,len(p)+1));k.curve_lib.stock_curve(nt,n,len(k.ends),0,arr(x),arr(dw),arr(dt),seg,arr(nn),arr(p),old)
                    np.testing.assert_allclose(actual,old[:,0],rtol=3e-14)


if __name__=='__main__':unittest.main()
