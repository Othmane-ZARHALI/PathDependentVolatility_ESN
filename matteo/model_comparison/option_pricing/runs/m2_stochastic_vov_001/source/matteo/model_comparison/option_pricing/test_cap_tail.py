"""Independent Python equations, quadrature and adapted stock-step correspondence."""
import ctypes as ct,json
import numpy as np
from .full_readout_fit import HERE,Kernel,arr
from .maturity_curve import CurveKernel
from esn_eq.pricing import _normal_rule
from esn_eq.states import symmetric_square_root


def reference(x,p,mode=0,a=0.):
    eta=p[0]*x[...,0]+p[1]*(x[...,1]-p[2])**2+p[3]*(x[...,2]-p[4])**2
    def h(z,L):return L-np.logaddexp(0,p[6]*(L-z))/p[6]
    if mode==1:
        gap=p[5]-eta;hh=p[5]-.5*(gap+np.hypot(gap,2*np.log(2)/p[6]))
    elif mode==2:hh=h(eta,p[5])+a*np.log1p(np.logaddexp(0,p[6]*(eta-p[5]))/p[6])
    elif mode==3:hh=p[0]*x[...,0]+h(p[1]*(x[...,1]-p[2])**2,3.)+h(p[3]*(x[...,2]-p[4])**2,2.5)
    elif mode==4:
        w=.5+.45*np.tanh(.5*x[...,0]);return (1-w)*np.exp(2*h(eta,p[5]))+w*np.exp(2*h(eta+a*p[3]*(x[...,2]-p[4])**2,p[5]))
    else:hh=h(eta,p[5])
    return np.exp(2*hh)


def test_experimental_equations_and_tensor_quadrature():
    out=HERE/'runs/m2_cap_tail_001';p=np.array(json.loads((out/'design.json').read_text())['baseline']);rng=np.random.default_rng(392)
    for lib,modes in [(out/'cap.dylib',[(0,0),(1,0),(2,.25),(2,.5),(3,0)]),(out/'cluster_vov/cluster.dylib',[(4,0.),(4,.5),(4,1.)])]:
        k=Kernel(lib);k.lib.configure.argtypes=[ct.c_int,ct.c_double]
        for mode,a in modes:
            k.lib.configure(mode if mode<4 else 0,a);x=np.r_[rng.normal(size=(100,3)),[[0,0,0],[0,20,20],[0,-20,-20],[1000,0,0],[-1000,0,0]]]
            np.testing.assert_allclose(k.multiplier(x,p)[:,0],reference(x,p,mode,a),rtol=3e-13,atol=1e-13)
            means=rng.normal(size=(34,3));cov=np.array([[.4,0,0],[0,.6,-.15],[0,-.15,1.]])
            z,w=_normal_rule(32);o=arr(z@symmetric_square_root(cov).T);pred=k.expectation(means,o,arr(w),p)[:,0]
            ids=[0,10,33];expected=np.array([reference(m+o,p,mode,a)@w for m in means[ids]])
            np.testing.assert_allclose(pred[ids],expected,rtol=1e-6,atol=1e-10)


def test_stock_step_and_no_extra_noise():
    out=HERE/'runs/m2_cap_tail_001';p=np.array(json.loads((out/'design.json').read_text())['baseline']);nc=16;rng=np.random.default_rng(131);nt=20;npth=17
    x=arr(rng.normal(size=(nt,npth,3)));dt=np.full(nt,.001);dw=arr(rng.normal(size=(nt,npth))*np.sqrt(dt[:,None]));seg=np.zeros(nt,dtype='int32');norm=np.zeros((nt,8));norm[:,0]=np.log(30.)
    for lib,mode,a in [(out/'cap.dylib',2,.25),(out/'cluster_vov/cluster.dylib',4,.5)]:
        k=CurveKernel(lib,lib,np.arange(1,17));k.lib.configure.argtypes=[ct.c_int,ct.c_double];k.lib.configure(mode if mode<4 else 0,a);result=np.zeros((npth,len(p)+1))
        k.curve_lib.stock_curve(nt,npth,nc,0,x,dw,arr(dt),seg,arr(norm),arr(p),result)
        v=p[-1]+(p[7]-p[-1])*reference(x,p,mode,a)/30
        expected=np.exp(np.sum(-.5*v*dt[:,None]+np.sqrt(v)*dw,axis=0));np.testing.assert_allclose(result[:,0],expected,rtol=1e-13)
