import subprocess
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scipy.special import ndtr
from .m2_spx import frozen_stage3
from .full_readout_fit import Kernel,arr,HERE,NAMES
from .parity_hedge_controls import FrozenHedgeKernel,correction_at_strikes
from .parity_stock import StockValues
from esn_eq.pricing import DatedM2Pricer,ForwardVarianceCurve
from esn_eq.config import ReadoutParameters

class FrozenHedgeTests(unittest.TestCase):
    def test_native_predictable_hedges_match_independent_python(self):
        with tempfile.TemporaryDirectory() as tmp:
            c=Path(tmp);lib=c/'hedge.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'parity_hedge_kernel.cpp'),'-o',str(lib)],check=True)
            a,pars,premium,_=frozen_stage3();p=DatedM2Pricer(a,premium);initial,_=p.initial_state([1,0,0])
            x=np.array([getattr(pars,n) for n in NAMES[:7]]+[.06,.07,.065,.0001]);x[4]=.3
            grid=np.array([0,.00001,.00002,.001,.002,.25,.25001,.4,.75,.9]);cache=p.prepare(grid,initial,7,185)
            for name,data in [('features',cache.projections),('equity',cache.equity_increments),('grid',grid)]:np.save(c/f'{name}.npy',data)
            k=Kernel(lib)
            def norm(times):
                from esn_eq.pricing import _normal_rule
                from esn_eq.states import symmetric_square_root
                result=[]
                for t in times:
                    m,v=p.gaussian_projection_law(float(t),initial);z,w=_normal_rule(12)
                    e=k.expectation(m,arr(z@symmetric_square_root(v).T),arr(w),x)[0]
                    result.append(np.r_[np.log(e[0]),e[1:]/e[0]])
                return np.array(result)
            strikes=np.array([.7,1.,1.3]);got=FrozenHedgeKernel(lib).evaluate(c,norm,x,strikes)
            dt=np.diff(grid);seg=np.minimum(np.searchsorted([.25,.75],grid[:-1],side='right'),2);xi=x[7+seg]
            normal=norm(grid[:-1]);mult=k.multiplier(cache.projections,x)[:,0].reshape(len(dt),7)
            var=x[10]+(xi[:,None]-x[10])*mult/np.exp(normal[:,0,None])
            total=np.sum(xi*dt);remain=total;logs=np.zeros(7);last=np.ones(7);delta=np.zeros((7,3,2));hedge=np.zeros_like(delta);next_time=0.
            for j,t in enumerate(grid[:-1]):
                if t+1e-14>=next_time:
                    current=np.exp(logs);hedge+=delta*(current-last)[:,None,None]
                    for axis in range(2):
                        v=np.full(7,remain) if axis==0 else remain*var[j]/xi[j]
                        delta[:,:,axis]=ndtr((logs[:,None]-np.log(strikes)+.5*v[:,None])/np.sqrt(v[:,None]))
                    last=current;next_time=t+1/(252*4)
                logs+=-.5*var[j]*dt[j]+np.sqrt(var[j])*cache.equity_increments[j];remain-=xi[j]*dt[j]
            stock=np.exp(logs);hedge+=delta*(stock-last)[:,None,None]
            np.testing.assert_allclose(got[:,0],stock,rtol=2e-13)
            price_lib=c/'price.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'parity_stock_kernel.cpp'),'-o',str(price_lib)],check=True)
            np.testing.assert_allclose(StockValues(price_lib)(c,normal,x),stock,rtol=2e-13)
            np.testing.assert_allclose(got[:,1:].reshape(7,3,2),hedge,rtol=2e-12,atol=2e-13)
    def test_controls_interpolate_across_forward_without_side_switch(self):
        rng=np.random.default_rng(39);c=rng.normal(size=(5,3,4));beta=rng.normal(size=(3,4));anchors=[.8,1.,1.2]
        left=correction_at_strikes(c,beta,anchors,[1-1e-9]);right=correction_at_strikes(c,beta,anchors,[1+1e-9])
        np.testing.assert_allclose(left,right,atol=5e-8)

if __name__=='__main__':unittest.main()
