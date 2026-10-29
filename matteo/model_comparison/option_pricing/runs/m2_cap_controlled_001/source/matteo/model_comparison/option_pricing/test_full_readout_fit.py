"""Independent production correspondence for compiled readout/stock derivatives."""
import tempfile
import subprocess
from pathlib import Path
import unittest
import numpy as np
from .full_readout_fit import Kernel,checks,NAMES,arr,HERE
from .m2_spx import frozen_stage3
from esn_eq.pricing import DatedM2Pricer,ForwardVarianceCurve
from esn_eq.config import ReadoutParameters

class KernelCorrespondence(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.temp=tempfile.TemporaryDirectory();cls.folder=Path(cls.temp.name)
  lib=cls.folder/'kernel.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'full_readout_kernel.cpp'),'-o',str(lib)],check=True)
  cls.k=Kernel(lib);a,p,r,_=frozen_stage3();cls.dated=DatedM2Pricer(a,r);cls.pars=ReadoutParameters(**(p.__dict__ if hasattr(p,'__dict__') else {n:getattr(p,n) for n in NAMES[:7]})|{'spike_shift':.2})
  cls.x=np.array([getattr(cls.pars,n) for n in NAMES[:7]]+[.06,.07,.065,.0001])
 @classmethod
 def tearDownClass(cls):cls.temp.cleanup()
 def test_multiplier_derivatives(self):self.assertTrue(checks(self.k,self.x)['multiplier_matches_production'])
 def test_gaussian_expectation_matches_production(self):
  from esn_eq.pricing import _normal_rule
  from esn_eq.states import symmetric_square_root
  initial,_=self.dated.initial_state([1,0,0]);m,c=self.dated.gaussian_projection_law(.17,initial);z,w=_normal_rule(12)
  got=self.k.expectation(m,arr(z@symmetric_square_root(c).T),arr(w),self.x)[0,0]
  expected=self.dated.gaussian_expectation(m,c,self.pars,12)
  np.testing.assert_allclose(got,expected,rtol=3e-14)
 def test_stock_and_all_eleven_derivatives(self):
  p=self.dated;initial,_=p.initial_state([1,0,0]);grid=np.array([0,.0001,.002,.02,.25,.4,.75,.9]);cache=p.prepare(grid,initial,8,113,checkpoints=[.4,.9]);c=self.folder/'cache';c.mkdir(exist_ok=True)
  for n,a in [('features',cache.projections),('equity',cache.equity_increments),('grid',grid),('events',np.array([5,7],dtype='int32'))]:np.save(c/(n+'.npy'),a)
  from esn_eq.pricing import _normal_rule
  from esn_eq.states import symmetric_square_root
  def run(x):
   normal=[]
   for t in grid[:-1]:
    m,cov=p.gaussian_projection_law(t,initial);z,w=_normal_rule(12);a=self.k.expectation(m,arr(z@symmetric_square_root(cov).T),arr(w),x)[0];normal.append(np.r_[np.log(a[0]),a[1:]/a[0]])
   return self.k.stock(c,arr(normal),x)
  result=run(self.x);prod=p.spx(cache,self.pars,ForwardVarianceCurve((.25,.75,1.),tuple(self.x[7:10])),[.4,.9],order=12,variance_floor=self.x[10]);np.testing.assert_allclose(result[:,:,0],np.array([prod.terminal_ratios[t] for t in [.4,.9]]),rtol=3e-14)
  for j in range(11):
   h=1e-5*max(abs(self.x[j]),.001);a=self.x.copy();b=a.copy();a[j]+=h;b[j]-=h
   fd=(run(a)[:,:,0]-run(b)[:,:,0])/(2*h);np.testing.assert_allclose(fd,result[:,:,1+j],rtol=2e-5,atol=2e-7,err_msg=NAMES[j])
if __name__=='__main__':unittest.main()
