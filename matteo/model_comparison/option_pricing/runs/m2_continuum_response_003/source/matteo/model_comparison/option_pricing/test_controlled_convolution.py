"""Same-GH32 value/derivative correspondence across actual conditional laws."""
from pathlib import Path
import ctypes
import json
import subprocess
import tempfile
import unittest
import numpy as np
from .full_readout_fit import Kernel,HERE,arr,NAMES
from .m2_spx import frozen_stage3
from esn_eq.pricing import DatedM2Pricer,_normal_rule
from esn_eq.states import symmetric_square_root

class ConvolutionTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.tmp=tempfile.TemporaryDirectory();cls.folder=Path(cls.tmp.name)
  for name,source in [('direct','controlled_readout_kernel.cpp'),('fast','controlled_readout_convolution_kernel.cpp')]:subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/source),'-o',str(cls.folder/(name+'.dylib'))],check=True)
  cls.direct=Kernel(cls.folder/'direct.dylib');cls.fast=Kernel(cls.folder/'fast.dylib');a,p,r,_=frozen_stage3();cls.dated=DatedM2Pricer(a,r);cls.base=np.array([getattr(p,n) for n in NAMES[:7]]+[.06,.07,.06,.0001]);cls.base[4]=.0001
 @classmethod
 def tearDownClass(cls):cls.tmp.cleanup()
 def law(self,days,lag):
  states=np.load(HERE/f'runs/m2_controlled_readout_cache_001/search/vix_d{days}.npz')['states'][::16];m,c=self.dated.gaussian_projection_law(lag,states);z,w=_normal_rule(32);return arr(m),arr(z@symmetric_square_root(c).T),arr(w)
 def test_values_and_derivatives_across_starts_and_lags(self):
  starts=[self.base.copy(),np.array(json.loads((HERE/'runs/m2_full_readout_fit_001/run.json').read_text())['selected_parameters']),self.base.copy()];starts[2][[1,2,3,4,5,6]]=[.12,2.,.35,.5,6.5,6.]
  worstvalue=worstder=0.
  for days in [27,55,118,153]:
   for lag in [.00002,.01,.08]:
    law=self.law(days,lag)
    for p in starts:
     direct=self.direct.expectation(*law,p);fast=self.fast.expectation(*law,p);ve=float(np.max(abs(fast[:,0]/direct[:,0]-1)));de=float(np.max(abs(fast[:,1:]-direct[:,1:])/(1+abs(direct[:,1:]))));worstvalue=max(worstvalue,ve);worstder=max(worstder,de);self.assertLess(ve,2e-8);self.assertLess(de,2e-5)
  print('Convolution correspondence max value-relative / derivative-scaled errors:',worstvalue,worstder,flush=True)
 def test_independence_failure_falls_back(self):
  m,o,w=self.law(55,.01);o=o.copy();o[:,0]+=.01*o[:,1];p=self.base;np.testing.assert_allclose(self.fast.expectation(m,o,w,p),self.direct.expectation(m,o,w,p),rtol=1e-14,atol=1e-12)
 def test_used_table_is_exercised(self):
  self.fast.expectation(*self.law(153,.01),self.base);counts=(ctypes.c_longlong*3)();self.fast.lib.convolution_counts(counts);self.assertGreater(counts[0],0);print('Convolution table/fallback/max-knots counters:',list(counts),flush=True)

if __name__=='__main__':unittest.main()
