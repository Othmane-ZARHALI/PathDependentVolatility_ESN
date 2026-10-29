"""Production/serial correspondence of parallel kernel and corrected oracle."""
from pathlib import Path
import subprocess
import unittest
import numpy as np
from . import test_full_readout_fit as legacy
from .full_readout_fit import Kernel,HERE,arr
from .controlled_readout_fit import ControlledEngine

class ControlledKernelCorrespondence(legacy.KernelCorrespondence):
 @classmethod
 def setUpClass(cls):
  super().setUpClass();lib=cls.folder/'parallel.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'controlled_readout_kernel.cpp'),'-o',str(lib)],check=True);cls.k=Kernel(lib)
 def test_parallel_expectation_matches_serial_for_many_states(self):
  serial=Kernel(self.folder/'kernel.dylib');rng=np.random.default_rng(6);m=arr(rng.normal(size=(128,3)));o=arr(rng.normal(size=(125,3)));w=arr(np.full(125,1/125));a=self.k.expectation(m,o,w,self.x);b=serial.expectation(m,o,w,self.x);np.testing.assert_allclose(a,b,rtol=2e-12,atol=1e-11)
 def test_adaptive_normalizer_matches_direct_production(self):
  from concurrent.futures import ThreadPoolExecutor
  e=object.__new__(ControlledEngine);e.k=self.k;e.dated=self.dated;e.initial,_=self.dated.initial_state([1,0,0]);e.tau=1/max(self.dated.layout.rates);e.horizon=.9;e.pool=ThreadPoolExecutor(max_workers=4);e.normal_laws={};e.last_readout=None;e.normal_cached=None;e.exposure_cached={}
  try:
   f,audit=e.normalizer(self.x);t=np.array([0,1e-5,.0007,.017,.13,.73]);got=np.exp(f(t)[:,0]);expected=self.dated.normalizer(t,e.initial,self.pars,48);np.testing.assert_allclose(got,expected,rtol=2e-7);self.assertLess(audit['independent_error'],2e-7)
  finally:e.pool.shutdown()

if __name__=='__main__':unittest.main()
