"""Stock/stress rotation used by the exploratory fixed-factor scan."""
import ctypes as ct
import numpy as np

def arr(x):
    return np.ascontiguousarray(x, dtype=np.float64)

class Kernel:
    def __init__(self,lib):
        self.k=ct.CDLL(str(lib));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.k.leverage_stats.argtypes=[ct.c_int]*4+[I]+[P]*7+[I]*4+[P]*3+[ct.c_uint64,P]
        self.k.leverage_reference.argtypes=[ct.c_int]*3+[P]*2+[I]*2+[P,ct.c_uint64,P]
    def stats(self,d,pairs,p,x,seed):
        out=np.zeros((2*pairs,int(d['expiries']),3,4));self.k.leverage_stats(len(d['dt']),pairs,d['projection'].shape[1],int(d['expiries']),d['transition'],d['decay'],d['root'],d['cross'],d['mean'],d['projection'],d['time'],d['dt'],d['seg'],d['coarse'],d['reference'],d['expiry'],d['lognorm'],arr(p),arr(x),int(seed),out);return out
    def reference(self,d,pairs,x,seed):
        idx=np.flatnonzero(d['reference']);grid=np.r_[d['time'][idx],d['time'][-1]+d['dt'][-1]];expiry=np.full(len(idx),-1,np.int32)
        for t,e in enumerate(d['expiry']):
            if e>=0:
                end=d['time'][t]+d['dt'][t];j=np.searchsorted(grid,end);assert abs(grid[j]-end)<1e-12;expiry[j-1]=e
        out=np.zeros((2*pairs,int(d['expiries']),4));self.k.leverage_reference(len(idx),pairs,int(d['expiries']),arr(grid[:-1]),arr(np.diff(grid)),d['seg'][idx].copy(),expiry,arr(x),int(seed),out);return out

def stock(s,theta,cs,cf):
    return np.exp(-.5*s[...,0]+np.sqrt(1-theta*theta)*s[...,1]-theta*(cs*s[...,2]+cf*s[...,3]))
