import ctypes as ct
import numpy as np
from .full_readout_fit import arr

class StockValues:
    def __init__(self,lib):
        self.lib=ct.CDLL(str(lib));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.lib.stock_values.argtypes=[ct.c_int,ct.c_int,P,P,P,I,P,P,P]
    def __call__(self,folder,norm,p):
        x=np.load(folder/'features.npy',mmap_mode='r');dw=np.load(folder/'equity.npy',mmap_mode='r');grid=np.load(folder/'grid.npy')
        seg=np.minimum(np.searchsorted([.25,.75],grid[:-1],side='right'),2).astype('int32');out=np.empty(x.shape[1])
        self.lib.stock_values(len(dw),x.shape[1],x,dw,arr(np.diff(grid)),seg,arr(norm),arr(p),out)
        return out
