"""Independent-pilot hedging controls, shared across the put/call boundary."""
import ctypes as ct
from pathlib import Path
import numpy as np
from scipy.special import ndtr
from .full_readout_fit import arr
from .parity_adaptive import reference_calls,fit_control_coefficients
from .continuous_strike_controls import interpolation_weights


def anchors_for(q):
    f=float(q.forward.iloc[0])
    k=q.K_over_F.to_numpy()
    return np.unique(np.r_[np.geomspace(k.min(),1.,5),np.geomspace(1.,k.max(),5)])


class FrozenHedgeKernel:
    def __init__(self,lib):
        self.lib=ct.CDLL(str(lib));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.lib.frozen_hedges.argtypes=[ct.c_int,ct.c_int,ct.c_int,P,P,P,I,I,P,P,P,P,P]
    def evaluate(self,folder,norm,p,anchors):
        grid=np.load(folder/'grid.npy');x=np.load(folder/'features.npy',mmap_mode='r');dw=np.load(folder/'equity.npy',mmap_mode='r');dt=np.diff(grid)
        seg=np.minimum(np.searchsorted([.25,.75],grid[:-1],side='right'),2).astype('int32')
        integrated=np.r_[0,np.cumsum(p[7+seg]*dt)]
        remaining=integrated[-1]-integrated[:-1]
        # About four hedge trades per trading day; time grid itself unchanged.
        rebalance=np.zeros(len(dt),dtype='int32');next_time=0.
        for j,t in enumerate(grid[:-1]):
            if t+1e-14>=next_time:rebalance[j]=1;next_time=t+1/(252*4)
        out=np.empty((x.shape[1],1+2*len(anchors)))
        self.lib.frozen_hedges(len(dt),x.shape[1],len(anchors),x,dw,arr(dt),seg,rebalance,arr(remaining),arr(norm(grid[:-1])),arr(p),arr(anchors),out)
        return out


def control_data(folder,output,anchors,p):
    w=np.load(folder/'weights.npy');a=np.load(folder/'reference.npz');b=a['brownian'];duration=a['duration'];total=float(duration@p[7:10])
    ref=np.exp(b@np.sqrt(p[7:10])-.5*total);d1=-np.log(anchors)/np.sqrt(total)+.5*np.sqrt(total)
    known=ndtr(d1)-anchors*ndtr(d1-np.sqrt(total))
    refpay=np.maximum(ref[:,None]-anchors,0)
    # One zero-mean family per anchor; all are calls. No put/call switching.
    controls=np.stack([np.broadcast_to((w-1)[:,None],refpay.shape),
                       np.broadcast_to((w*(output[:,0]-1))[:,None],refpay.shape),
                       w[:,None]*(refpay-known),w[:,None]*output[:,1::2],w[:,None]*output[:,2::2]],axis=-1)
    y,_=reference_calls(output[:,0],np.zeros((len(w),1)),ref,np.zeros((len(w),1)),w,anchors,total,np.zeros(1))
    return controls,y


def fit_anchors(controls,base_draws):
    return np.stack([fit_control_coefficients(controls[:,j],base_draws[:,j,None])[:,0] for j in range(base_draws.shape[1])])


def correction_at_strikes(controls,beta,anchors,strikes):
    anchor_draws=np.einsum('nac,ac->na',controls,beta)
    return anchor_draws@interpolation_weights(np.asarray(strikes),np.asarray(anchors)).T
