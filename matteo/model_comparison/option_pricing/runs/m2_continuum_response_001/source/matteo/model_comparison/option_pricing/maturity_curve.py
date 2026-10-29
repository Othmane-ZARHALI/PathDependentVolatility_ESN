"""Maturity-aligned curve policy and generic reusable-OU SPX/VIX engine."""
import ctypes as ct
import numpy as np
from .calibration_policy import LOWER as OLD_LOWER, UPPER as OLD_UPPER
from .forward_variance_constraints import validate_forward_variance_levels
from .full_readout_fit import Kernel, arr, SCALE
from .parity_adaptive_engine import ParityAdaptiveEngine
from .parity_parameter_core import FitEngine
from .maturity_curve_predict import GenericPredict
from .all_maturity_core import InSampleLoss

ENDS_DAYS=np.array([27,29,55,57,85,91.25,118,119,148,153,183,211,273.75,302,666,1030],float)


def vix_only_loss(quotes,forwards):
    from esn_eq.calibration_loss import JointPriceObjective,LossConfig
    JointPriceObjective(quotes,forwards,config=LossConfig(spx_weight=0.),split='in_sample_fit')
    loss=object.__new__(InSampleLoss);loss.quotes=quotes.copy().reset_index(drop=True);loss.forwards=forwards.copy().reset_index(drop=True)
    loss.weights=np.zeros(len(quotes))
    for _,g in loss.quotes.groupby('expiry'):loss.weights[g.index]=1/(loss.quotes.expiry.nunique()*len(g))
    loss.sqrtweights=np.sqrt(loss.weights)
    return loss


def bounds(n):
    lower=np.r_[OLD_LOWER[:7],np.full(n,OLD_LOWER[-1]),OLD_LOWER[-1]]
    upper=np.r_[OLD_UPPER[:7],np.full(n,.2),OLD_UPPER[-1]]
    scale=np.r_[SCALE[:7],np.full(n,.06),SCALE[-1]]
    return lower,upper,scale


def validate(p):
    p=np.asarray(p,float);n=len(p)-8
    if n<1 or p.ndim!=1 or not np.isfinite(p).all():raise ValueError('Invalid curve parameters')
    lo,hi,_=bounds(n)
    if np.any(p<lo) or np.any(p>hi):raise ValueError('Outside maturity-curve bounds')
    validate_forward_variance_levels(p[7:-1],p[-1])
    return p


def embed_three_piece(p,ends_days=ENDS_DAYS):
    """Exactly represent the old curve; retain its 91.25/273.75-day knots."""
    mids=(np.r_[0.,ends_days[:-1]]+ends_days)/2
    segments=np.searchsorted([91.25,273.75],mids,side='right')
    return np.r_[p[:7],np.asarray(p[7:10])[segments],p[-1]]


def window_nodes(expiry,ends,tau,order=8):
    delta=30/365;lags=[0.,delta];lag=tau
    while lag<delta:lags.append(lag);lag*=2
    events=np.unique([expiry+x for x in lags]+[x for x in ends if expiry<x<expiry+delta])
    z,w=np.polynomial.legendre.leggauss(order);half=np.diff(events)/2
    return ((events[:-1]+half)[:,None]+half[:,None]*z).ravel(),(half[:,None]*w).ravel()


class CurveKernel(Kernel):
    def __init__(self,convolution,stock,ends):
        super().__init__(convolution);self.curve_lib=ct.CDLL(str(stock));self.ends=np.array(ends);self.grids={}
        P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.curve_lib.stock_curve.argtypes=[ct.c_int,ct.c_int,ct.c_int,ct.c_int,P,P,P,I,P,P,P]

    def stock(self,folder,norm,p,derivatives=True):
        x=np.load(folder/'features.npy',mmap_mode='r');dw=np.load(folder/'equity.npy',mmap_mode='r')
        key=folder.parent
        if key not in self.grids:
            grid=np.load(folder/'grid.npy');segment=np.minimum(np.searchsorted(self.ends,grid[:-1],side='right'),len(self.ends)-1).astype('int32')
            self.grids[key]=(arr(np.diff(grid)),segment)
        dt,segment=self.grids[key]
        out=np.zeros((x.shape[1],len(p)+1))
        self.curve_lib.stock_curve(len(dw),x.shape[1],len(self.ends),int(derivatives),x,dw,dt,segment,arr(norm),arr(p),out)
        return out[None,:,:]


class CurveEngine(GenericPredict,ParityAdaptiveEngine):
    vix_prediction=FitEngine.vix_prediction

    def __init__(self,*args,ends,reference_levels,**kwargs):
        super().__init__(*args,**kwargs)
        self.ends=np.asarray(ends);self.ncurve=len(ends);self.npar=8+self.ncurve
        self.reference_levels=np.asarray(reference_levels);self.readout_exposures={};self.exposure_readout_key=None

    def vix_values(self,days,p,norm,role='production'):
        readout=tuple(p[:7])
        if readout!=self.exposure_readout_key:self.readout_exposures={};self.exposure_readout_key=readout
        key=(days,role)
        if key not in self.vixlaws:
            with np.load(self.cache/role/f'vix_d{days}.npz') as saved:a={k:saved[k] for k in saved.files}
            times,weights=window_nodes(days/365,self.ends,self.tau)
            a['times']=times;a['time_weights']=weights
            self.vixlaws[key]=(a,[self.law(float(t-days/365),a['states'],32) for t in times])
        a,laws=self.vixlaws[key]
        if key not in self.readout_exposures:
            expo=np.zeros((len(a['weights']),self.ncurve));de=np.zeros((len(expo),self.ncurve,7))
            for t,wt,law in zip(a['times'],a['time_weights'],laws):
                e=self.k.expectation(*law,p);nn=norm(t);m=np.exp(nn[0]);ratio=e[:,0]/m
                j=min(np.searchsorted(self.ends,t,side='right'),self.ncurve-1)
                expo[:,j]+=wt*ratio;de[:,j]+=wt*(e[:,1:]/m-ratio[:,None]*nn[1:])
            self.readout_exposures[key]=(expo,de)
        expo,de=self.readout_exposures[key];delta=30/365
        integral=p[-1]*delta+expo@(p[7:-1]-p[-1]);di=np.zeros((len(expo),self.npar))
        di[:,:7]=np.einsum('nsj,s->nj',de,p[7:-1]-p[-1]);di[:,7:-1]=expo;di[:,-1]=delta-expo.sum(1)
        return 100*np.sqrt(integral/delta),50/np.sqrt(delta*integral)[:,None]*di,a
