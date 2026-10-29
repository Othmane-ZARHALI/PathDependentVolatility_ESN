"""Independent OU spike-sensitivity layer; experimental, value-only pricing.

The original 3D Gaussian rule is kept. A checked Chebyshev expansion in the
bounded spike multiplier lets the added Gaussian integral share old work across
the six strength/persistence choices. It changes no model equation.
"""
import ctypes as ct
import hashlib
import json
import time
from pathlib import Path
from functools import lru_cache

import numpy as np
from scipy.fft import dct
from scipy.interpolate import CubicSpline
from scipy.special import ndtri
from scipy.stats import qmc

from .full_readout_fit import arr
from .maturity_curve import CurveEngine, CurveKernel, window_nodes

HALF_LIVES = (7, 30, 90)
FACTORS = (1.5, 2.)
SPECS = [(f'z{factor:g}_h{half}', float(np.log(factor)),
          float(np.log(2)*365/half), half, factor)
         for factor in FACTORS for half in HALF_LIVES]
SEED = 2026102607


@lru_cache(maxsize=16)
def normal_rule(order):
    z, w = np.polynomial.hermite.hermgauss(order)
    return arr(z*np.sqrt(2)), arr(w/np.sqrt(np.pi))


def coefficients(values):
    """Values on descending first-kind Chebyshev nodes, output row per state."""
    c = dct(values, type=2, axis=0)/len(values)
    c[0] *= .5
    return arr(c.T)


def seed_for(*parts):
    return int(np.random.SeedSequence([SEED, *parts]).generate_state(1, dtype=np.uint64)[0])


def outer_z(cache, day):
    """Independent scrambled normals, aligned with eight original 512-point nets."""
    normals=[]; records=[]
    for net in range(8):
        seed=seed_for(1, day, net)
        u=qmc.Sobol(1,scramble=True,seed=seed).random_base2(9)[:,0]
        if np.any((u<=0)|(u>=1)):raise ValueError('Nonfinite Gaussian draw')
        normals.append(ndtri(u));records.append(dict(net=net,seed=seed))
    return arr(np.concatenate(normals)),records


class LayerKernel(CurveKernel):
    def __init__(self, library, ends):
        super().__init__(library,library,ends)
        P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.lib.configure.argtypes=[ct.c_int,ct.c_double]
        self.lib.direct.argtypes=self.lib.expectation.argtypes
        self.lib.integrate_response.argtypes=[ct.c_int,ct.c_int,ct.c_int,P,P,ct.c_double,P,P,ct.c_double,P]
        self.lib.innovations.argtypes=[ct.c_uint64,ct.c_int,P]
        self.lib.zstock_curve.argtypes=[ct.c_int,ct.c_int,ct.c_int,P,P,P,I,P,P,ct.c_double,ct.c_double,ct.c_uint64,P,P]
        self.spec=SPECS[0];self.stock_audit=[]

    def configure(self, continuation):
        self.lib.configure(0 if continuation==0 else 2,continuation)

    def response(self,c,mean,sd,sensitivity,order):
        z,w=normal_rule(order);mean=arr(np.broadcast_to(mean,len(c)))
        out=np.empty(len(c));self.lib.integrate_response(len(c),c.shape[1],len(z),c,mean,float(sd),z,w,float(sensitivity),out)
        if np.any(out<=0) or not np.isfinite(out).all():raise ValueError('Nonpositive/nonfinite conditional multiplier')
        return out

    def stock(self,folder,norm,p,derivatives=False):
        if derivatives:raise ValueError('Value-only layer experiment')
        name,s,rate,half,factor=self.spec
        day=int(folder.parent.name[1:]);batch=int(folder.name[5:]);seed=seed_for(2,day,batch)
        x=np.load(folder/'features.npy',mmap_mode='r');dw=np.load(folder/'equity.npy',mmap_mode='r')
        key=folder.parent
        if key not in self.grids:
            grid=np.load(folder/'grid.npy');segment=np.minimum(np.searchsorted(self.ends,grid[:-1],side='right'),len(self.ends)-1).astype('int32')
            self.grids[key]=(arr(np.diff(grid)),segment)
        dt,segment=self.grids[key];out=np.zeros((x.shape[1],len(p)+1));tz=np.empty(x.shape[1])
        if batch==0:self.begin=time.perf_counter();print('SPX START',name,day,flush=True)
        self.lib.zstock_curve(len(dw),x.shape[1],len(self.ends),x,dw,dt,segment,arr(norm),arr(p),rate,s,seed,out,tz)
        self.stock_audit.append(dict(candidate=name,day=day,batch=batch,seed=seed,terminal_z=tz.tolist()))
        if batch==31:print('SPX DONE',name,day,'seconds',round(time.perf_counter()-self.begin,2),flush=True)
        return out[None,:,:]


def response_table(k,law,p,checks,label,tolerance=2e-6):
    """Adapt degree using off-node checks including extreme projected states.

    Off-node checks use direct original finite GH quadrature (<=12 rows), rather
    than trusting the accelerated scalar convolution on both sides.
    """
    mean,offset,w=law;mean=arr(mean).reshape(-1,3);n=len(mean)
    ids=np.unique(np.r_[np.arange(min(3,n)),np.argmax(abs(mean),axis=0),np.argmin(mean,axis=0),np.argmax(mean,axis=0)])
    probes=np.array([.503,.71,1.,1.237,1.731,1.997])
    actual=[]
    for A in probes:
        pp=p.copy();pp[3]*=A
        actual.append(k.expectation(arr(mean[ids]),offset,w,pp)[:,0])
    actual=np.array(actual)
    for degree in [32,64,128]:
        nodes=1.25+.75*np.cos(np.pi*(np.arange(degree)+.5)/degree)
        values=[]
        for A in nodes:
            pp=p.copy();pp[3]*=A;values.append(k.expectation(mean,offset,w,pp)[:,0])
        c=coefficients(np.array(values))
        pred=np.polynomial.chebyshev.chebval((2*probes-2.5)/1.5,c[ids].T).T
        err=float(np.max(abs(pred/actual-1)))
        if err<tolerance:break
    else:raise RuntimeError(f'Response interpolation failed {label}: {err}')
    checks.append(dict(label=label,degree=degree,relative_error=err,states=ids.tolist(),factors=probes.tolist()))
    return c


class LayerEngine(CurveEngine):
    """Install precomputed candidate-specific normalizer and VIX exposures."""
    def normalizer(self,p):
        return self.layer_normalizer,self.layer_normalizer_audit

    def vix_values(self,days,p,norm,role='production'):
        a=self.layer_draws[days];expo=self.layer_exposures[days];delta=30/365
        integral=p[-1]*delta+expo@(p[7:-1]-p[-1]);di=np.zeros((len(expo),self.npar))
        di[:,7:-1]=expo;di[:,-1]=delta-expo.sum(1)
        return 100*np.sqrt(integral/delta),50/np.sqrt(delta*integral)[:,None]*di,a


def make_normalizers(e,k,p,checks,out):
    """All six candidates share one original-state integration per response node."""
    direct_values={};rules={};tau=e.tau
    def direct(s):
        s=float(s)
        if s not in direct_values:
            t=float(tau*np.expm1(s));law=e.law(t,e.initial,48)
            c=response_table(k,law,p,checks,f'normalizer_t{t:.12g}',3e-7)
            values=[]
            for _,strength,rate,_,_ in SPECS:
                sd=np.sqrt(-np.expm1(-2*rate*t));values.append(k.response(c,0.,sd,strength,192)[0])
            direct_values[s]=np.log(values)
        return direct_values[s]
    knots=np.linspace(0,np.log1p(e.horizon/tau),257)
    for iteration in range(10):
        vals=np.array([direct(x) for x in knots]);spline=CubicSpline(knots,vals,axis=0,extrapolate=False)
        mid=(knots[:-1]+knots[1:])/2;actual=np.array([direct(s) for s in mid]);err=np.max(abs(np.expm1(spline(mid)-actual)),axis=1)
        if np.max(err)<1e-7:break
        knots=np.unique(np.r_[knots,mid[err>=1e-7]])
        if len(knots)>2048:raise RuntimeError('Normalizer knot budget exhausted')
        print('NORMALIZER refine',iteration,'knots',len(knots),flush=True)
    else:raise RuntimeError('Normalizer failed')
    probes=np.linspace(knots[0],knots[-1],131)[1:-1];actual=np.array([direct(x) for x in probes]);error=float(np.max(abs(np.expm1(spline(probes)-actual))))
    if error>2e-7:raise RuntimeError(f'Independent normalizer probes failed {error}')
    audit=dict(knots=len(knots),midpoint_error=float(np.max(err)),independent_error=error,original_GH=48,new_Z_GH=192,interpolation='Checked Chebyshev response expansion and cubic log-time normalizer')
    np.savez_compressed(out,knots=knots,log_values=vals)
    def single(j):
        def norm(t):
            val=spline(np.log1p(np.asarray(t)/tau))[...,j];a=np.zeros(np.shape(val)+(8,));a[...,0]=val;return a
        return norm
    return {spec[0]:(single(j),audit) for j,spec in enumerate(SPECS)}


def make_exposures(e,k,p,normalizers,cache,checks,out,time_order=8,z_order=128,days=(27,55,118,153)):
    exposures={spec[0]:{} for spec in SPECS};draws={};seeds={}
    for day in days:
        tick=time.perf_counter()
        with np.load(cache/'production'/f'vix_d{day}.npz') as f:a={key:f[key] for key in f.files}
        if 'z_normals' in a:eps=a['z_normals'];records=json.loads((cache/'production'/f'vix_d{day}_seeds.json').read_text())
        else:eps,records=outer_z(cache,day)
        a['z_normals']=eps;seeds[str(day)]=records
        times,tw=window_nodes(day/365,e.ends,e.tau,time_order);draws[day]=a
        for name,*_ in SPECS:exposures[name][day]=np.zeros((4096,len(e.ends)))
        for index,(t,weight) in enumerate(zip(times,tw)):
            lag=float(t-day/365);law=e.law(lag,a['states'],32)
            c=response_table(k,law,p,checks,f'vix_d{day}_node{index}')
            j=min(np.searchsorted(e.ends,t,side='right'),len(e.ends)-1)
            for name,strength,rate,_,_ in SPECS:
                zT=np.sqrt(-np.expm1(-2*rate*day/365))*eps
                mu=np.exp(-rate*lag)*zT;sd=np.sqrt(-np.expm1(-2*rate*lag))
                value=k.response(c,mu,sd,strength,z_order);m=np.exp(normalizers[name][0](t)[0])
                exposures[name][day][:,j]+=weight*value/m
            if index%20==0:print('VIX TABLE',day,index,'/',len(times),'seconds',round(time.perf_counter()-tick,1),flush=True)
        print('VIX TABLE DONE',day,round(time.perf_counter()-tick,1),flush=True)
    np.savez_compressed(out,**{f'{name}_d{day}':value for name,ee in exposures.items() for day,value in ee.items()})
    return exposures,draws,seeds
