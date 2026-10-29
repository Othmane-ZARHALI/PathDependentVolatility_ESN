"""Locally scaled protected corrections: same OU law and positive floor f/2."""
import ctypes as ct
import hashlib,json,shutil,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import quad
from scipy.special import erf
from .full_readout_fit import HERE,arr,dump
from .protected_curve_check import context
from .asymmetric_readout import PARENT,LEVERAGE,Kernel as OldKernel,evaluate_stocks
from .forward_leverage import relabel
from .leverage_scan import combine
from .m2_protected_variance import ENDS,MODES,envelope,Kernel as PreviousKernel

OUT=HERE/'runs/m2_rich_protected_variance_001'
CONFIG=[dict(mode=2,width=1.,kind=-1,name='previous protected 27-day tanh')]+[
    dict(mode=j,width=float(w),kind=kind,name=f'local {half}d '+('tanh' if kind==0 else 'erf')+f' width{w}')
    for j,half in enumerate([3.3,9.5,27.]) for w in [1,2] for kind in [0,1]]

def time_basis(t,k,width,kind):
    t=np.asarray(t)
    if kind<0:
        a,ap=envelope(t);K=365/10+k*(.5+4/(3*np.sqrt(3)));return a/K,ap/K
    # Piecewise C1: q and q_t vanish at every quoted window endpoint.
    q=np.zeros_like(t);qd=np.zeros_like(t)
    C=.5+4/(3*np.sqrt(3)*width**2) if kind==0 else np.sqrt(2/(np.pi*np.e))*(1+width**-2)
    def put(mask,a,ap,A):
        B=k*C;K=B/2+np.sqrt(A*A+B*B/4);q[mask]=a/K;qd[mask]=ap/K
    mask=t<ENDS[0];T=ENDS[0];v=np.pi*t[mask]/(2*T);A=np.pi/(2*T);put(mask,np.cos(v)**2,-A*np.sin(2*v),A)
    for left,right in zip(ENDS[:-1],ENDS[1:]):
        mask=(t>=left)&(t<right);A=np.pi/(right-left);v=A*(t[mask]-left);put(mask,np.sin(v)**2,A*np.sin(2*v),A)
    mask=t>=ENDS[-1];T=10/365;u=t[mask]-ENDS[-1];A=np.pi/(2*T);v=A*np.minimum(u,T);put(mask,np.where(u>=T,1.,np.sin(v)**2),np.where(u>=T,0.,A*np.sin(2*v)),A)
    return q,qd

def hvalue(t,u,k,width,kind):
    q,qd=time_basis(t,k,width,kind);u=np.asarray(u)
    if kind<=0:
        z=np.tanh(u/width);return qd*z-q*k*(u/width+2*z/width**2)*(1-z*z)
    derivative=np.sqrt(2/np.pi)/width*np.exp(-u*u/(2*width**2))
    return qd*erf(u/(np.sqrt(2)*width))-q*k*u*(1+width**-2)*derivative

class Kernel:
    def __init__(self,path):
        self.old=OldKernel(path);self.k=self.old.k;self.k.rich_stock.argtypes=self.k.asymmetric_stock.argtypes
        P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS');self.k.rich_probe.argtypes=[ct.c_int]+[P]*4+[I,P]
    def probe(self,t,u,k,width,kind):
        t,u,k,width,kind=np.broadcast_arrays(t,u,k,width,kind);out=np.zeros((t.size,3));self.k.rich_probe(t.size,arr(t.ravel()),arr(u.ravel()),arr(k.ravel()),arr(width.ravel()),np.ascontiguousarray(kind.ravel(),dtype=np.int32),out);return out.reshape(t.shape+(3,))
    def stocks(self,d,pairs,p,x,coeff,rates,seed):
        coeff=arr(coeff);assert coeff.shape[1]==13 and np.max(abs(coeff).sum(1))<=1+1e-14
        models=arr(np.tile(x,(len(coeff),1)));out=np.zeros((2*pairs,int(d['expiries']),len(coeff)+1,2));self.k.rich_stock(len(d['dt']),pairs,d['projection'].shape[1],int(d['expiries']),len(coeff),d['transition'],d['decay'],d['root'],d['cross'],d['mean'],d['projection'],d['time'],d['dt'],d['seg'],d['coarse'],d['reference'],d['expiry'],d['lognorm'],arr(p),models,coeff,arr(rates),1.,int(seed),out);assert np.isfinite(out).all() and (out>0).all();return out

def audit(k,d,p,x,rates):
    rng=np.random.default_rng(20261029101);records=[];maxerr=0.;maxh=0.;boundary=0.
    for cfg in CONFIG:
        rate=rates[cfg['mode']];width=cfg['width'];kind=cfg['kind'];t=np.r_[rng.uniform(0,1030/365,200),ENDS,ENDS-1e-10,ENDS+1e-10];u=rng.normal(0,3,len(t));got=k.probe(t,u,rate,width,kind);q,qd=time_basis(t,rate,width,kind);expected=np.column_stack([q,qd,hvalue(t,u,rate,width,kind)]);maxerr=max(maxerr,float(abs(got-expected).max()));boundary=max(boundary,float(abs(k.probe(ENDS,np.ones(8),rate,width,kind)[:,:2]).max()))
        tt=np.linspace(0,1030/365,1500);uu=np.linspace(-10,10,301);maxh=max(maxh,float(abs(k.probe(tt[:,None],uu,rate,width,kind)[...,2]).max()))
        for order in [64,96]:
            z,w=hermgauss(order);z=np.sqrt(2)*z;w=w/np.sqrt(np.pi)
            for day in [27,55,118,153]:
                T=day/365
                for state in [-4.,-1.,1.,4.]:
                    def integrand(s):
                        lag=max(0,s-T);u=np.exp(-rate*lag)*state+np.sqrt(-np.expm1(-2*rate*lag))*z
                        return float(k.probe(s,u,rate,width,kind)[:,2]@w)
                    value,error=quad(integrand,T,T+30/365,epsabs=2e-11,points=ENDS[(ENDS>T+1e-13)&(ENDS<T+30/365-1e-13)].tolist(),limit=100)
                    records.append(dict(**cfg,GH=order,day=day,state=state,integral=value,estimated_quad_error=error))
    residual=max(abs(r['integral']) for r in records);assert maxerr<1e-11 and boundary<1e-10 and maxh<=1 and residual<1e-8
    small=dict(d);n=97
    for key in ['time','dt','mean','seg','coarse','reference','expiry','transition','lognorm']:small[key]=d[key][:n].copy()
    small['expiry'][:]=-1;small['expiry'][-1]=0;small['expiries']=1
    coeff=arr([np.zeros(13),np.eye(13)[0]]);got=k.stocks(small,4,p,x,coeff,rates,778982)
    old=PreviousKernel(OUT/'rich.dylib').stocks(small,4,p,x,arr([[0,0,0],[0,0,1]]),rates,778982)
    replay=float(abs(got-old).max());assert replay<1e-12
    dump(OUT/'implementation_audit.json',dict(compiled_formula_max_error=maxerr,endpoint_q_and_dtq_max_error=boundary,observed_max_abs_h=maxh,analytic_global_h_bound=1.,max_conditional_integral_error=residual,baseline_and_previous_correction_replay_error=replay,conditional_records=records,positivity_floor=float((.05+.07*x[19])**2/2),no_clipping=True,Lean_status='no new Lean proof'))
    print('RICH AUDIT PASS',maxerr,residual,maxh,replay,flush=True)

def stage(k,d,p,x,rates,coeff,name,pairs,seedtag):
    folder=OUT/name;folder.mkdir();pieces=[];seeds=[]
    for batch in range(8):
        seed=int(np.random.SeedSequence([20261029109,seedtag,batch]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed);start=time.perf_counter();pieces.append(k.stocks(d,pairs,p,x,coeff,rates,seed));print(name,'batch',batch,'seconds',round(time.perf_counter()-start,2),flush=True)
    stocks=combine(pieces);np.savez_compressed(folder/'stocks.npz',stocks=stocks,seeds=np.array(seeds,dtype=np.uint64),coeff=coeff)
    evaluate_stocks(stocks,np.load(LEVERAGE/'reference.npz')['stats'],arr(np.tile(x,(len(coeff),1))),np.arange(len(coeff),dtype=float),pd.read_csv(PARENT/'spx_membership.csv'),folder/'SPX');return relabel(folder/'SPX')

def run():
    OUT.mkdir(exist_ok=False);(OUT/'source').mkdir();sources=[Path(__file__),HERE/'m2_rich_protected_variance.cpp',HERE/'m2_protected_variance.cpp',HERE/'m2_protected_variance.py',HERE/'asymmetric_readout.cpp',HERE/'leverage_scan.cpp',HERE/'two_speed_precision_kernel.cpp',PARENT/'selection.json']
    for p in sources:shutil.copy2(p,OUT/'source'/p.name)
    dump(OUT/'source_hashes.json',{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    lib=OUT/'rich.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'m2_rich_protected_variance.cpp'),'-o',str(lib)],check=True)
    prior,dated,_,normalizer=context();x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);p=arr(prior['selected']);rates=arr(dated.layout.rates[MODES]);assert np.allclose((dated.layout.driver_loadings[MODES]**2).sum(1),1,atol=1e-14);d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));k=Kernel(lib)
    coeff=arr([np.zeros(13),np.eye(13)[0]]+[sign*np.eye(13)[j] for j in range(1,13) for sign in [-1,1]])
    dump(OUT/'design.json',dict(config=CONFIG,pilot_coefficients=coeff.tolist(),mode_indices=MODES.tolist(),half_life_days=(np.log(2)*365/rates).tolist(),baseline=x.tolist(),normalizer=normalizer,OU_dynamics_unchanged=True,stock_coupling_unchanged=True,positivity_floor=float((.05+.07*x[19])**2/2),amplitude='f/2 times L1-unit coefficient vector',quoted_VIX_days=[27,55,118,153],curve_refit=False,adopted=False))
    audit(k,d,p,x,rates)
    scores=stage(k,d,p,x,rates,coeff,'pilot',32,111)
    best=min([r for r in scores if r['steps']==256 and r['candidate']>=2],key=lambda r:r['loss'])['candidate'];selected=arr([coeff[0],coeff[1],coeff[best]])
    dump(OUT/'selection.json',dict(pilot_new_best=best,coefficients=selected.tolist(),criterion='minimum pilot fine-grid primary SPX interval loss among new shapes',frozen_before_fresh=True,new_shape=CONFIG[int(np.flatnonzero(coeff[best])[0])]))
    fresh=stage(k,d,p,x,rates,selected,'fresh',64,222)
    dump(OUT/'run.json',dict(status='pricing complete; independent recurrence and reporting pending',pilot_paths=512,fresh_paths=1024,pilot_scores=scores,fresh_scores=fresh,adopted=False,production_baseline_unchanged=True))
    print('RICH PRICING COMPLETE',flush=True)

if __name__=='__main__':run()
