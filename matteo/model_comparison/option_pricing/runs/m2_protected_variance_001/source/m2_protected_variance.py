"""VIX-protected additive variance experiment, fixed OU bank and stress model."""
import ctypes as ct
import hashlib,json,shutil,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import quad
from .full_readout_fit import HERE,arr,dump
from .protected_curve_check import context
from .asymmetric_readout import PARENT,LEVERAGE,Kernel as OldKernel,evaluate_stocks
from .forward_leverage import relabel
from .leverage_scan import combine

OUT=HERE/'runs/m2_protected_variance_001'
ENDS=np.array([27,55,57,85,118,148,153,183])/365
MODES=np.array([7,6,5])

def envelope(t):
    z=np.tanh((np.asarray(t)[...,None]-ENDS)/(10/365));f=z*z
    return np.prod(f,axis=-1),sum(2*z[...,j]*(1-f[...,j])/(10/365)*np.prod(np.delete(f,j,axis=-1),axis=-1) for j in range(8))

def correction(t,u,k):
    a,ap=envelope(t);z=np.tanh(u);d=1-z*z
    return (ap*z-a*k*(u+2*z)*d)/(365/10+k*(.5+4/(3*np.sqrt(3))))

class Kernel:
    def __init__(self,path):
        self.old=OldKernel(path);self.k=self.old.k
        self.k.protected_stock.argtypes=self.k.asymmetric_stock.argtypes
        P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        self.k.protected_probe.argtypes=[ct.c_int]+[P]*4
    def probe(self,t,u,k):
        t,u,k=np.broadcast_arrays(t,u,k);out=np.zeros((t.size,3))
        self.k.protected_probe(t.size,arr(t.ravel()),arr(u.ravel()),arr(k.ravel()),out)
        return out.reshape(t.shape+(3,))
    def stocks(self,d,pairs,p,x,coeff,rates,seed):
        assert np.max(np.sum(abs(coeff),axis=1))<=1+1e-14
        models=arr(np.tile(x,(len(coeff),1)));out=np.zeros((2*pairs,int(d['expiries']),len(coeff)+1,2))
        self.k.protected_stock(len(d['dt']),pairs,d['projection'].shape[1],int(d['expiries']),len(coeff),d['transition'],d['decay'],d['root'],d['cross'],d['mean'],d['projection'],d['time'],d['dt'],d['seg'],d['coarse'],d['reference'],d['expiry'],d['lognorm'],arr(p),models,arr(coeff),arr(rates),1.,int(seed),out)
        assert np.isfinite(out).all() and (out>0).all()
        return out

def audit(k,d,p,x,rates):
    rng=np.random.default_rng(2026102901);t=rng.uniform(0,1030/365,1000);u=rng.normal(0,3,1000);r=rng.choice(rates,1000)
    probe=k.probe(t,u,r);a,ap=envelope(t);err=float(np.max(abs(probe-np.column_stack([a,ap,correction(t,u,r)]))))
    assert err<1e-11
    rec=[]
    for rate in rates:
        for order in [64,96]:
            z,w=hermgauss(order);z=z*np.sqrt(2);w=w/np.sqrt(np.pi)
            for day in [27,55,118,153]:
                T=day/365
                for state in [-4.,-1.,0.,1.,4.]:
                    def integrand(s):
                        lag=max(0,s-T);v=np.exp(-rate*lag)*state+np.sqrt(-np.expm1(-2*rate*lag))*z
                        return float(k.probe(np.full(order,s),v,np.full(order,rate))[:,2]@w)
                    val,e=quad(integrand,T,T+30/365,epsabs=1e-11,points=ENDS[(ENDS>T)&(ENDS<T+30/365)].tolist(),limit=100)
                    rec.append(dict(rate=float(rate),GH=order,day=day,state=state,integral=val,quad_error=e))
    residual=max(abs(r['integral']) for r in rec);assert residual<1e-8
    # Zero-amplitude compiled stock replay, on a prefix with a forced endpoint.
    dd={key:value for key,value in d.items()};n=97
    for key in ['transition','mean','time','dt','seg','coarse','reference','expiry','lognorm']:dd[key]=d[key][:n].copy()
    dd['expiry'][:]=-1;dd['expiry'][-1]=0;dd['expiries']=np.array(1)
    new=k.stocks(dd,8,p,x,arr([[0,0,0]]),rates,918276)
    old=k.old.stocks(dd,8,p,arr([x]),arr([0]),np.ones((n,1)),1.,918276)
    replay=float(abs(new-old).max());assert replay<1e-13
    result=dict(compiled_formula_max_error=err,max_conditional_integral_error=residual,zero_amplitude_stock_replay_max_error=replay,conditional_records=rec,global_h_bound=1.,L1_amplitude_budget=float((.05+.07*x[19])**2/2),variance_lower_bound=float((.05+.07*x[19])**2/2),protection='four quoted expiry VIX values, same full Gaussian filtration; continuous time only; inherits baseline quote and forward errors',no_clipping=True,no_new_Lean_claim=True)
    dump(OUT/'implementation_audit.json',result);print('AUDIT PASS',err,residual,replay,flush=True)

def stage(k,d,p,x,rates,coeff,name,pairs,seedtag):
    folder=OUT/name;folder.mkdir();pieces=[];seeds=[]
    for b in range(8):
        seed=int(np.random.SeedSequence([2026102907,seedtag,b]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed);tick=time.perf_counter()
        pieces.append(k.stocks(d,pairs,p,x,coeff,rates,seed));print(name,'batch',b,'seconds',round(time.perf_counter()-tick,2),flush=True)
    stocks=combine(pieces);np.savez_compressed(folder/'stocks.npz',stocks=stocks,seeds=np.array(seeds,dtype=np.uint64),coeff=coeff)
    models=arr(np.tile(x,(len(coeff),1)));q=pd.read_csv(PARENT/'spx_membership.csv')
    evaluate_stocks(stocks,np.load(LEVERAGE/'reference.npz')['stats'],models,np.arange(len(coeff),dtype=float),q,folder/'SPX');scores=relabel(folder/'SPX');return scores

def run():
    OUT.mkdir(exist_ok=False);(OUT/'source').mkdir()
    sources=[Path(__file__),HERE/'m2_protected_variance.cpp',HERE/'asymmetric_readout.cpp',HERE/'leverage_scan.cpp',HERE/'two_speed_precision_kernel.cpp',PARENT/'selection.json']
    for s in sources:shutil.copy2(s,OUT/'source'/s.name)
    dump(OUT/'source_hashes.json',{str(s):hashlib.sha256(s.read_bytes()).hexdigest() for s in sources})
    lib=OUT/'protected.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'m2_protected_variance.cpp'),'-o',str(lib)],check=True)
    k=Kernel(lib);prior,dated,norm,normalizer_audit=context();p=arr(prior['selected']);x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);rates=arr(dated.layout.rates[MODES]);d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'))
    coeff=arr([[0,0,0]]+[list(sign*np.eye(3)[j]*strength) for j in range(3) for strength in [.5,1.] for sign in [-1,1]])
    dump(OUT/'design.json',dict(mode_indices=MODES.tolist(),half_life_days=(np.log(2)*365/rates).tolist(),coefficients=coeff.tolist(),baseline=x.tolist(),normalizer=normalizer_audit,OU_dynamics_unchanged=True,quoted_VIX_expiries=[27,55,118,153],window_days=30,envelope_width_days=10,positivity_budget='sum(abs(epsilon))<=floor/2',curve_refit=False,adopted=False))
    audit(k,d,p,x,rates)
    scores=stage(k,d,p,x,rates,coeff,'pilot',32,101)
    fine=[s for s in scores if s['steps']==256];best=min(fine,key=lambda s:s['loss'])['candidate']
    dump(OUT/'selection.json',dict(pilot_best=best,coefficients=coeff[best].tolist(),criterion='minimum all-primary SPX interval loss on pilot fine grid',selection_frozen_before_fresh=True))
    selected=arr([coeff[0],coeff[best]])
    fresh=stage(k,d,p,x,rates,selected,'fresh',64,202)
    dump(OUT/'run.json',dict(status='pricing completed; result analysis follows',pilot_paths=512,fresh_paths=1024,selection=best,pilot_scores=scores,fresh_scores=fresh,production_baseline_unchanged=True,adopted=False))
    print('PROTECTED VARIANCE PRICING COMPLETE',flush=True)

if __name__=='__main__':run()
