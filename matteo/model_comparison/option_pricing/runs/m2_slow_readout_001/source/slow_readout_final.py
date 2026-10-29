"""Slow downside signals from existing OU modes: paired VIX/SPX research pilot."""
import ctypes as ct
import hashlib,json,shutil,subprocess,time
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from numpy.polynomial.hermite import hermgauss
from .targeted_skew import Kernel as ParentKernel
from . import asymmetric_readout as old
from .asymmetric_readout import PARENT,LEVERAGE,VIX,fit_vix,save_vix,evaluate_stocks,load_outer,stress_cond
from .full_readout_fit import HERE,arr,dump
from .protected_curve_check import context
from .leverage_scan import combine
from .vix_only_curve_fit import PRIOR,ACTIVE,DAYS
from .maturity_curve import window_nodes
from .forward_leverage import relabel
from esn_eq.states import symmetric_square_root
OUT=HERE/'runs/m2_slow_readout_001'
LEVELS=arr([.5,1])
ALPHA=arr(np.r_[0,np.tile(LEVELS,3)])

class Kernel(ParentKernel):
    def __init__(self,path):
        super().__init__(path);P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.k.slow_moments.argtypes=[ct.c_int]*4+[P]*6+[ct.c_double]*2+[P]*4+[ct.c_int,P]
        self.k.slow_stock.argtypes=[ct.c_int]*5+[I]+[P]*7+[I]*4+[P]*5+[I,ct.c_uint64,P]
    def moments4(self,mean,cov,p,center,alpha,order):
        means=arr(mean).reshape(-1,4);assert abs(cov[0,1:]).max()<1e-10
        z,w=hermgauss(order);z=arr(z*np.sqrt(2));w=arr(w/np.sqrt(np.pi));r=symmetric_square_root(cov[1:3,1:3]);yy,ff=np.meshgrid(z,z,indexing='ij');off=arr(np.column_stack([yy.ravel(),ff.ravel()])@r.T)
        beta=cov[3,1:3]@np.linalg.pinv(cov[1:3,1:3],rcond=1e-12);res=max(0,float(cov[3,3]-beta@cov[1:3,3]));out=np.zeros((len(means),len(alpha),2))
        self.k.slow_moments(len(means),order,len(off),len(alpha),means,arr(z*np.sqrt(max(0,cov[0,0]))),w,off,arr(np.outer(w,w).ravel()),arr(p),float(center),1.,arr(alpha),arr(beta),arr(z*np.sqrt(res)),w,order,out)
        # Direct marginal expectation avoids a second approximation for E[D].
        out[:,:,0]=(np.exp(np.logaddexp(0,-(means[:,3,None]-center+z[None,:]*np.sqrt(max(0,cov[3,3]))))[:,:,None]*arr(alpha)[None,None,:])*w[None,:,None]).sum(1)
        return out
    def stocks(self,d,pairs,p,models,alpha,moments,indices,seed):
        out=np.zeros((2*pairs,int(d['expiries']),len(alpha)+1,2));self.k.slow_stock(len(d['dt']),pairs,d['projection'].shape[1],int(d['expiries']),len(alpha),d['transition'],d['decay'],d['root'],d['cross'],d['mean'],d['projection'],d['time'],d['dt'],d['seg'],d['coarse'],d['reference'],d['expiry'],d['lognorm'],arr(p),arr(models),arr(alpha),arr(moments),np.ascontiguousarray(indices,dtype=np.int32),int(seed),out);return out

class Law:
    def __init__(self,dated,index):
        self.d=dated;self.P=np.vstack([dated.projection,np.eye(len(dated.layout.rates))[index]])
    def __call__(self,t,state):
        d=self.d;r=d.layout.rates;dec=np.exp(-r*t);shift=d.q_drift*(-np.expm1(-r*t))/r;cov=d.stationary*(-np.expm1(-(r[:,None]+r[None,:])*t));return (np.asarray(state)*dec+shift)@self.P.T,self.P@cov@self.P.T

class Means:
    def __init__(self,k,laws,initial,p,norm):self.k,self.laws,self.initial,self.p,self.norm=k,laws,initial,p,norm;self.tau=1/max(laws[0].d.layout.rates)
    def direct(self,t,order=48):
        result=np.ones((7,2))
        for j,law in enumerate(self.laws):
            m,c=law(float(t),self.initial);v=self.k.moments4(m,c,self.p,m[3],LEVELS,order)[0];v[:,1]/=np.exp(self.norm(t)[0]);result[1+2*j:3+2*j]=v
        return result
    def build(self):
        knots=np.linspace(0,np.log1p((1060/365)/self.tau),129)
        for _ in range(7):
            values=np.array([np.log(self.direct(self.tau*np.expm1(t))) for t in knots]);s=CubicSpline(knots,values,axis=0);mid=(knots[:-1]+knots[1:])/2;actual=np.array([self.direct(self.tau*np.expm1(t)) for t in mid]);err=np.max(abs(np.expm1(s(mid)-np.log(actual))),axis=(1,2));bad=err>2e-6
            if not bad.any():break
            knots=np.unique(np.r_[knots,mid[bad]])
        assert err.max()<2e-6;self.spline=s;self.audit=dict(knots=len(knots),GH=48,midpoint_relative_error=float(err.max()));return self
    def __call__(self,t):return np.exp(self.spline(np.log1p(np.asarray(t)/self.tau)))

def exposures(k,laws,initial,p,norm,M,x,draws,order,timeorder,folder):
    folder.mkdir(exist_ok=False);ends=np.array(json.loads((PRIOR/'run.json').read_text())['ends_years']);L=np.log(365);hf=np.exp(L*x[21]);hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf]);b=np.array([5*x[16],5*x[20]]);H={}
    for day,a in draws.items():
        tick=time.perf_counter();T=day/365;times,tw=window_nodes(T,ends,1/max(laws[0].d.layout.rates),timeorder);out=np.zeros((7,len(a['weights']),16));stateZ=np.sqrt(-np.expm1(-2*rates*T))[None,:]*np.column_stack([a['z_normals'],a['fast_normals']])
        for t,wt in zip(times,tw):
            lag=t-T;EE=np.ones((len(a['weights']),7,2));oldmean,oldcov=laws[0].d.gaussian_projection_law(float(lag),a['states']);oldE=k.moments(oldmean,oldcov,p,0,1,arr([0]),order)[:,0,1]/np.exp(norm(t)[0]);EE[:,0,1]=oldE
            for j,law in enumerate(laws):
                mean,cov=law(float(lag),a['states']);uncond,_=law(float(t),initial);E=k.moments4(mean,cov,p,uncond[3],LEVELS,order);E[:,:,1]/=np.exp(norm(t)[0]);EE[:,1+2*j:3+2*j]=E
            stress=k.stress(stateZ*np.exp(-rates*lag),np.sqrt(-np.expm1(-2*rates*lag)),b,-np.expm1(-2*rates*t),x[22],x[17],M(t),24 if order<32 else 32);response=EE[:,:,1]*stress[:,:,0]+EE[:,:,0]*stress[:,:,1];response[:,0]=(1-x[17])*oldE+x[17]*stress_cond(a,T,t,x)
            seg=min(np.searchsorted(ends,t,side='right'),15);out[:,:,seg]+=wt*response.T
        H[day]=out;np.savez_compressed(folder/f'd{day}.npz',H=out,weights=a['weights'],nets=a['nets'],times=times,time_weights=tw);print('SLOW VIX',folder.name,day,'seconds',round(time.perf_counter()-tick,2),flush=True)
    return H

def run():
    begin=time.perf_counter();OUT.mkdir(exist_ok=False);(OUT/'source').mkdir();sources=[__file__,HERE/'slow_readout.cpp',HERE/'targeted_skew.cpp',HERE/'asymmetric_readout.cpp',HERE/'leverage_scan.cpp',HERE/'two_speed_precision_kernel.cpp',PARENT/'selection.json',PRIOR/'run.json'];hashes={str(s):hashlib.sha256(open(s,'rb').read()).hexdigest() for s in sources}
    for s in sources:shutil.copy2(s,OUT/'source'/__import__('pathlib').Path(s).name)
    lib=OUT/'slow.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'slow_readout.cpp'),'-o',str(lib)],check=True);k=Kernel(lib);prior,dated,norm,na=context();initial=arr(prior['initial_state']);p=arr(prior['selected']);x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);half=np.log(2)*365/dated.layout.rates[:12];indices=np.array([np.argmin(abs(np.log(half/t))) for t in [2,8,32]]);assert len(set(indices))==3;laws=[Law(dated,i) for i in indices];M=Means(k,laws,initial,p,norm).build();dump(OUT/'design.json',dict(target_half_lives_days=[2,8,32],mode_indices=indices.tolist(),actual_half_lives_days=half[indices].tolist(),alpha=ALPHA.tolist(),baseline=x.tolist(),source_hashes=hashes,moment_spline=M.audit,normalizer=na,stationary_slow_scale=1,OU_dynamics_unchanged=True,fast_M2_readout_unchanged=True,normalization='E[Hnew | stress]=(1-w)+w A',forward_scope='four training-sample equalities, independent checking required',new_bank=False,adopted=False));print('DESIGN',indices,half[indices],flush=True)
    probes=np.geomspace(1e-6,1030/365,11);aa=np.array([M.direct(t,48) for t in probes]);bb=np.array([M.direct(t,64) for t in probes]);dump(OUT/'moment_order_audit.json',dict(max_GH48_64_relative_error=float(abs(aa/bb-1).max()),times=probes.tolist()))
    qv=pd.read_csv(old.OUT/'VIX_membership.csv');f=pd.read_csv(old.OUT/'forward_targets.csv');draws=load_outer(64,20261028613);H=exposures(k,laws,initial,p,norm,M,x,draws,24,4,OUT/'training_exposures');v=VIX(H,draws,x,qv,f);ledger=[];old.ALPHA=ALPHA;(OUT/'training_VIX').mkdir()
    for i in range(7):
        y,rec=fit_vix(v,i,x[ACTIVE]);rec.update(candidate=i,half_life_days=0 if i==0 else float(half[indices[(i-1)//2]]));ledger.append(rec);save_vix(v,i,y,f'a{i}',OUT/'training_VIX');dump(OUT/'VIX_fit_ledger.json',ledger);print('SLOW FIT',i,rec['loss'],rec['max_forward_gap'],flush=True)
    models=[x.copy() for _ in range(7)];names=['baseline']+[f'slow_{half[idx]:.2f}d_alpha{a:g}_fixed' for idx in indices for a in LEVELS];mi=list(range(7));alph=list(ALPHA);modes=[int(indices[0])]+[int(idx) for idx in indices for _ in range(2)]
    for i in range(1,7):
        xx=x.copy();xx[ACTIVE]=ledger[i]['y'];models.append(xx);names.append(names[i].replace('fixed','refit'));mi.append(i);alph.append(ALPHA[i]);modes.append(modes[i])
    models=arr(models);alph=arr(alph);dump(OUT/'candidates.json',dict(models=models.tolist(),alpha=alph.tolist(),names=names,mean_indices=mi,mode_indices=modes));d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));moments=arr(M(d['time'])[:,mi]);np.savez_compressed(OUT/'component_means.npz',time=d['time'],moments=moments)
    # All 13 candidates share exactly the same factors and stock driver.
    pieces=[];seeds=[]
    for b in range(8):
        seed=int(np.random.SeedSequence([20261028613,501,b]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed);tick=time.perf_counter();pieces.append(k.stocks(d,32,p,models,alph,moments,modes,seed));print('SLOW STOCK',b,'seconds',round(time.perf_counter()-tick,2),flush=True)
    stocks=combine(pieces);np.savez_compressed(OUT/'stocks.npz',stocks=stocks,seeds=np.array(seeds,dtype=np.uint64));q=pd.read_csv(PARENT/'spx_membership.csv');scores=evaluate_stocks(stocks,np.load(LEVERAGE/'reference.npz')['stats'],models,alph,q,OUT/'SPX');scores=relabel(OUT/'SPX')
    checking={day:dict(np.load(old.OUT/f'fresh_VIX_draws/d{day}.npz')) for day in DAYS};HH=exposures(k,laws,initial,p,norm,M,x,checking,32,8,OUT/'checking_exposures');cv=VIX(HH,checking,x,qv,f);(OUT/'checking_VIX').mkdir();vr=[]
    for c in range(13):
        rec=save_vix(cv,mi[c],models[c,ACTIVE],f'c{c}',OUT/'checking_VIX');rec.update(candidate=c,name=names[c]);vr.append(rec)
    dump(OUT/'checking_VIX/scores.json',vr);dump(OUT/'run.json',dict(status='completed pilot; audits follow',wall_seconds=time.perf_counter()-begin,SPX_scores=scores,VIX_scores=vr,SPX_paths=len(stocks),fine_steps=256,coarse_steps=128,VIX_check_draws_previously_examined=True,production_baseline_unchanged=True));print('SLOW COMPLETE',flush=True)
if __name__=='__main__':
    import sys
    if '--weak' in sys.argv:
        OUT=HERE/'runs/m2_slow_readout_002';LEVELS=arr([.125,.25]);ALPHA=arr(np.r_[0,np.tile(LEVELS,3)])
    run()
