"""Linear memory of normalized nonlinear variance signals on frozen OU factors."""
import ctypes as ct
import hashlib,json,shutil,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
from .asymmetric_readout import Kernel as OldKernel,PARENT,LEVERAGE,VIX,fit_vix,save_vix,evaluate_stocks,load_outer,stress_cond,OUT as OLDOUT
from . import asymmetric_readout as old
from .full_readout_fit import HERE,arr,dump
from .protected_curve_check import context
from .leverage_scan import combine
from .vix_only_curve_fit import PRIOR,ACTIVE,DAYS
from .maturity_curve import window_nodes
from .forward_leverage import relabel
from esn_eq.states import symmetric_square_root
OUT=HERE/'runs/m2_memory_filter_001'
HALF=arr([3,10,30]);RATES=arr(np.log(2)*365/HALF)
CONFIG=[(0,0.,0)]+[(kind,share,j) for kind in [1,2] for j in range(3) for share in [.25,.5]]

def stress_parameters(x):
    hf=np.exp(np.log(365)*x[21]);hs=hf*(365/hf)**x[18];return np.log(2)*365/np.array([hs,hf]),5*x[[16,20]]

class Kernel(OldKernel):
    def __init__(self,path):
        super().__init__(path);P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.k.memory_stock.argtypes=[ct.c_int]*5+[I]+[P]*7+[I]*4+[P]*4+[I]*2+[P,ct.c_uint64,P]
        self.k.memory_history.argtypes=[ct.c_int]*4+[ct.c_double]+[P]*14+[ct.c_uint64,P]
    def stocks(self,d,pairs,p,models,config,seed):
        out=np.zeros((2*pairs,int(d['expiries']),len(config)+1,2));kind,share,index=zip(*config);self.k.memory_stock(len(d['dt']),pairs,d['projection'].shape[1],int(d['expiries']),len(config),d['transition'],d['decay'],d['root'],d['cross'],d['mean'],d['projection'],d['time'],d['dt'],d['seg'],d['coarse'],d['reference'],d['expiry'],d['lognorm'],arr(p),arr(models),arr(share),np.array(kind,np.int32),np.array(index,np.int32),RATES,int(seed),out);return out
    def history(self,b,terminal,zterminal,p,x,seed):
        n=len(terminal);out=np.zeros((n,3,3,2));self.k.memory_history(b['nt'],n,24,3,float(b['dt']),b['decay'],b['root'],b['projection'],b['mean'],b['bridge'],arr(terminal),b['lognorm'],b['stressbridge'],arr(zterminal),b['stress_q'],arr(p),arr(x),RATES,int(seed),out);return out

def history_basis(dated,initial,norm,x,day,steps=128):
    T=day/365;nt=int(np.ceil(T*252*steps/2))*2;h=T/nt;t=np.arange(nt)*h;r=dated.layout.rates[:24];P=dated.projection[:,:24];S=dated.stationary[:24,:24];CT=S*(-np.expm1(-(r[:,None]+r[None,:])*T));ev,Q=np.linalg.eigh(CT);keep=ev>ev.max()*1e-11;inverse=(Q[:,keep]/ev[keep])@Q[:,keep].T
    C=S[None,:,:]*(-np.expm1(-(r[None,:,None]+r[None,None,:])*t[:,None,None]));cross=C*np.exp(-r[None,None,:]*(T-t[:,None,None]));bridge=arr(P[None,:,:]@cross@inverse);dec=np.exp(-r*h);cov=S*(-np.expm1(-(r[:,None]+r[None,:])*h));root=symmetric_square_root(cov);qmean=np.exp(-t[:,None]*r)*initial[:24]+dated.q_drift[:24]*(-np.expm1(-t[:,None]*r))/r;mean=arr(qmean@P.T);terminalmean=np.exp(-T*r)*initial[:24]+dated.q_drift[:24]*(-np.expm1(-T*r))/r;a,_=stress_parameters(x);q=-np.expm1(-2*a[None,:]*t[:,None]);qb=-np.expm1(-2*a*T);sb=arr(q*np.exp(-a[None,:]*(T-t[:,None]))/qb)
    errors=[]
    for j in np.linspace(0,nt-1,13,dtype=int):
        expected=P@cross[j];got=bridge[j]@CT;errors.append(float(abs(got-expected).max()))
    audit=dict(day=int(day),steps_per_trading_day=steps,nt=nt,dt=h,pseudoinverse_relative_cutoff=1e-11,retained_terminal_rank=int(keep.sum()),discarded_terminal_variance=float(ev[~keep].sum()),max_projected_bridge_crosscov_error=max(errors),terminal_projection_error=float(abs(P@CT@inverse-P).max()),innovation_covariance_error=float(abs(root@root.T-cov).max()))
    assert max(errors)<2e-7,audit
    return dict(nt=nt,dt=h,decay=arr(dec),root=arr(root),projection=arr(P),mean=mean,bridge=bridge,lognorm=arr(norm(t)[:,0]),stressbridge=sb,stress_q=arr(q),terminalmean=arr(terminalmean)),audit

def histories(k,dated,initial,norm,p,x,draws,folder,seedbase,steps=128):
    folder.mkdir(exist_ok=False);out={};audits=[];a,_=stress_parameters(x)
    for day,z in draws.items():
        tick=time.perf_counter();b,audit=history_basis(dated,initial,norm,x,day,steps);zt=np.sqrt(-np.expm1(-2*a*day/365))[None,:]*np.column_stack([z['z_normals'],z['fast_normals']]);out[day]=k.history(b,z['states'][:,:24]-b['terminalmean'],zt,p,x,seedbase+int(day));np.savez_compressed(folder/f'd{day}.npz',memory=out[day],weights=z['weights'],nets=z['nets']);audits.append(audit);print('FILTER HISTORY',folder.name,day,round(time.perf_counter()-tick,2),flush=True)
    dump(folder/'bridge_audit.json',audits);return out

def exposures(k,dated,initial,p,norm,x,draws,memory,folder,order=32,timeorder=8,lagcount=32,grid=1):
    folder.mkdir(exist_ok=False);ends=np.array(json.loads((PRIOR/'run.json').read_text())['ends_years']);H={}
    for day,a in draws.items():
        tick=time.perf_counter();T=day/365;times,tw=window_nodes(T,ends,1/max(dated.layout.rates),timeorder);times=np.asarray(times);lags=np.unique(np.r_[0,np.geomspace(1e-7,30/365,lagcount),times-T]);n=len(a['weights']);E=[]
        for lag in lags:
            t=T+lag;mean,cov=dated.gaussian_projection_law(float(lag),a['states']);ER=k.moments(mean,cov,p,0,1,arr([0]),order)[:,0,1]/np.exp(norm(t)[0]);rates,b=stress_parameters(x);zz=np.column_stack([a['z_normals'],a['fast_normals']]);loading=b*np.sqrt(-np.expm1(-2*rates*T))*np.exp(-rates*lag);EA=np.exp(zz*loading-.5*loading**2);E.append(np.column_stack([ER,EA]))
        E=arr(E);mem=memory[day][:,:,:,grid].copy();MM=[mem.copy()]
        for j,h in enumerate(np.diff(lags)):
            z=RATES*h;e=np.exp(-z);u=-np.expm1(-z);linear=1-u/z;mem=e[None,:,None]*mem+u[None,:,None]*E[j,:,None,:]+linear[None,:,None]*(E[j+1]-E[j])[:,None,:];MM.append(mem.copy())
        MM=arr(MM);out=np.zeros((len(CONFIG),n,16))
        for t,wt in zip(times,tw):
            j=np.searchsorted(lags,t-T);ER,AS,AF=E[j].T;A=(1-x[22])*AS+x[22]*AF;baseline=(1-x[17])*ER+x[17]*A;seg=min(np.searchsorted(ends,t,side='right'),15)
            for c,(kind,share,node) in enumerate(CONFIG):
                mm=MM[j,:,node];filtered=(1-x[17])*mm[:,0]+x[17]*(A if kind==1 else (1-x[22])*mm[:,1]+x[22]*mm[:,2]);out[c,:,seg]+=wt*((1-share)*baseline+share*filtered)
        H[day]=out;np.savez_compressed(folder/f'd{day}.npz',H=out,times=times,time_weights=tw,lags=lags,weights=a['weights'],nets=a['nets']);print('FILTER EXPOSURE',folder.name,day,round(time.perf_counter()-tick,2),flush=True)
    return H

def run():
    OUT.mkdir(exist_ok=False);(OUT/'source').mkdir();sources=[Path(__file__),HERE/'memory_filter.cpp',HERE/'MEMORY_FILTER_PLAN.md',HERE/'asymmetric_readout.cpp',HERE/'leverage_scan.cpp',HERE/'two_speed_precision_kernel.cpp',PARENT/'selection.json',PRIOR/'run.json'];hashes={str(s):hashlib.sha256(s.read_bytes()).hexdigest() for s in sources}
    for s in sources:shutil.copy2(s,OUT/'source'/s.name)
    subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'memory_filter.cpp'),'-o',str(OUT/'memory.dylib')],check=True);k=Kernel(OUT/'memory.dylib');prior,dated,norm,na=context();initial=arr(prior['initial_state']);p=arr(prior['selected']);x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);dump(OUT/'design.json',dict(source_hashes=hashes,baseline=x.tolist(),half_lives_days=HALF.tolist(),config=CONFIG,OU_dynamics_unchanged=True,filter_initial_value=1,scope='normalized signal filtering; not literal variance-level feedback',population_forward_variance_mean_preserved=True,exact_population_VIX_forwards_not_claimed=True,adopted=False,old_normalizer_audit=na));draws=load_outer(64,20261028821);mem=histories(k,dated,initial,norm,p,x,draws,OUT/'training_histories',20261028821);H=exposures(k,dated,initial,p,norm,x,draws,mem,OUT/'training_exposures');qv=pd.read_csv(OLDOUT/'VIX_membership.csv');f=pd.read_csv(OLDOUT/'forward_targets.csv');v=VIX(H,draws,x,qv,f);old.ALPHA=arr([c[1] for c in CONFIG]);ledger=[];(OUT/'training_VIX').mkdir()
    for i in range(len(CONFIG)):
        y,rec=fit_vix(v,i,x[ACTIVE]);rec['configuration']=CONFIG[i];ledger.append(rec);save_vix(v,i,y,f'c{i}',OUT/'training_VIX');dump(OUT/'VIX_fit_ledger.json',ledger);print('FILTER FIT',i,rec['loss'],rec['max_forward_gap'],flush=True)
    models=[x.copy() for _ in CONFIG];cfg=list(CONFIG);indices=list(range(len(CONFIG)));names=['baseline']+[f'{"M2" if kind==1 else "whole"}_{HALF[j]:g}d_share{share:g}_fixed' for kind,share,j in CONFIG[1:]]
    for i,c in enumerate(CONFIG):
        xx=x.copy();xx[ACTIVE]=ledger[i]['y'];models.append(xx);cfg.append(c);indices.append(i);names.append('baseline_refit' if i==0 else names[i].replace('fixed','refit'))
    models=arr(models);dump(OUT/'candidates.json',dict(models=models.tolist(),config=cfg,exposure_indices=indices,names=names));d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));pieces=[];seeds=[]
    for b in range(8):
        seed=int(np.random.SeedSequence([20261028821,501,b]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed);tick=time.perf_counter();pieces.append(k.stocks(d,32,p,models,cfg,seed));print('FILTER STOCK',b,round(time.perf_counter()-tick,2),flush=True)
    stocks=combine(pieces);np.savez_compressed(OUT/'stocks.npz',stocks=stocks,seeds=np.array(seeds,dtype=np.uint64));q=pd.read_csv(PARENT/'spx_membership.csv');evaluate_stocks(stocks,np.load(LEVERAGE/'reference.npz')['stats'],models,arr([c[1] for c in cfg]),q,OUT/'SPX');scores=relabel(OUT/'SPX');feasible=[s for s in scores if s['steps']==256 and s['candidate']>=13 and ledger[indices[s['candidate']]]['eligible'] and ledger[indices[s['candidate']]]['loss']<=ledger[0]['loss']+.1];chosen=min(feasible,key=lambda s:s['IV_MAE_pp'])['candidate'];dump(OUT/'selection.json',dict(candidate=chosen,name=names[chosen],selected_before_fresh_stock_and_VIX=True,criterion='four training equalities; VIX training loss within alpha0-refit+.1; smallest pilot fine-grid IV MAE among refits',adopted=False));dump(OUT/'run.json',dict(status='pilot complete; independent checking follows',SPX_paths=512,SPX_scores=scores,adopted=False));print('FILTER PILOT COMPLETE',chosen,names[chosen],flush=True)
if __name__=='__main__':run()
