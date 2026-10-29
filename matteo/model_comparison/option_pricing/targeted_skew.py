"""State-dependent stock loading and stress-regime-normalized readout pilots."""
import ctypes as ct
import hashlib,json,shutil,subprocess,time
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from numpy.polynomial.hermite import hermgauss
from .asymmetric_readout import Kernel as OldKernel,OUT as OLD,ALPHA as OLD_ALPHA,PARENT,LEVERAGE,VBASE,VIX,fit_vix,save_vix,evaluate_stocks,load_outer,stress_cond
from .full_readout_fit import HERE,arr,dump
from .protected_curve_check import context
from .leverage_scan import combine,stock
from .vix_only_curve_fit import PRIOR,ACTIVE,DAYS
from .maturity_curve import window_nodes
from esn_eq.states import symmetric_square_root

OUT=HERE/'runs/m2_targeted_skew_001'
ALPHA=arr([0,.5,1,1.5,2])

class Kernel(OldKernel):
    def __init__(self,lib):
        super().__init__(lib);P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.k.regime_stress.argtypes=[ct.c_int]*3+[P]*4+[ct.c_double]*2+[P]*4
        self.k.targeted_stock.argtypes=[ct.c_int]*6+[I]+[P]*10+[I]*4+[P]*4+[I,P,ct.c_double,ct.c_uint64,P]
    def stress(self,means,sd,b,q,r,w,m,order):
        z,weights=hermgauss(order);n=len(means);na=len(m);out=np.zeros((n,na,2));self.k.regime_stress(n,na,order,arr(means),arr(sd),arr(b),arr(q),float(r),float(w),arr(m),arr(np.sqrt(2)*z),arr(weights/np.sqrt(np.pi)),out);return out
    def stocks(self,d,pairs,p,models,alpha,mode,moments,scale,seed):
        nb=d['diffusion'].shape[1]+2;out=np.zeros((2*pairs,int(d['expiries']),len(alpha)+1,2));self.k.targeted_stock(len(d['dt']),pairs,d['projection'].shape[1],nb,int(d['expiries']),len(alpha),d['transition'],d['decay'],d['root'],d['cross'],d['bmean'],d['broot'],d['diffusion'],d['mean'],d['projection'],d['time'],d['dt'],d['seg'],d['coarse'],d['reference'],d['expiry'],d['lognorm'],arr(p),arr(models),arr(alpha),np.ascontiguousarray(mode,dtype='int32'),arr(moments),float(scale),int(seed),out);return out

def primitive_bridge(d,dated):
    keep=np.flatnonzero(np.any(dated.projection!=0,axis=0));load=dated.layout.driver_loadings[keep];used=np.flatnonzero(np.any(load!=0,axis=0));assert used[0]==0;L=load[:,used];rates=dated.layout.rates[keep];a=[];roots=[];errors=[];eigenvalues=[]
    for tr in range(len(d['root'])):
        h=float(d['dt'][np.flatnonzero(d['transition']==tr)[0]]);C=L[:,1:].T*np.sqrt(2*rates)*(-np.expm1(-rates*h))/rates;A=np.linalg.lstsq(d['root'][tr],C.T,rcond=1e-6)[0].T;rem=h*np.eye(len(used)-1)-A@A.T;residual=symmetric_square_root(rem);actual=A@d['root'][tr].T;scale=np.sqrt(h*np.diag(d['root'][tr]@d['root'][tr].T));error=float(np.max(abs(actual-C)/scale));assert error<1e-8,error;assert np.max(abs(A@A.T+residual@residual.T-h*np.eye(len(used)-1)))/h<1e-9;a.append(A);roots.append(residual);errors.append(error);eigenvalues.append(float(np.linalg.eigvalsh(rem).min()/h))
    d.update(bmean=arr(a),broot=arr(roots),diffusion=arr(d['projection']@(np.sqrt(2*rates)[:,None]*L)))
    return dict(old_primitive_names=[dated.layout.primitive_names[i] for i in used],maximum_standardized_Brownian_OU_cross_error=max(errors),minimum_relative_residual_eigenvalue=min(eigenvalues),OU_stream_unchanged=True,conditional_inverse_relative_cutoff=1e-6,scope='finite-precision augmentation of existing root; reconstructed covariance checked, not symbolic exact arithmetic')

class ComponentMeans:
    def __init__(self,k,dated,initial,p,norm,alpha,scale):self.k,self.dated,self.initial,self.p,self.norm,self.alpha,self.scale=k,dated,initial,p,norm,alpha,scale;self.tau=1/max(dated.layout.rates)
    def direct(self,t,order=48):
        mean,cov=self.dated.gaussian_projection_law(float(t),self.initial);value=self.k.moments(mean,cov,self.p,mean[1],self.scale,self.alpha,order)[0];value[:,1]/=np.exp(self.norm(t)[0]);value[self.alpha==0]=1.;return value
    def build(self):
        knots=np.linspace(0,np.log1p((1060/365)/self.tau),257)
        for _ in range(6):
            values=np.array([np.log(self.direct(self.tau*np.expm1(t))) for t in knots]);s=CubicSpline(knots,values,axis=0);mid=(knots[:-1]+knots[1:])/2;actual=np.array([self.direct(self.tau*np.expm1(t)) for t in mid]);error=np.max(abs(np.expm1(s(mid)-np.log(actual))),axis=(1,2));bad=error>2e-7
            if not bad.any():break
            knots=np.unique(np.r_[knots,mid[bad]])
        assert error.max()<2e-7;self.spline=s;self.audit=dict(knots=len(knots),GH=48,maximum_midpoint_relative_error=float(error.max()));return self
    def __call__(self,t):return np.exp(self.spline(np.log1p(np.asarray(t)/self.tau)))

def exposures(k,dated,initial,p,norm,M,x,draws,alpha,scale,oldorder,stressorder,timeorder,folder):
    folder.mkdir(exist_ok=False);ends=np.array(json.loads((PRIOR/'run.json').read_text())['ends_years']);L=np.log(365);hf=np.exp(L*x[21]);hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf]);b=np.array([5*x[16],5*x[20]]);H={}
    for day,a in draws.items():
        tick=time.perf_counter();T=day/365;times,tw=window_nodes(T,ends,1/max(dated.layout.rates),timeorder);out=np.zeros((len(alpha),len(a['weights']),16));stateZ=np.sqrt(-np.expm1(-2*rates*T))[None,:]*np.column_stack([a['z_normals'],a['fast_normals']])
        for t,wt in zip(times,tw):
            lag=t-T;mean,cov=dated.gaussian_projection_law(float(lag),a['states']);unconditional,_=dated.gaussian_projection_law(float(t),initial);E=k.moments(mean,cov,p,unconditional[1],scale,alpha,oldorder);E[:,:,1]/=np.exp(norm(t)[0]);futuremeans=stateZ*np.exp(-rates*lag);sd=np.sqrt(-np.expm1(-2*rates*lag));q=-np.expm1(-2*rates*t);stress=k.stress(futuremeans,sd,b,q,x[22],x[17],M(t),stressorder);response=E[:,:,1]*stress[:,:,0]+E[:,:,0]*stress[:,:,1];zero=np.flatnonzero(alpha==0)
            if len(zero):response[:,zero]=(1-x[17])*E[:,zero,1]+x[17]*stress_cond(a,T,t,x)[:,None]
            seg=min(np.searchsorted(ends,t,side='right'),15);out[:,:,seg]+=wt*response.T
        H[day]=out;np.savez_compressed(folder/f'd{day}.npz',H=out,weights=a['weights'],nets=a['nets'],times=times,time_weights=tw);print('REGIME VIX',folder.name,day,'seconds',round(time.perf_counter()-tick,2),flush=True)
    return H

def make_config(x,ledger):
    models=[x.copy(),x.copy()];alpha=[0.,0.];mode=[0,1];names=['original','adaptive_stock'];indices=[0,0]
    for idx,a in enumerate(ALPHA[1:],1):models.append(x.copy());alpha.append(a);mode.append(2);names.append(f'regime_{a:g}_fixed');indices.append(idx)
    for idx,a in enumerate(ALPHA):
        xx=x.copy();xx[ACTIVE]=ledger[idx]['y'];models.append(xx);alpha.append(a);mode.append(2);names.append(f'regime_{a:g}_refit');indices.append(idx)
    return arr(models),arr(alpha),np.array(mode,np.int32),names,indices

def run():
    begin=time.perf_counter();OUT.mkdir(exist_ok=False);(OUT/'source').mkdir();sources=[__file__,HERE/'targeted_skew.cpp',HERE/'asymmetric_readout.cpp',HERE/'leverage_scan.cpp',HERE/'two_speed_precision_kernel.cpp',HERE/'TARGETED_SKEW_PLAN.md',PARENT/'selection.json',PRIOR/'run.json'];hashes={str(p):hashlib.sha256(open(p,'rb').read()).hexdigest() for p in sources}
    for p in sources:shutil.copy2(p,OUT/'source'/__import__('pathlib').Path(p).name)
    lib=OUT/'targeted.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'targeted_skew.cpp'),'-o',str(lib)],check=True);k=Kernel(lib);prior,dated,norm,oldnormaudit=context();initial=arr(prior['initial_state']);p=arr(prior['selected']);x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);d=dict(np.load(PARENT/'fresh_check/grid_and_transitions.npz'));bridge=primitive_bridge(d,dated);dump(OUT/'primitive_bridge_audit.json',bridge);np.savez_compressed(OUT/'augmented_grid.npz',**d);scale=float(json.loads((OLD/'design.json').read_text())['scale']);M=ComponentMeans(k,dated,initial,p,norm,ALPHA,scale).build();m=M(d['time']);np.savez_compressed(OUT/'component_means.npz',times=d['time'],moments=m,alpha=ALPHA);dump(OUT/'design.json',dict(alpha=ALPHA.tolist(),source_hashes=hashes,baseline=x.tolist(),scale=scale,moment_audit=M.audit,old_normalizer_audit=oldnormaudit,VIX_adaptive_law_unchanged=True,gaussian_factor_dynamics_unchanged=True,regime_identity='E[Hnew|stress state]=(1-w)+wA',no_population_forward_certificate=True,stock_primary_quotes=1294,stock_eligible_quotes=1751,VIX_quotes=217))
    probes=np.geomspace(1e-7,1030/365,13);a=np.array([M.direct(t,48) for t in probes]);b=np.array([M.direct(t,64) for t in probes]);dump(OUT/'moment_order_audit.json',dict(maximum_GH48_vs64_relative_error=float(np.max(abs(a/b-1))),probe_times=probes.tolist()))
    qv=pd.read_csv(OLD/'VIX_membership.csv');f=pd.read_csv(OLD/'forward_targets.csv');qv.to_csv(OUT/'VIX_membership.csv',index=False);f.to_csv(OUT/'forward_targets.csv',index=False);draws=load_outer(128,20261028491);H=exposures(k,dated,initial,p,norm,M,x,draws,ALPHA,scale,24,20,4,OUT/'training_exposures');v=VIX(H,draws,x,qv,f);ledger=[];vf=[];(OUT/'training_VIX').mkdir()
    for idx,a in enumerate(ALPHA):
        y,rec=fit_vix(v,idx,x[ACTIVE]);rec['alpha']=float(a);ledger.append(rec);vf.append(save_vix(v,idx,y,f'a{a:g}',OUT/'training_VIX'));dump(OUT/'VIX_fit_ledger.json',ledger);print('REGIME FIT',rec,flush=True)
    models,alpha,mode,names,indices=make_config(x,ledger);dump(OUT/'candidates.json',dict(models=models.tolist(),alpha=alpha.tolist(),mode=mode.tolist(),names=names,mean_indices=indices,VIX_training=vf));moments=arr(m[:,indices]);spx=pd.read_csv(PARENT/'spx_membership.csv');spx.to_csv(OUT/'SPX_membership.csv',index=False);pieces=[];seeds=np.load(PARENT/'fresh_check/stocks.npz')['seeds']
    for b,seed in enumerate(seeds):
        tick=time.perf_counter();pieces.append(k.stocks(d,32,p,models,alpha,mode,moments,scale,int(seed)));print('TARGETED REPLAY',b,'seconds',round(time.perf_counter()-tick,2),flush=True)
    replay=combine(pieces);old=np.load(LEVERAGE/'replay_stats.npz')['stats'];expected=stock(old[:,:,:2],.8,np.sqrt(1-x[22]),np.sqrt(x[22]));error=float(np.max(abs(replay[:,:,0]-expected)));assert error<2e-10;dump(OUT/'baseline_replay.json',dict(maximum_stock_error=error));np.savez_compressed(OUT/'replay_stocks.npz',stocks=replay,seeds=seeds);ref=np.load(LEVERAGE/'reference.npz')['stats'];scores=evaluate_stocks(replay,ref,models,alpha,spx,OUT/'replay');feasible=[]
    for c in range(6,11):
        s=next(z for z in scores if z['candidate']==c and z['steps']==128);idx=c-6
        if ledger[idx]['eligible'] and ledger[idx]['loss']<=ledger[0]['loss']+.1:feasible.append(s)
    chosen=min(feasible,key=lambda s:s['IV_MAE_pp'])['candidate'] if feasible else 6;dump(OUT/'selection.json',dict(regime_candidate=int(chosen),name=names[chosen],selected_before_fresh_stock=True,adaptive_candidate=1,adopted=False,criterion='VIX forward feasibility and loss <= alpha0-refit+.1; minimize replay native-IV MAE among regime-refit candidates'))
    pieces=[];freshseeds=[]
    for b in range(8):
        seed=int(np.random.SeedSequence([20261028491,500,b]).generate_state(1,dtype=np.uint64)[0]);freshseeds.append(seed);pieces.append(k.stocks(d,64,p,models,alpha,mode,moments,scale,seed));print('TARGETED FRESH',b,flush=True)
    fresh=combine(pieces);np.savez_compressed(OUT/'fresh_stocks.npz',stocks=fresh,seeds=np.array(freshseeds,dtype=np.uint64));fs=evaluate_stocks(fresh,ref,models,alpha,spx,OUT/'fresh');checking={day:dict(np.load(OLD/f'fresh_VIX_draws/d{day}.npz')) for day in DAYS};HH=exposures(k,dated,initial,p,norm,M,x,checking,ALPHA,scale,32,32,8,OUT/'checking_exposures');cv=VIX(HH,checking,x,qv,f);(OUT/'checking_VIX').mkdir();vr=[]
    for c in [0]+list(range(2,11)):
        idx=indices[c];rec=save_vix(cv,idx,models[c,ACTIVE],f'c{c}',OUT/'checking_VIX');rec.update(candidate=c,name=names[c],alpha=float(alpha[c]));vr.append(rec)
    dump(OUT/'checking_VIX/scores.json',vr);dump(OUT/'run.json',dict(status='completed two targeted experiments; audit and plots follow',wall_seconds=time.perf_counter()-begin,regime_selected=int(chosen),SPX_scores=fs,VIX_scores=vr,adaptive_VIX_prices_identical_to_original=True,checking_stock_paths=len(fresh),stock_auxiliary_reference_reused=True,VIX_check_draws_previously_used=True,baseline_adopted=False));print('TARGETED COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':run()
