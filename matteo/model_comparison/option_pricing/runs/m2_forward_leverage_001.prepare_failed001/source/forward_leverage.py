"""Fixed-variance stock coupling to conditional next30day integrated variance."""
import ctypes as ct
import hashlib,json,shutil,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
from .full_readout_fit import HERE,arr,dump
from .protected_curve_check import context
from .asymmetric_readout import PARENT,PRIOR,LEVERAGE,evaluate_stocks
from .leverage_scan import combine
from esn_eq.states import symmetric_square_root

OUT=HERE/'runs/m2_forward_leverage_001'
TARGET=HERE/'runs/m2_targeted_skew_001'
L=30/365

class Kernel:
    def __init__(self,path):
        self.k=ct.CDLL(str(path));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.k.forward_table.argtypes=[ct.c_int]*2+[P]*7
        self.k.forward_probe.argtypes=[ct.c_int]*5+[P]*16+[ct.c_double,P]
        self.k.forward_stock.argtypes=[ct.c_int]*8+[I]+[P]*10+[I]*5+[P]*13+[ct.c_uint64,P]
    def tables(self,sd,p):
        z,w=hermgauss(48);meta=np.zeros((len(sd),2));table=np.zeros((len(sd),8193,2));self.k.forward_table(len(sd),48,arr(sd),arr(np.sqrt(2)*z),arr(w/np.sqrt(np.pi)),arr(p),meta,table);return meta,table
    def probe(self,states,stress,signal,b,p,x,t):
        states=arr(states);out=np.zeros((len(states),b['diffusion_lag'].shape[2]+2));self.k.forward_probe(len(states),len(b['lags']),len(b['weight']),states.shape[1],out.shape[1]-2,states,arr(stress),arr(signal['futuremean']),arr(signal['inverse_norm']),arr(signal['hat']),b['projection_lag'],b['diffusion_lag'],b['offset'],b['weight'],b['meta'],b['table'],b['lags'],arr(p),arr(x),float(t),out);return out
    def stocks(self,d,b,signals,p,x,pairs,stride,seed):
        out=np.zeros((2*pairs,int(d['expiries']),5,2));self.k.forward_stock(len(d['dt']),pairs,stride,d['projection'].shape[1],d['diffusion'].shape[1],int(d['expiries']),len(b['lags']),len(b['weight']),d['transition'],d['decay'],d['root'],d['cross'],d['bmean'],d['broot'],d['diffusion'],d['mean'],d['projection'],d['time'],d['dt'],d['seg'],d['coarse'],d['reference'],d['expiry'],signals['update'],d['lognorm'],arr(p),arr(x),signals['futuremean'],signals['inverse_norm'],signals['hat'],b['projection_lag'],b['diffusion_lag'],b['offset'],b['weight'],b['meta'],b['table'],b['lags'],int(seed),out);return out

def hat_weights(times,lags,ends,coeff):
    result=np.zeros((len(times),len(lags)))
    for i,t in enumerate(times):
        for j,(a,b) in enumerate(zip(lags[:-1],lags[1:])):
            breaks=np.r_[a,ends[(ends>t+a)&(ends<t+b)]-t,b]
            for u,v in zip(breaks[:-1],breaks[1:]):
                c=coeff[min(np.searchsorted(ends,t+(u+v)/2,side='right'),len(coeff)-1)];result[i,j]+=c*(b*(v-u)-.5*(v*v-u*u))/(b-a);result[i,j+1]+=c*(.5*(v*v-u*u)-a*(v-u))/(b-a)
    assert result.min()>-1e-12;return arr(result)

def basis(k,dated,d,p,lag_count=22,order=12):
    lags=arr(np.r_[0,np.geomspace(1e-5,L,lag_count)]);nd=d['projection'].shape[1];rates=dated.layout.rates[:nd];keep=np.flatnonzero(np.any(dated.projection!=0,axis=0));assert np.array_equal(keep,np.arange(nd));load=dated.layout.driver_loadings[:nd];used=np.flatnonzero(np.any(load!=0,axis=0));B=np.sqrt(2*rates)[:,None]*load[:,used];P=d['projection'];decay=np.exp(-lags[:,None]*rates);proj=arr(P[None,:,:]*decay[:,None,:]);diff=arr(proj@B);z,w=hermgauss(order);z=np.sqrt(2)*z;w=w/np.sqrt(np.pi);u,v=np.meshgrid(z,z,indexing='ij');points=np.column_stack([u.ravel(),v.ravel()]);offset=[];sd=[]
    for lag in lags:
        _,cov=dated.gaussian_projection_law(float(lag),np.zeros(len(dated.layout.rates)));root=symmetric_square_root(cov);assert abs(root[0,1:]).max()<1e-10;offset.append(points@root[1:,1:].T);sd.append(np.sqrt(cov[0,0]))
    meta,table=k.tables(arr(sd),p);return dict(lags=lags,projection_lag=proj,diffusion_lag=diff,offset=arr(offset),weight=arr(np.outer(w,w).ravel()),meta=meta,table=table,conditional_C_sd=arr(sd))

def signals(d,dated,prior,norm,x,b,times=None):
    if times is None:
        ids=np.arange(0,len(d['dt']),32);times=d['time'][ids];update=np.full(len(d['dt']),-1,np.int32);update[ids]=np.arange(len(ids),dtype='int32')
    else:update=np.array([],np.int32);times=arr(times)
    u=times[:,None]+b['lags'];rates=dated.layout.rates[:24];initial=arr(prior['initial_state'])[:24];drift=dated.q_drift[:24];states=np.exp(-u[:,:,None]*rates)*initial+(-np.expm1(-u[:,:,None]*rates))*drift/rates;mean=arr(states@d['projection'].T);inverse=arr(np.exp(-norm(u.ravel())[:,0]).reshape(u.shape));f=(.05+.07*x[19])**2;coeff=(.2-f)*x[:16];ends=arr(prior['ends_years']);hat=hat_weights(times,b['lags'],ends,coeff)
    return dict(update=update,times=arr(times),futuremean=mean,inverse_norm=inverse,hat=hat)

def relabel(folder):
    # Input evaluator's64/128 grid labels become the actual128/256 labels.
    for path in list(folder.glob('SPX_prices_*.csv')):
        q=pd.read_csv(path);old=int(q.steps.iloc[0]);q['steps']=old*2;q.to_csv(path.with_name(path.name.replace(f'_{old}.csv',f'_{old*2}.csv')+'.new'),index=False);path.unlink()
    for path in folder.glob('*.new'):path.rename(path.with_name(path.name[:-4]))
    scores=json.loads((folder/'SPX_scores.json').read_text());scores=[{**s,'steps':2*s['steps']} for s in scores];dump(folder/'SPX_scores.json',scores);pd.DataFrame(scores).to_csv(folder/'SPX_scores.csv',index=False);B=dict(np.load(folder/'SPX_blocks.npz'));np.savez_compressed(folder/'SPX_blocks.npz',**{key if key=='cheap' else key.rsplit('_',1)[0]+'_'+str(2*int(key.rsplit('_',1)[1])):v for key,v in B.items()});return scores

def setup(out):
    prior,dated,norm,na=context();p=arr(prior['selected']);x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);d=dict(np.load(TARGET/'grid_refinement/grid.npz'));k=Kernel(out/'forward.dylib');return prior,dated,norm,na,p,x,d,k

def run():
    begin=time.perf_counter();OUT.mkdir(exist_ok=False);(OUT/'source').mkdir();sources=[Path(__file__),HERE/'forward_leverage.cpp',HERE/'FORWARD_LEVERAGE_PLAN.md',HERE/'asymmetric_readout.cpp',HERE/'leverage_scan.cpp',HERE/'two_speed_precision_kernel.cpp',PARENT/'selection.json',PRIOR/'run.json',TARGET/'grid_refinement/grid.npz'];hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    for p in sources:
        if p.suffix!='.npz':shutil.copy2(p,OUT/'source'/p.name)
    subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'forward_leverage.cpp'),'-o',str(OUT/'forward.dylib')],check=True);prior,dated,norm,na,p,x,d,k=setup(OUT);b=basis(k,dated,d,p);ss=signals(d,dated,prior,norm,x,b);np.savez_compressed(OUT/'basis.npz',**b);np.savez_compressed(OUT/'signals.npz',**ss);dump(OUT/'design.json',dict(status='prepared; integration audit must pass before stock run',source_hashes=hashes,baseline=x.tolist(),horizon_days=30,GH_YF=12,GH_C=48,lag_count=len(b['lags']),normalizer=na,all_factor_dynamics_and_VIX_law_unchanged=True,stock_grids=[128,256],direction_updates_per_trading_day=[4,8],pricing_counts=dict(primary=1294,eligible=1751,VIX=217,SPX_maturities=7),planned_fresh_paths=1024,old_auxiliary_reference_paths=65536,stock_selection='none: all four couplings fixed before pricing'))
    print('FORWARD PREPARED',time.perf_counter()-begin,'seconds; signal updates',len(ss['times']),flush=True)

if __name__=='__main__':run()
