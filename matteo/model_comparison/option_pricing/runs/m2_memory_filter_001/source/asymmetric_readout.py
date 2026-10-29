"""Asymmetric existing-feedback readout, forward-constrained VIX curve fits."""
import ctypes as ct
import hashlib,json,shutil,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares,minimize
from scipy.special import ndtr
from numpy.polynomial.hermite import hermgauss
from .full_readout_fit import HERE,arr,dump
from .protected_curve_check import context,loss
from .leverage_scan import combine,blocks,stock
from .two_speed_joint_core import physical
from .vix_only_curve_fit import PRIOR,ACTIVE,DAYS,DELTA
from .maturity_curve import window_nodes
from .plot_joint_curve_fit import implied_vol
from esn_eq.states import symmetric_square_root

OUT=HERE/'runs/m2_asymmetric_readout_001'
PARENT=HERE/'runs/m2_protected_curve_001'
LEVERAGE=HERE/'runs/m2_leverage_scan_001'
VBASE=HERE/'runs/m2_two_speed_001'
ALPHA=arr([0,.25,.5,1,1.5,2])

class Kernel:
    def __init__(self,lib):
        self.k=ct.CDLL(str(lib));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.k.asymmetric_moments.argtypes=[ct.c_int]*4+[P]*6+[ct.c_double]*2+[P]*2
        self.k.asymmetric_stock.argtypes=[ct.c_int]*5+[I]+[P]*7+[I]*4+[P]*5+[ct.c_double,ct.c_uint64,P]
    def moments(self,mean,cov,p,center,scale,alpha,order):
        root=symmetric_square_root(cov);assert np.max(abs(root[0,1:]))<1e-10 and np.max(abs(root[1:,0]))<1e-10
        z,w=hermgauss(order);z=np.sqrt(2)*z;w=w/np.sqrt(np.pi);a,b=np.meshgrid(z,z,indexing='ij');off=np.column_stack([a.ravel(),b.ravel()])@root[1:,1:].T;yw=np.outer(w,w).ravel();means=arr(mean).reshape(-1,3);out=np.zeros((len(means),len(alpha),2));self.k.asymmetric_moments(len(means),order,len(off),len(alpha),means,arr(z*root[0,0]),arr(w),arr(off),arr(yw),arr(p),float(center),float(scale),arr(alpha),out);return out
    def stocks(self,d,pairs,p,models,alpha,Z,scale,seed):
        out=np.zeros((2*pairs,int(d['expiries']),len(alpha)+1,2));self.k.asymmetric_stock(len(d['dt']),pairs,d['projection'].shape[1],int(d['expiries']),len(alpha),d['transition'],d['decay'],d['root'],d['cross'],d['mean'],d['projection'],d['time'],d['dt'],d['seg'],d['coarse'],d['reference'],d['expiry'],d['lognorm'],arr(p),arr(models),arr(alpha),arr(Z),float(scale),int(seed),out);return out

class Normalizer:
    def __init__(self,k,dated,initial,p,norm,alpha,scale):
        self.k,self.dated,self.initial,self.p,self.norm,self.alpha,self.scale=k,dated,initial,p,norm,alpha,scale;self.tau=1/max(dated.layout.rates)
    def direct(self,t,order=48):
        m,c=self.dated.gaussian_projection_law(float(t),self.initial);a=self.k.moments(m,c,self.p,m[1],self.scale,self.alpha,order)[0];z=(1-self.w)*a[:,1]/np.exp(self.norm(t)[0])+self.w*a[:,0];z[self.alpha==0]=1.;return z
    def build(self,w):
        self.w=w;knots=np.linspace(0,np.log1p((1060/365)/self.tau),257);maxerr=1.
        for iteration in range(6):
            values=np.array([np.log(self.direct(self.tau*np.expm1(t))) for t in knots]);s=CubicSpline(knots,values,axis=0);mid=(knots[:-1]+knots[1:])/2;actual=np.array([self.direct(self.tau*np.expm1(t)) for t in mid]);error=np.max(abs(np.expm1(s(mid)-np.log(actual))),axis=1);maxerr=float(error.max());bad=error>2e-7
            if not bad.any():break
            knots=np.unique(np.r_[knots,mid[bad]])
        assert maxerr<2e-7,maxerr
        self.spline=s;self.audit=dict(knots=len(knots),midpoint_relative_error=maxerr,GH=48);return self
    def __call__(self,t):return np.exp(self.spline(np.log1p(np.asarray(t)/self.tau)))

def stress_cond(a,T,t,x):
    L=np.log(365);hf=np.exp(L*x[21]);hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf]);b=np.array([5*x[16],5*x[20]]);z=np.array([a['z_normals'],a['fast_normals']]);loading=np.exp(-rates*(t-T))*np.sqrt(-np.expm1(-2*rates*T))*b
    v=np.exp(loading[:,None]*z-.5*loading[:,None]**2);return (1-x[22])*v[0]+x[22]*v[1]

def load_outer(count,seed):
    result={}
    for day in DAYS:
        old=dict(np.load(HERE/f'runs/m2_continuum_response_003/fresh_base_cache/production/vix_d{day}.npz'));fast=dict(np.load(VBASE/f'checking_fast/d{day}.npz'));rng=np.random.default_rng(seed+day);idx=np.concatenate([rng.choice(np.flatnonzero(old['nets']==j),count,replace=False) for j in range(8)]);result[day]={key:value[idx].copy() for key,value in old.items()};result[day]['fast_normals']=fast['fast_normals'][idx];result[day]['weights']*=fast['fast_weight'][idx];result[day]['source_indices']=idx
    return result

def exposures(k,dated,initial,p,norm,Z,x,draws,alpha,scale,order,timeorder,folder):
    folder.mkdir(exist_ok=False);ends=np.array(json.loads((PRIOR/'run.json').read_text())['ends_years']);H={}
    for day,a in draws.items():
        tick=time.perf_counter();T=day/365;times,tw=window_nodes(T,ends,1/max(dated.layout.rates),timeorder);out=np.zeros((len(alpha),len(a['weights']),16))
        for j,(t,wt) in enumerate(zip(times,tw)):
            means,cov=dated.gaussian_projection_law(float(t-T),a['states']);m,_=dated.gaussian_projection_law(float(t),initial);E=k.moments(means,cov,p,m[1],scale,alpha,order);stress=stress_cond(a,T,t,x);response=((1-x[17])*E[:,:,1]/np.exp(norm(t)[0])+x[17]*stress[:,None]*E[:,:,0])/Z(t)[None,:];seg=min(np.searchsorted(ends,t,side='right'),15);out[:,:,seg]+=wt*response.T
        H[day]=out;np.savez_compressed(folder/f'd{day}.npz',H=out,times=times,time_weights=tw,weights=a['weights'],nets=a['nets'],source_indices=a['source_indices']);print('VIX EXPOSURES',folder.name,day,'order',order,'points',len(a['weights']),'seconds',round(time.perf_counter()-tick,2),flush=True)
    return H

class VIX:
    def __init__(self,H,draws,x,q,f):
        self.H,self.draws,self.x,self.q,self.f=H,draws,x,q,f;self.levels,self.floor=physical(x);self.anchor=pd.read_csv(VBASE/'fresh_check/prices_otm_two_speed.csv').model_price_U.to_numpy();self.anchorF=pd.read_csv(VBASE/'fresh_check/forwards_two_speed.csv').model_forward.to_numpy();self.base={}
        for day in DAYS:self.base[day]=100*np.sqrt((self.floor*DELTA+H[day][0]@(self.levels-self.floor))/DELTA)
    def predict(self,idx,y):
        levels=self.levels.copy();levels[ACTIVE]=self.floor+(.2-self.floor)*y;prices=[];jac=[];f=[];fj=[];B=[];BF=[];raw={};offset=0
        for j,day in enumerate(DAYS):
            g=self.q[self.q.expiry_days==day];a=self.draws[day];h=self.H[day][idx];I=self.floor*DELTA+h@(levels-self.floor);V=100*np.sqrt(I/DELTA);dV=50/np.sqrt(DELTA*I)[:,None]*(.2-self.floor)*h[:,ACTIVE];w=a['weights'];K=g.strike.to_numpy();sign=np.where(g.option_type.eq('put'),-1,1);v0=self.base[day];values=w[:,None]*(np.maximum(sign*(V[:,None]-K),0)-np.maximum(sign*(v0[:,None]-K),0));block=values.reshape(8,-1,len(g)).mean(1);price=self.anchor[offset:offset+len(g)]+block.mean(0);J=np.einsum('nq,nj->qj',w[:,None]*(sign*(V[:,None]-K)>0)*sign,dV)/len(w);fb=(w*(V-v0)).reshape(8,-1).mean(1);prices.extend(price);jac.extend(J);f.append(self.anchorF[j]+fb.mean());fj.append((w[:,None]*dV).mean(0));B.append(block);BF.append(fb);raw[day]=dict(V=V,weights=w,nets=a['nets']);offset+=len(g)
        price=np.array(prices);J=np.array(jac);f=np.array(f);fj=np.array(fj);q=self.q;dist=np.where(price<q.price_bid_U,price-q.price_bid_U,np.where(price>q.price_ask_U,price-q.price_ask_U,0))/q.price_scale_U;weights=np.concatenate([np.full(len(g),1/(4*len(g))) for _,g in q.groupby('expiry_days',sort=True)]);r=np.sqrt(weights)*dist;outside=(price<q.price_bid_U)|(price>q.price_ask_U);jr=np.sqrt(weights)[:,None]*outside.to_numpy()[:,None]*J/q.price_scale_U.to_numpy()[:,None]
        return dict(prices=price,jac=jr,residual=np.asarray(r),loss=float(r@r),forwards=f,fjac=fj,blocks=np.concatenate(B,axis=1),forward_blocks=np.column_stack(BF),raw=raw,levels=levels)

def fit_vix(e,idx,start):
    target=e.f.forward_vix_points.to_numpy();memo={}
    def p(y):
        key=tuple(y)
        if key not in memo:memo.clear();memo[key]=e.predict(idx,y)
        return memo[key]
    pre=least_squares(lambda y:p(y)['forwards']-target,start,jac=lambda y:p(y)['fjac'],bounds=(np.zeros(7),np.ones(7)),max_nfev=45,ftol=1e-11,xtol=1e-11,gtol=1e-10)
    opt=minimize(lambda y:p(y)['loss'],pre.x,jac=lambda y:2*p(y)['residual']@p(y)['jac'],method='SLSQP',bounds=[(0,1)]*7,constraints=[dict(type='eq',fun=lambda y:p(y)['forwards']-target,jac=lambda y:p(y)['fjac'])],options=dict(maxiter=80,ftol=1e-9));r=p(opt.x);return opt.x,dict(alpha=float(ALPHA[idx]),loss=r['loss'],max_forward_gap=float(max(abs(r['forwards']-target))),eligible=bool(max(abs(r['forwards']-target))<1e-5),success=bool(opt.success),message=opt.message,levels=r['levels'].tolist(),y=opt.x.tolist())

def save_vix(e,idx,y,label,folder):
    r=e.predict(idx,y);q=e.q.copy();q['model_price_U']=r['prices'];anchorblocks=np.load(VBASE/'fresh_check/blocks_two_speed.npz');q['difference_SE_U']=r['blocks'].std(0,ddof=1)/np.sqrt(8);q['anchor_SE_U']=anchorblocks['native'].std(0,ddof=1)/np.sqrt(8);q['interval_residual']=np.where(q.model_price_U<q.price_bid_U,q.model_price_U-q.price_bid_U,np.where(q.model_price_U>q.price_ask_U,q.model_price_U-q.price_ask_U,0))/q.price_scale_U;q.to_csv(folder/f'VIX_prices_{label}.csv',index=False);f=e.f.copy();f['model_forward']=r['forwards'];f['gap']=f.model_forward-f.forward_vix_points;f['difference_SE']=r['forward_blocks'].std(0,ddof=1)/np.sqrt(8);f['anchor_SE']=anchorblocks['forward'].std(0,ddof=1)/np.sqrt(8);f.to_csv(folder/f'VIX_forwards_{label}.csv',index=False);np.savez_compressed(folder/f'VIX_blocks_{label}.npz',native=r['blocks'],forward=r['forward_blocks']);return dict(label=label,option_loss=r['loss'],forwards=r['forwards'].tolist(),max_forward_gap=float(abs(f.gap).max()),negative=int((q.model_price_U<0).sum()),inside=int(q.interval_residual.eq(0).sum()))

def evaluate_stocks(stocks,ref,models,alpha,q,folder):
    folder.mkdir(exist_ok=False);days=sorted(q.expiry_days.unique());ends=np.array(json.loads((PRIOR/'run.json').read_text())['ends_years']);levels,f=physical(models[0]);cs=np.sqrt(1-models[0,22]);cf=np.sqrt(models[0,22]);C=[];deltas={}
    for j,day in enumerate(days):
        g=q[q.expiry_days==day];K=g.K_over_F.to_numpy();a=ref[:,j];I=a[:,0];v=.36*I;m=-.5*I-.8*(cs*a[:,2]+cf*a[:,3]);sd=np.sqrt(v)[:,None];d1=(m[:,None]+v[:,None]-np.log(K)[None,:])/sd;B=np.exp(m+v/2)[:,None]*ndtr(d1)-K*ndtr(d1-sd);C.append(blocks(B)*g.forward.to_numpy());T=day/365;detI=np.maximum(0,np.minimum(ends,T)-np.minimum(np.r_[0,ends[:-1]],T))@levels;deltas[day]=ndtr(-np.log(K)/np.sqrt(detI)+.5*np.sqrt(detI))
    C=np.concatenate(C,axis=1);scores=[];saved={}
    for a in range(len(alpha)):
        for grid,steps in enumerate([64,128]):
            pieces=[];eb=[]
            for j,day in enumerate(days):
                g=q[q.expiry_days==day].copy();K=g.K_over_F.to_numpy();F=g.forward.to_numpy();S=stocks[:,j,a,grid];R=stocks[:,j,-1,grid];y=np.maximum(S[:,None]-K,0)-np.maximum(R[:,None]-K,0)-deltas[day]*(S-R)[:,None];E=blocks(y)*F;eb.append(E);g['stock_mean']=float(S.mean());g['stock_mean_SE']=float(blocks(S[:,None]).std(ddof=1)/np.sqrt(8));pieces.append(g)
            frame=pd.concat(pieces,ignore_index=True);E=np.concatenate(eb,axis=1);frame['model_price_U']=E.mean(0)+C.mean(0)-frame.option_type.eq('put').to_numpy()*(frame.forward-frame.strike).to_numpy();frame['mc_se_U']=np.sqrt(E.var(0,ddof=1)/8+C.var(0,ddof=1)/8);frame['model_IV']=[implied_vol(r.model_price_U,r.forward,r.strike,r.T_years,r.option_type) for r in frame.itertuples()];frame['market_IV']=[implied_vol(r.price_mid_U,r.forward,r.strike,r.T_years,r.option_type) for r in frame.itertuples()];frame['alpha']=alpha[a];frame['candidate']=a;frame['steps']=steps;primary=frame.priority_selected.to_numpy();p=frame[primary].reset_index(drop=True);ivmae=float(np.mean([np.nanmean(abs(g.model_IV-g.market_IV)) for _,g in p.groupby('expiry_days')]));inside=(p.model_price_U>=p.price_bid_U)&(p.model_price_U<=p.price_ask_U);row=dict(candidate=a,alpha=float(alpha[a]),steps=steps,loss=loss(p,p.model_price_U.to_numpy()),IV_MAE_pp=ivmae,inside=int(inside.sum()),invalid_IV=int(p.model_IV.isna().sum()),negative=int((p.model_price_U<0).sum()),median_SE_scales=float(np.median(p.mc_se_U/p.price_scale_U)));scores.append(row);frame.to_csv(folder/f'SPX_prices_c{a}_{steps}.csv',index=False);saved[(a,steps)]=(E,frame)
            if steps==128:print('SPX',folder.name,row,flush=True)
    # Paired candidate-minus-original uncertainty, including the same auxiliary.
    for a in range(len(alpha)):
        E,b=saved[(a,128)];E0,b0=saved[(0,128)];paired=b[['expiry_days','strike','option_type','priority_selected','price_scale_U']].copy();paired['candidate_minus_original_U']=b.model_price_U-b0.model_price_U;paired['paired_SE_U']=(E-E0).std(0,ddof=1)/np.sqrt(8);paired.to_csv(folder/f'paired_c{a}.csv',index=False)
    dump(folder/'SPX_scores.json',scores);pd.DataFrame(scores).to_csv(folder/'SPX_scores.csv',index=False);np.savez_compressed(folder/'SPX_blocks.npz',cheap=C,**{f'c{a}_{steps}':E for (a,steps),(E,_) in saved.items()});return scores

def run():
    begin=time.perf_counter();OUT.mkdir(exist_ok=False);(OUT/'source').mkdir();sources=[Path(__file__),HERE/'asymmetric_readout.cpp',HERE/'leverage_scan.cpp',HERE/'two_speed_precision_kernel.cpp',HERE/'ASYMMETRIC_READOUT_PLAN.md',PARENT/'selection.json',PRIOR/'run.json'];hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    for p in sources:shutil.copy2(p,OUT/'source'/p.name)
    lib=OUT/'asymmetric.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'asymmetric_readout.cpp'),'-o',str(lib)],check=True);k=Kernel(lib);prior,dated,norm,oldnormaudit=context();initial=arr(prior['initial_state']);p=arr(prior['selected']);x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);d=dict(np.load(PARENT/'fresh_check/grid_and_transitions.npz'));_,cov=dated.gaussian_projection_law(1000,np.zeros(len(initial)));scale=float(np.sqrt(cov[1,1]));covsign=float(d['projection'][1]@d['cross'][0]);assert covsign>0
    Z=Normalizer(k,dated,initial,p,norm,ALPHA,scale).build(x[17]);dump(OUT/'design.json',dict(alpha=ALPHA.tolist(),theta=.8,X='(Y-E_Q[Y])/stationary_sd_Y; existing feedback projection',scale=scale,feedback_equity_covariance_loading=covsign,normalizer=Z.audit,old_normalizer=oldnormaudit,source_hashes=hashes,baseline=x.tolist(),frozen_gaussian_dynamics=True,forward_constraints='four training-estimator equalities, not population certificate',SPX_eligible=1751,SPX_primary=1294,VIX_quotes=217))
    probes=np.geomspace(1e-7,1030/365,13);z48=np.array([Z.direct(t,48) for t in probes]);z64=np.array([Z.direct(t,64) for t in probes]);dump(OUT/'normalizer_order_audit.json',dict(times=probes.tolist(),maximum_relative_GH48_vs64_by_alpha=np.max(abs(z48/z64-1),axis=0).tolist()));np.savez_compressed(OUT/'normalizers.npz',times=d['time'],Z=Z(d['time']),alpha=ALPHA)
    draws=load_outer(128,20261028203);H=exposures(k,dated,initial,p,norm,Z,x,draws,ALPHA,scale,24,4,OUT/'training_exposures');q=pd.read_csv(VBASE/'fresh_check/prices_otm_two_speed.csv').drop(columns=['model_price_U','mc_se_U','interval_residual']);f=pd.read_csv(VBASE/'fresh_check/forwards_two_speed.csv').drop(columns=['model_forward','forward_error','mc_se']);q.to_csv(OUT/'VIX_membership.csv',index=False);f.to_csv(OUT/'forward_targets.csv',index=False);v=VIX(H,draws,x,q,f);models=[x.copy()];alphas=[0.];ledger=[];vrecs=[];start=x[ACTIVE].copy();training=OUT/'training';training.mkdir()
    vrecs.append(save_vix(v,0,start,'c0',training))
    for idx,alpha in enumerate(ALPHA):
        y,rec=fit_vix(v,idx,start);ledger.append(rec);xx=x.copy();xx[ACTIVE]=y;models.append(xx);alphas.append(alpha);vrecs.append(save_vix(v,idx,y,f'c{idx+1}',training));dump(OUT/'VIX_fit_ledger.json',ledger);print('VIX FIT',rec,flush=True)
    models=arr(models);alphas=arr(alphas);Zstock=arr(np.column_stack([np.ones(len(d['time'])),Z(d['time'])]));dump(OUT/'candidates.json',dict(models=models.tolist(),alpha=alphas.tolist(),VIX_training=vrecs));spx=pd.read_csv(PARENT/'spx_membership.csv');spx.to_csv(OUT/'SPX_membership.csv',index=False);replay=[];seeds=np.load(PARENT/'fresh_check/stocks.npz')['seeds']
    for b,seed in enumerate(seeds):
        tick=time.perf_counter();replay.append(k.stocks(d,32,p,models,alphas,Zstock,scale,int(seed)));print('STOCK REPLAY',b,'seconds',round(time.perf_counter()-tick,2),flush=True)
    replay=combine(replay);np.savez_compressed(OUT/'replay_stocks.npz',stocks=replay,seeds=seeds);stats=np.load(LEVERAGE/'replay_stats.npz')['stats'];expect=stock(stats[:,:,:2],.8,np.sqrt(1-x[22]),np.sqrt(x[22]));error=float(np.max(abs(replay[:,:,0]-expect)));assert error<2e-10,error;dump(OUT/'zero_multiplier_stock_replay.json',dict(max_absolute_stock_error=error));ref=np.load(LEVERAGE/'reference.npz')['stats'];scores=evaluate_stocks(replay,ref,models,alphas,spx,OUT/'replay');eligible=[s for s in scores if s['steps']==128 and s['candidate']>0 and ledger[s['candidate']-1]['eligible']]
    # Predeclared selection: retain VIX option loss within baseline+.1, then
    # minimize SPX IV error (not price loss alone). If none, record tradeoff.
    feasible=[s for s in eligible if ledger[s['candidate']-1]['loss']<=ledger[0]['loss']+.1];pool=feasible or eligible;best=min(pool,key=lambda s:s['IV_MAE_pp']);chosen=int(best['candidate']);dump(OUT/'selection.json',dict(candidate=chosen,alpha=float(alphas[chosen]),selected_before_fresh_stock=True,VIX_preservation_filter_satisfied=bool(feasible),selection_metric='minimum SPX IV MAE among forward-feasible strengths with VIX loss <= refitted alpha0 loss + .1; if none retain best tradeoff explicitly',SPX=best,adopted=False));print('SELECTION',chosen,alphas[chosen],flush=True)
    # Fresh stock streams, jointly simulated once for all frozen candidates.
    fresh=[];freshseeds=[]
    for b in range(8):
        seed=int(np.random.SeedSequence([20261028203,400,b]).generate_state(1,dtype=np.uint64)[0]);freshseeds.append(seed);fresh.append(k.stocks(d,64,p,models,alphas,Zstock,scale,seed));print('FRESH STOCK',b,flush=True)
    fresh=combine(fresh);np.savez_compressed(OUT/'fresh_stocks.npz',stocks=fresh,seeds=np.array(freshseeds,dtype=np.uint64));fscores=evaluate_stocks(fresh,ref,models,alphas,spx,OUT/'fresh');checking=load_outer(512,20261028291);selected_indices=sorted(set([0,chosen-1]));aa=ALPHA[selected_indices];checkH=exposures(k,dated,initial,p,norm,lambda t:Z(t)[selected_indices],x,checking,aa,scale,32,8,OUT/'checking_exposures');cv=VIX(checkH,checking,x,q,f);cvrec=[]
    for idx in [0,chosen]:
        local=0 if idx==0 else selected_indices.index(chosen-1);y=models[idx,ACTIVE];cvrec.append(save_vix(cv,local,y,f'c{idx}',OUT/'fresh'))
    dump(OUT/'run.json',dict(status='completed exploratory asymmetric-readout trial',wall_seconds=time.perf_counter()-begin,selection=chosen,alpha=float(alphas[chosen]),SPX_fresh=fscores,VIX_check=cvrec,checking_VIX_states='randomized subsets of existing large state/stress stream, not wholly new',checking_stock_paths=len(fresh),reference_paths=len(ref),reference_reused=True,baseline_adopted=False))
    print('ASYMMETRIC COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':run()
