"""Eleven-parameter M2 calibration on an immutable exogenous reservoir."""
import argparse
import ctypes as ct
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares
from scipy.special import expit
from .m2_spx import ROOT,M2,frozen_stage3
from .all_maturity_core import populations,InSampleLoss
from .spx_volatility_sampling import proposal
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from esn_eq.pricing import DatedM2Pricer,ForwardVarianceCurve,dated_grid,_normal_rule
from esn_eq.states import symmetric_square_root
from esn_eq.market_targets import MarketTargets
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.vix_pricing import ConditionalVIXPricer
from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths

HERE=Path(__file__).resolve().parent
NAMES=['cluster_loading','feedback_curvature','feedback_shift','spike_curvature','spike_shift','cap_level','cap_sharpness','xi_1','xi_2','xi_3','variance_floor']
LOW=np.array([.05,.005,.05,.002,.01,2.,2.,.005,.005,.005,1e-6])
HIGH=np.array([.8,.4,5.,.6,2.,8.,15.,.20,.20,.20,.0025])
SCALE=np.array([.25,.1,2.,.25,.5,6.,8.,.06,.06,.06,.001])

def dump(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def arr(x):return np.ascontiguousarray(x,dtype=np.float64)

class Kernel:
 def __init__(self,path):
  self.lib=ct.CDLL(str(path));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
  self.lib.multiplier.argtypes=[ct.c_int,P,P,P]
  self.lib.expectation.argtypes=[ct.c_int,ct.c_int,P,P,P,P,P]
  self.lib.stock.argtypes=[ct.c_int,ct.c_int,ct.c_int,P,P,P,I,I,P,P,P]
 def multiplier(self,x,p):
  x=arr(x).reshape(-1,3);o=np.empty((len(x),8));self.lib.multiplier(len(x),x,arr(p),o);return o
 def expectation(self,mean,offset,w,p):
  mean=arr(mean).reshape(-1,3);o=np.empty((len(mean),8));self.lib.expectation(len(mean),len(w),mean,offset,w,arr(p),o);return o
 def stock(self,c,norm,p):
  x=np.load(c/'features.npy',mmap_mode='r');dw=np.load(c/'equity.npy',mmap_mode='r');grid=np.load(c/'grid.npy');events=np.load(c/'events.npy');seg=np.minimum(np.searchsorted([.25,.75],grid[:-1],side='right'),2).astype('int32')
  out=np.empty((len(events),x.shape[1],12));self.lib.stock(len(dw),x.shape[1],len(events),x,dw,arr(np.diff(grid)),seg,events,arr(norm),arr(p),out);return out

class Engine:
 def __init__(self,out,role,kernel,dated,initial,spx,vix,forwards,ends):
  self.out,self.role,self.k,self.dated,self.initial,self.spx,self.vix,self.forwards,self.ends=out,role,kernel,dated,initial,spx,vix,forwards,ends
  self.quotes=pd.concat([spx,vix],ignore_index=True);self.loss=InSampleLoss(self.quotes,forwards)
  self.folder=out/role
  self.caches=sorted(self.folder.glob('stock_batch*'))
  self.grid=np.load(self.caches[0]/'grid.npy');self.times=self.grid[:-1]
  tau=1/max(dated.layout.rates);self.tau=tau
  self.knots=np.linspace(0,np.log1p((self.grid[-1]+30/365)/tau),257)
  self.normlaws=[self.law(float(tau*np.expm1(s)),initial,48) for s in self.knots]
  self.vixlaws={}
  self.seed_records=json.loads((self.folder/'seeds.json').read_text())
 def law(self,t,state,order):
  m,c=self.dated.gaussian_projection_law(t,state);points,w=_normal_rule(order)
  return arr(m),arr(points@symmetric_square_root(c).T),arr(w)
 def normalizer(self,p,probes=False):
  values=np.array([self.k.expectation(m,o,w,p)[0] for m,o,w in self.normlaws]);lg=np.log(values[:,0]);d=values[:,1:]/values[:,:1]
  spline=CubicSpline(self.knots,np.column_stack([lg,d]),axis=0,extrapolate=False)
  error=None
  if probes:
   probes=np.linspace(self.knots[0],self.knots[-1],130)[1:-1]
   direct=np.array([self.k.expectation(*self.law(float(self.tau*np.expm1(s)),self.initial,48),p)[0,0] for s in probes])
   error=float(np.max(np.abs(np.expm1(spline(probes)[:,0]-np.log(direct)))))
  return lambda t:spline(np.log1p(np.asarray(t)/self.tau)),error
 def vix_laws(self,days,order):
  key=(days,order)
  if key not in self.vixlaws:
   a=np.load(self.folder/f'vix_d{days}.npz');states=a['states'];T=days/365
   # Window nodes/weights are candidate-independent; its normalizer is unused here.
   times=a['times'];weights=a['time_weights'];laws=[self.law(float(t-T),states,order) for t in times]
   self.vixlaws[key]=(times,weights,laws,a['weights'],a['nets'])
  return self.vixlaws[key]
 def predict(self,p,order=12,detail=False):
  norm,error=self.normalizer(p,probes=detail);ng=norm(self.times)
  terminal=np.concatenate([self.k.stock(c,ng,p) for c in self.caches],axis=1)
  weights=np.concatenate([np.load(c/'weights.npy') for c in self.caches]);prices=[];errors=[];jac=[];raw={};fwd=[];fse=[];fj=[]
  for j,(days,q) in enumerate(self.spx.groupby('expiry_days',sort=True)):
   R=terminal[j,:,0];dR=terminal[j,:,1:];sign=np.where(q.option_type.to_numpy()=='call',1.,-1.);distance=sign*(R[:,None]-q.K_over_F.to_numpy());pay=np.maximum(distance,0)*q.forward.to_numpy()
   y=weights[:,None]*pay;prices.extend(y.mean(0));errors.extend(y.std(0,ddof=1)/np.sqrt(len(y)))
   jac.extend(np.einsum('nk,nj->kj',weights[:,None]*(distance>0)*sign*q.forward.to_numpy(),dR)/len(y))
   raw[f'spx_{days}']=y
  for days,q in self.vix.groupby('expiry_days',sort=True):
   times,tw,laws,w,nets=self.vix_laws(int(days),order);integral=np.zeros(len(w));der=np.zeros((len(w),11));delta=30/365
   for t,wt,law in zip(times,tw,laws):
    E=self.k.expectation(*law,p);nn=norm(t);M=np.exp(nn[0]);ratio=E[:,0]/M;segment=min(np.searchsorted(self.ends,t,side='right'),2);xi=p[7+segment]
    integral+=wt*(p[10]+(xi-p[10])*ratio)
    der[:,:7]+=wt*(xi-p[10])*(E[:,1:]/M-ratio[:,None]*nn[1:])
    der[:,7+segment]+=wt*ratio;der[:,10]+=wt*(1-ratio)
   V=100*np.sqrt(integral/delta);dV=50/np.sqrt(delta*integral)[:,None]*der
   sign=np.where(q.option_type.to_numpy()=='call',1.,-1.);distance=sign*(V[:,None]-q.strike.to_numpy());pay=np.maximum(distance,0);y=w[:,None]*pay
   b=np.array([y[nets==i].mean(0) for i in range(8)]);prices.extend(b.mean(0));errors.extend(b.std(0,ddof=1)/np.sqrt(8))
   jac.extend(np.einsum('nk,nj->kj',w[:,None]*(distance>0)*sign,dV)/len(w))
   bf=np.array([(w*V)[nets==i].mean() for i in range(8)]);fwd.append(bf.mean());fse.append(bf.std(ddof=1)/np.sqrt(8));fj.append(np.mean(w[:,None]*dV,axis=0));raw[f'vix_{days}']=V;raw[f'vix_weight_{days}']=w;raw[f'vix_blocks_{days}']=b
  prices,errors,jac,fwd,fse,fj=map(np.asarray,[prices,errors,jac,fwd,fse,fj]);r=self.loss.evaluate(prices,fwd,p[7:10]);outside=(prices<self.quotes.price_bid_U.to_numpy())|(prices>self.quotes.price_ask_U.to_numpy())
  pj=self.loss.sqrtweights[:,None]*outside[:,None]*jac/self.quotes.price_scale_U.to_numpy()[:,None];pen=np.zeros((2,11));pen[:,7:10]=np.sqrt(.1)*np.array([[-1,1,0],[0,-1,1]])
  r['jacobian']=np.concatenate([pj,fj/(.5*2),pen]);r['total']+=.1/3;r['components']['state_penalty']=.1/3
  # Constant state penalty included as a residual to keep sum squares identical.
  r['residuals']=np.r_[r['residuals'],np.sqrt(.1/3)];r['jacobian']=np.vstack([r['jacobian'],np.zeros(11)])
  r.update(prices=prices,errors=errors,forwards=fwd,forward_errors=fse,raw=raw,normalizer_probe_error=error)
  return r

def prepare(out,role,p,initial,parameters,levels,spx,vix,paths,seed):
 folder=out/role;folder.mkdir();ends=[.25,.75,2129/365];maturities=np.sort(spx.T_years.unique());events=np.r_[maturities,7/365,21/365,.25,.75]
 grid=dated_grid(events,1/(252*512));basis=volatility_path_basis(p,grid,7/365);records=[]
 for batch,start in enumerate(range(0,paths,128)):
  n=min(128,paths-start);s=int(np.random.SeedSequence([seed,0 if role=='search' else 1,batch]).generate_state(1)[0]);tick=time.perf_counter()
  r=prepare_volatility_paths(p,grid,initial,n,s,basis,proposal(1));cache=r.cache;c=folder/f'stock_batch{batch}';c.mkdir()
  for name,a in [('features',cache.projections),('equity',cache.equity_increments),('grid',grid),('weights',r.weights),('events',np.searchsorted(grid,maturities).astype('int32'))]:np.save(c/(name+'.npy'),a)
  records.append(dict(kind='stock',batch=batch,seed=s,paths=n));print(role,'stock batch',batch,'seconds',time.perf_counter()-tick,flush=True);del r,cache
 del basis
 declaration=json.loads((HERE/'runs/m2_vix_targeted_check_001/proposals.json').read_text());curve=ForwardVarianceCurve(tuple(ends),tuple(levels));engine=ConditionalVIXPricer(p,parameters,curve,initial,normalizer_order=48)
 for days in [27,55,118,153]:
  T=days/365;sampler=GaussianOuterSampler(p,T,initial);d=declaration[str(days)];mixture=DefensiveGaussianMixture(d['shifts'],d['probabilities']);states=[];weights=[];nets=[];power=4 if role=='search' else 5
  for i in range(8):
   s=int(np.random.SeedSequence([seed,2 if role=='search' else 3,days,i]).generate_state(1)[0]);z,w,_=mixture.draw(power,s);states.append(sampler.mean+z@sampler.root.T);weights.append(w);nets.extend([i]*len(w));records.append(dict(kind='vix',days=days,net=i,seed=s))
  window=engine.gaussian_window(T,time_order=8);np.savez_compressed(folder/f'vix_d{days}.npz',states=np.concatenate(states),weights=np.concatenate(weights),nets=nets,times=window.times,time_weights=window.weights)
 dump(folder/'seeds.json',records)

def checks(k,p):
 from esn_eq.config import ReadoutParameters
 from esn_eq.features import FeaturePaths
 from esn_eq.readout import TemperedExponentialQuadratic
 rng=np.random.default_rng(24);x=rng.normal(size=(17,3));x[0]=[20,20,-20];x[1]=[-20,0,0]
 y=k.multiplier(x,p);pars=ReadoutParameters(**dict(zip(NAMES[:7],p[:7])))
 f=FeaturePaths(x[:,1],x[:,2],x[:,0],np.zeros(len(x)),np.zeros(len(x)));expected=TemperedExponentialQuadratic().multiplier(f,pars);np.testing.assert_allclose(y[:,0],expected,rtol=2e-14)
 errs=[]
 for j in range(7):
  h=1e-5*max(abs(p[j]),.1);a=p.copy();b=p.copy();a[j]+=h;b[j]-=h;fd=(k.multiplier(x,a)[:,0]-k.multiplier(x,b)[:,0])/(2*h)
  err=float(np.max(np.abs(fd-y[:,j+1])/(1+np.abs(fd))));assert err<2e-6,(j,err);errs.append(err)
 return dict(multiplier_matches_production=True,maximum_scaled_derivative_error=max(errs))

def run(args):
 out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();a,pars,premium,hashes=frozen_stage3();p=DatedM2Pricer(a,premium);initial,_=p.initial_state([1,0,0]);targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');spx,vix,membership,atm=populations(targets.quotes);forwards=targets.vix_forwards.sort_values('expiry').reset_index(drop=True)
 assert len(spx)==1294 and len(vix)==217
 for q in [spx,vix]:assert np.all(np.where(q.option_type=='put',q.K_over_F<=1,q.K_over_F>=1))
 old=json.loads((HERE/'runs/m2_all_maturity_fit_001/run.json').read_text());baseline=np.array([getattr(pars,n) for n in NAMES[:7]]+old['selected_levels']+[1e-4]);start=baseline.copy();start[4]=.15
 sources=[Path(__file__),HERE/'full_readout_kernel.cpp',HERE/'FULL_READOUT_FIT_PLAN.md',HERE/'m2_spx.py',HERE/'all_maturity_core.py',HERE/'joint_curve_core.py',HERE/'spx_volatility_sampling.py',HERE/'runs/m2_vix_targeted_check_001/proposals.json',HERE/'runs/m2_all_maturity_fit_001/run.json',HERE/'active_scope.json',*sorted((M2/'src/esn_eq').glob('*.py'))];captured=capture_sources(sources);save_sources(out,captured,ROOT)
 lib=out/'kernel.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'full_readout_kernel.cpp'),'-o',str(lib)],check=True);k=Kernel(lib)
 design=dict(status='preparing immutable reservoir',scope='all market observations in sample',parameter_names=NAMES,bounds_lower=LOW.tolist(),bounds_upper=HIGH.tolist(),baseline=baseline.tolist(),start=start.tolist(),baseline_levels=baseline[7:10].tolist(),ends_years=[.25,.75,2129/365],architecture=asdict(a),premium=asdict(premium),initial_state=initial.tolist(),free_parameters=11,fixed='OU grids/counts/correlations, projection powers, Q prices, initial state, disabled optional channels',implementation_sha256=source_hashes(captured,ROOT),target_sha256=targets.artifact_sha256,frozen_candidate_sha256=hashes,seed=args.seed,budget=dict(stock_steps_per_trading_day=512,search_stock_paths=128,check_stock_paths=256,vix_nets=8,search_vix_points_per_net=16,check_vix_points_per_net=32,search_conditional_GH=12,check_conditional_GH=24,normalizer_GH=48,time_order=8,max_nfev_per_start=35,starts=2,optimization_seconds=1800),estimator='positive raw frozen-proposal importance payoffs; no strike-dependent controls',market_out_of_sample=False)
 dump(out/'design.json',design);dump(out/'kernel_checks.json',checks(k,start))
 for name,t in [('spx_contracts',spx),('vix_contracts',vix),('spx_membership',membership),('atm_proxies',atm),('forward_targets',forwards)]:t.to_csv(out/(name+'.csv'),index=False)
 prepare(out,'search',p,initial,pars,baseline[7:10],spx,vix,128,args.seed)
 engine=Engine(out,'search',k,p,initial,spx,vix,forwards,design['ends_years']);records=[];evals=[]
 def record(engine,name,x,role,order):
  r=engine.predict(x,order,True);q,f=engine.loss.tables(r['prices'],r['errors'],r['forwards'],r['forward_errors']);q.to_csv(out/f'prices_{role}_{name}.csv',index=False);f.to_csv(out/f'forwards_{role}_{name}.csv',index=False);np.savez_compressed(out/f'raw_{role}_{name}.npz',**r['raw']);rec=dict(role=role,candidate=name,conditional_GH=order,total=r['total'],components=r['components'],inside_bidask=int(q.interval_residual.eq(0).sum()),negative_estimates=int(q.negative_estimate.sum()),nominal_se_passes=int(q.nominal_se_screen.sum()),normalizer_probe_error=r['normalizer_probe_error']);records.append(rec);print('record',json.dumps(rec),flush=True);return r
 record(engine,'baseline',baseline,'search',12)
 # End-to-end Jacobian at the positive start on three informative directions.
 r0=engine.predict(start);rng=np.random.default_rng(17);errs=[]
 for _ in range(2):
  direction=rng.normal(size=11)*SCALE;h=2e-5
  fd=(engine.predict(start+h*direction)['residuals']-engine.predict(start-h*direction)['residuals'])/(2*h);an=r0['jacobian']@direction;err=float(np.linalg.norm(fd-an)/max(1,np.linalg.norm(fd)));assert err<.01,err;errs.append(err)
 dump(out/'jacobian_checks.json',dict(directional_relative_errors=errs,warning='Payoff/interval kinks can affect finite differences.'))
 best=[r0['total'],start.copy()];ledger=[];optbegin=time.perf_counter();budget=1800
 def evaluate(x):
  if time.perf_counter()-optbegin>budget:raise TimeoutError('Declared 30-minute optimization budget exhausted')
  tick=time.perf_counter();r=engine.predict(x);evals.append(dict(index=len(evals),parameters=x.tolist(),total=r['total'],seconds=time.perf_counter()-tick))
  if r['total']<best[0]:best[:]=[r['total'],x.copy()];dump(out/'best_so_far.json',dict(parameters=x.tolist(),total=r['total']))
  dump(out/'evaluations.json',evals);print('evaluation',len(evals),'loss',r['total'],'seconds',time.perf_counter()-tick,flush=True);return r
 starts=[start,start.copy()];starts[1][4]=.5;starts[1][2]=1.5;starts[1][5]=5.5
 for s in starts:
  cache={}
  def current(x):
   key=tuple(x)
   if key not in cache:cache.clear();cache[key]=evaluate(x)
   return cache[key]
  try:
   fit=least_squares(lambda x:current(x)['residuals'],s,jac=lambda x:current(x)['jacobian'],bounds=(LOW,HIGH),x_scale=SCALE,max_nfev=35,ftol=1e-5,xtol=1e-5,gtol=1e-5)
   ledger.append(dict(start=s.tolist(),parameters=fit.x.tolist(),success=fit.success,nfev=fit.nfev,message=fit.message,optimality=fit.optimality))
  except TimeoutError as e:ledger.append(dict(status='budget exhausted',message=str(e)));break
 selected=best[1];assert selected[4]>0;dump(out/'selected.json',dict(parameters=selected.tolist(),levels=selected[7:10].tolist(),checking_not_used=True,optimizer_seconds=time.perf_counter()-optbegin));dump(out/'optimizer_ledger.json',ledger)
 record(engine,'selected',selected,'search',12);del engine
 prepare(out,'check',p,initial,pars,baseline[7:10],spx,vix,256,args.seed)
 check=Engine(out,'check',k,p,initial,spx,vix,forwards,design['ends_years'])
 record(check,'baseline',baseline,'check',24);record(check,'selected',selected,'check',24)
 record(check,'baseline',baseline,'check_gh12',12);record(check,'selected',selected,'check_gh12',12)
 check_sources(captured)
 final=design|dict(status='completed bounded eleven-parameter calibration and separate numerical check',selected_parameters=selected.tolist(),selected_levels=selected[7:10].tolist(),scores=records,optimizer_ledger=ledger,wall_seconds=time.perf_counter()-begin,limitations=['Finite-path and GH/time/stock-grid errors unresolved.','Fixed normalizer interpolation probes are sampled checks, not uniform bounds.','P compatibility has not been re-established for the new readout.','Search conditional GH12 and checking GH24 differ; matching-state GH12 diagnostics retained.','No formal certificate or market out-of-sample claim.'])
 dump(out/'run.json',final)
 rows=[]
 for role in ['search','check']:
  for name in ['baseline','selected']:
   q=pd.read_csv(out/f'prices_{role}_{name}.csv')
   for (instrument,days),g in q.groupby(['instrument','expiry_days']):rows.append(dict(role=role,candidate=name,instrument=instrument,expiry_days=int(days),quotes=len(g),interval_MSE=float(np.mean(g.interval_residual**2)),inside_bidask=int(g.interval_residual.eq(0).sum())))
 pd.DataFrame(rows).to_csv(out/'expiry_scores.csv',index=False);print('COMPLETED',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);parser.add_argument('--seed',type=int,default=2026102409);args=parser.parse_args()
 try:run(args)
 except Exception as e:
  if args.output.exists():dump(args.output/'failure.json',dict(type=type(e).__name__,message=str(e)))
  raise
