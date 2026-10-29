"""Broader 11-parameter exploration with restored precision and continuous controls."""
import argparse
from concurrent.futures import ThreadPoolExecutor
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
from .full_readout_fit import ROOT,M2,HERE,NAMES,LOW,HIGH,SCALE,Kernel,arr,dump,checks
from .m2_spx import frozen_stage3
from .all_maturity_core import InSampleLoss
from .audit_continuous_strike_controls import shape
from .continuous_strike_controls import warm_start
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.pricing import DatedM2Pricer,_normal_rule
from esn_eq.states import symmetric_square_root

LOWER=LOW.copy();UPPER=HIGH.copy();LOWER[4]=1e-4;LOWER[10]=1e-8;UPPER[6]=30.

class ControlledEngine:
 def __init__(self,cache,role,kernel,dated,initial,spx,vix,forwards):
  self.cache,self.role,self.k,self.dated,self.initial,self.spx,self.vix,self.forwards=cache,role,kernel,dated,initial,spx,vix,forwards
  self.quotes=pd.concat([spx,vix],ignore_index=True);self.loss=InSampleLoss(self.quotes,forwards);self.folder=cache/role;self.ends=[.25,.75,2129/365];self.tau=1/max(dated.layout.rates);self.horizon=float(spx.T_years.max())+30/365
  self.pool=ThreadPoolExecutor(max_workers=4);self.normal_laws={};self.vixlaws={};self.last_readout=None;self.normal_cached=None;self.exposure_cached={}
  self.stock={}
  for days in sorted(spx.expiry_days.unique()):
   c=self.folder/f'd{days}';batches=sorted(c.glob('batch*'));self.stock[int(days)]=(batches,np.load(c/'weights.npy'),np.load(c/'correction.npy'))
 def close(self):self.pool.shutdown()
 def law(self,t,state,order):
  m,c=self.dated.gaussian_projection_law(t,state);z,w=_normal_rule(order);return arr(m),arr(z@symmetric_square_root(c).T),arr(w)
 def normalizer(self,p):
  key=tuple(p[:7])
  if key==self.last_readout and self.normal_cached is not None:return self.normal_cached
  self.last_readout=key;self.exposure_cached={};direct_values={}
  def direct(s):
   s=float(s)
   if s not in direct_values:
    if s not in self.normal_laws:
     t=float(self.tau*np.expm1(s));m,c=self.dated.gaussian_projection_law(t,self.initial);self.normal_laws[s]=(arr(m),symmetric_square_root(c))
    m,root=self.normal_laws[s];z,w=_normal_rule(48);a=self.k.expectation(m,arr(z@root.T),arr(w),p)[0]
    direct_values[s]=np.r_[np.log(a[0]),a[1:]/a[0]]
   return direct_values[s]
  knots=np.linspace(0,np.log1p(self.horizon/self.tau),257)
  for _ in range(10):
   values=np.array(list(self.pool.map(direct,knots)));spline=CubicSpline(knots,values,axis=0,extrapolate=False);mid=(knots[:-1]+knots[1:])/2;actual=np.array(list(self.pool.map(direct,mid)));errors=abs(np.expm1(spline(mid)[:,0]-actual[:,0]));bad=errors>1e-7
   probes=np.linspace(knots[0],knots[-1],130)[1:-1];probe_actual=np.array(list(self.pool.map(direct,probes)));probe_errors=abs(np.expm1(spline(probes)[:,0]-probe_actual[:,0]));probe_bad=probe_errors>1e-7
   if not bad.any() and not probe_bad.any():break
   knots=np.unique(np.r_[knots,mid[bad],probes[probe_bad]])
   if len(knots)>2048:raise RuntimeError('Candidate normalizer knot budget exhausted')
  else:raise RuntimeError('Candidate normalizer midpoint refinement failed')
  probes=np.linspace(knots[0],knots[-1],131)[1:-1];actual=np.array(list(self.pool.map(direct,probes)));err=float(np.max(abs(np.expm1(spline(probes)[:,0]-actual[:,0]))))
  if err>2e-7:raise RuntimeError(f'Candidate independent normalizer probe failed: {err}')
  audit=dict(knots=len(knots),midpoint_error=float(errors.max()),independent_error=err,GH=48,probe_note='Declared probes may enter refinement; additional final independent probes checked separately.')
  self.normal_cached=(lambda t:spline(np.log1p(np.asarray(t)/self.tau)),audit);return self.normal_cached
 def exposures(self,days,p,norm):
  if days in self.exposure_cached:return self.exposure_cached[days]
  if days not in self.vixlaws:
   a=np.load(self.folder/f'vix_d{days}.npz');T=days/365;times=a['times'];laws=[self.law(float(t-T),a['states'],32) for t in times]
   self.vixlaws[days]=(times,a['time_weights'],laws,a['weights'],a['nets'],a['correction'],a['forward_correction'])
  times,tw,laws,w,nets,correction,fc=self.vixlaws[days];expo=np.zeros((len(w),3));der=np.zeros((len(w),3,7))
  for t,wt,law in zip(times,tw,laws):
   E=self.k.expectation(*law,p);nn=norm(t);M=np.exp(nn[0]);ratio=E[:,0]/M;s=min(np.searchsorted(self.ends,t,side='right'),2);expo[:,s]+=wt*ratio;der[:,s]+=wt*(E[:,1:]/M-ratio[:,None]*nn[1:])
  result=(expo,der,w,nets,correction,fc);self.exposure_cached[days]=result;return result
 def predict(self,p,detail=False):
  norm,na=self.normalizer(p);prices=[];errors=[];jac=[];raw={};fwd=[];fse=[];fj=[];shapes=[]
  for days,q0 in self.spx.groupby('expiry_days',sort=True):
   q=q0.reset_index(drop=True);batches,w,c=self.stock[int(days)];term=[]
   for folder in batches:
    grid=np.load(folder/'grid.npy');term.append(self.k.stock(folder,norm(grid[:-1]),p)[0])
   term=np.concatenate(term);R=term[:,0];dR=term[:,1:];sign=np.where(q.option_type=='call',1.,-1.);distance=sign*(R[:,None]-q.K_over_F.to_numpy());pay=np.maximum(distance,0)*q.forward.to_numpy();uncontrolled=w[:,None]*pay;y=uncontrolled-c
   prices.extend(y.mean(0));errors.extend(y.std(0,ddof=1)/np.sqrt(len(y)));jac.extend(np.einsum('nk,nj->kj',w[:,None]*(distance>0)*sign*q.forward.to_numpy(),dR)/len(y));shapes.append(dict(instrument='SPX',days=int(days),**shape(q,y.mean(0))))
   if detail:raw[f'spx_{days}']=y;raw[f'spx_raw_{days}']=uncontrolled;raw[f'stock_{days}']=R
  for days,q0 in self.vix.groupby('expiry_days',sort=True):
   q=q0.reset_index(drop=True);expo,de,w,nets,c,fc=self.exposures(int(days),p,norm);delta=30/365;integral=p[10]*delta+expo@(p[7:10]-p[10]);di=np.zeros((len(w),11));di[:,:7]=np.einsum('nsj,s->nj',de,p[7:10]-p[10]);di[:,7:10]=expo;di[:,10]=delta-expo.sum(1);V=100*np.sqrt(integral/delta);dv=50/np.sqrt(delta*integral)[:,None]*di
   sign=np.where(q.option_type=='call',1.,-1.);distance=sign*(V[:,None]-q.strike.to_numpy());y=w[:,None]*np.maximum(distance,0)-c;b=y.reshape(8,64,len(q)).mean(1);prices.extend(b.mean(0));errors.extend(b.std(0,ddof=1)/np.sqrt(8));jac.extend(np.einsum('nk,nj->kj',w[:,None]*(distance>0)*sign,dv)/len(w));bf=(w*V-fc).reshape(8,64).mean(1);fwd.append(bf.mean());fse.append(bf.std(ddof=1)/np.sqrt(8));fj.append(np.mean(w[:,None]*dv,axis=0));shapes.append(dict(instrument='VIX',days=int(days),**shape(q,y.mean(0))))
   if detail:raw[f'vix_{days}']=V;raw[f'vix_values_{days}']=y;raw[f'vix_forward_values_{days}']=w*V-fc;raw[f'vix_blocks_{days}']=b
  prices,errors,jac,fwd,fse,fj=map(np.asarray,[prices,errors,jac,fwd,fse,fj]);r=self.loss.evaluate(prices,fwd,p[7:10]);outside=(prices<self.quotes.price_bid_U.to_numpy())|(prices>self.quotes.price_ask_U.to_numpy());pj=self.loss.sqrtweights[:,None]*outside[:,None]*jac/self.quotes.price_scale_U.to_numpy()[:,None];pen=np.zeros((2,11));pen[:,7:10]=np.sqrt(.1)*np.array([[-1,1,0],[0,-1,1]]);r['jacobian']=np.concatenate([pj,fj/(.5*2),pen,np.zeros((1,11))]);r['residuals']=np.r_[r['residuals'],np.sqrt(.1/3)];r['total']+=.1/3;r['components']['state_penalty']=.1/3;r.update(prices=prices,errors=errors,forwards=fwd,forward_errors=fse,raw=raw,shapes=shapes,normalizer=na);return r

def run(args):
 out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();cache=args.cache.resolve();crecord=json.loads((cache/'design.json').read_text());a,pars,premium,hashes=frozen_stage3();p=DatedM2Pricer(a,premium);initial,_=p.initial_state([1,0,0]);spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');forwards=pd.read_csv(cache/'forward_targets.csv');baseline=np.array([getattr(pars,n) for n in NAMES[:7]]+crecord['baseline_levels']+[1e-4]);warm=warm_start(baseline,LOWER,UPPER,NAMES);previous=np.array(json.loads((HERE/'runs/m2_full_readout_fit_001/run.json').read_text())['selected_parameters']);third=baseline.copy();third[[1,2,3,4,5,6]]=[.12,2.,.35,.5,6.5,6.];starts=[warm,previous,third]
 sources=[Path(__file__),HERE/'controlled_readout_kernel.cpp',HERE/'controlled_readout_cache.py',HERE/'continuous_strike_controls.py',HERE/'CONTROLLED_READOUT_EXPLORATION_PLAN.md',HERE/'full_readout_fit.py',HERE/'m2_spx.py',HERE/'all_maturity_core.py',HERE/'joint_curve_core.py',HERE/'audit_continuous_strike_controls.py',cache/'design.json',cache/'anchors.csv',*sorted((M2/'src/esn_eq').glob('*.py'))];captured=capture_sources(sources);save_sources(out,captured,ROOT);lib=out/'kernel.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'controlled_readout_kernel.cpp'),'-o',str(lib)],check=True);k=Kernel(lib)
 design=dict(status='starting restored-precision broader calibration',cache=str(cache),architecture=asdict(a),premium=asdict(premium),initial_state=initial.tolist(),parameter_names=NAMES,baseline=baseline.tolist(),baseline_levels=baseline[7:10].tolist(),starts=[s.tolist() for s in starts],bounds_lower=LOWER.tolist(),bounds_upper=UPPER.tolist(),ends_years=crecord['ends_years'],implementation_sha256=source_hashes(captured,ROOT),target_sha256=crecord['target_sha256'],budget=dict(spx_search_paths=512,spx_check_paths=512,vix_nets=8,vix_points_per_net=64,conditional_GH=32,normalizer_GH=48,time_order=8,stock_steps_per_trading_day=512,starts=3,curve_prefit_nfev=6,all_parameter_nfev=12,total_optimization_seconds=2700,per_start_seconds=900),estimator='Frozen zero-mean controls interpolated continuously path by path; signed estimates/shape failures retained. No nearest-donor switching or post-price smoothing.',market_out_of_sample=False,checking='Existing separate numerical streams; not new independent complete research validation.',free_parameters=11)
 dump(out/'design.json',design);dump(out/'kernel_checks.json',checks(k,warm));records=[];evaluations=[];ledger=[]
 for name in ['spx_contracts','vix_contracts','atm_proxies','spx_membership','forward_targets']:(out/(name+'.csv')).write_bytes((cache/(name+'.csv')).read_bytes())
 search=ControlledEngine(cache,'search',k,p,initial,spx,vix,forwards)
 def record(engine,name,x,role):
  tick=time.perf_counter();r=engine.predict(x,True);q,f=engine.loss.tables(r['prices'],r['errors'],r['forwards'],r['forward_errors']);q.to_csv(out/f'prices_{role}_{name}.csv',index=False);f.to_csv(out/f'forwards_{role}_{name}.csv',index=False);np.savez_compressed(out/f'paths_{role}_{name}.npz',**r['raw']);rec=dict(role=role,candidate=name,total=r['total'],components=r['components'],inside_bidask=int(q.interval_residual.eq(0).sum()),negatives=int(q.negative_estimate.sum()),nominal_se_passes=int(q.nominal_se_screen.sum()),normalizer=r['normalizer'],shapes=r['shapes'],seconds=time.perf_counter()-tick);records.append(rec);dump(out/'scores_so_far.json',records);print('RECORD',json.dumps(rec),flush=True);return r
 base_record=record(search,'baseline',baseline,'search')
 prior=np.load(HERE/'runs/m2_continuous_strike_controls_003/saved_estimator_draws.npz');replay=[]
 for days,g in spx.groupby('expiry_days',sort=True):
  previous=prior[f'search_raw_d{days}'].mean(0);current=base_record['raw'][f'spx_raw_{days}'].mean(0);change=float(np.max(abs(current-previous)/g.price_scale_U.to_numpy()));replay.append(dict(days=int(days),maximum_raw_price_change_in_quote_scales=change));assert change<.05,(days,change)
 dump(out/'baseline_replay.json',dict(note='Same original primitives/parameters; candidate GH48 interpolator differs within its declared tolerance.',SPX=replay))
 w=record(search,'warm',warm,'search')
 # Directional check across all eleven parameters before optimization.
 direction=np.array([.05,.02,.2,.03,.00002,.2,.5,.01,-.01,.01,.00001]);h=1e-4;fd=(search.predict(warm+h*direction)['residuals']-search.predict(warm-h*direction)['residuals'])/(2*h);err=float(np.linalg.norm(fd-w['jacobian']@direction)/max(1.,np.linalg.norm(fd)));dump(out/'jacobian_check.json',dict(relative_error=err,tolerance=.01));assert err<.01,err
 best=[w['total'],warm.copy()];optbegin=time.perf_counter();startbegin=optbegin
 def evaluate(x,stage,si):
  if time.perf_counter()-optbegin>2700 or time.perf_counter()-startbegin>900:raise TimeoutError('Declared optimization time budget exhausted')
  tick=time.perf_counter();r=search.predict(x);entry=dict(index=len(evaluations),start=si,stage=stage,parameters=x.tolist(),total=r['total'],components=r['components'],negative_count=r['negative_count'],shapes=r['shapes'],normalizer=r['normalizer'],seconds=time.perf_counter()-tick);evaluations.append(entry);dump(out/'evaluations.json',evaluations)
  if r['total']<best[0]:best[:]=[r['total'],x.copy()];dump(out/'best_so_far.json',entry)
  print('EVAL',len(evaluations),'start',si,stage,'loss',r['total'],'seconds',entry['seconds'],flush=True);return r
 for si,s in enumerate(starts):
  if time.perf_counter()-optbegin>2700:break
  startbegin=time.perf_counter();memo={};fixed=s.copy()
  def curvecurrent(levels):
   x=fixed.copy();x[7:10]=levels;key=tuple(x)
   if key not in memo:memo.clear();memo[key]=evaluate(x,'curve_prefit',si)
   return memo[key]
  try:
   cf=least_squares(lambda x:curvecurrent(x)['residuals'],fixed[7:10],jac=lambda x:curvecurrent(x)['jacobian'][:,7:10],bounds=(LOWER[7:10],UPPER[7:10]),x_scale=SCALE[7:10],max_nfev=6,ftol=1e-6,xtol=1e-6,gtol=1e-6);s=s.copy();s[7:10]=cf.x;memo={}
   def current(x):
    key=tuple(x)
    if key not in memo:memo.clear();memo[key]=evaluate(x,'all_11',si)
    return memo[key]
   fit=least_squares(lambda x:current(x)['residuals'],s,jac=lambda x:current(x)['jacobian'],bounds=(LOWER,UPPER),x_scale=SCALE,max_nfev=12,ftol=1e-6,xtol=1e-6,gtol=1e-6);ledger.append(dict(start=si,parameters=fit.x.tolist(),curve_prefit_nfev=cf.nfev,all_parameter_nfev=fit.nfev,success=fit.success,message=fit.message,optimality=fit.optimality))
  except TimeoutError as e:ledger.append(dict(start=si,status='time budget',message=str(e)))
  dump(out/'optimizer_ledger.json',ledger)
 selected=best[1];assert selected[4]>0;dump(out/'selected.json',dict(parameters=selected.tolist(),levels=selected[7:10].tolist(),search_loss=best[0],checking_not_used=True,optimizer_seconds=time.perf_counter()-optbegin));record(search,'selected',selected,'search');search.close();del search
 assert json.loads((cache/'run.json').read_text())['status']=='complete immutable per-expiry 512-path cache'
 check=ControlledEngine(cache,'check',k,p,initial,spx,vix,forwards);record(check,'baseline',baseline,'check');record(check,'warm',warm,'check');record(check,'selected',selected,'check');check.close();check_sources(captured)
 dump(out/'run.json',design|dict(status='completed bounded broader exploration with continuous controls',selected_parameters=selected.tolist(),selected_levels=selected[7:10].tolist(),scores=records,optimizer_ledger=ledger,wall_seconds=time.perf_counter()-begin,limitations=['Exploratory signed corrected prices; any shape/negative failures are retained.','No population price-convergence or global optimum certificate.','P compatibility unverified.','Existing checking streams reused; not independent validation of entire research history.']))
 rows=[]
 for role in ['search','check']:
  for name in ['baseline','warm','selected']:
   q=pd.read_csv(out/f'prices_{role}_{name}.csv')
   for (instrument,days),g in q.groupby(['instrument','expiry_days']):rows.append(dict(role=role,candidate=name,instrument=instrument,expiry_days=int(days),quotes=len(g),interval_MSE=float(np.mean(g.interval_residual**2)),inside_bidask=int(g.interval_residual.eq(0).sum())))
 pd.DataFrame(rows).to_csv(out/'expiry_scores.csv',index=False);print('COMPLETED',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 try:run(a)
 except Exception as e:
  if a.output.exists():dump(a.output/'failure.json',dict(type=type(e).__name__,message=str(e)))
  raise
