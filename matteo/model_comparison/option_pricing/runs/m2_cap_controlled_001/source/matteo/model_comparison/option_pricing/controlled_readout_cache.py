"""Replay original per-expiry primitives into a readout-independent cache."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .spx_conditional_control_study import setup
from .all_maturity_core import nearest_same_side
from .continuous_strike_controls import correction_draws
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from .spx_volatility_sampling import proposal
from .joint_state_core import projected_initial_shift
from esn_eq.pricing import DatedM2Pricer,ForwardVarianceCurve
from esn_eq.vix_pricing import ConditionalVIXPricer
from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths

HERE=Path(__file__).resolve().parent
OLD=HERE/'runs/m2_all_maturity_fit_001'
STATE=HERE/'runs/m2_joint_state_fit_003'
CURVE=HERE/'runs/m2_joint_curve_fit_002'

def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def interpolate_controls(q,draws,spx=True):
 out=np.empty_like(draws);records=[]
 for side,g in q.groupby('option_type'):
  if spx:idx=sorted(set(int((g.scaled_moneyness-z).abs().idxmin()) for z in [.5,1.5,2.5]),key=lambda j:q.strike.iloc[j])
  else:
   ix=g.sort_values('strike').index.to_numpy();idx=ix[np.round(np.linspace(0,len(ix)-1,3)).astype(int)]
  k=q.strike.iloc[idx].to_numpy();out[:,g.index]=correction_draws(g.strike,k,draws[:,idx])
  records.extend(dict(side=side,strike=float(z)) for z in k)
 return out,records

def write_features(folder,projections,equity,grid,shift):
 folder.mkdir(parents=True,exist_ok=False);out=np.lib.format.open_memmap(folder/'features.npy',mode='w+',dtype='float64',shape=projections.shape)
 for i in range(0,len(projections),2048):out[i:i+2048]=projections[i:i+2048]+shift[i:i+2048,None,:]
 out.flush();del out
 if isinstance(equity,Path):(folder/'equity.npy').symlink_to(equity.resolve())
 else:np.save(folder/'equity.npy',equity)
 np.save(folder/'grid.npy',grid);np.save(folder/'events.npy',np.array([len(grid)-1],dtype='int32'))

def new_worker(out,role,days):
 out=Path(out);a,pars,premium,_=frozen_stage3();p=DatedM2Pricer(a,premium);initial,_=p.initial_state([1,0,0]);zero=np.zeros_like(initial);source=OLD/f'new_spx_{role}_d{days}';grid=np.load(source/'grid.npy');basis=volatility_path_basis(p,grid,7/365);shift=projected_initial_shift(p,grid[:-1],initial);record=json.loads((source/'completed.json').read_text());weights=np.load(source/'controls.npz')['weights'];destination=out/role/f'd{days}';destination.mkdir(parents=True,exist_ok=True)
 for batch,seed in enumerate(record['seeds']):
  tick=time.perf_counter();r=prepare_volatility_paths(p,grid,zero,128,seed,basis,proposal(1));np.testing.assert_array_equal(r.weights,weights[batch*128:(batch+1)*128]);equitypath=Path(json.loads((OLD/'run.json').read_text())['scratch'])/f'{role}_d{days}_b{batch}_equity.npy';np.testing.assert_array_equal(r.cache.equity_increments,np.load(equitypath,mmap_mode='r'))
  write_features(destination/f'batch{batch}',r.cache.projections,equitypath,grid,shift);print(role,days,'batch',batch,'exact primitive replay',time.perf_counter()-tick,flush=True);del r
 dump(destination/'completed.json',dict(role=role,days=days,seeds=record['seeds'],exact_equity_and_weight_replay=True))

def run(out):
 out.mkdir(parents=True,exist_ok=False);tick=time.perf_counter();targets,p,parameters,zero,anchor,panel,master,_,hashes=setup();initial,_=p.initial_state([1,0,0]);old=json.loads((OLD/'run.json').read_text());parent=json.loads((STATE/'run.json').read_text());spx=pd.read_csv(OLD/'spx_contracts.csv');vix=pd.read_csv(OLD/'vix_contracts.csv');oldspx=pd.read_csv(CURVE/'spx_contracts.csv');oldvix=pd.read_csv(CURVE/'vix_contracts.csv');sources=[Path(__file__),HERE/'continuous_strike_controls.py',HERE/'m2_spx.py',HERE/'all_maturity_core.py',HERE/'spx_conditional_control_study.py',HERE/'joint_state_core.py',HERE/'spx_volatility_sampling.py',OLD/'run.json',OLD/'spx_contracts.csv',OLD/'vix_contracts.csv',STATE/'run.json',CURVE/'spx_contracts.csv',CURVE/'vix_contracts.csv',*sorted((M2/'src/esn_eq').glob('*.py'))]
 for role,oldrole in [('search','search'),('check','fresh')]:
  for days in [29,57,211,302]:sources.append(STATE/f'controls_{oldrole}_t0_d{days}.npz')
  for days in [119,666,1030]:sources += [OLD/f'new_spx_{role}_d{days}/controls.npz',OLD/f'new_spx_{role}_d{days}/grid.npy',OLD/f'new_spx_{role}_d{days}/completed.json']
  for days in [27,55,118,153]:sources.append(STATE/f'vix_state2_{oldrole}_t0_d{days}.npz')
 captured=capture_sources(sources);save_sources(out,captured,ROOT);design=dict(status='preparing fixed per-expiry cache',architecture=asdict(p.architecture),premium=asdict(p.premium),initial_state=initial.tolist(),spx_paths=512,vix_points_per_net=64,vix_nets=8,stock_steps_per_trading_day=512,source_sha256=source_hashes(captured,ROOT),target_sha256=targets.artifact_sha256,baseline_levels=old['selected_levels'],ends_years=old['ends_years'],parent=str(OLD),checking='Existing separate checking streams, not a new independent validation of research history.');dump(out/'design.json',design)
 for name in ['spx_contracts','vix_contracts','atm_proxies','spx_membership','forward_targets']:(out/(name+'.csv')).write_bytes((OLD/(name+'.csv')).read_bytes())
 curve=ForwardVarianceCurve(tuple(old['ends_years']),tuple(old['selected_levels']));vengine=ConditionalVIXPricer(p,parameters,curve,initial,normalizer_order=48);anchors=[]
 for role,oldrole in [('search','search'),('check','fresh')]:
  (out/role).mkdir()
  for days,g in spx.groupby('expiry_days'):
   q=g.reset_index(drop=True);dest=out/role/f'd{days}';dest.mkdir()
   if days in [29,57,211,302]:
    donor=oldspx[oldspx.expiry_days==days].reset_index(drop=True);ix=nearest_same_side(q,donor);saved=np.load(STATE/f'controls_{oldrole}_t0_d{days}.npz');c=saved['controls'][:,ix];beta=saved['beta'][ix]
   else:saved=np.load(OLD/f'new_spx_{role}_d{days}/controls.npz');c=saved['controls'];beta=saved['beta']
   correction,aa=interpolate_controls(q,np.einsum('njk,jk->nj',c,beta));np.save(dest/'correction.npy',correction);np.save(dest/'weights.npy',saved['weights']);anchors.extend(dict(role=role,instrument='SPX',days=int(days),**z) for z in aa)
   if days in [29,57,211,302]:
    grid=master[:np.searchsorted(master,q.T_years.iloc[0])+1];shift=projected_initial_shift(p,grid[:-1],initial)
    for b in range(4):
     stem=Path(parent['scratch'])/f'{oldrole}_t0_d{days}_b{b}';proj=np.load(str(stem)+'_projections.npy',mmap_mode='r');write_features(dest/f'batch{b}',proj,Path(str(stem)+'_equity.npy'),grid,shift)
    print(role,days,'copied original projections',flush=True)
  for days,g in vix.groupby('expiry_days'):
   q=g.reset_index(drop=True);donor=oldvix[oldvix.expiry_days==days].reset_index(drop=True);ix=nearest_same_side(q,donor);saved=np.load(STATE/f'vix_state2_{oldrole}_t0_d{days}.npz');z=saved['controls']@saved['beta'][:,1+ix];correction,aa=interpolate_controls(q,z,False);window=vengine.gaussian_window(float(q.T_years.iloc[0]),time_order=8)
   np.savez_compressed(out/role/f'vix_d{days}.npz',states=saved['states'],weights=saved['weights'],nets=np.repeat(np.arange(8),64),times=window.times,time_weights=window.weights,correction=correction,forward_correction=saved['controls']@saved['beta'][:,0]);anchors.extend(dict(role=role,instrument='VIX',days=int(days),**z) for z in aa)
 pd.DataFrame(anchors).to_csv(out/'anchors.csv',index=False)
 with ProcessPoolExecutor(max_workers=2) as pool:
  futures=[pool.submit(new_worker,str(out),role,days) for role in ['search','check'] for days in [119,666,1030]]
  for f in futures:f.result()
 check_sources(captured);dump(out/'run.json',design|dict(status='complete immutable per-expiry 512-path cache',wall_seconds=time.perf_counter()-tick));print('CACHE COMPLETE',time.perf_counter()-tick,flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 try:run(a.output.resolve())
 except Exception as e:
  if a.output.exists():dump(a.output/'failure.json',dict(type=type(e).__name__,message=str(e)))
  raise
