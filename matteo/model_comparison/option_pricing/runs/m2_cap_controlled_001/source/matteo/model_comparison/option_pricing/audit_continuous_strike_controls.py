"""Saved-draw experiment: continuous, known-zero control interpolation.

No repricing, new simulation, coefficient fit, optimization or price smoothing.
Market-only anchors at scaled moneyness .5/1.5/2.5 on each side are prespecified.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .continuous_strike_controls import correction_draws,affine_correction_draws,warm_start
from .m2_spx import ROOT
from .all_maturity_core import nearest_same_side
from .full_readout_fit import HERE,ROOT,NAMES,LOW,HIGH,dump
from .provenance import capture_sources,source_hashes,save_sources,check_sources


def shape(q,prices):
 violations=convex=0
 for side,g in q.assign(price=prices).groupby('option_type'):
  g=g.sort_values('strike');d=np.diff(g.price);violations+=int(np.sum(d< -1e-9) if side=='put' else np.sum(d>1e-9))
  slopes=d/np.diff(g.strike);convex+=int(np.sum(np.diff(slopes)< -1e-9))
 return dict(monotonicity_violations=violations,convexity_violations=convex,negative_estimates=int(np.sum(prices<0)))


def run(out):
 out.mkdir(parents=True,exist_ok=False);parent=HERE/'runs/m2_all_maturity_fit_001';record=json.loads((parent/'run.json').read_text());stateparent=Path(record['parent_run']);old=HERE/'runs/m2_joint_curve_fit_002';oldq=pd.read_csv(old/'spx_contracts.csv');quotes=pd.read_csv(parent/'spx_contracts.csv');saved=np.load(parent/'comparison_paths.npz')
 sources=[Path(__file__),HERE/'continuous_strike_controls.py',HERE/'test_continuous_strike_controls.py',parent/'run.json',parent/'comparison_paths.npz',parent/'spx_contracts.csv',old/'spx_contracts.csv',HERE/'full_readout_fit.py',HERE/'all_maturity_core.py']
 for role,oldrole in [('search','search'),('check','fresh')]:
  for days in [29,57,119,211,302,666,1030]:sources.append(stateparent/f'controls_{oldrole}_t0_d{days}.npz' if days in [29,57,211,302] else parent/f'new_spx_{role}_d{days}/controls.npz')
 captured=capture_sources(sources);save_sources(out,captured,ROOT);rows=[];archives={};anchorrows=[]
 for role,oldrole in [('search','search'),('check','fresh')]:
  for days,g in quotes.groupby('expiry_days',sort=True):
   q=g.reset_index(drop=True)
   if days in [29,57,211,302]:
    donors=oldq[oldq.expiry_days==days].reset_index(drop=True);ix=nearest_same_side(q,donors);a=np.load(stateparent/f'controls_{oldrole}_t0_d{days}.npz');c=a['controls'][:,ix];beta=a['beta'][ix]
   else:a=np.load(parent/f'new_spx_{role}_d{days}/controls.npz');c=a['controls'];beta=a['beta']
   frozen=np.einsum('njk,jk->nj',c,beta);oldvalues=saved[f'{role}_selected_spx_{days}'];raw=oldvalues+frozen;continuous=np.empty_like(frozen);affine=np.empty_like(frozen)
   for side,sub in q.groupby('option_type'):
    indices=sorted(set(int((sub.scaled_moneyness-z).abs().idxmin()) for z in [.5,1.5,2.5]),key=lambda j:q.strike.iloc[j]);anchors=q.strike.iloc[indices].to_numpy()
    continuous[:,sub.index]=correction_draws(sub.strike.to_numpy(),anchors,frozen[:,indices])
    affine[:,sub.index]=affine_correction_draws(sub.strike.to_numpy(),anchors,frozen[:,indices])
    for j in indices:anchorrows.append(dict(role=role,days=int(days),side=side,strike=float(q.strike.iloc[j])))
   newvalues=raw-continuous;archives[f'{role}_continuous_d{days}']=newvalues;archives[f'{role}_raw_d{days}']=raw
   archives[f'{role}_affine_d{days}']=raw-affine;archives[f'{role}_old_d{days}']=oldvalues
   for method,values in [('old_nearest',oldvalues),('continuous_controls',newvalues),('affine_controls',raw-affine),('raw',raw)]:
    mean=values.mean(0);se=values.std(0,ddof=1)/np.sqrt(len(values));row=dict(role=role,days=int(days),method=method,median_SE=float(np.median(se)),**shape(q,mean));rows.append(row)
   outtable=q[['expiry_days','strike','option_type','price_bid_U','price_ask_U','price_scale_U']].copy();outtable['old_price']=oldvalues.mean(0);outtable['continuous_price']=newvalues.mean(0);outtable['raw_price']=raw.mean(0);outtable.to_csv(out/f'{role}_d{days}.csv',index=False)
 df=pd.DataFrame(rows);df.to_csv(out/'shape_and_precision.csv',index=False);pd.DataFrame(anchorrows).to_csv(out/'control_anchors.csv',index=False);np.savez_compressed(out/'saved_estimator_draws.npz',**archives)
 full=json.loads((HERE/'runs/m2_full_readout_fit_001/run.json').read_text());x=warm_start(full['baseline'],LOW,HIGH,NAMES);dump(out/'warm_start.json',dict(parameter_names=NAMES,parameters=x.tolist(),only_change_from_previous='spike_shift 0 -> 0.01; all ten other entries unchanged',algorithm='scipy.optimize.least_squares, method=trf, analytic Jacobian; same objective at start and subsequent candidates'))
 check_sources(captured);summary=df.groupby(['role','method'])[['monotonicity_violations','convexity_violations','negative_estimates']].sum().reset_index().to_dict('records');dump(out/'run.json',dict(status='completed saved-draw continuity experiment; not a new calibration',implementation_sha256=source_hashes(captured,ROOT),summary=summary,construction='Piecewise-linear deterministic interpolation of frozen zero-mean control DRAWS at three prespecified market-only anchors per option side; endpoint-constant extrapolation.',guarantee='Removes discontinuous nearest-donor switching, preserves known-zero expectation under the original control assumptions; does not guarantee finite-sample shape.',adoption='Do not adopt until positivity, monotonicity, convexity and precision checks are assessed; no price clipping or projection.',new_simulation=False,optimization=False))
 print(df.to_string(index=False));print(json.dumps(summary,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.output)
