"""Replay archived means/SEs, membership, monotonicity and source integrity."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .full_readout_fit import ROOT,HERE,LOW,HIGH,NAMES,dump

def audit(run):
 r=json.loads((run/'run.json').read_text());x=np.array(r['selected_parameters']);assert len(x)==11 and np.all(x>=LOW) and np.all(x<=HIGH) and x[4]>0
 hashes={key:hashlib.sha256((ROOT/key).read_bytes()).hexdigest()==value for key,value in r['implementation_sha256'].items()};assert all(hashes.values())
 from .all_maturity_core import populations,InSampleLoss
 from esn_eq.market_targets import MarketTargets
 targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');spx,vix,_,_=populations(targets.quotes);quotes=pd.concat([spx,vix],ignore_index=True);f=targets.vix_forwards.sort_values('expiry').reset_index(drop=True);loss=InSampleLoss(quotes,f)
 assert len(spx)==1294 and len(vix)==217
 checks=[]
 for role in ['search','check','check_gh12']:
  for name in ['baseline','selected']:
   q=pd.read_csv(run/f'prices_{role}_{name}.csv');fw=pd.read_csv(run/f'forwards_{role}_{name}.csv');raw=np.load(run/f'raw_{role}_{name}.npz');maxerr=0.;violations=0
   pd.testing.assert_frame_equal(q[['instrument','expiry','strike','option_type']],quotes[['instrument','expiry','strike','option_type']],check_dtype=False)
   for (instrument,days),g in q.groupby(['instrument','expiry_days'],sort=False):
    if instrument=='SPX':
     y=raw[f'spx_{days}'];mean=y.mean(0);se=y.std(0,ddof=1)/np.sqrt(len(y));assert np.all(y>=0)
    else:
     blocks=raw[f'vix_blocks_{days}'];mean=blocks.mean(0);se=blocks.std(0,ddof=1)/np.sqrt(8)
    np.testing.assert_allclose(mean,g.model_price_U,rtol=1e-12,atol=1e-10);np.testing.assert_allclose(se,g.mc_se_U,rtol=1e-12,atol=1e-10);maxerr=max(maxerr,float(np.max(abs(mean-g.model_price_U))))
    for side,h in g.groupby('option_type'):
     d=np.diff(h.sort_values('strike').model_price_U);violations+=int(np.sum(d< -1e-10) if side=='put' else np.sum(d>1e-10))
   for row in fw.itertuples():
    V=raw[f'vix_{row.expiry_days}'];w=raw[f'vix_weight_{row.expiry_days}'];b=(w*V).reshape(8,-1).mean(1);np.testing.assert_allclose([b.mean(),b.std(ddof=1)/np.sqrt(8)],[row.model_forward,row.mc_se],rtol=1e-12)
   assert violations==0
   p=np.array(r['baseline'] if name=='baseline' else r['selected_parameters']);sc=loss.evaluate(q.model_price_U,fw.model_forward,p[7:10]);saved=next(t for t in r['scores'] if t['role']==role and t['candidate']==name);np.testing.assert_allclose(sc['total']+.1/3,saved['total'],rtol=1e-12)
   checks.append(dict(role=role,candidate=name,maximum_mean_replay_error=maxerr,strike_monotonicity_violations=violations,total=saved['total']))
 seedrecords=[]
 for role in ['search','check']:seedrecords.extend(json.loads((run/role/'seeds.json').read_text()))
 seeds=[z['seed'] for z in seedrecords];assert len(set(seeds))==len(seeds)
 # Paired error diagnostics on the SAME checking states, no new pricing.
 low=pd.read_csv(run/'prices_check_gh12_selected.csv');high=pd.read_csv(run/'prices_check_selected.csv');gh=(high.model_price_U-low.model_price_U).abs()/high.price_scale_U
 pointwise=high[['instrument','expiry_days','strike','option_type','price_scale_U']].copy();pointwise['conditional_GH24_minus_GH12_U']=high.model_price_U-low.model_price_U;pointwise.to_csv(run/'conditional_order_changes.csv',index=False)
 base=np.load(run/'raw_check_baseline.npz');sel=np.load(run/'raw_check_selected.npz');bprices=[];sprices=[]
 for days in sorted(spx.expiry_days.unique()):
  bprices.append(base[f'spx_{days}'].reshape(8,-1,len(spx[spx.expiry_days==days])).mean(1));sprices.append(sel[f'spx_{days}'].reshape(8,-1,len(spx[spx.expiry_days==days])).mean(1))
 for days in sorted(vix.expiry_days.unique()):bprices.append(base[f'vix_blocks_{days}']);sprices.append(sel[f'vix_blocks_{days}'])
 bprices=np.concatenate(bprices,axis=1);sprices=np.concatenate(sprices,axis=1)
 bfor=np.array([(base[f'vix_{d}']*base[f'vix_weight_{d}']).reshape(8,-1).mean(1) for d in [27,55,118,153]]).T
 sfor=np.array([(sel[f'vix_{d}']*sel[f'vix_weight_{d}']).reshape(8,-1).mean(1) for d in [27,55,118,153]]).T
 changes=[]
 for i in range(8):
  mask=np.arange(8)!=i
  bl=loss.evaluate(bprices[mask].mean(0),bfor[mask].mean(0),np.array(r['baseline'])[7:10])['total'];sl=loss.evaluate(sprices[mask].mean(0),sfor[mask].mean(0),x[7:10])['total'];changes.append(sl-bl)
 se=float(np.sqrt(7/8*np.sum((changes-np.mean(changes))**2)))
 cachehash={str(p.relative_to(run)):hashlib.sha256(p.read_bytes()).hexdigest() for role in ['search','check'] for p in sorted((run/role).rglob('*.npy'))}
 for role in ['search','check']:
  for p in sorted((run/role).glob('*.npz')):cachehash[str(p.relative_to(run))]=hashlib.sha256(p.read_bytes()).hexdigest()
 dump(run/'cache_manifest.json',dict(sha256=cachehash,note='Post-run integrity capture; fixed cache is read-only in candidate evaluation.'))
 result=dict(status='passed saved-price replay, scope, seeds, monotonicity and source checks',records=checks,source_count=len(hashes),unique_seed_count=len(seeds),parameters_in_bounds=True,beta_J_positive=True,GH24_vs_GH12_selected_max_price_scales=float(gh.max()),GH24_vs_GH12_selected_median_price_scales=float(gh.median()),GH24_vs_GH12_selected_VIX_median_price_scales=float(gh[high.instrument=='VIX'].median()),paired_eight_block_jackknife_loss_change_SE=se,leave_one_block_loss_changes=changes,limitations='Replay does not regenerate paths or establish price convergence. Stock blocks and VIX scramble blocks are paired by index for descriptive joint jackknife, not temporal observations.')
 dump(run/'score_audit.json',result);return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();print(json.dumps(audit(a.run),indent=2))
