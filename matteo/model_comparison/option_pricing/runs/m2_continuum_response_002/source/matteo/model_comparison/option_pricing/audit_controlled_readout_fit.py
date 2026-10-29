"""Independent saved-price, objective, shape and source replay."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT,frozen_stage3
from dataclasses import asdict
from .all_maturity_core import InSampleLoss,populations
from .full_readout_fit import dump
from .audit_continuous_strike_controls import shape
from esn_eq.market_targets import MarketTargets


def audit(run):
 r=json.loads((run/'run.json').read_text());LOWER=np.array(r['bounds_lower']);UPPER=np.array(r['bounds_upper']);cache=Path(r['cache']);sources={key:hashlib.sha256((ROOT/key).read_bytes()).hexdigest()==v for key,v in r['implementation_sha256'].items()};assert all(sources.values());cr=json.loads((cache/'run.json').read_text());assert cr['status']=='complete immutable per-expiry 512-path cache'
 cs={key:hashlib.sha256((ROOT/key).read_bytes()).hexdigest()==v for key,v in cr['source_sha256'].items()};assert all(cs.values())
 assert r['architecture']==cr['architecture'] and r['premium']==cr['premium'];np.testing.assert_array_equal(r['initial_state'],cr['initial_state'])
 architecture,parameters,premium,frozen_hashes=frozen_stage3();assert json.loads(json.dumps(asdict(architecture)))==r['architecture'];assert json.loads(json.dumps(asdict(premium)))==r['premium']
 targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');assert r['target_sha256']==cr['target_sha256']==targets.artifact_sha256;spx,vix,_,_=populations(targets.quotes);quotes=pd.concat([spx,vix],ignore_index=True);f=targets.vix_forwards.sort_values('expiry').reset_index(drop=True);loss=InSampleLoss(quotes,f);assert (len(spx),len(vix))==(1294,217);selected=np.array(r['selected_parameters']);assert np.all(selected>=LOWER) and np.all(selected<=UPPER) and selected[4]>0;baseline=np.array(r['baseline']);warm=np.array(r['starts'][0]);np.testing.assert_array_equal(np.delete(warm,4),np.delete(baseline,4));assert warm[4]==LOWER[4]
 records=[]
 for role in ['search','check']:
  for name in ['baseline','warm','selected']:
   q=pd.read_csv(run/f'prices_{role}_{name}.csv');fw=pd.read_csv(run/f'forwards_{role}_{name}.csv');a=np.load(run/f'paths_{role}_{name}.npz');pd.testing.assert_frame_equal(q[['instrument','expiry','strike','option_type']],quotes[['instrument','expiry','strike','option_type']],check_dtype=False);shapes=[];worst=0.
   for (instrument,days),g in q.groupby(['instrument','expiry_days'],sort=False):
    if instrument=='SPX':values=a[f'spx_{days}'];mean=values.mean(0);se=values.std(0,ddof=1)/np.sqrt(len(values));assert len(values)==512
    else:values=a[f'vix_values_{days}'];blocks=values.reshape(8,64,len(g)).mean(1);mean=blocks.mean(0);se=blocks.std(0,ddof=1)/np.sqrt(8)
    np.testing.assert_allclose(mean,g.model_price_U,rtol=1e-12,atol=1e-10);np.testing.assert_allclose(se,g.mc_se_U,rtol=1e-12,atol=1e-10);worst=max(worst,float(max(abs(mean-g.model_price_U))));shapes.append(dict(instrument=instrument,days=int(days),**shape(g,mean)))
   for row in fw.itertuples():
    b=a[f'vix_forward_values_{row.expiry_days}'].reshape(8,64).mean(1);np.testing.assert_allclose([b.mean(),b.std(ddof=1)/np.sqrt(8)],[row.model_forward,row.mc_se],rtol=1e-12,atol=1e-12)
   p=baseline if name=='baseline' else warm if name=='warm' else selected;computed=loss.evaluate(q.model_price_U,fw.model_forward,p[7:10]);saved=next(z for z in r['scores'] if z['role']==role and z['candidate']==name);np.testing.assert_allclose(computed['total']+.1/3,saved['total'],rtol=1e-12)
   records.append(dict(role=role,candidate=name,total=saved['total'],maximum_saved_mean_error=worst,shapes=shapes,normalizer=saved['normalizer']))
 searchbase=next(z for z in r['scores'] if z['role']=='search' and z['candidate']=='warm');searchsel=next(z for z in r['scores'] if z['role']=='search' and z['candidate']=='selected');assert searchsel['total']<=searchbase['total']+1e-10
 # Eight paired independent-block comparisons, descriptive only.
 blockprices={};blockforwards={}
 for name in ['baseline','warm','selected']:
  a=np.load(run/f'paths_check_{name}.npz');blocks=[];fb=[]
  for days,g in spx.groupby('expiry_days'):blocks.append(a[f'spx_{days}'].reshape(8,64,len(g)).mean(1))
  for days,g in vix.groupby('expiry_days'):blocks.append(a[f'vix_values_{days}'].reshape(8,64,len(g)).mean(1));fb.append(a[f'vix_forward_values_{days}'].reshape(8,64).mean(1))
  blockprices[name]=np.concatenate(blocks,axis=1);blockforwards[name]=np.array(fb).T
 paired={}
 for reference,p in [('baseline',baseline),('warm',warm)]:
  differences=[]
  for i in range(8):
   mask=np.arange(8)!=i;bl=loss.evaluate(blockprices[reference][mask].mean(0),blockforwards[reference][mask].mean(0),p[7:10])['total'];sl=loss.evaluate(blockprices['selected'][mask].mean(0),blockforwards['selected'][mask].mean(0),selected[7:10])['total'];differences.append(sl-bl)
  se=float(np.sqrt(7/8*np.sum((differences-np.mean(differences))**2)));paired[reference]=dict(descriptive_loss_change_SE=se,leave_one_block_changes=differences)
 result=dict(status='passed independent saved-price/SE/objective/scope/source/warm-start replay',records=records,source_count=len(sources),cache_source_count=len(cs),quotes=1511,parameters_in_bounds=True,beta_J_positive=True,starting_incumbent_retained=True,paired_checking_comparisons=paired,limitations='Saved-array replay does not regenerate paths or certify accuracy. Shape failures are reported rather than asserted away. Existing checking streams are not new independent full-history validation.')
 # Final selected checking point independently evaluated with direct tensor GH32.
 fast=pd.read_csv(run/'prices_check_selected.csv');direct=pd.read_csv(run/'prices_check_direct_selected.csv');scale=fast.price_scale_U.to_numpy();changes=direct.model_price_U.to_numpy()-fast.model_price_U.to_numpy();maximum=float(np.max(abs(changes)/scale));assert maximum<1e-3,maximum
 dr=np.load(run/'paths_check_direct_selected.npz')
 for (instrument,days),g in direct.groupby(['instrument','expiry_days'],sort=False):
  values=dr[f'spx_{days}'] if instrument=='SPX' else dr[f'vix_values_{days}'];mean=values.mean(0);np.testing.assert_allclose(mean,g.model_price_U,rtol=1e-12,atol=1e-10)
  se=values.std(0,ddof=1)/np.sqrt(512) if instrument=='SPX' else values.reshape(8,64,len(g)).mean(1).std(0,ddof=1)/np.sqrt(8);np.testing.assert_allclose(se,g.mc_se_U,rtol=1e-12,atol=1e-10)
 dfor=pd.read_csv(run/'forwards_check_direct_selected.csv');ffor=pd.read_csv(run/'forwards_check_selected.csv');dscore=loss.evaluate(direct.model_price_U,dfor.model_forward,selected[7:10])['total']+.1/3;stored=next(z for z in r['scores'] if z['role']=='check_direct');np.testing.assert_allclose(dscore,stored['total'],rtol=1e-12)
 for row in dfor.itertuples():
  b=dr[f'vix_forward_values_{row.expiry_days}'].reshape(8,64).mean(1);np.testing.assert_allclose([b.mean(),b.std(ddof=1)/np.sqrt(8)],[row.model_forward,row.mc_se],rtol=1e-12,atol=1e-12)
 result['direct_GH32_replay']=dict(maximum_price_change_in_quote_scales=maximum,maximum_forward_change=float(np.max(abs(dfor.model_forward-ffor.model_forward))),loss=dscore)
 manifest=json.loads((cache/'immutable_manifest.json').read_text());assert manifest['before_first_saved_optimization_evaluation']
 for key,entry in manifest['files'].items():
  h=hashlib.sha256()
  with (cache/key).open('rb') as stream:
   for block in iter(lambda:stream.read(16*1024*1024),b''):h.update(block)
  assert h.hexdigest()==entry['sha256'],key
 result['immutable_cache_arrays_verified']=len(manifest['files'])
 result['frozen_stage3_endpoint_sha256']=frozen_hashes
 moments=[]
 for name,p in [('baseline',baseline),('warm',warm),('selected',selected)]:
  arrays=np.load(run/f'paths_check_{name}.npz')
  for days in [27,55,118,153]:
   data=np.load(cache/'check'/f'vix_d{days}.npz');v=arrays[f'vix_{days}'];weights=data['weights'];contribution=weights*v*v;blocks=contribution.reshape(8,64).mean(1);t=days/365;end=t+30/365;left=np.array([0,.25,.75]);right=np.array(r['ends_years']);length=np.maximum(0,np.minimum(end,right)-np.maximum(t,left));target=float(1e4*length@p[7:10]/(30/365));estimate=float(blocks.mean());se=float(blocks.std(ddof=1)/np.sqrt(8));moments.append(dict(candidate=name,days=days,target=target,estimate=estimate,SE=se,residual=estimate-target,largest_contribution_share=float(contribution.max()/contribution.sum())))
 result['VIX_second_moment_checking']=moments
 stocks=[]
 for name in ['baseline','warm','selected']:
  arrays=np.load(run/f'paths_check_{name}.npz')
  for days in [29,57,119,211,302,666,1030]:
   w=np.load(cache/'check'/f'd{days}'/'weights.npy');y=w*arrays[f'stock_{days}'];stocks.append(dict(candidate=name,days=days,weighted_stock_mean=float(y.mean()),target=1.,iid_SE=float(y.std(ddof=1)/np.sqrt(len(y))),mean_importance_mass=float(w.mean()),largest_weighted_stock_share=float(y.max()/y.sum())))
 result['SPX_forward_moment_checking']=stocks
 tail=[];arrays=np.load(run/'paths_check_selected.npz')
 for (instrument,days),g in quotes.groupby(['instrument','expiry_days'],sort=False):
  if instrument=='SPX':positive=arrays[f'spx_raw_{days}']>0
  else:
   v=arrays[f'vix_{days}'];sign=np.where(g.option_type=='call',1.,-1.);positive=sign*(v[:,None]-g.strike.to_numpy())>0
  tail.append(dict(instrument=instrument,days=int(days),contracts=len(g),zero_raw_payoff_contracts=int((positive.sum(0)==0).sum()),minimum_raw_hit_count=int(positive.sum(0).min())))
 result['selected_checking_tail_hits']=tail
 dump(run/'score_audit.json',result);return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();print(json.dumps(audit(a.run),indent=2))
