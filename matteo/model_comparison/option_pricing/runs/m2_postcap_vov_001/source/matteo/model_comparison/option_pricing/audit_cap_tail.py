"""Independent saved-draw price/parity/SE/loss replay and paired sensitivity."""
import json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from .full_readout_fit import HERE,dump
from .audit_maturity_curve import replay,loss_factory
from .cap_tail_study import VARIANTS


def blocks(folder,name,quotes):
    bs=[];fs=[];rawstats=[]
    with np.load(folder/f'paths_{name}.npz') as a:
        for (ins,day),g in quotes.groupby(['instrument','expiry_days'],sort=False):
            if ins=='VIX':
                bs.append(a[f'vix_native_blocks_{day}']);fs.append(a[f'vix_forward_values_{day}'].reshape(8,512).mean(1))
            else:
                call=a[f'spx_{day}'].reshape(8,512,len(g)).mean(1);adj=np.where(g.option_type.eq('put'),g.forward-g.strike,0);bs.append(call-adj)
                stock=a[f'stock_{day}'];w=a[f'weights_spx_{day}'];weighted=w*stock;rawstats.append(dict(days=int(day),maximum_stock_ratio=float(stock.max()),weighted_forward_ratio=float(weighted.mean()),raw_forward_SE=float(weighted.std(ddof=1)/64),largest_weighted_stock_share=float(weighted.max()/weighted.sum())))
    return np.concatenate(bs,axis=1),np.stack(fs,axis=1),rawstats


def paired(folder,name,baseline_folder,baseline_name,quotes,forwards,levels,baselevels):
    b,f,raw=blocks(folder,name,quotes);bb,bf,baseraw=blocks(baseline_folder,baseline_name,quotes);objective=loss_factory(quotes,forwards);values=[]
    for omit in [None,*range(8)]:
        keep=np.ones(8,dtype=bool) if omit is None else np.arange(8)!=omit
        r=objective.evaluate(b[keep].mean(0),f[keep].mean(0),levels);base=objective.evaluate(bb[keep].mean(0),bf[keep].mean(0),baselevels)
        values.append([r['total']-base['total'],r['components']['SPX']-base['components']['SPX'],r['components']['VIX']+r['components']['forwards']-base['components']['VIX']-base['components']['forwards']])
    values=np.array(values);se=np.sqrt(7/8*np.sum((values[1:]-values[1:].mean(0))**2,axis=0))
    return dict(order=['native_joint','SPX','native_VIX_plus_forwards'],change=values[0].tolist(),jackknife_SE=se.tolist(),leave_one_block_changes=values[1:].tolist(),raw_stock=raw,baseline_raw_stock=baseraw,limitation='Descriptive eight-block paired sensitivity, same reused primitives and fixed controls; not a confidence or convergence certificate.')


def run():
    out=HERE/'runs/m2_cap_tail_001';prior=HERE/'runs/m2_maturity_curve_fit_001';d=json.loads((prior/'run.json').read_text());p=np.array(d['selected']);cache=Path(d['cache']);spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');quotes=pd.concat([spx,vix],ignore_index=True);f=pd.read_csv(cache/'forward_targets.csv');records=[]
    for name,*_ in VARIANTS:
        score=json.loads((out/f'score_{name}.json').read_text());x=np.array(score['parameters']);r=replay(out,name,x[7:-1],vix);records.append(dict(candidate=name,role='VIX fixed-curve profile',native_loss=r['original_otm_total']))
    for name in ['cluster050','cluster100']:
        folder=out/'cluster_vov';r=replay(folder,name,p[7:-1],vix);records.append(dict(candidate=name,role='VIX cluster profile',native_loss=r['original_otm_total']))
    if (out/'followup/run.json').exists():
        for name in ['higher_L','logtail025']:
            folder=out/'followup';x=np.array(json.loads((folder/f'parameters_{name}.json').read_text())['fitted']);assert x[4]>0 and x[-1]>=.0025 and np.all(x[7:-1]>=x[-1]) and np.all(x[7:-1]<=.2)
            r=replay(folder,name+'_fit',x[7:-1],quotes);records.append(dict(candidate=name+'_fit',role='full SPX/VIX VIX-focused fit',native_loss=r['original_otm_total']))
            dump(folder/f'paired_{name}.json',paired(folder,name+'_fit',prior,'selected',quotes,f,x[7:-1],p[7:-1]))
            for label in [name,name+'_fit']:
                score=json.loads((folder/'checking_vix'/f'score_{label}.json').read_text());xx=np.array(score['parameters']);r=replay(folder/'checking_vix',label,xx[7:-1],vix);records.append(dict(candidate=label,role='reused separate checking VIX',native_loss=r['original_otm_total']))
                dump(folder/f'checking_paired_{label}.json',paired(folder/'checking_vix',label,prior/'checking_vix','selected',vix,f,xx[7:-1],p[7:-1]))
    checked=[]
    for folder in [out,out/'followup',out/'cluster_vov']:
        design=json.loads((folder/'design.json').read_text())
        for rel,sha in design['source_sha256'].items():
            file=folder/'source'/rel;assert hashlib.sha256(file.read_bytes()).hexdigest()==sha,rel;checked.append(rel)
    dump(out/'audit.json',dict(status='PASS saved-draw price/parity/SE/objective and frozen source replay',records=records,frozen_source_checks=len(checked),unique_frozen_sources=len(set(checked)),limitation='Does not independently regenerate primitives or certify sampler, continuous stock, GH/time errors, fitted-model performance or the whole simulator.'))
    print('PASS',len(records),'records',len(checked),'source checks')
if __name__=='__main__':run()
