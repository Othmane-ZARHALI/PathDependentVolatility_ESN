"""Paired eight-net VIX loss sensitivity and model second-moment targets."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .full_readout_fit import dump
from .audit_maturity_curve import loss_factory


def check(run):
    r=json.loads((run/'run.json').read_text());folder=run/'checking_vix';scores={};moments=[]
    for name in ['baseline','selected']:
        p=np.array(r[name]);q=pd.read_csv(folder/f'prices_canonical_{name}.csv');o=pd.read_csv(folder/f'prices_otm_{name}.csv');f=pd.read_csv(folder/f'forwards_{name}.csv')
        with np.load(folder/f'paths_{name}.npz') as a:
            blocks=[];native=[];forward=[]
            for day,g in q.groupby('expiry_days',sort=False):
                blocks.append(a[f'vix_blocks_{day}']);native.append(a[f'vix_native_blocks_{day}']);forward.append(a[f'vix_forward_values_{day}'].reshape(8,512).mean(1))
                second=(a[f'weights_vix_{day}']*a[f'vix_{day}']**2).reshape(8,512).mean(1)
                ends=np.array(r['ends_days']);left=np.r_[0,ends[:-1]];length=np.maximum(0,np.minimum(ends,day+30)-np.maximum(left,day))
                target=float(10000*length@p[7:-1]/30)
                moments.append(dict(candidate=name,days=int(day),target_second_moment=target,estimated_second_moment=float(second.mean()),net_SE=float(second.std(ddof=1)/np.sqrt(8))))
            b=np.concatenate(blocks,axis=1);nb=np.concatenate(native,axis=1);fb=np.stack(forward,axis=1)
            scores[name]={}
            for convention,quotes,values in [('canonical',q,b),('native',o,nb)]:
                objective=loss_factory(quotes,f);all_values=[]
                for omit in [None,*range(8)]:
                    keep=np.arange(8)!=omit if omit is not None else np.ones(8,dtype=bool)
                    x=objective.evaluate(values[keep].mean(0),fb[keep].mean(0),p[7:-1])
                    all_values.append(x['components']['VIX']+x['components']['forwards'])
                scores[name][convention]=all_values
    result=dict(moments=moments,paired={})
    for convention in ['canonical','native']:
        a=np.array(scores['baseline'][convention]);b=np.array(scores['selected'][convention]);delta=b-a
        se=float(np.sqrt(7/8*np.sum((delta[1:]-delta[1:].mean())**2)))
        result['paired'][convention]=dict(before=float(a[0]),after=float(b[0]),change=float(delta[0]),leave_one_net_changes=delta[1:].tolist(),jackknife_SE=se)
    result['limitation']='Descriptive paired eight-net sensitivity on previously used separate streams; not market OOS, untouched holdout or convergence certificate.'
    dump(run/'checking_vix/uncertainty.json',result);pd.DataFrame(moments).to_csv(run/'checking_vix/second_moments.csv',index=False);print(json.dumps(result['paired'],indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();check(a.run.resolve())
