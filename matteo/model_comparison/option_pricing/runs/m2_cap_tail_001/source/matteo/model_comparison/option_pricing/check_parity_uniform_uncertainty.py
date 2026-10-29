"""Paired eight-block loss sensitivity and observed rare-stock contribution."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT
from .all_maturity_core import InSampleLoss
from .full_readout_fit import dump


def check(run):
    design=json.loads((run/'run.json').read_text());out={}
    results={}
    for name in ['baseline','selected']:
        q=pd.read_csv(run/f'prices_canonical_{name}.csv');f=pd.read_csv(run/f'forwards_{name}.csv');a=np.load(run/f'paths_{name}.npz')
        blocks=[]
        for (ins,day),g in q.groupby(['instrument','expiry_days'],sort=False):
            y=a[f'spx_{day}' if ins=='SPX' else f'vix_values_{day}'];blocks.append(y.reshape(8,512,len(g)).mean(1))
        b=np.concatenate(blocks,axis=1);bf=np.stack([a[f'vix_forward_values_{d}'].reshape(8,512).mean(1) for d in f.expiry_days],axis=1)
        loss=InSampleLoss(q,f);levels=design[name][7:10]
        full=loss.evaluate(b.mean(0),bf.mean(0),levels)['total']+.1/3
        loo=[loss.evaluate(np.delete(b,j,axis=0).mean(0),np.delete(bf,j,axis=0).mean(0),levels)['total']+.1/3 for j in range(8)]
        results[name]=np.array(loo)
        out[name]=dict(full_loss=full,leave_one_block_losses=loo)
    delta=results['selected']-results['baseline'];se=np.sqrt(7/8*np.sum((delta-delta.mean())**2))
    out['paired_change']=dict(full=out['selected']['full_loss']-out['baseline']['full_loss'],leave_one_block_changes=delta.tolist(),jackknife_SE=float(se),
                              interpretation='Descriptive sensitivity of nonlinear loss to eight independent path/net blocks, conditional on frozen controls. Not a global price-convergence certificate.')
    arrays={n:np.load(run/f'paths_{n}.npz') for n in ['baseline','selected']};rare=[]
    for day in [29,57,119,211,302,666,1030]:
        for name,a in arrays.items():
            r=a[f'stock_{day}'];w=a[f'weights_spx_{day}'];j=int(np.argmax(w*r));weighted=w*r
            rare.append(dict(candidate=name,days=day,path_index=j,stock_ratio=float(r[j]),weighted_stock_ratio=float(weighted[j]),
                             share_of_weighted_stock_sum=float(weighted[j]/weighted.sum()),raw_forward_ratio=float(weighted.mean()),
                             mass=float(w.mean()),mass_SE=float(w.std(ddof=1)/np.sqrt(len(w)))))
    out['largest_stock_contributions']=rare
    out['specific_caution']='Old reference SPX302 contains a stock ratio about 56.6 at path1222; this single sample drives a large wing-price/loss change. Controls are frozen near selected parameters, so baseline variance can be much larger. Do not describe the large loss reduction as a precise population improvement.'
    dump(run/'validation/block_uncertainty.json',out);print(json.dumps(out['paired_change'],indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--run',type=Path,required=True);args=a.parse_args();check(args.run.resolve())
