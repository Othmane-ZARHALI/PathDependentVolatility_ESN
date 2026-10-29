"""Export same-sample prefix precision diagnostics from saved checking draws."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from .run_parity_adaptive import ivs


def export(run):
    for name in ['baseline','selected']:
        q=pd.read_csv(run/f'prices_canonical_{name}.csv');a=np.load(run/f'paths_{name}.npz');rows=[]
        for (ins,day),g in q.groupby(['instrument','expiry_days'],sort=False):
            y=a[f'spx_{day}' if ins=='SPX' else f'vix_values_{day}']
            for n in [512,1024,2048,4096]:
                if ins=='SPX':price=y[:n].mean(0);se=y[:n].std(0,ddof=1)/np.sqrt(n)
                else:
                    b=y.reshape(8,512,len(g))[:,:n//8].mean(1);price=b.mean(0);se=b.std(0,ddof=1)/np.sqrt(8)
                t=g.copy();t['uniform_paths']=n;t['model_price_U']=price;t['mc_se_U']=se;t['model_iv_percent']=ivs(t,price);rows.append(t)
        pd.concat(rows,ignore_index=True).to_csv(run/f'uniform_budget_prices_{name}.csv',index=False)
    print('PREFIX TABLES SAVED',run)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--run',type=Path,required=True);args=a.parse_args();export(args.run.resolve())
