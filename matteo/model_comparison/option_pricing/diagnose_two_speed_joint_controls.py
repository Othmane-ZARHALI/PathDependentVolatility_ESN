"""Read-only variance comparison with existing frozen SPX hedge controls.

No estimator selection, parameter fit, or price substitution is performed.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from .prepare_two_speed_joint import CACHE


def diagnose(run):
    rows=[]
    for phase,role in [('training','search'),('checking','check')]:
        for name in ['vix_only','joint']:
            q=pd.read_csv(run/phase/f'prices_spx_{name}.csv')
            with np.load(run/phase/f'stocks_{name}.npz') as a:
                for day,g in q.groupby('expiry_days',sort=True):
                    S=a[f'd{day}_stock'];w=a[f'd{day}_weights'];side=np.where(g.option_type.eq('call'),1.,-1.);raw=w[:,None]*np.maximum(side[None,:]*(S[:,None]-g.K_over_F.to_numpy()[None,:]),0)*g.forward.to_numpy()[None,:];c=np.load(CACHE/role/f'd{day}'/'correction.npy');assert c.shape==raw.shape
                    for estimator,y in [('current_reference',a[f'd{day}_values']),('legacy_frozen_hedges',raw-c)]:
                        blocks=y.reshape(8,-1,len(g)).mean(1);se=blocks.std(0,ddof=1)/np.sqrt(8);rows.append(dict(phase=phase,candidate=name,day=int(day),estimator=estimator,median_SE_scales=float(np.median(se/g.price_scale_U)),SE_gate_passes=int(np.sum(2*se<=g.price_scale_U.to_numpy()/4)),quotes=len(g)))
    pd.DataFrame(rows).to_csv(run/'control_variance_diagnostic.csv',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();diagnose(a.run)
