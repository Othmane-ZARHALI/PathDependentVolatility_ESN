"""Read-only first-order mean-one M2 variance control; not an adopted estimator.

Its zero mean assumes the exact Q normalizer. A numerical-normalizer error gate
would be necessary before using it for qualified prices. It fails to improve
the pilot residual SE consistently, so no fit or price replacement is performed.
"""
import json
import numpy as np
import pandas as pd
from .full_readout_fit import HERE
from .two_speed_joint_study import load_spx


def run():
    root=HERE/'runs/m2_two_speed_precision_conditional_001';parent=HERE/'runs/m2_two_speed_precision_002';x=np.array(json.loads((root/'design.json').read_text())['baseline']);f=(.05+.07*x[19])**2;rows=[]
    for role in ['search','check']:
        data,_=load_spx(role);q=pd.read_csv(root/role/'prices_conditional.csv')
        for day,g in q.groupby('expiry_days'):
            pieces=[]
            for d in data[day]:
                coefficient=(.2-f)*(1-x[17])*x[d['seg']]*d['dt'];pieces.append(coefficient@(d['old']-1))
            deltaI=np.concatenate(pieces);t=np.load(parent/role/f'd{day}_target.npz');iv=t['out'][:,2];K=g.K_over_F.to_numpy();d1=-np.log(K)[None,:]/np.sqrt(iv[:,None])+.5*np.sqrt(iv[:,None]);vega=np.exp(-.5*d1*d1)/np.sqrt(2*np.pi)/(2*np.sqrt(iv[:,None]));control=t['weights'][:,None]*vega*deltaI[:,None];r=np.load(root/role/f'd{day}_correction.npz')['residual'];se0=r.reshape(8,-1,len(K)).mean(1).std(0,ddof=1)/np.sqrt(8);se1=(r-control).reshape(8,-1,len(K)).mean(1).std(0,ddof=1)/np.sqrt(8)
            rows.append(dict(role=role,day=int(day),median_residual_SE_ratio=float(np.median(se1/se0)),weighted_mean_deltaI=float(np.mean(t['weights']*deltaI)),adopted=False))
    pd.DataFrame(rows).to_csv(root/'mean_control_diagnostic.csv',index=False)


if __name__=='__main__':run()
