"""Direct tensor GH32/GH48 probes at selected low/median/high VIX states."""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from .m2_spx import frozen_stage3
from .full_readout_fit import HERE,Kernel,dump
from .parity_parameter_core import FitEngine
from esn_eq.pricing import DatedM2Pricer


def probe(run):
    r=json.loads((run/'run.json').read_text());p=np.array(r['selected']);cache=Path(r['check_cache']);saved=np.load(run/'checking/paths_selected.npz')
    spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv')
    a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial,_=dated.initial_state([1,0,0])
    e=FitEngine(cache,Kernel(run/'convolution.dylib'),dated,initial,spx,vix,f,controls=json.loads((run/'vix_controls.json').read_text()))
    directlib=run/'validation/direct_tensor.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'controlled_readout_kernel.cpp'),'-o',str(directlib)],check=True);direct=Kernel(directlib)
    records=[]
    try:
        norm,na=e.normalizer(p)
        for days in sorted(vix.expiry_days.unique()):
            v=saved[f'vix_{days}'];order=np.argsort(v);indices=[int(order[0]),int(order[len(order)//2]),int(order[-1])]
            a=np.load(cache/'production'/f'vix_d{days}.npz');states=a['states'][indices];values={}
            for gh in [32,48]:
                integral=np.zeros(len(indices))
                for t,wt in zip(a['times'],a['time_weights']):
                    m,offset,w=e.law(float(t-days/365),states,gh);expect=direct.expectation(m,offset,w,p)[:,0];M=np.exp(norm(t)[0]);segment=min(np.searchsorted([.25,.75,2129/365],t,side='right'),2)
                    integral+=wt*(p[10]+(p[7+segment]-p[10])*expect/M)
                values[gh]=100*np.sqrt(integral/(30/365))
            records.append(dict(days=int(days),state_indices=indices,production_GH32=v[indices].tolist(),direct_GH32=values[32].tolist(),direct_GH48=values[48].tolist(),
                                maximum_production_vs_direct32=float(np.max(abs(v[indices]-values[32]))),maximum_32_vs48_VIX_points=float(np.max(abs(values[32]-values[48])))))
    finally:e.close()
    result=dict(status='representative selected-state quadrature probes complete',records=records,normalizer=na,
                maximum_GH32_vs48_VIX_points=max(x['maximum_32_vs48_VIX_points'] for x in records),
                limitation='Twelve selected states only; not whole-surface convergence or a uniform error bound.')
    dump(run/'validation/selected_quadrature.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--run',type=Path,required=True);args=a.parse_args();probe(args.run.resolve())
