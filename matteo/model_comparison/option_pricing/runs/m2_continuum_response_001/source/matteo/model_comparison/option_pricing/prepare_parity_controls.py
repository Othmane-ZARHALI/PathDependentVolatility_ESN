"""Freeze independently piloted controls while the exogenous cache is built."""
import argparse
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .full_readout_fit import HERE,Kernel,dump
from .controlled_readout_fast_fit import ControlledEngine
from .parity_adaptive_cache import stock_expiry,PARENT
from .parity_hedge_controls import FrozenHedgeKernel,anchors_for,control_data,fit_anchors
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.pricing import DatedM2Pricer


def run(cache,out):
    out.mkdir(parents=True,exist_ok=True);begin=time.perf_counter()
    parent=HERE/'runs/m2_controlled_readout_fit_002';selected=np.array(json.loads((parent/'selected.json').read_text())['parameters'])
    cache_design=json.loads((cache/'design.json').read_text());seed=cache_design['seed']+1000
    sources=[Path(__file__),HERE/'parity_hedge_kernel.cpp',HERE/'parity_hedge_controls.py',HERE/'controlled_readout_kernel.cpp',
             HERE/'controlled_readout_convolution_kernel.cpp',HERE/'parity_adaptive.py',HERE/'parity_adaptive_cache.py',
             HERE/'PARITY_ADAPTIVE_PLAN_V2.md',parent/'selected.json',cache/'design.json',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(sources);save_sources(out,captured,ROOT)
    lib=out/'hedges.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'parity_hedge_kernel.cpp'),'-o',str(lib)],check=True)
    convolution=out/'convolution.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'controlled_readout_convolution_kernel.cpp'),'-o',str(convolution)],check=True)
    k=Kernel(convolution);hk=FrozenHedgeKernel(lib);a,_,premium,_=frozen_stage3();p=DatedM2Pricer(a,premium);initial,_=p.initial_state([1,0,0])
    spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');forward=pd.read_csv(cache/'forward_targets.csv')
    engine=ControlledEngine(PARENT,'check',k,p,initial,spx,vix,forward);norm,na=engine.normalizer(selected)
    design=dict(status='preparing independent pilot/frozen controls',parameters=selected.tolist(),pilot_seed=seed,pilot_paths_per_expiry=256,
                source_sha256=source_hashes(captured,ROOT),normalizer=na,hedge_trades_per_trading_day=4)
    dump(out/'design.json',design);coefficients={}
    for days,q in spx.groupby('expiry_days',sort=True):
        days=int(days);stock_expiry(out/'spx_pilot',days,256,seed);anchors=anchors_for(q);cs=[];ys=[]
        for folder in sorted((out/'spx_pilot/production'/f'd{days}').glob('batch*')):
            output=hk.evaluate(folder,norm,selected,anchors);c,y=control_data(folder,output,anchors,selected)
            cs.append(c);ys.append(y);np.savez_compressed(folder/'hedge_controls.npz',output=output,controls=c,base_draws=y)
        beta=fit_anchors(np.concatenate(cs),np.concatenate(ys));coefficients[str(days)]=dict(anchors=anchors.tolist(),beta=beta.tolist())
        dump(out/'coefficients.json',coefficients);print('PILOT',days,'coefficients frozen',flush=True)
    completed=set()
    while len(completed)<7*32:
        progress=False
        for days,q in spx.groupby('expiry_days',sort=True):
            days=int(days);anchors=np.array(coefficients[str(days)]['anchors'])
            for folder in sorted((cache/'production'/f'd{days}').glob('batch*')):
                key=(days,folder.name);dest=out/'production'/f'd{days}';dest.mkdir(parents=True,exist_ok=True);target=dest/f'{folder.name}.npz'
                if key in completed:continue
                if target.exists():completed.add(key);continue
                if not (folder/'completed.json').exists():continue
                tick=time.perf_counter();output=hk.evaluate(folder,norm,selected,anchors);c,y=control_data(folder,output,anchors,selected)
                np.savez_compressed(target,output=output,controls=c,base_anchor_draws=y)
                completed.add(key);progress=True;print('FROZEN',days,folder.name,len(completed),'/224 seconds',time.perf_counter()-tick,flush=True)
        if not progress:time.sleep(3)
    engine.close();check_sources(captured)
    dump(out/'run.json',design|dict(status='complete independently fitted frozen controls',wall_seconds=time.perf_counter()-begin,production_batches=len(completed)))
    print('CONTROLS COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',type=Path,required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();run(args.cache.resolve(),args.output.resolve())
