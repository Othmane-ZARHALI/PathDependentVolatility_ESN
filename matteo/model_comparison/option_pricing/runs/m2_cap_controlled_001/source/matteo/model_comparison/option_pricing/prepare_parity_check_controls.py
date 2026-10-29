"""Evaluate fixed training-pilot controls on independent checking primitives."""
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
from .parity_hedge_controls import FrozenHedgeKernel,control_data
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.pricing import DatedM2Pricer


def run(cache,out):
    out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();parent=HERE/'runs/m2_parity_adaptive_controls_001'
    design=json.loads((parent/'design.json').read_text());parameters=np.array(design['parameters']);coefficients=json.loads((parent/'coefficients.json').read_text())
    sources=[Path(__file__),HERE/'parity_hedge_kernel.cpp',HERE/'parity_hedge_controls.py',HERE/'controlled_readout_kernel.cpp',HERE/'controlled_readout_convolution_kernel.cpp',
             HERE/'controlled_readout_fast_fit.py',HERE/'PARITY_PARAMETER_FIT_PLAN.md',cache/'design.json',parent/'design.json',parent/'coefficients.json',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(sources);save_sources(out,captured,ROOT);(out/'coefficients.json').write_bytes((parent/'coefficients.json').read_bytes())
    lib=out/'hedges.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'parity_hedge_kernel.cpp'),'-o',str(lib)],check=True)
    conv=out/'convolution.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'controlled_readout_convolution_kernel.cpp'),'-o',str(conv)],check=True)
    a,_,premium,_=frozen_stage3();p=DatedM2Pricer(a,premium);initial,_=p.initial_state([1,0,0]);spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');forwards=pd.read_csv(cache/'forward_targets.csv')
    engine=ControlledEngine(HERE/'runs/m2_controlled_readout_cache_001','check',Kernel(conv),p,initial,spx,vix,forwards);norm,na=engine.normalizer(parameters)
    record=dict(status='evaluating frozen training-pilot controls on checking draws',parameters=parameters.tolist(),coefficient_parent=str(parent),
                coefficients_fitted_on_checking=False,source_sha256=source_hashes(captured,ROOT),normalizer=na)
    dump(out/'design.json',record);hk=FrozenHedgeKernel(lib);completed=set()
    while len(completed)<224:
        progress=False
        for day in sorted(spx.expiry_days.unique()):
            anchors=np.array(coefficients[str(day)]['anchors'])
            for folder in sorted((cache/'production'/f'd{day}').glob('batch*')):
                key=(int(day),folder.name);dest=out/'production'/f'd{day}';dest.mkdir(parents=True,exist_ok=True);target=dest/f'{folder.name}.npz'
                if key in completed:continue
                if not (folder/'completed.json').exists():continue
                tick=time.perf_counter();output=hk.evaluate(folder,norm,parameters,anchors);c,y=control_data(folder,output,anchors,parameters)
                np.savez_compressed(target,output=output,controls=c,base_anchor_draws=y);completed.add(key);progress=True
                print('CHECK CONTROL',day,folder.name,len(completed),'/224 seconds',time.perf_counter()-tick,flush=True)
        if not progress:time.sleep(3)
        if time.perf_counter()-begin>14400:raise TimeoutError('Checking cache/control preparation exceeded four hours')
    engine.close();check_sources(captured);dump(out/'run.json',record|dict(status='complete fixed controls on independent checking samples',batches=len(completed),wall_seconds=time.perf_counter()-begin))
    print('CHECK CONTROLS COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',type=Path,required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();run(args.cache.resolve(),args.output.resolve())
