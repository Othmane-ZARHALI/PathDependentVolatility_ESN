"""Fresh reusable exogenous paths for the single-payoff estimator experiment."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .full_readout_fit import HERE,dump
from .spx_volatility_sampling import proposal
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths

PARENT=HERE/'runs/m2_controlled_readout_cache_001'

def stock_expiry(out,days,paths,seed):
    out=Path(out);a,_,premium,_=frozen_stage3();p=DatedM2Pricer(a,premium);initial,_=p.initial_state([1,0,0])
    grid=np.load(PARENT/f'check/d{days}/batch0/grid.npy');basis=volatility_path_basis(p,grid,7/365)
    dest=out/'production'/f'd{days}';dest.mkdir(parents=True,exist_ok=True)
    records=[]
    for batch,start in enumerate(range(0,paths,128)):
        n=min(128,paths-start);s=int(np.random.SeedSequence([seed,days,batch]).generate_state(1)[0]);c=dest/f'batch{batch:03d}'
        if (c/'completed.json').exists():
            record=json.loads((c/'completed.json').read_text());assert record['seed']==s;records.append(record);continue
        c.mkdir(exist_ok=True);tick=time.perf_counter()
        r=prepare_volatility_paths(p,grid,initial,n,s,basis,proposal(1));cache=r.cache
        # Initial state already included: no additional translation.
        np.save(c/'features.npy',cache.projections);np.save(c/'equity.npy',cache.equity_increments)
        np.save(c/'grid.npy',grid);np.save(c/'weights.npy',r.weights);np.save(c/'events.npy',np.array([len(grid)-1],dtype='int32'))
        seg=np.minimum(np.searchsorted([.25,.75],grid[:-1],side='right'),2);dt=np.diff(grid)
        brownian=np.stack([cache.equity_increments[seg==j].sum(0) for j in range(3)],axis=1)
        duration=np.array([dt[seg==j].sum() for j in range(3)])
        np.savez(c/'reference.npz',brownian=brownian,duration=duration)
        record=dict(days=days,batch=batch,paths=n,seed=s,seconds=time.perf_counter()-tick)
        dump(c/'completed.json',record);records.append(record)
        print('CACHE',days,'paths',start+n,'/',paths,'seconds',record['seconds'],flush=True);del cache,r
    dump(dest/'completed.json',dict(days=days,paths=paths,records=records))
    return dict(days=days,seconds=sum(r['seconds'] for r in records),paths=paths)

def vix_cache(out,days,p,initial,seed):
    saved=np.load(PARENT/'check'/f'vix_d{days}.npz');T=days/365
    declaration=json.loads((HERE/'runs/m2_vix_targeted_check_001/proposals.json').read_text())[str(days)]
    sampler=GaussianOuterSampler(p,T,initial);mixture=DefensiveGaussianMixture(declaration['shifts'],declaration['probabilities'])
    for role,power in [('pilot',7),('production',9)]:
        states=[];weights=[];normals=[];records=[]
        for net in range(8):
            s=int(np.random.SeedSequence([seed,1 if role=='pilot' else 2,days,net]).generate_state(1)[0]);z,w,_=mixture.draw(power,s)
            states.append(sampler.mean+z@sampler.root.T);weights.append(w);normals.append(z);records.append(dict(net=net,seed=s,power=power))
        folder=out/role;folder.mkdir(exist_ok=True)
        np.savez_compressed(folder/f'vix_d{days}.npz',states=np.concatenate(states),weights=np.concatenate(weights),normals=np.concatenate(normals),
                            nets=np.repeat(np.arange(8),2**power),times=saved['times'],time_weights=saved['time_weights'])
        dump(folder/f'vix_d{days}_seeds.json',records)

def run(out,paths,seed):
    out.mkdir(parents=True,exist_ok=True);tick=time.perf_counter()
    a,_,premium,_=frozen_stage3();p=DatedM2Pricer(a,premium);initial,_=p.initial_state([1,0,0])
    sources=[Path(__file__),HERE/'parity_adaptive.py',HERE/'PARITY_ADAPTIVE_PLAN.md',PARENT/'design.json',PARENT/'spx_contracts.csv',PARENT/'vix_contracts.csv',
             HERE/'spx_volatility_sampling.py',HERE/'runs/m2_vix_targeted_check_001/proposals.json',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(sources);save_sources(out,captured,ROOT)
    design=dict(status='preparing fresh fixed reservoir',paths=paths,seed=seed,architecture=asdict(a),premium=asdict(premium),initial_state=initial.tolist(),
                source_sha256=source_hashes(captured,ROOT),stock_steps_per_trading_day=512,proposal='unchanged vol7 per expiry',
                vix_nets=8,vix_points_per_net=512,independent_vix_pilot_points_per_net=128,market_out_of_sample=False)
    if (out/'design.json').exists():
        old=json.loads((out/'design.json').read_text());assert old==design,'Resume design changed'
    else:dump(out/'design.json',design)
    for name in ['spx_contracts','vix_contracts','forward_targets','atm_proxies']:
        (out/f'{name}.csv').write_bytes((PARENT/f'{name}.csv').read_bytes())
    for days in [27,55,118,153]:vix_cache(out,days,p,initial,seed)
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(stock_expiry,str(out),int(days),paths,seed) for days in [1030,666,302,211,119,57,29]]
        results=[f.result() for f in futures]
    check_sources(captured);dump(out/'run.json',design|dict(status='complete reusable fresh reservoir',wall_seconds=time.perf_counter()-tick,expiries=results))
    print('CACHE COMPLETE',time.perf_counter()-tick,flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);a.add_argument('--paths',type=int,default=4096);a.add_argument('--seed',type=int,default=2026102409);args=a.parse_args()
    run(args.output.resolve(),args.paths,args.seed)
