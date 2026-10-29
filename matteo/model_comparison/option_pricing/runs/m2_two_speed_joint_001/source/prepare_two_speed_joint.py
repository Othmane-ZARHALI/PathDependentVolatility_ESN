"""Coarsen reusable SPX primitives and rebuild the frozen M2 ratio."""
import argparse,json,hashlib,time,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from .full_readout_fit import HERE,dump,Kernel
from .controlled_readout_fast_fit import ControlledEngine
from .m2_spx import frozen_stage3
from .vix_only_curve_fit import PRIOR
from .maturity_curve import ENDS_DAYS
from esn_eq.pricing import DatedM2Pricer

CACHE=HERE/'runs/m2_controlled_readout_cache_001'
DAYS=[29,57,119,211,302,666,1030]


def indices(grid,stride):
    idx=np.r_[np.arange(0,len(grid)-1,stride),len(grid)-1];j=np.searchsorted(grid,ENDS_DAYS/365);j=j[(j>0)&(j<len(grid))]
    return np.unique(np.r_[idx,j,j-1])


def prepare(out):
    out.mkdir(parents=True,exist_ok=False);start=time.perf_counter();prior=json.loads((PRIOR/'run.json').read_text());cd=json.loads((CACHE/'design.json').read_text());assert cd['architecture']==prior['architecture'] and cd['premium']==prior['premium'] and cd['initial_state']==prior['initial_state']
    p=np.array(prior['selected']);ends=np.array(prior['ends_years']);a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);ctx=object.__new__(ControlledEngine);ctx.dated=dated;ctx.initial=np.array(prior['initial_state']);ctx.k=Kernel(PRIOR/'convolution.dylib');ctx.ends=ends;ctx.tau=1/max(dated.layout.rates);ctx.horizon=1060/365;ctx.pool=ThreadPoolExecutor(max_workers=4);ctx.normal_laws={};ctx.last_readout=None;ctx.normal_cached=None
    try:norm,na=ctx.normalizer(p)
    finally:ctx.close()
    dump(out/'normalizer_audit.json',na);hashes={};records=[]
    for role in ['search','check']:
        for day in DAYS:
            dest=out/role/f'd{day}';dest.mkdir(parents=True);source=CACHE/role/f'd{day}';weights=np.load(source/'weights.npy')
            for batch,folder in enumerate(sorted(source.glob('batch*'))):
                tick=time.perf_counter();grid=np.load(folder/'grid.npy');idx=indices(grid,8);times=grid[idx];feat=np.load(folder/'features.npy',mmap_mode='r');eq=np.load(folder/'equity.npy',mmap_mode='r');X=np.asarray(feat[idx[:-1]])
                eta=p[0]*X[:,:,0]+p[1]*(X[:,:,1]-p[2])**2+p[3]*(X[:,:,2]-p[4])**2;lg=2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6]);old=np.exp(lg-norm(times[:-1])[:,0,None]);dw=np.add.reduceat(eq,idx[:-1],axis=0)
                seg=np.minimum(np.searchsorted(ends,times[:-1],side='right'),15).astype('int32');refseg=np.minimum(np.searchsorted([.25,.75],times[:-1],side='right'),2)
                refbrownian=np.stack([dw[refseg==j].sum(0) for j in range(3)],axis=1);duration=np.array([np.diff(times)[refseg==j].sum() for j in range(3)])
                seed=int(np.random.SeedSequence([2026102801,0 if role=='search' else 1,day,batch]).generate_state(1,dtype=np.uint64)[0]);eps=np.random.default_rng(seed).standard_normal((len(dw),dw.shape[1],2));c=dest/f'batch{batch}';c.mkdir()
                for name,val in dict(old=old,equity=dw,eps=eps,time=times[:-1],dt=np.diff(times),seg=seg,grid=times,refbrownian=refbrownian,duration=duration,weights=weights[batch*128:(batch+1)*128]).items():np.save(c/f'{name}.npy',np.ascontiguousarray(val))
                for f in [folder/'features.npy',folder/'equity.npy',folder/'grid.npy',source/'weights.npy']:hashes[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
                records.append(dict(role=role,day=day,batch=batch,seed=seed,steps=len(dw),original_steps=len(eq),paths=dw.shape[1]));print('PREPARED',role,day,batch,'steps',len(dw),'seconds',round(time.perf_counter()-tick,2),flush=True)
    dump(out/'design.json',dict(source=str(CACHE),input_hashes=hashes,records=records,old_readout=p[:7].tolist(),initial_state=prior['initial_state'],ends_years=ends.tolist(),steps_per_trading_day=64,knot_treatment='both bracketing original fine-grid points retained; any knot delay is bounded by one512-step interval',stress_innovations='independent iid standard normals, separate from old path importance sampler',market_out_of_sample=False,wall_seconds=time.perf_counter()-start))
    print('JOINT PREPARATION COMPLETE',time.perf_counter()-start,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();prepare(args.output.resolve())
