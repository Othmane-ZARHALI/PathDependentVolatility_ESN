"""Coupled stock timestep check and complete eligible-SPX diagnostics."""
import argparse,json,ctypes as ct,time,subprocess,hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from .full_readout_fit import HERE,dump,arr,Kernel
from .controlled_readout_fast_fit import ControlledEngine
from .prepare_two_speed_joint import CACHE,DAYS,indices
from .vix_only_curve_fit import PRIOR
from .m2_spx import frozen_stage3
from .parity_adaptive import reference_calls
from .plot_joint_curve_fit import implied_vol
from esn_eq.pricing import DatedM2Pricer


def check(run):
    out=run/'timestep_check';out.mkdir(exist_ok=False);start=time.perf_counter();selection=json.loads((run/'selection.json').read_text());design=json.loads((run/'design.json').read_text());prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);baseline=np.array(design['baseline']);chosen=np.array(selection['parameters']);models=[('vix_only',baseline),('selected',chosen)]
    a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);ctx=object.__new__(ControlledEngine);ctx.dated=dated;ctx.initial=np.array(prior['initial_state']);ctx.k=Kernel(PRIOR/'convolution.dylib');ctx.ends=ends;ctx.tau=1/max(dated.layout.rates);ctx.horizon=1060/365;ctx.pool=ThreadPoolExecutor(max_workers=4);ctx.normal_laws={};ctx.last_readout=None;ctx.normal_cached=None
    try:norm,na=ctx.normalizer(p)
    finally:ctx.close()
    libpath=out/'coupled.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'two_speed_joint_refine.cpp'),'-o',str(libpath)],check=True);lib=ct.CDLL(str(libpath));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS');lib.coupled_stock.argtypes=[ct.c_int,ct.c_int,P,P,P,P,I,I,I,P,P,P,ct.c_uint64,P]
    quotes=pd.read_csv(run/'spx_membership.csv');reflevels=np.array(json.loads((PRIOR/'design.json').read_text())['reference_levels']);summary=[];seedrecords=[];momentrows=[]
    rows_by_model={name:[] for name,_ in models};paired_by_model={name:[] for name,_ in models}
    for day in DAYS:
        q=quotes[quotes.expiry_days.eq(day)].copy();g=[];weights=[];refstocks=[];refvariance=None
        for batch,folder in enumerate(sorted((CACHE/'check'/f'd{day}').glob('batch*'))):
            tick=time.perf_counter();grid=np.load(folder/'grid.npy');times=arr(grid[:-1]);dt=arr(np.diff(grid));feat=np.load(folder/'features.npy',mmap_mode='r');dw=np.load(folder/'equity.npy',mmap_mode='r');seg=np.minimum(np.searchsorted(ends,times,side='right'),15).astype('int32');mask64=np.zeros(len(dt),dtype='int32');mask128=mask64.copy();mask64[indices(grid,8)[:-1]]=1;mask128[indices(grid,4)[:-1]]=1;lg=arr(norm(times)[:,0]);seed=int(np.random.SeedSequence([2026102809,day,batch]).generate_state(1,dtype=np.uint64)[0]);seedrecords.append(dict(day=day,batch=batch,seed=seed))
            terms=[]
            for name,x in models:
                outS=np.zeros((3,dw.shape[1]));lib.coupled_stock(len(dt),dw.shape[1],feat,dw,times,dt,seg,mask64,mask128,lg,arr(p),arr(x),seed,outS);terms.append(outS)
            g.append(np.array(terms));ww=np.load(CACHE/'check'/f'd{day}'/'weights.npy')[batch*128:(batch+1)*128];weights.append(ww)
            rs=np.minimum(np.searchsorted([.25,.75],times,side='right'),2);brown=np.stack([np.asarray(dw[rs==j]).sum(0) for j in range(3)],axis=1);dur=np.array([dt[rs==j].sum() for j in range(3)]);refvariance=float(dur@reflevels);refstocks.append(np.exp(brown@np.sqrt(reflevels)-.5*refvariance));print('TIMESTEP',day,batch,'seconds',round(time.perf_counter()-tick,1),flush=True)
        stocks=np.concatenate(g,axis=2);w=np.concatenate(weights);ref=np.concatenate(refstocks)
        for j,(name,x) in enumerate(models):
            values=[]
            for level,steps in enumerate([64,128,512]):
                S=stocks[j,level];y,_=reference_calls(S,np.zeros((len(S),1)),ref,np.zeros((len(S),1)),w,q.K_over_F.to_numpy(),refvariance,np.zeros(1));y=y*q.forward.to_numpy();y-=q.option_type.eq('put').to_numpy()[None,:]*(q.forward-q.strike).to_numpy()[None,:];values.append(y);b=y.reshape(8,-1,len(q)).mean(1)
                qq=q.copy();qq['steps_per_day']=steps;qq['model_price_U']=b.mean(0);qq['mc_se_U']=b.std(0,ddof=1)/np.sqrt(8);qq['interval_residual']=np.where(qq.model_price_U<qq.price_bid_U,qq.model_price_U-qq.price_bid_U,np.where(qq.model_price_U>qq.price_ask_U,qq.model_price_U-qq.price_ask_U,0))/qq.price_scale_U;qq['model_native_IV']=[implied_vol(a.model_price_U,a.forward,a.strike,a.T_years,a.option_type) for a in qq.itertuples()];rows_by_model[name].append(qq)
                momentrows.append(dict(model=name,day=day,steps=steps,mass=float(w.mean()),stock_mean=float(np.mean(w*S)),stock_mean_se=float(np.std(w*S,ddof=1)/np.sqrt(len(S)))))
            for lo,hi in [(0,1),(1,2),(0,2)]:
                dy=values[hi]-values[lo];b=dy.reshape(8,-1,len(q)).mean(1);pp=q[['expiry_days','strike','option_type','priority_selected','price_scale_U']].copy();pp['coarse_steps']=[64,128,512][lo];pp['fine_steps']=[64,128,512][hi];pp['paired_change_U']=b.mean(0);pp['paired_SE_U']=b.std(0,ddof=1)/np.sqrt(8);paired_by_model[name].append(pp)
        np.savez_compressed(out/f'stock_d{day}.npz',stock=stocks,weights=w,reference=ref)
    for name,_ in models:
        table=pd.concat(rows_by_model[name],ignore_index=True);table.to_csv(out/f'prices_all_{name}.csv',index=False);paired=pd.concat(paired_by_model[name],ignore_index=True);paired.to_csv(out/f'paired_timestep_{name}.csv',index=False)
        for steps,z in table.groupby('steps_per_day'):
            primary=z[z.priority_selected];loss=float(np.mean([np.mean(g.interval_residual**2) for _,g in primary.groupby('expiry_days')]));summary.append(dict(model=name,steps=int(steps),SPX_primary_loss=loss,inside=int(primary.interval_residual.eq(0).sum()),quotes=len(primary),negative=int((primary.model_price_U<0).sum()),median_SE_scales=float(np.median(primary.mc_se_U/primary.price_scale_U))))
    pd.DataFrame(momentrows).to_csv(out/'stock_moments.csv',index=False);dump(out/'run.json',dict(status='completed coupled stock timestep diagnostic; accuracy gates assessed separately',normalizer=na,seeds=seedrecords,models=[dict(name=name,parameters=x.tolist()) for name,x in models],summaries=summary,source_hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),HERE/'two_speed_joint_refine.cpp']},wall_seconds=time.perf_counter()-start,limitations=['Same512old checking paths; independent new stress innovations coupled across grids and candidates.','All1751eligible SPX quotes retained; cutoff3defines1294primary quotes.','No market OOS, uniform timestep bound, or sample precision certificate.']))
    print('TIMESTEP CHECK COMPLETE',time.perf_counter()-start,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);args=ap.parse_args();check(args.run.resolve())
