"""Fresh target-Q streaming paired validation of the protected SPX curve fit."""
import argparse
import ctypes as ct
from dataclasses import asdict
import hashlib
import json
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr

from .full_readout_fit import HERE,arr,dump,Kernel
from .m2_spx import frozen_stage3
from .controlled_readout_fast_fit import ControlledEngine
from .vix_only_curve_fit import PRIOR,ACTIVE
from .two_speed_joint_core import physical,vix_x
from .two_speed_precision_conditional import ConditionalKernel
from .two_speed_precision import black
from .protected_curve_fit import black_parts
from esn_eq.pricing import DatedM2Pricer,dated_grid
from esn_eq.states import joint_transition_covariance,symmetric_square_root


class FreshKernel:
    def __init__(self,lib):
        self.lib=ct.CDLL(str(lib));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.lib.protected_fresh.argtypes=[ct.c_int]*5+[I]+[P]*7+[I]*4+[P]*3+[ct.c_uint64,P]
    def simulate(self,d,pairs,readout,models,seed):
        result=np.zeros((2*pairs,d['expiries'],2,4));self.lib.protected_fresh(len(d['dt']),pairs,d['projection'].shape[1],d['expiries'],len(d['decay']),d['transition'],d['decay'],d['root'],d['cross'],d['mean'],d['projection'],d['time'],d['dt'],d['seg'],d['coarse'],d['reference'],d['expiry'],d['lognorm'],arr(readout),arr(models),seed,result);return result


def context():
    prior=json.loads((PRIOR/'run.json').read_text());a,_,premium,_=frozen_stage3();assert json.loads(json.dumps(asdict(a)))==prior['architecture'] and json.loads(json.dumps(asdict(premium)))==prior['premium'];dated=DatedM2Pricer(a,premium);ctx=object.__new__(ControlledEngine);ctx.dated=dated;ctx.initial=arr(prior['initial_state']);ctx.k=Kernel(PRIOR/'convolution.dylib');ctx.ends=arr(prior['ends_years']);ctx.tau=1/max(dated.layout.rates);ctx.horizon=1060/365;ctx.pool=ThreadPoolExecutor(max_workers=4);ctx.normal_laws={};ctx.last_readout=None;ctx.normal_cached=None
    try:norm,audit=ctx.normalizer(arr(prior['selected']))
    finally:ctx.close()
    return prior,dated,norm,audit


def transition_arrays(dated,initial,times):
    keep=np.flatnonzero(np.any(dated.projection!=0,axis=0));drop=np.setdiff1d(np.arange(len(dated.layout.rates)),keep);assert len(keep)==24 and not dated.architecture.echo.enabled
    assert np.all(dated.projection[:,drop]==0) # retained OU drift is diagonal
    dt=np.diff(times);unique,inv=np.unique(dt,return_inverse=True);roots=[];cross=[];decay=[];errors=[]
    for step in unique:
        full=joint_transition_covariance(dated.layout,float(step));c=full[1:,0];legacy=symmetric_square_root(full[1:,1:]-np.outer(c,c)/step);conditional=(legacy@legacy.T)[np.ix_(keep,keep)];root=symmetric_square_root(conditional);implemented=conditional+np.outer(c[keep],c[keep])/step;actual=root@root.T+np.outer(c[keep],c[keep])/step;scale=np.sqrt(np.diag(implemented));error=float(np.max(abs(actual-implemented)/np.outer(scale,scale)));assert error<2e-9
        roots.append(root);cross.append(c[keep]/step);decay.append(np.exp(-dated.layout.rates[keep]*step));errors.append(error)
    rates=dated.layout.rates[keep];dec=np.exp(-times[:-1,None]*rates);states=dec*initial[keep]+dated.q_drift[keep]*(-np.expm1(-times[:-1,None]*rates))/rates;projection=arr(dated.projection[:,keep]);mean=arr(states@projection.T)
    return dict(transition=inv.astype('int32'),root=arr(roots),cross=arr(cross),decay=arr(decay),projection=projection,mean=mean),dict(retained_modes=keep.tolist(),marginalized_modes=drop.tolist(),unique_transition_sizes=len(unique),maximum_standardized_covariance_reconstruction_error=max(errors),marginalization_reason='Discarded bank has zero readout projections, diagonal exogenous OU drift and no return/state feedback. Retained full joint OU/equity marginal covariance is preserved.')


def build_data(prior,dated,norm,quotes):
    days=np.sort(quotes.expiry_days.unique());ends=arr(prior['ends_years']);events=np.unique(np.r_[days/365,ends]);coarse=dated_grid(events,1/(252*64));fine=arr(np.r_[0,((coarse[:-1,None]+np.diff(coarse)[:,None]*np.array([.5,1.]))).ravel()]);nt=len(fine)-1
    assert np.array_equal(fine[::2],coarse)
    d,audit=transition_arrays(dated,arr(prior['initial_state']),fine);d.update(time=arr(fine[:-1]),dt=arr(np.diff(fine)),seg=np.minimum(np.searchsorted(ends,fine[:-1],side='right'),15).astype('int32'),coarse=(np.arange(nt)%2==0).astype('int32'),expiries=len(days))
    idx=np.unique(np.r_[np.arange(0,len(coarse)-1,64),np.searchsorted(coarse,ends[ends<coarse[-1]])]);refindices=np.r_[2*idx,nt];refmask=np.zeros(nt,np.int32);refmask[refindices[:-1]]=1;expiry=np.full(nt,-1,np.int32)
    for j,day in enumerate(days):
        pos=np.flatnonzero(fine==day/365);assert len(pos)==1;expiry[pos[0]-1]=j
    d.update(reference=refmask,expiry=expiry,lognorm=arr(norm(fine[:-1])[:,0]));grid=(fine[refindices],np.diff(fine[refindices]),d['seg'][refindices[:-1]])
    return d,grid,days,audit


def terminal_residual(stocks,K,level):
    S=stocks[:,level];R=stocks[:,2];iv=stocks[:,3];delta=ndtr(-np.log(K)[None,:]/np.sqrt(iv[:,None])+.5*np.sqrt(iv[:,None]));return np.maximum(S[:,None]-K,0)-np.maximum(R[:,None]-K,0)-delta*(S-R)[:,None]


def loss(q,prices):
    dist=np.where(prices<q.price_bid_U,prices-q.price_bid_U,np.where(prices>q.price_ask_U,prices-q.price_ask_U,0))/q.price_scale_U
    return float(np.mean([np.mean(np.asarray(dist)[g.index]**2) for _,g in q.groupby('expiry_days')]))


def evaluate(run,stocks,H,days,label):
    quotes=pd.read_csv(run/'spx_membership.csv');models=np.array([json.loads((run/'design.json').read_text())['baseline'],json.loads((run/'selection.json').read_text())['parameters']]);out=run/'fresh_check';frames=[];expensive={};cheap={};_,f=physical(models[0]);N=len(stocks);pairs=N//2
    for model,name in enumerate(['baseline','selected']):
        for level,steps in enumerate([64,128]):
            pieces=[];eb=[];cb=[]
            for j,day in enumerate(days):
                q=quotes[quotes.expiry_days==day].copy();K=q.K_over_F.to_numpy();F=q.forward.to_numpy();r=terminal_residual(stocks[:,j,model,:],K,level);pair=(r[:pairs]+r[pairs:])/2;rb=pair.reshape(8,-1,len(K)).mean(1);I=f*(day/365)+(.2-f)*(H[int(day)]@models[model,:16]);B=black(I.reshape(-1),K).reshape(8,-1,len(K)).mean(1);price=(rb.mean(0)+B.mean(0))*F-q.option_type.eq('put').to_numpy()*(F-q.strike.to_numpy());se=np.sqrt(rb.var(0,ddof=1)/8+B.var(0,ddof=1)/8)*F;q['model']=name;q['steps_per_day']=steps;q['model_price_U']=price;q['mc_se_U']=se;q['interval_residual']=np.where(price<q.price_bid_U,price-q.price_bid_U,np.where(price>q.price_ask_U,price-q.price_ask_U,0))/q.price_scale_U;pieces.append(q);eb.append(rb*F);cb.append(B*F)
            frame=pd.concat(pieces,ignore_index=True);frames.append(frame);frame.to_csv(out/f'prices_{name}_{steps}.csv',index=False);expensive[(name,steps)]=np.concatenate(eb,axis=1);cheap[(name,steps)]=np.concatenate(cb,axis=1)
    paired=[];summary=[];rng=np.random.default_rng(2026102843)
    for frame in frames:
        p=frame[frame.priority_selected];summary.append(dict(model=frame.model.iloc[0],steps=int(frame.steps_per_day.iloc[0]),loss=float(np.mean([np.mean(g.interval_residual**2) for _,g in p.groupby('expiry_days')])),inside=int(p.interval_residual.eq(0).sum()),negative=int((p.model_price_U<0).sum()),median_SE_scales=float(np.median(p.mc_se_U/p.price_scale_U))))
    for steps in [64,128]:
        a=next(z for z in frames if z.model.iloc[0]=='baseline' and z.steps_per_day.iloc[0]==steps);b=next(z for z in frames if z.model.iloc[0]=='selected' and z.steps_per_day.iloc[0]==steps);de=expensive[('selected',steps)]-expensive[('baseline',steps)];dc=cheap[('selected',steps)]-cheap[('baseline',steps)];pp=a[['expiry_days','strike','option_type','priority_selected','price_scale_U']].copy();pp['candidate_minus_baseline_U']=b.model_price_U.to_numpy()-a.model_price_U.to_numpy();pp['paired_SE_U']=np.sqrt(de.var(0,ddof=1)/8+dc.var(0,ddof=1)/8);pp.to_csv(out/f'paired_models_{steps}.csv',index=False)
        # Resample independent expensive and cheap batch means separately, using
        # identical batch resamples for the two models. Fixed selected model only.
        deltas=[]
        for _ in range(1024):
            ie=rng.integers(0,8,8);ic=rng.integers(0,8,8);ps=[]
            for name in ['baseline','selected']:
                p0=(expensive[(name,steps)][ie].mean(0)+cheap[(name,steps)][ic].mean(0))-a.option_type.eq('put').to_numpy()*(a.forward-a.strike).to_numpy();ps.append(loss(a[a.priority_selected].reset_index(drop=True),p0[a.priority_selected.to_numpy()]))
            deltas.append(ps[1]-ps[0])
        paired.append(dict(steps=steps,loss_delta_selected_minus_baseline=summary[[ (z['model'],z['steps']) for z in summary].index(('selected',steps))]['loss']-summary[[ (z['model'],z['steps']) for z in summary].index(('baseline',steps))]['loss'],paired_batch_bootstrap_loss_delta_interval=np.quantile(deltas,[.025,.975]).tolist(),note='1024paired bootstrap resamples of8expensive and8cheap batches; approximate fixed-candidate sampling interval, not a convergence bound.'))
    for name in ['baseline','selected']:
        a=next(z for z in frames if z.model.iloc[0]==name and z.steps_per_day.iloc[0]==64);b=next(z for z in frames if z.model.iloc[0]==name and z.steps_per_day.iloc[0]==128);de=expensive[(name,128)]-expensive[(name,64)];pp=a[['expiry_days','strike','option_type','priority_selected','price_scale_U']].copy();pp['fine_minus_coarse_U']=b.model_price_U.to_numpy()-a.model_price_U.to_numpy();pp['paired_SE_U']=de.std(0,ddof=1)/np.sqrt(8);pp.to_csv(out/f'paired_grid_{name}.csv',index=False)
    np.savez_compressed(out/'blocks.npz',**{f'{name}_{step}_expensive':v for (name,step),v in expensive.items()},**{f'{name}_{step}_cheap':v for (name,step),v in cheap.items()});dump(out/'scores.json',dict(label=label,summaries=summary,paired=paired));return summary,paired


def check(args):
    run=args.run.resolve();out=run/'fresh_check';out.mkdir(exist_ok=False);begin=time.perf_counter();base=np.array(json.loads((run/'design.json').read_text())['baseline']);selected=np.array(json.loads((run/'selection.json').read_text())['parameters']);np.testing.assert_array_equal(vix_x(base),vix_x(selected));models=arr([base,selected]);prior,dated,norm,na=context();quotes=pd.read_csv(run/'spx_membership.csv') if (run/'spx_membership.csv').exists() else pd.read_csv(HERE/'runs/m2_two_speed_joint_001/spx_membership.csv');quotes.to_csv(run/'spx_membership.csv',index=False);d,grid,days,ta=build_data(prior,dated,norm,quotes)
    lib=out/'fresh.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'protected_curve_fresh.cpp'),'-o',str(lib)],check=True);k=FreshKernel(lib);stocks=[];seeds=[];np.savez_compressed(out/'grid_and_transitions.npz',**d,reference_grid=grid[0],reference_segments=grid[2]);dump(out/'transition_audit.json',ta)
    sources=[Path(__file__),HERE/'protected_curve_fresh.cpp',HERE/'two_speed_precision_kernel.cpp',HERE/'two_speed_precision_conditional.cpp',HERE/'two_speed_precision_conditional.py',PRIOR/'run.json',PRIOR/'convolution.dylib',run/'selection.json'];hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources};(out/'source').mkdir()
    for f in sources:shutil.copy2(f,out/'source'/f.name)
    dump(out/'design.json',dict(source_hashes=hashes,seed_base=args.seed,paths=16*args.pairs,pairs_per_batch=args.pairs,batches=8,antithetic_sampling_unit='pair',normalizer=na,normalizer_order=48,steps_per_day=[64,128],fine_steps=len(d['dt']),same_path_all7maturities=True,all_old_equity_and_stress_draws_new=True,target_law='direct Q; no importance sampling; exact retained Gaussian OU/equity transitions',VIX_map_identical=True,reference_points_per_scramble=2**args.power,reference_scrambles=8))
    for batch in range(8):
        seed=int(np.random.SeedSequence([args.seed,batch]).generate_state(1,dtype=np.uint64)[0]);tick=time.perf_counter();s=k.simulate(d,args.pairs,arr(prior['selected']),models,seed);assert np.isfinite(s).all() and np.all(s>0);stocks.append(s);seeds.append(seed);print('FRESH BATCH',batch,'paths',2*args.pairs,'seconds',round(time.perf_counter()-tick,1),flush=True)
    # Put all plus partners before all minus partners for pair-level statistics.
    stocks=np.concatenate([np.concatenate([s[:args.pairs] for s in stocks]),np.concatenate([s[args.pairs:] for s in stocks])]);np.savez_compressed(out/'stocks.npz',stock=stocks,seeds=np.array(seeds,dtype=np.uint64))
    # Same reference path, evaluated at each expiry prefix; no new equity paths.
    bridge=ConditionalKernel(HERE/'runs/m2_two_speed_precision_conditional_002/conditional.dylib');H={};hseeds={}
    for day in days:
        end=np.searchsorted(grid[0],day/365);assert grid[0][end]==day/365;g=(grid[0][:end+1],np.diff(grid[0][:end+1]),grid[2][:end]);pieces=[];hseeds[int(day)]=[]
        for batch in range(8):
            seed=int(np.random.SeedSequence([args.seed,1000,int(day),batch]).generate_state(1)[0]);pieces.append(bridge.integrals(g,base,args.power,seed));hseeds[int(day)].append(seed)
        H[int(day)]=np.array(pieces);np.savez_compressed(out/f'reference_d{day}.npz',H=H[int(day)],grid=g[0],segments=g[2],seeds=np.array(hseeds[int(day)],dtype=np.uint64));print('FRESH REFERENCE',day,flush=True)
    scores,paired=evaluate(run,stocks,H,days,'fresh512' if args.pairs==32 else 'fresh_extended');dump(out/'run.json',dict(status='fresh paired Q check completed; final accuracy gates separate',stock_seeds=seeds,reference_seeds=hseeds,summaries=scores,paired=paired,wall_seconds=time.perf_counter()-begin,market_OOS=False,VIX_parameters_changed=False))
    print('FRESH CHECK COMPLETE',json.dumps(scores),'seconds',time.perf_counter()-begin,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--pairs',type=int,default=32);ap.add_argument('--power',type=int,default=12);ap.add_argument('--seed',type=int,default=2026102847);args=ap.parse_args();check(args)
