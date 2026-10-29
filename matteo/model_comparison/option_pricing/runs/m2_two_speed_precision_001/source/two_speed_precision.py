"""Frozen two-speed SPX precision pilot with an independent-stress reference.

No optimizer is run unless a later experiment explicitly clears pricing gates.
Importance weights multiply the actual-minus-reference correction only. The
independent stress-only Black expectation is estimated on fresh cheap paths.
"""
import argparse
import ctypes as ct
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr

from .full_readout_fit import HERE, arr, dump
from .two_speed_joint_study import load_spx
from .two_speed_joint_core import StockKernel, physical, vix_x
from .vix_only_curve_fit import ACTIVE
from esn_eq.conditional_spx import lognormal_payoffs

PARENT = HERE/'runs/m2_two_speed_joint_001'
SELECTED = HERE/'runs/m2_two_speed_001'
FREE = np.array([j for j in range(16) if j not in ACTIVE])


def reference_grid(d, stride=64):
    """Daily-ish reference only; target stock remains on its original fine grid."""
    n=len(d['dt'])
    starts=np.unique(np.r_[np.arange(0,n,stride),np.flatnonzero(np.diff(d['seg']))+1])
    indices=np.r_[starts,n];mask=np.zeros(n,np.int32);mask[starts]=1
    times=np.asarray(d['grid'][indices]);return indices,mask,times,np.diff(times),np.asarray(d['seg'][starts])


class PrecisionKernel:
    def __init__(self,lib):
        self.lib=ct.CDLL(str(lib));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.lib.precision_stock.argtypes=[ct.c_int]*3+[P]*5+[I,I,P,P,P]
        self.lib.precision_integrals.argtypes=[ct.c_int,ct.c_int,P,P,I,P,ct.c_uint64,P]

    def stock(self,d,x,stride=64):
        indices,mask,t,dt,seg=reference_grid(d,stride);n=d['old'].shape[1];out=np.zeros((n,51));path=np.zeros((len(indices),n,3))
        self.lib.precision_stock(len(d['dt']),n,len(indices),d['old'],d['equity'],d['eps'],d['time'],d['dt'],d['seg'],mask,arr(x),out,path)
        return out,path,(t,dt,seg)

    def integrals(self,grid,x,pairs,seed):
        t,dt,seg=grid;out=np.zeros((pairs,32));self.lib.precision_integrals(len(dt),pairs,arr(t[:-1]),arr(dt),np.ascontiguousarray(seg,dtype=np.int32),arr(x),seed,out)
        return out.reshape(pairs,2,16)


def black(iv,k):
    return lognormal_payoffs(-.5*iv[:,None],iv[:,None],k[None,:],np.full((1,len(k)),'call'))


def residuals(out,path,k,weights):
    S,R,iv=out[:,:3].T;sd=np.sqrt(iv[:,None]);delta=ndtr(-np.log(k)[None,:]/sd+.5*sd)
    payoff=np.maximum(S[:,None]-k,0)-np.maximum(R[:,None]-k,0)
    terminal=payoff-delta*(S-R)[:,None]
    # Full future stress is independent of all old equity/OU drivers. Conditional
    # on it, every stock increment has mean zero under Q. These are predictable
    # daily-ish hedges in that enlarged filtration, with no continuous-time claim.
    hedge=np.zeros_like(terminal)
    for j in range(len(path)-1):
        remaining=iv-path[j,:,2]
        if np.any(remaining<=0):raise FloatingPointError('Nonpositive remaining reference variance')
        sd=np.sqrt(remaining)[:,None];Sa,Ra=path[j,:,:2].T
        da=ndtr(np.log(Sa[:,None]/k[None,:])/sd+.5*sd)
        dr=ndtr(np.log(Ra[:,None]/k[None,:])/sd+.5*sd)
        dSa,dRa=(path[j+1,:,:2]-path[j,:,:2]).T
        hedge+=da*dSa[:,None]-dr*dRa[:,None]
    return dict(unhedged=weights[:,None]*payoff,terminal=weights[:,None]*terminal,
                dynamic=weights[:,None]*(payoff-hedge))


def summarize(q,prices,se,resid_se,cheap_se,method,role):
    frame=q.copy();frame['model_price_U']=prices;frame['mc_se_U']=se;frame['residual_SE_U']=resid_se;frame['cheap_SE_U']=cheap_se
    frame['interval_residual']=np.where(prices<q.price_bid_U,prices-q.price_bid_U,np.where(prices>q.price_ask_U,prices-q.price_ask_U,0))/q.price_scale_U
    frame['SE_gate_pass']=2*se<=q.price_scale_U/4
    rows=[]
    for day,g in frame.groupby('expiry_days',sort=True):
        h=g[g.priority_selected]
        rows.append(dict(role=role,method=method,day=int(day),quotes=len(h),loss=float(np.mean(h.interval_residual**2)),inside=int(h.interval_residual.eq(0).sum()),negative=int((h.model_price_U<0).sum()),SE_gate_passes=int(h.SE_gate_pass.sum()),median_SE_scales=float(np.median(h.mc_se_U/h.price_scale_U)),median_residual_SE_scales=float(np.median(h.residual_SE_U/h.price_scale_U)),median_cheap_SE_scales=float(np.median(h.cheap_SE_U/h.price_scale_U))))
    return frame,rows


def kernel_checks(k,legacy,data,x,out):
    d0=data[119][0];ids=np.array([0,7,19,47]);d={name:arr(d0[name][:,ids]) for name in ['old','equity','eps']};d.update({name:np.ascontiguousarray(d0[name]) for name in ['time','dt','seg','grid']})
    actual,path,grid=k.stock(d,x);old=legacy.stock(d,x,False)[:,0];np.testing.assert_allclose(actual[:,0],old,rtol=2e-12)
    # Scalar Python reference of the held-stress control and its integrated variance.
    hf=np.exp(np.log(365)*x[21]);hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf]);b=5*x[[16,20]];levels,floor=physical(x)
    z=np.zeros((4,2));logR=np.zeros(4);iv=np.zeros(4);indices,mask,*_=reference_grid(d)
    qr=np.ones(4)
    for t,dt in enumerate(d['dt']):
        if mask[t]:
            A=np.exp(z*b-.5*b*b*(-np.expm1(-2*rates*d['time'][t])));qr=(1-x[17])+x[17]*((1-x[22])*A[:,0]+x[22]*A[:,1])
        v=floor+(levels[d['seg'][t]]-floor)*qr;logR+=-.5*v*dt+np.sqrt(v)*d['equity'][t];iv+=v*dt
        z=np.exp(-rates*dt)*z+np.sqrt(-np.expm1(-2*rates*dt))*d['eps'][t]
    np.testing.assert_allclose(actual[:,1],np.exp(logR),rtol=2e-12);np.testing.assert_allclose(actual[:,2],iv,rtol=2e-12)
    errors=[]
    for j in range(16):
        dx=np.zeros(23);dx[j]=1e-6;a=k.stock(d,x+dx)[0];c=k.stock(d,x-dx)[0];fd=(a[:,:3]-c[:,:3])/2e-6;analytic=actual[:,[3+j,19+j,35+j]];error=float(np.max(abs(fd-analytic)/(1+abs(fd))));assert error<3e-6,(j,error);errors.append(error)
    # If old response is identically one the target equals the control pathwise.
    unit=d.copy();unit['old']=np.ones_like(d['old']);a,p,_=k.stock(unit,x)
    zero=residuals(a,p,np.array([.8,1.,1.2]),np.array([.4,1.1,.8,2.]))
    assert max(np.max(abs(v)) for v in zero.values())==0
    # A probe of every free coordinate leaves the exact VIX parameter map intact.
    probe=x.copy();probe[FREE]=np.linspace(.1,.9,len(FREE));np.testing.assert_array_equal(vix_x(probe),vix_x(x))
    audit=dict(target_replay_max_error=float(np.max(abs(actual[:,0]-old))),reference_replay_max_error=float(np.max(abs(actual[:,1]-np.exp(logR)))),integral_replay_max_error=float(np.max(abs(actual[:,2]-iv))),curve_derivative_errors=errors,identical_model_residuals_exact_zero=True,VIX_parameter_map_unchanged=True,free_indices=FREE.tolist())
    dump(out/'kernel_checks.json',audit)


def run(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    sources=[Path(__file__),HERE/'two_speed_precision_kernel.cpp',HERE/'TWO_SPEED_PRECISION_PLAN.md',HERE/'two_speed_joint_core.py',HERE/'two_speed_joint_study.py',HERE/'parity_adaptive.py',PARENT/'design.json',PARENT/'spx_membership.csv',SELECTED/'selection.json',SELECTED/'experimental_model_two_speed.json',HERE/'active_scope.json']
    hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources};(out/'source').mkdir()
    for f in sources:shutil.copy2(f,out/'source'/f.name)
    x=np.array(json.loads((PARENT/'design.json').read_text())['baseline']);q=pd.read_csv(PARENT/'spx_membership.csv');assert len(q)==1751 and q.priority_selected.sum()==1294
    q.to_csv(out/'contracts.csv',index=False)
    lib=out/'precision.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'two_speed_precision_kernel.cpp'),'-o',str(lib)],check=True)
    k=PrecisionKernel(lib);legacy=StockKernel(PARENT/'stock.dylib');rows=[];archive={};inputs={};costs=[]
    dump(out/'design.json',dict(baseline=x.tolist(),source_hashes=hashes,cheap_pairs_per_batch=args.pairs,batches=8,cheap_seed_base=args.seed,target_stock_steps_per_trading_day=64,control_reference_stride=64,model_changed=False,optimization_enabled=False,free_curve_indices=FREE.tolist(),protection='All stress parameters, floor, original readout and seven VIX-sensitive xi levels frozen. No refit until precision gate passes.',precision_gate='2SE<=price_scale/4 for every one of1294 primary SPX quotes, then stock timestep gate separately.',estimator='Fresh stress-only Black mean plus likelihood-weighted actual-minus-reference residual; optional terminal or daily dynamic difference hedges.',uncertainty='8 old-path blocks for residual plus8 independently seeded cheap antithetic-pair batches; independent sums of variance estimates; excludes normalizer, timestep and fitting error.'))
    for role in ['search','check']:
        data,_=load_spx(role);folder=out/role;folder.mkdir()
        if role=='search':kernel_checks(k,legacy,data,x,out)
        frames={name:[] for name in ['unhedged','terminal','dynamic']}
        for day,g in q.groupby('expiry_days',sort=True):
            tick=time.perf_counter();outputs=[];paths=[];weights=[];grid=None
            for d in data[day]:
                a,p,gr=k.stock(d,x);outputs.append(a);paths.append(p);weights.append(d['weights'])
                if grid is not None:
                    for u,v in zip(gr,grid):np.testing.assert_array_equal(u,v)
                grid=gr
                for name in ['old','equity','eps','grid','weights']:
                    f=d['folder']/f'{name}.npy';inputs[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
            a=np.concatenate(outputs);p=np.concatenate(paths,axis=1);w=np.concatenate(weights);K=g.K_over_F.to_numpy();F=g.forward.to_numpy();rs=residuals(a,p,K,w)
            cheap=[];Hs=[];seeds=[]
            for batch in range(8):
                seed=int(np.random.SeedSequence([args.seed,int(role=='check'),int(day),batch]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed)
                H=k.integrals(grid,x,args.pairs,seed);_,floor=physical(x);iv=floor*grid[0][-1]+(.2-floor)*(H@x[:16]);B=black(iv.reshape(-1),K).reshape(args.pairs,2,len(K)).mean(1);cheap.append(B.mean(0));Hs.append(H)
            cheap=np.array(cheap);cheap_mean=cheap.mean(0);cheap_se=cheap.std(0,ddof=1)/np.sqrt(8)
            # Persist small curve-linear stress integrals for future protected fits.
            np.savez_compressed(folder/f'd{day}_reference.npz',H=np.concatenate(Hs),grid=grid[0],segments=grid[2],seeds=np.array(seeds,dtype=np.uint64))
            np.savez_compressed(folder/f'd{day}_target.npz',out=a,weights=w)
            for name,r in rs.items():
                rb=r.reshape(8,-1,len(K)).mean(1);resid_se=rb.std(0,ddof=1)/np.sqrt(8);prices=(cheap_mean+r.mean(0))*F
                prices-=g.option_type.eq('put').to_numpy()*(F-g.strike.to_numpy())
                se=np.sqrt(resid_se**2+cheap_se**2)*F;frame,record=summarize(g,prices,se,resid_se*F,cheap_se*F,name,role);frames[name].append(frame);rows.extend(record)
                archive[f'{role}_d{day}_{name}_residual']=rb;archive[f'{role}_d{day}_{name}_cheap']=cheap
            costs.append(dict(role=role,day=int(day),seconds=time.perf_counter()-tick));print('PRECISION',role,day,'seconds',round(costs[-1]['seconds'],1),'dynamic medianSEscales',rows[-1]['median_SE_scales'],'passes',rows[-1]['SE_gate_passes'],flush=True)
        for name,fs in frames.items():pd.concat(fs).to_csv(folder/f'prices_{name}.csv',index=False)
        pd.DataFrame(rows).to_csv(out/'summary.csv',index=False)
    np.savez_compressed(out/'blocks.npz',**archive);pd.DataFrame(costs).to_csv(out/'costs.csv',index=False)
    for f,h in hashes.items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h
    dump(out/'input_hashes.json',inputs)
    check=[r for r in rows if r['role']=='check' and r['method']=='dynamic'];passed=sum(r['SE_gate_passes'] for r in check)
    dump(out/'run.json',dict(status='frozen-model precision pilot completed',primary_sampling_passes=passed,primary_quotes=1294,SPX_calibration_allowed=passed==1294,timestep_gate_cleared=False,model_changed=False,VIX_parameters_changed=False,baseline_replaced=False,wall_seconds=time.perf_counter()-begin))
    manifest={str(f.relative_to(out)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(out.rglob('*')) if f.is_file()};dump(out/'artifact_manifest.json',manifest)
    print('COMPLETE',passed,'/1294 sampling passes; seconds',time.perf_counter()-begin,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--pairs',type=int,default=4096);ap.add_argument('--seed',type=int,default=2026102811);args=ap.parse_args()
    if args.pairs<64:ap.error('At least64 antithetic pairs per batch')
    run(args)
