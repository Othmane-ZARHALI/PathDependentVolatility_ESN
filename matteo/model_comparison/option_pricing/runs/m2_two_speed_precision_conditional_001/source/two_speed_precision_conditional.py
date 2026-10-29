"""One-step conditional correction and scrambled OU-bridge reference pilot."""
import argparse
import ctypes as ct
import hashlib
import heapq
import json
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtri
from scipy.stats import qmc

from .full_readout_fit import HERE,arr,dump
from .two_speed_precision import reference_grid,black,summarize,PARENT,SELECTED
from .two_speed_joint_study import load_spx
from .two_speed_joint_core import physical


def bridge(grid,x):
    times=grid[0][:-1];nt=len(times);hf=365**x[21];hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf])
    nodes=[nt-1];left=[0];right=[0];coeff=[np.ravel(np.c_[np.zeros(2),np.zeros(2),np.sqrt(-np.expm1(-2*rates*times[-1]))])]
    heap=[(-(times[-1]-times[0]),0,nt-1)]
    while heap:
        _,l,r=heapq.heappop(heap)
        if r-l<2:continue
        mid=(times[l]+times[r])/2;j=int(np.clip(np.searchsorted(times,mid),l+1,r-1))
        if j>l+1 and abs(times[j-1]-mid)<abs(times[j]-mid):j-=1
        a=times[j]-times[l];b=times[r]-times[j];qa=-np.expm1(-2*rates*a);qt=-np.expm1(-2*rates*(a+b));c=np.exp(-rates*b)*qa/qt
        variance=qa-c*c*qt
        if np.any(variance<=0):raise FloatingPointError('Invalid OU bridge variance')
        nodes.append(j);left.append(l);right.append(r);coeff.append(np.ravel(np.c_[np.exp(-rates*a)-c*np.exp(-rates*(a+b)),c,np.sqrt(variance)]))
        heapq.heappush(heap,(-(times[j]-times[l]),l,j));heapq.heappush(heap,(-(times[r]-times[j]),j,r))
    assert sorted(nodes)==list(range(1,nt))
    return np.array(nodes,np.int32),np.array(left,np.int32),np.array(right,np.int32),arr(coeff)


class ConditionalKernel:
    def __init__(self,path):
        self.lib=ct.CDLL(str(path));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.lib.precision_conditional.argtypes=[ct.c_int]*3+[P]*5+[I,I]+[P]*4
        self.lib.precision_bridge_integrals.argtypes=[ct.c_int]*2+[P]*3+[I]*4+[P]*3

    def correction(self,d,x,iv,K):
        _,mask,*_=reference_grid(d);out=np.zeros((len(iv),len(K)+1))
        self.lib.precision_conditional(len(d['dt']),len(iv),len(K),d['old'],d['equity'],d['eps'],d['time'],d['dt'],d['seg'],mask,arr(x),arr(iv),arr(K),out)
        return out

    def integrals(self,grid,x,m,seed):
        nodes,left,right,c=bridge(grid,x);engine=qmc.Sobol(d=2*len(nodes),scramble=True,seed=seed);u=engine.random_base2(m)
        if np.any((u<=0)|(u>=1)):raise FloatingPointError('Sobol normal endpoint')
        normals=ndtri(u);out=np.zeros((len(u),16));t,dt,seg=grid
        self.lib.precision_bridge_integrals(len(dt),len(u),normals,arr(t[:-1]),arr(dt),np.ascontiguousarray(seg,dtype=np.int32),nodes,left,right,c,arr(x),out)
        return out


def run(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();parent=args.parent.resolve();x=np.array(json.loads((parent/'design.json').read_text())['baseline']);contracts=pd.read_csv(parent/'contracts.csv');sources=[Path(__file__),HERE/'two_speed_precision_conditional.cpp',HERE/'two_speed_precision_kernel.cpp',HERE/'two_speed_precision.py',HERE/'TWO_SPEED_PRECISION_PLAN.md',parent/'design.json',parent/'run.json',parent/'artifact_manifest.json',SELECTED/'selection.json',SELECTED/'experimental_model_two_speed.json']
    hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources};(out/'source').mkdir()
    for f in sources:shutil.copy2(f,out/'source'/f.name)
    lib=out/'conditional.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'two_speed_precision_conditional.cpp'),'-o',str(lib)],check=True);k=ConditionalKernel(lib)
    dump(out/'design.json',dict(parent=str(parent),source_hashes=hashes,baseline=x.tolist(),reference_points_per_scramble=2**args.power,reference_scrambles=8,reference_seed_base=args.seed,representatives_only=not args.all_quotes,estimator='Exact discrete one-step Black-value martingale correction plus independent scrambled Sobol OU-bridge stress reference.',model_changed=False,calibration_enabled=False,expensive_paths_per_expiry=512,checking_streams_reused=True))
    blocks={};rows=[];costs=[]
    for role in ['search','check']:
        folder=out/role;folder.mkdir();data,_=load_spx(role);frames=[]
        for day,g0 in contracts.groupby('expiry_days',sort=True):
            tick=time.perf_counter();g0=g0.reset_index(drop=True)
            if args.all_quotes:ids=np.arange(len(g0))
            else:
                h=g0[g0.priority_selected];ids=np.unique([int(h.iloc[np.argmin(abs(np.log(h.K_over_F.to_numpy())-z*h.atm_proxy.to_numpy()*np.sqrt(h.T_years.to_numpy())))].name) for z in [-2,-1,0,1,2]])
            g=g0.iloc[ids].copy();K=g.K_over_F.to_numpy();F=g.forward.to_numpy();target=dict(np.load(parent/role/f'd{day}_target.npz'));residual=[]
            for batch,d in enumerate(data[day]):
                ii=slice(batch*128,(batch+1)*128);a=k.correction(d,x,target['out'][ii,2],K);np.testing.assert_allclose(a[:,-1],target['out'][ii,0],rtol=2e-12);residual.append(a[:,:-1]*d['weights'][:,None]);_,_,t,dt,seg=reference_grid(d);grid=(t,dt,seg)
            r=np.concatenate(residual);rb=r.reshape(8,-1,len(K)).mean(1);rs=rb.std(0,ddof=1)/np.sqrt(8);cheap=[];Hs=[];seeds=[]
            for batch in range(8):
                seed=int(np.random.SeedSequence([args.seed,int(role=='check'),int(day),batch]).generate_state(1)[0]);seeds.append(seed);H=k.integrals(grid,x,args.power,seed);_,f=physical(x);iv=f*grid[0][-1]+(.2-f)*(H@x[:16]);cheap.append(black(iv,K).mean(0));Hs.append(H)
            cheap=np.array(cheap);cs=cheap.std(0,ddof=1)/np.sqrt(8);prices=(cheap.mean(0)+r.mean(0))*F-g.option_type.eq('put').to_numpy()*(F-g.strike.to_numpy());se=np.sqrt(rs**2+cs**2)*F
            frame,record=summarize(g,prices,se,rs*F,cs*F,'conditional',role);frames.append(frame);rows.extend(record);blocks[f'{role}_d{day}_residual']=rb;blocks[f'{role}_d{day}_reference']=cheap
            np.savez_compressed(folder/f'd{day}_reference.npz',H=np.stack(Hs),grid=grid[0],segments=grid[2],seeds=np.array(seeds,dtype=np.uint64));np.savez_compressed(folder/f'd{day}_correction.npz',residual=r,quote_indices=ids)
            costs.append(dict(role=role,day=int(day),seconds=time.perf_counter()-tick));pd.concat(frames).to_csv(folder/'prices_conditional.csv',index=False);pd.DataFrame(rows).to_csv(out/'summary.csv',index=False)
            print('CONDITIONAL',role,day,'seconds',round(costs[-1]['seconds'],1),'medianSE',record[0]['median_SE_scales'],'resid',record[0]['median_residual_SE_scales'],'cheap',record[0]['median_cheap_SE_scales'],'passes',record[0]['SE_gate_passes'],'/',len(g),flush=True)
    np.savez_compressed(out/'blocks.npz',**blocks);pd.DataFrame(costs).to_csv(out/'costs.csv',index=False)
    for f,h in hashes.items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h
    selected=[r for r in rows if r['role']=='check'];passes=sum(r['SE_gate_passes'] for r in selected);n=sum(r['quotes'] for r in selected)
    dump(out/'run.json',dict(status='conditional precision pilot completed',sampling_passes=passes,primary_quotes_evaluated=n,all_primary_quotes_evaluated=bool(args.all_quotes),SPX_calibration_allowed=False,timestep_gate_cleared=False,model_changed=False,VIX_parameters_changed=False,wall_seconds=time.perf_counter()-begin))
    dump(out/'artifact_manifest.json',{str(f.relative_to(out)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(out.rglob('*')) if f.is_file()})


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--parent',type=Path,default=HERE/'runs/m2_two_speed_precision_002');ap.add_argument('--output',type=Path,required=True);ap.add_argument('--power',type=int,default=12);ap.add_argument('--seed',type=int,default=2026102817);ap.add_argument('--all-quotes',action='store_true');args=ap.parse_args();run(args)
