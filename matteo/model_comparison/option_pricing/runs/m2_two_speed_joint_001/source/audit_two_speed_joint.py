"""Source, cache and independently reconstructed coupled-stock checks."""
import argparse,json,hashlib,ctypes as ct,subprocess
from pathlib import Path
import numpy as np
from .full_readout_fit import HERE,dump,arr
from .vix_only_curve_fit import PRIOR
from .two_speed_joint_core import physical
from .two_speed_joint_study import CACHE


def audit(run):
    d=json.loads((run/'design.json').read_text());checked=0
    for name,sha in d['source_hashes'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha,name;checked+=1
    cd=json.loads((CACHE/'design.json').read_text())
    for name,sha in cd['input_hashes'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha,name;checked+=1
    paths=[CACHE/rec['role']/f"d{rec['day']}"/f"batch{rec['batch']}" for rec in cd['records']]
    for rec,folder in zip(cd['records'],paths):
        eps=np.load(folder/'eps.npy',mmap_mode='r');expected=np.random.default_rng(rec['seed']).standard_normal(eps.shape);assert np.array_equal(eps,expected),str(folder)
    prior=json.loads((PRIOR/'run.json').read_text());p=arr(prior['selected']);assert p[4]>0
    check=json.loads((run/'stock_kernel_checks.json').read_text());assert len(check['derivatives'])==23 and max(a['error'] for a in check['derivatives'])<2e-5
    libfile=run/'coupled_audit.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'two_speed_joint_refine.cpp'),'-o',str(libfile)],check=True);lib=ct.CDLL(str(libfile));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS');lib.coupled_stock.argtypes=[ct.c_int,ct.c_int,P,P,P,P,I,I,I,P,P,P,ct.c_uint64,P]
    rng=np.random.default_rng(201);nt=17;n=5;dt=arr(rng.uniform(.0005,.003,nt));time=arr(np.r_[0,np.cumsum(dt)[:-1]]);features=arr(rng.normal(size=(nt,n,3)));dw=arr(rng.normal(size=(nt,n))*np.sqrt(dt)[:,None]);seg=np.minimum(np.arange(nt),15).astype('int32');m64=(np.arange(nt)%4==0).astype('int32');m128=(np.arange(nt)%2==0).astype('int32');lg=arr(np.full(nt,1.5));x=np.full(23,.4);x[:16]=.2;seed=2026102849;actual=np.zeros((3,n));lib.coupled_stock(nt,n,features,dw,time,dt,seg,m64,m128,lg,p,arr(x),seed,actual)
    mask=(1<<64)-1
    def mix(v):
        v&=mask;v=((v^(v>>30))*0xbf58476d1ce4e5b9)&mask;v=((v^(v>>27))*0x94d049bb133111eb)&mask;return v^(v>>31)
    def normal(i):
        a=mix(seed+2*i*0x9e3779b97f4a7c15);b=mix(seed+(2*i+1)*0x9e3779b97f4a7c15);u=((a>>11)+.5)*2**-53;v=((b>>11)+.5)*2**-53;return np.sqrt(-2*np.log(u))*np.cos(2*np.pi*v)
    hf=365**x[21];hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf]);b=5*x[[16,20]];z=np.zeros((n,2));logs=np.zeros((3,n));held=np.zeros((2,n));levels,floor=physical(x)
    for t in range(nt):
        eta=p[0]*features[t,:,0]+p[1]*(features[t,:,1]-p[2])**2+p[3]*(features[t,:,2]-p[4])**2;old=np.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6])-lg[t]);q=-np.expm1(-2*rates*time[t]);A=np.exp(z*b-.5*b*b*q);Q=(1-x[17])*old+x[17]*((1-x[22])*A[:,0]+x[22]*A[:,1]);v=floor+(levels[seg[t]]-floor)*Q
        if m64[t]:held[0]=v
        if m128[t]:held[1]=v
        vv=np.r_[held,v[None,:]];logs+=-.5*vv*dt[t]+np.sqrt(vv)*dw[t];eps=np.array([[normal(2*(t*n+i)+j) for j in range(2)] for i in range(n)]);z=z*np.exp(-rates*dt[t])+np.sqrt(-np.expm1(-2*rates*dt[t]))*eps
    oracle=np.exp(logs);err=float(np.max(abs(actual-oracle)));np.testing.assert_allclose(actual,oracle,rtol=2e-13,atol=2e-13)
    quality=[]
    r=json.loads((run/'run.json').read_text())
    for sc in r['checking_scores']:
        ss=sc['instruments'][0];quality.append(dict(candidate=sc['name'],SPX_loss=sc['SPX_loss'],VIX_loss=sc['VIX_loss'],SPX_SE_gate_passes=sum(q['SE_gate_passes'] for q in ss['per_maturity']),SPX_quotes=ss['quotes'],maximum_VIX_forward_gap=float(max(abs(np.array(sc['forward_errors']))))))
    result=dict(status='implementation_checks_PASS; fit_quality_separate',verified_source_and_input_files=checked,all56_stress_seed_replays=True,stock_derivative_directions=23,numpy_primary_stock_replay=check['numpy_replay_max_error'],independent_coupled_stock_replay=err,positive_spike_shift=float(p[4]),quality=quality,limitations=['Existing old SPX/VIX streams reused; not market OOS.','Sampled replays do not establish uniform convergence.','Finite bounded fit with original M2 readout fixed; no structural infeasibility theorem.','SPX sampling, coupled timestep and VIX forward accuracy are separate gates.'])
    dump(run/'audit.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);args=ap.parse_args();audit(args.run.resolve())
