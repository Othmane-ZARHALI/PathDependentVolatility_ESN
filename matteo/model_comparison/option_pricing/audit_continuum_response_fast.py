"""Same independent audit with an algebraically batched NumPy power oracle.

For uniform u grid with denominator K, compute B=exp(log(G)/K) once at
all original tensor-Gaussian nodes. Recurrence B^j gives all G^(j/K).
This preserves the exact finite GH48 rule and uses no interpolation in G.
It accelerates the independent refinement, without changing calibration code.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from . import audit_continuum_response as audit
from .maturity_curve import window_nodes
from .full_readout_fit import dump


def powers(law,p,grid_size):
    m,o,w=law;m=np.asarray(m).reshape(-1,3);xx=m[:,None,:]+o[None,:,:]
    eta=p[0]*xx[:,:,0]+p[1]*(xx[:,:,1]-p[2])**2+p[3]*(xx[:,:,2]-p[4])**2
    h=p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6]
    step=np.exp(2*h/(grid_size-1));current=np.ones_like(h);values=[current@w]
    for j in range(1,grid_size):
        current*=step;values.append(current@w)
    return np.array(values)[::-1]


def tables(ctx,k,p,initial,draws,folder,grid_size=65,time_order=16,order=48):
    folder.mkdir(parents=True,exist_ok=True);u=np.linspace(0,1,grid_size);nodes={'u_grid':u};checks=[]
    for day,a in draws.items():
        begin=time.perf_counter();times,tw=window_nodes(day/365,ctx.ends,ctx.tau,time_order);ratios=[]
        for it,t in enumerate(times):
            law=ctx.law(float(t-day/365),a['states'],order);nlaw=ctx.law(float(t),initial,48);v=powers(law,p,grid_size);m=powers(nlaw,p,grid_size);assert np.all(v>0) and np.all(m>0);ratios.append(np.log(v)-np.log(m))
            if it in [0,len(times)//2,len(times)-1]:
                for j in [0,grid_size//2,grid_size-1]:
                    k.lib.configure(0,float(1-u[j]));cpp=k.expectation(*law,p)[:,0];nm=k.expectation(*nlaw,p)[0,0];err=float(max(abs(cpp/v[j]-1)));me=float(abs(nm/m[j,0]-1));assert max(err,me)<1e-11;checks.append(dict(day=day,time_index=it,u=float(u[j]),relative_error=err,normalizer_error=me))
            if it%32==0:print('BATCHED INDEPENDENT GH48',day,it,'/',len(times),'seconds',round(time.perf_counter()-begin,1),flush=True)
        lr=np.array(ratios).transpose(1,0,2);np.savez_compressed(folder/f'd{day}.npz',times=times,weights=tw,log_ratios=lr,ratios=np.exp(lr[0]))
        for key,v in dict(times=times,weights=tw,log_ratios=lr,ratios=np.exp(lr[0])).items():nodes[f'd{day}_{key}']=v
    dump(folder.parent/'batched_power_scalar_checks.json',checks);return nodes


def run(folder):
    audit.tables=tables
    audit.run(folder)
    path=folder/'audit/audit.json';a=json.loads(path.read_text());a['refinement_backend']='independent NumPy tensor-GH48 with exact uniform-power recurrence';a['refinement_backend_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();a['refinement_scalar_crosschecks']='batched_power_scalar_checks.json';dump(path,a)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);run(ap.parse_args().run.resolve())
