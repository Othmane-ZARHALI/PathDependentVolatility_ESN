"""Independent block-level reconstruction and exact large-base redraw."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from .vix_only_curve_fit import PRIOR,ARCHIVE,DAYS,DELTA,load_draws
from .full_readout_fit import dump,HERE
from .m2_spx import ROOT,frozen_stage3
from .forward_floor_core import unpack
from .postcap_vov_core import conditional_multiplier
from .parity_adaptive import gaussian_controls
from .mixed_response_study import seed_for,ZSHIFTS,ZPROB
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture

INITIAL=HERE/'runs/m2_continuum_response_002'


def boundary_V(a,nodes,day,x,p,ends):
    levels,floor,theta=unpack(x,p);b,m,lh=theta;T=day/365;rate=np.log(2)*365/np.exp(lh);times=nodes[f'd{day}_times'];tw=nodes[f'd{day}_weights'];old=nodes[f'd{day}_ratios'];z=np.sqrt(-np.expm1(-2*rate*T))*a['z_normals'];integ=np.full(len(z),floor*DELTA)
    for it,(t,wt) in enumerate(zip(times,tw)):
        seg=min(np.searchsorted(ends,t,side='right'),len(ends)-1);L=conditional_multiplier(z,T,t-T,b,rate);integ+=wt*(levels[seg]-floor)*((1-m)*old[it]+m*L)
    return 100*np.sqrt(integ/DELTA)


def run(run):
    out=run/'multilevel_audit';out.mkdir(exist_ok=False);prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);sel=json.loads((run/'selection.json').read_text());big=load_draws(run/'fresh_base_cache');small=load_draws(INITIAL/'fresh_cache');bn=dict(np.load(run/'fresh_base_nodes.npz'));sn={}
    for day in DAYS:
        with np.load(INITIAL/'fresh_power_tables'/f'd{day}.npz') as aa:
            for k,v in aa.items():
                if k!='log_ratios':sn[f'd{day}_{k}']=v
    market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv');controls=json.loads((PRIOR/'vix_controls.json').read_text());rows=[];moments=[]
    models=[('continuous',sel['parameters'])]+([('interior',sel['secondary']['parameters'])] if sel['secondary'] else [])
    for name,xx in models:
        x=np.array(xx);paths=np.load(run/'fresh_check'/f'paths_{name}.npz');npblock=[];cpblock=[];fpblock=[];e2block=[]
        for day in DAYS:
            q=market[market.expiry_days.eq(day)];K=q.strike.to_numpy();put=q.option_type.eq('put').to_numpy();components=[];second=[]
            for a,V in [(big[day],boundary_V(big[day],bn,day,x,p,ends)),(small[day],paths[f'd{day}_V']),(small[day],boundary_V(small[day],sn,day,x,p,ends))]:
                cc=gaussian_controls(a['normals'],a['weights']);fixed=controls[str(day)];beta=np.stack([np.interp(np.log(K),np.log(fixed['anchors']),row) for row in fixed['beta']]);f=a['weights']*V-cc@np.array(fixed['forward_beta']);call=a['weights'][:,None]*np.maximum(V[:,None]-K,0)-cc@beta;native=call-put[None,:]*(f[:,None]-K)
                block=lambda v:np.array([v[a['nets']==j].mean(0) for j in range(8)])
                components.append((block(native),block(call),block(f)));second.append(block(a['weights']*V**2))
            n,c,f=[components[0][j]+components[1][j]-components[2][j] for j in range(3)];npblock.append(n);cpblock.append(c);fpblock.append(f);e2block.append(second[0]+second[1]-second[2]);levels,_,_=unpack(x,p);lengths=np.maximum(0,np.minimum(ends,day/365+DELTA)-np.maximum(np.r_[0,ends[:-1]],day/365));exact=10000/DELTA*(lengths@levels);mc=e2block[-1];moments.append(dict(candidate=name,day=day,exact_second_moment=float(exact),multilevel_second_moment=float(mc.mean()),SE=float(mc.std(ddof=1)/np.sqrt(8))))
        independent=dict(native=np.concatenate(npblock,axis=1),call=np.concatenate(cpblock,axis=1),forward=np.column_stack(fpblock));saved=np.load(run/'fresh_check'/f'blocks_{name}.npz')
        for key,v in independent.items():
            error=float(max(abs(v-saved[key]).ravel()));assert error<2e-8,(name,key,error);rows.append(dict(candidate=name,field=key,maximum_block_replay_error=error))
    dump(out/'independent_block_replay.json',rows);dump(out/'second_moments.json',moments)
    a,_,pr,_=frozen_stage3();dated=DatedM2Pricer(a,pr);initial=np.array(prior['initial_state']);design=json.loads((run/'fresh_base_cache/design.json').read_text());redraw=[]
    for day,aa in big.items():
        dec=design['old_proposals'][str(day)];mix=DefensiveGaussianMixture(np.array([np.r_[old,z] for old in dec['shifts'] for z in ZSHIFTS]),np.array([p0*q for p0 in dec['probabilities'] for q in ZPROB]));sampler=GaussianOuterSampler(dated,day/365,initial)
        for net in range(8):
            ids=np.flatnonzero(aa['nets']==net);z,w,c=mix.draw(13,seed_for(design['seed_base'],day,net));assert np.array_equal(z[:,:-1],aa['normals'][ids]);assert np.array_equal(z[:,-1],aa['z_normals'][ids]);assert np.array_equal(w,aa['weights'][ids]);assert np.array_equal(c,aa['components'][ids]);assert np.array_equal(sampler.mean+z[:,:-1]@sampler.root.T,aa['states'][ids]);redraw.append(dict(day=day,net=net,pass_exact_redraw=True))
    dump(out/'new_base_redraw.json',redraw);d=json.loads((run/'design.json').read_text());hashes=[]
    for path,h in (d['source_sha256']|d['input_sha256']).items():
        good=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h;assert good,path;hashes.append(dict(path=path,pass_hash=good))
    dump(out/'audit.json',dict(status='PASS',source_input_checks=hashes,independent_price_block_replays=rows,new_large_scrambles_exact=32,estimator='base mean plus paired continuum-minus-boundary mean;8-net covariance retained',limitations=['Difference checking stream deliberately reused from002.','Finite-sample forward constraints are not population certificates.','Direct model-VIX quadrature audit is separate.']))
    print('MULTILEVEL AUDIT PASS',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);run(ap.parse_args().run.resolve())
