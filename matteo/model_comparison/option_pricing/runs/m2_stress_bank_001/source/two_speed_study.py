"""Joint four-expiry, forward-constrained two-speed stress pilot."""
import argparse,json,time,hashlib,shutil
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares,minimize
from .full_readout_fit import HERE,dump
from .vix_only_curve_fit import PRIOR,ARCHIVE,ACTIVE,DAYS,load_draws,paired
from .forward_floor_core import ForwardFloorFit
from .two_speed_core import TwoSpeedFit,metadata,boundary_x,LOG_HALF_MAX
from .two_speed_multilevel_core import TwoSpeedMultilevelFit
from .continuum_response_study import fit,checks,PARENT
from .mixed_response_study import ZSHIFTS,ZPROB,seed_for
from .plot_joint_curve_fit import implied_vol
from esn_eq.importance_sampling import DefensiveGaussianMixture

PREVIOUS=HERE/'runs/m2_continuum_response_003'


def augmented(draws,seed,out):
    """Independent randomized fast net, independently permuted within each block.

    Extra coordinate has target N(0,1); f/g multiplies existing likelihood.
    Independent scrambles and permutations avoid treating coincident 1D net
    ranks as a joint low-discrepancy design. SE units remain whole blocks.
    """
    result={};design=[]
    mix=DefensiveGaussianMixture(ZSHIFTS[:,None],ZPROB)
    for day,a in draws.items():
        b={k:v.copy() for k,v in a.items()};z=np.zeros(len(a['weights']));w=z.copy()
        for net in range(8):
            ids=np.flatnonzero(a['nets']==net);n=len(ids);assert n&(n-1)==0
            ss=seed_for(seed,day,net);zz,ww,_=mix.draw(int(np.log2(n)),ss)
            perm=np.random.default_rng(ss+12345).permutation(n)
            z[ids]=zz[perm,0];w[ids]=ww[perm]
            design.append(dict(day=int(day),net=net,seed=int(ss),permutation_seed=int(ss+12345),count=n))
        b['fast_normals']=z;b['weights']=a['weights']*w;result[day]=b
    out.mkdir(parents=True);dump(out/'design.json',dict(seeds=design,shifts=ZSHIFTS.tolist(),probabilities=ZPROB.tolist(),independent_old_fast=True,within_net_permutation=True))
    for day,a in result.items():np.savez_compressed(out/f'd{day}.npz',fast_normals=a['fast_normals'],fast_weight=a['weights']/draws[day]['weights'])
    return result


def engine(nodes,draws,p,ends,market,forwards,controls,seed,out):
    extra=augmented(draws,seed,out)
    args=(nodes,extra,market,forwards,controls,p,ends)
    return TwoSpeedMultilevelFit(TwoSpeedFit(*args),ForwardFloorFit(*args),ForwardFloorFit(nodes,draws,market,forwards,controls,p,ends))


def record(e,x,folder,name):
    folder.mkdir(parents=True,exist_ok=True);r=e.predict(x,True);blocks=r['blocks']
    for typ,q,key,values in [('otm',e.market,'native',r['prices']),('canonical',e.calls,'call',r['calls'])]:
        q=q.copy();q['model_price_U']=values;q['mc_se_U']=blocks[key].std(0,ddof=1)/np.sqrt(8)
        q['interval_residual']=np.where(values<q.price_bid_U,values-q.price_bid_U,np.where(values>q.price_ask_U,values-q.price_ask_U,0))/q.price_scale_U
        q.to_csv(folder/f'prices_{typ}_{name}.csv',index=False)
    f=e.forwards.copy();f['model_forward']=r['forwards'];f['forward_error']=r['forwards']-f.forward_vix_points.to_numpy();f['mc_se']=blocks['forward'].std(0,ddof=1)/np.sqrt(8);f.to_csv(folder/f'forwards_{name}.csv',index=False)
    np.savez_compressed(folder/f'blocks_{name}.npz',**blocks)
    np.savez_compressed(folder/f'paths_{name}.npz',**{f'd{day}_{k}':v for day,a in r['raw'].items() for k,v in a.items()})
    q=e.market.copy();q['model_price_U']=r['prices'];q['model_forward']=q.expiry_days.map(dict(zip(f.expiry_days,f.model_forward)))
    q['market_native_IV']=[implied_vol(a.price_mid_U,a.forward,a.strike,a.T_years,a.option_type) for a in q.itertuples()]
    q['model_native_IV']=[implied_vol(a.model_price_U,a.model_forward,a.strike,a.T_years,a.option_type) for a in q.itertuples()];q.to_csv(folder/f'native_iv_{name}.csv',index=False)
    score=dict(name=name,parameters=x.tolist(),model=metadata(x,e.p),option_loss=r['option_loss'],loss=r['loss'],inside_bidask=int(np.sum((r['prices']>=q.price_bid_U)&(r['prices']<=q.price_ask_U))),negative_native_prices=int(np.sum(r['prices']<0)),forward_errors=f.forward_error.tolist(),forward_SEs=f.mc_se.tolist(),moments=r['moments'])
    per=[]
    for day,g in q.groupby('expiry_days'):
        ids=g.index;wing=g.strike/g.forward>2
        per.append(dict(day=int(day),option_loss=float(np.mean((np.where(r['prices'][ids]<g.price_bid_U,r['prices'][ids]-g.price_bid_U,np.where(r['prices'][ids]>g.price_ask_U,r['prices'][ids]-g.price_ask_U,0))/g.price_scale_U)**2)/4),wing_IV_MAE=float(np.nanmean(abs(g.loc[wing,'model_native_IV']-g.loc[wing,'market_native_IV'])))))
    score['per_maturity']=per;dump(folder/f'score_{name}.json',score);return score


def run(out):
    out.mkdir(parents=True);start=time.perf_counter();prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);old=json.loads((PREVIOUS/'selection.json').read_text());base=np.array(old['boundary_parameters'])
    market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv');forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv');controls=json.loads((PRIOR/'vix_controls.json').read_text())
    assert len(market)==217 and p[6]>0
    sources=[HERE/'two_speed_core.py',HERE/'two_speed_multilevel_core.py',Path(__file__),HERE/'TWO_SPEED_PLAN.md',HERE/'forward_floor_core.py',HERE/'mixed_response_core.py',HERE/'vix_only_curve_fit.py',HERE/'continuum_response_study.py',PREVIOUS/'selection.json',PRIOR/'run.json',PRIOR/'vix_controls.json',ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv',ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv']
    src=out/'source';src.mkdir();hashes={}
    for f in sources:hashes[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest();shutil.copy2(f,src/f.name)
    train_nodes=PARENT/'fresh_conditional_nodes.npz';check_nodes=PREVIOUS/'fresh_base_nodes.npz'
    inputs=[train_nodes,check_nodes,*[PARENT/'fresh_cache/production'/f'vix_d{day}.npz' for day in DAYS],*[PREVIOUS/'fresh_base_cache/production'/f'vix_d{day}.npz' for day in DAYS]]
    dump(out/'design.json',dict(source_hashes=hashes,input_hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},training='32768 reused old-state/slow draws plus new independent permuted fast normals',checking='65536 previously used checking old-state/slow draws plus new fast normals; not fully fresh',SPX_priced=False,baseline_replaced=False,maturities=list(DAYS),native_quote_count=len(market),active_indices=ACTIVE.tolist(),forward_constraints='four training-estimator equalities; not exact population means',bounds=dict(amplitudes=[0,5],half_lives_days=[1,365],ordering='fast <= slow',floor_vol=[.05,.12]),selection='training loss among forward-feasible candidates only'))
    nodes=dict(np.load(train_nodes));e=engine(nodes,load_draws(PARENT/'fresh_cache'),p,ends,market,forwards,controls,20261027201,out/'training_fast')
    x=np.r_[base,.7,np.log(30)/LOG_HALF_MAX,.2];x[9]=(base[9]-x[12])/(1-x[12])
    probe=x.copy();probe[:7]=np.clip(probe[:7],.03,.8);probe[7:]=[.6,.65,.5,.5,.7,.5,.3]
    checks(e,probe,out)
    # Nesting check: fast fraction zero recovers the old estimator exactly.
    nested=x.copy();nested[13]=0;rr=e.predict(nested);bb=e.bl.predict(boundary_x(nested));err=float(np.max(abs(rr['prices']-bb['prices'])));assert err<1e-10
    dump(out/'boundary_replay.json',dict(maximum_price_error=err))
    ledger=[];starts=[]
    for bf,hf,r in [(3.5,30,.2),(5,10,.35),(4,50,.5),(3,3,.15)]:
        xx=x.copy();xx[11:]=[bf/5,np.log(hf)/LOG_HALF_MAX,r];xx[9]=(base[9]-xx[12])/(1-xx[12]);starts.append(xx)
    for j,xx in enumerate(starts):fit(e,xx,f'two_speed_start{j}',out,ledger,100)
    eligible=[z for z in ledger if z['eligible']];assert eligible
    winner=min(eligible,key=lambda z:z['option_loss']);selected=np.array(winner['parameters'])
    # Boundary is always an eligible candidate to retain a genuine nested comparison.
    ref=nested.copy();bl=[];(out/'single_fit').mkdir();ref[:11]=fit(e.bl,base,'single_speed_refit',out/'single_fit',bl,100);ref[12]=0;ref[13]=0
    dump(out/'single_optimizer.json',bl)
    if e.predict(ref)['option_loss']<winner['option_loss']:selected=ref.copy()
    dump(out/'selection.json',dict(parameters=selected.tolist(),model=metadata(selected,p),single_parameters=ref.tolist(),selected_before_checking=True,adopted=False))
    scores=[record(e,xx,out/'training',name) for name,xx in [('two_speed',selected),('single_speed',ref)]];dump(out/'training/scores.json',scores)
    del e,nodes
    ce=engine(dict(np.load(check_nodes)),load_draws(PREVIOUS/'fresh_base_cache'),p,ends,market,forwards,controls,20261027237,out/'checking_fast')
    scores=[record(ce,xx,out/'fresh_check',name) for name,xx in [('two_speed',selected),('single_speed',ref)]];dump(out/'fresh_check/scores.json',scores)
    pairs=paired(ce,out/'fresh_check','two_speed','single_speed');dump(out/'fresh_check/paired_losses.json',pairs)
    levels=p[7:-1].copy();from .forward_floor_core import unpack
    for name,xx in [('two_speed',selected),('single_speed',ref)]:
        lev,floor,_=unpack(xx,p);pp=p.copy();pp[7:-1]=lev;pp[-1]=floor
        dump(out/f'experimental_model_{name}.json',dict(parameters=pp.tolist(),ends_years=ends.tolist(),readout=metadata(xx,p),architecture=prior['architecture'],premium=prior['premium'],initial_state=prior['initial_state'],adopted=False))
    for f in sources:assert hashlib.sha256(f.read_bytes()).hexdigest()==hashes[str(f)]
    dump(out/'run.json',dict(status='completed pilot; numerical audit follows',wall_seconds=time.perf_counter()-start,checking_scores=scores,paired_losses=pairs,SPX_priced=False,baseline_replaced=False));print('TWO SPEED COMPLETE',time.perf_counter()-start,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();run(args.output.resolve())
