"""Joint calibration of a four-node geometric stress bank with cached states."""
import argparse,json,time,hashlib,shutil
from pathlib import Path
import numpy as np
import pandas as pd
from .full_readout_fit import HERE,dump
from .vix_only_curve_fit import PRIOR,ARCHIVE,ACTIVE,DAYS,load_draws,paired
from .stress_bank_core import StressBankFit,metadata,unpack,HALF_LIVES,weights,stick_units
from .stress_bank_multilevel_core import StressBankMultilevelFit
from .continuum_response_study import fit,PARENT
from .mixed_response_study import seed_for
from .plot_joint_curve_fit import implied_vol
from esn_eq.importance_sampling import DefensiveGaussianMixture

PREVIOUS=HERE/'runs/m2_two_speed_001'
OLD_CHECK=HERE/'runs/m2_continuum_response_003'
EXTRA_SHIFTS=np.array([[0,0],[2,0],[4,0],[6,0],[0,2],[0,4],[0,6]],float)
EXTRA_PROB=np.array([.5,.1,.1,.05,.1,.1,.05])


def augmented(draws,fast_folder,seed,out):
    result={};design=[];mix=DefensiveGaussianMixture(EXTRA_SHIFTS,EXTRA_PROB)
    out.mkdir(parents=True)
    for day,a in draws.items():
        fast=dict(np.load(fast_folder/f'd{day}.npz'));assert len(fast['fast_normals'])==len(a['weights'])
        b=dict(a);z=np.zeros((len(a['weights']),4));z[:,2]=fast['fast_normals'];z[:,3]=a['z_normals'];ew=np.zeros(len(z))
        for net in range(8):
            ids=np.flatnonzero(a['nets']==net);n=len(ids);ss=seed_for(seed,day,net);zz,ww,_=mix.draw(int(np.log2(n)),ss);perm=np.random.default_rng(ss+12345).permutation(n)
            z[ids,:2]=zz[perm];ew[ids]=ww[perm];design.append(dict(day=int(day),net=net,count=n,seed=int(ss),permutation_seed=int(ss+12345)))
        b['bank_normals']=z;b['weights']=a['weights']*fast['fast_weight']*ew;result[day]=b
        np.savez_compressed(out/f'd{day}.npz',bank_normals=z,extra_weights=ew,combined_weights=b['weights'])
    dump(out/'design.json',dict(seeds=design,extra_shifts=EXTRA_SHIFTS.tolist(),extra_probabilities=EXTRA_PROB.tolist(),previous_fast_folder=str(fast_folder),coordinate_order=['7d','28d','112d','448d'],target_independent=True,within_net_permutation=True))
    return result


def engine(nodes,draws,fast_folder,p,ends,market,forwards,controls,seed,out):
    bank=augmented(draws,fast_folder,seed,out);base={day:dict(a,bank_normals=np.column_stack([np.zeros((len(a['weights']),3)),a['z_normals']])) for day,a in draws.items()}
    args=(nodes,bank,market,forwards,controls,p,ends)
    return StressBankMultilevelFit(StressBankFit(*args),StressBankFit(*args,reference=True),StressBankFit(nodes,base,market,forwards,controls,p,ends,reference=True))


def checks(e,x,out):
    rows=[];r=e.predict(x)
    for j in range(16):
        dx=np.zeros(16);dx[j]=1e-6;rp=e.predict(x+dx);rm=e.predict(x-dx)
        fp=(rp['prices']-rm['prices'])/2e-6;ff=(rp['forwards']-rm['forwards'])/2e-6
        pe=float(np.linalg.norm(fp-r['native_jacobian'][:,j])/max(1,np.linalg.norm(fp)));fe=float(np.linalg.norm(ff-r['forward_jacobian'][:,j])/max(1,np.linalg.norm(ff)))
        rows.append(dict(parameter=j,price_error=pe,forward_error=fe));assert pe<.003 and fe<1e-5,rows[-1]
    xx=x.copy();xx[12:15]=0;xx[15]=.65;a=e.predict(xx);b=e.bl.predict(xx);replay=float(max(abs(a['prices']-b['prices'])));assert replay<1e-10
    dump(out/'jacobian_checks.json',rows);dump(out/'reference_replay.json',dict(maximum_price_error=replay));print('ALL 16 JACOBIANS AND REFERENCE REPLAY PASS',flush=True)


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
    out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);old=json.loads((PREVIOUS/'selection.json').read_text());priorx=np.array(old['parameters'])
    market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv');forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv');controls=json.loads((PRIOR/'vix_controls.json').read_text());assert len(market)==217 and p[4]>0
    sources=[Path(__file__),HERE/'stress_bank_core.py',HERE/'stress_bank_multilevel_core.py',HERE/'STRESS_BANK_PLAN.md',HERE/'two_speed_core.py',HERE/'two_speed_study.py',HERE/'continuum_response_study.py',HERE/'vix_only_curve_fit.py',HERE/'mixed_response_core.py',PRIOR/'run.json',PRIOR/'vix_controls.json',PREVIOUS/'selection.json',ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv',ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv',HERE.parents[1]/'exogenous_reservoir_eq_m2/src/esn_eq/importance_sampling.py']
    src=out/'source';src.mkdir();hashes={}
    for f in sources:hashes[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest();shutil.copy2(f,src/f.name)
    train_nodes=PARENT/'fresh_conditional_nodes.npz';check_nodes=OLD_CHECK/'fresh_base_nodes.npz';inputs=[train_nodes,check_nodes]
    for phase,cache in [('training',PARENT/'fresh_cache'),('checking',OLD_CHECK/'fresh_base_cache')]:
        inputs.extend([cache/'production'/f'vix_d{day}.npz' for day in DAYS]);inputs.extend([PREVIOUS/f'{phase}_fast'/f'd{day}.npz' for day in DAYS])
    dump(out/'design.json',dict(source_hashes=hashes,input_hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},half_lives_days=HALF_LIVES.tolist(),grid_selected_before_fit=True,training='32768 reused old/slow/fast normals plus new joint two-dimensional extra nodes',checking='65536 reused old/slow/fast normals plus new extra nodes; not wholly fresh',SPX_priced=False,baseline_replaced=False,parameter_order=['xi_unit'+str(j) for j in ACTIVE]+['floor_vol_unit']+['b_unit'+str(j) for j in range(4)]+['stick_unit'+str(j) for j in range(4)],native_quote_count=len(market),forward_constraints='four training-estimator equalities; eligibility max gap<=1e-5',bounds=dict(amplitudes=[0,5],floor_volatility=[.05,.12],active_xi_upper=.2),selection='minimum training option loss among forward-feasible starts',two_speed_comparison='stored selected model and checking blocks; shared old-state streams; not a nested family comparison'))
    e=engine(dict(np.load(train_nodes)),load_draws(PARENT/'fresh_cache'),PREVIOUS/'training_fast',p,ends,market,forwards,controls,20261027301,out/'training_bank')
    x=np.r_[priorx[:7],priorx[10],[.6,.6,.7,.46],stick_units([.05,.10,.47,.30])]
    probe=x.copy();probe[:7]=np.clip(probe[:7],.03,.8);probe[7]=.5;probe[8:12]=[.5,.55,.6,.65];probe[12:]=stick_units([.1,.15,.2,.3]);checks(e,probe,out)
    starts=[]
    for b,alpha in [([3,3,3.5,2.3],[.05,.10,.47,.30]),([3.5]*4,[.225]*4),([2,4,3,2],[.10,.25,.35,.20])]:
        xx=x.copy();xx[8:12]=np.array(b)/5;xx[12:]=stick_units(alpha);starts.append(xx)
    ledger=[]
    for j,xx in enumerate(starts):fit(e,xx,f'bank_start{j}',out,ledger,100)
    eligible=[z for z in ledger if z['eligible']];assert eligible,'No forward-feasible candidate; retain diagnostics'
    winner=min(eligible,key=lambda z:z['option_loss']);selected=np.array(winner['parameters']);dump(out/'selection.json',dict(parameters=selected.tolist(),model=metadata(selected,p),selected_before_checking=True,adopted=False))
    train=record(e,selected,out/'training','bank');del e
    ce=engine(dict(np.load(check_nodes)),load_draws(OLD_CHECK/'fresh_base_cache'),PREVIOUS/'checking_fast',p,ends,market,forwards,controls,20261027337,out/'checking_bank');check=record(ce,selected,out/'fresh_check','bank')
    for prefix in ['prices_otm','prices_canonical','forwards','blocks','native_iv','score']:
        ext='npz' if prefix=='blocks' else ('json' if prefix=='score' else 'csv');shutil.copy2(PREVIOUS/'fresh_check'/f'{prefix}_two_speed.{ext}',out/'fresh_check'/f'{prefix}_two_speed.{ext}')
    comparison=paired(ce,out/'fresh_check','bank','two_speed');dump(out/'fresh_check/paired_losses.json',comparison)
    levels,floor,_=unpack(selected,p);pp=p.copy();pp[7:-1]=levels;pp[-1]=floor
    dump(out/'experimental_model_bank.json',dict(parameters=pp.tolist(),ends_years=ends.tolist(),readout=metadata(selected,p),architecture=prior['architecture'],premium=prior['premium'],initial_state=prior['initial_state'],adopted=False))
    table=pd.DataFrame(dict(start_day=np.r_[0,ends[:-1]]*365,end_day=ends*365,active_in_fit=np.isin(np.arange(16),ACTIVE),bank=levels));table.to_csv(out/'curve_levels.csv',index=False)
    for f in sources:assert hashlib.sha256(f.read_bytes()).hexdigest()==hashes[str(f)]
    dump(out/'run.json',dict(status='completed pilot; audit follows',training_score=train,checking_score=check,paired_losses=comparison,SPX_priced=False,baseline_replaced=False,wall_seconds=time.perf_counter()-begin));print('STRESS BANK COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();run(args.output.resolve())
