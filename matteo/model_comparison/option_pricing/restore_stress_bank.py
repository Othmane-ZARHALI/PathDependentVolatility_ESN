"""Restore all four forwards after bounded joint SLSQP, without changing the model.

Only seven free variance levels move in this finishing step. Each original and
restored candidate is retained. Selection uses training loss and feasibility;
new additional checking coordinates are generated only after selection.
"""
import argparse,json,time,hashlib,shutil
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from .full_readout_fit import HERE,dump
from .vix_only_curve_fit import PRIOR,ARCHIVE,ACTIVE,DAYS,load_draws,paired
from .stress_bank_study import PREVIOUS,OLD_CHECK,PARENT,engine,record
from .stress_bank_core import StressBankFit,metadata,unpack
from .stress_bank_multilevel_core import StressBankMultilevelFit


def cached_engine(parent,nodes,draws,p,ends,market,forwards,controls):
    bank={}
    for day,a in draws.items():
        saved=dict(np.load(parent/'training_bank'/f'd{day}.npz'));bank[day]=dict(a,bank_normals=saved['bank_normals'],weights=saved['combined_weights'])
    base={day:dict(a,bank_normals=np.column_stack([np.zeros((len(a['weights']),3)),a['z_normals']])) for day,a in draws.items()}
    args=(nodes,bank,market,forwards,controls,p,ends)
    return StressBankMultilevelFit(StressBankFit(*args),StressBankFit(*args,reference=True),StressBankFit(nodes,base,market,forwards,controls,p,ends,reference=True))


def run(parent,out):
    out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv');forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv');controls=json.loads((PRIOR/'vix_controls.json').read_text());target=forwards.forward_vix_points.to_numpy()
    src=out/'source';src.mkdir();sources=[Path(__file__),HERE/'stress_bank_study.py',HERE/'stress_bank_core.py',HERE/'stress_bank_multilevel_core.py',HERE/'STRESS_BANK_PLAN.md',HERE/'STRESS_BANK_REFINEMENT.md',HERE/'vix_only_curve_fit.py',HERE/'mixed_response_core.py',HERE/'continuum_response_study.py',PRIOR/'run.json',PRIOR/'vix_controls.json',ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv',ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv']
    hashes={}
    for f in sources:hashes[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest();shutil.copy2(f,src/f.name)
    inputs=[PARENT/'fresh_conditional_nodes.npz',OLD_CHECK/'fresh_base_nodes.npz',parent/'design.json',parent/'jacobian_checks.json',parent/'reference_replay.json']
    for day in DAYS:inputs.extend([PARENT/'fresh_cache/production'/f'vix_d{day}.npz',OLD_CHECK/'fresh_base_cache/production'/f'vix_d{day}.npz',parent/'training_bank'/f'd{day}.npz',PREVIOUS/'checking_fast'/f'd{day}.npz'])
    design=dict(status='forward restoration running',parent=str(parent),source_hashes=hashes,input_hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},model_unchanged=True,grid_unchanged=True,algorithm='bounded least-squares on seven xi units after each of three joint SLSQP starts; preserve originals; select training-feasible minimum option loss',forward_constraints='same four training-estimator equalities, max gap<=1e-5',SPX_priced=False,baseline_replaced=False,checking='same65536 old/node2/node3 draws plus new additional node0/1 seed20261027379; not wholly fresh')
    dump(out/'design.json',design)
    e=cached_engine(parent,dict(np.load(PARENT/'fresh_conditional_nodes.npz')),load_draws(PARENT/'fresh_cache'),p,ends,market,forwards,controls)
    processed=0;rows=[];candidates=[]
    while processed<3:
        try:ledger=json.loads((parent/'optimizer_ledger.json').read_text())
        except (FileNotFoundError,json.JSONDecodeError):ledger=[]
        if len(ledger)<=processed:
            if time.perf_counter()-begin>7200:raise TimeoutError('Parent multistart not complete')
            time.sleep(5);continue
        rec=ledger[processed];x=np.array(rec['parameters']);before=e.predict(x);memo={}
        def pred(y):
            key=tuple(y)
            if key not in memo:
                xx=x.copy();xx[:7]=y;memo.clear();memo[key]=e.predict(xx)
            return memo[key]
        opt=least_squares(lambda y:pred(y)['forwards']-target,np.clip(x[:7],0,1),jac=lambda y:pred(y)['forward_jacobian'][:,:7],bounds=(np.zeros(7),np.ones(7)),max_nfev=60,ftol=1e-12,xtol=1e-12,gtol=1e-11)
        xx=x.copy();xx[:7]=opt.x;after=e.predict(xx);gap=float(max(abs(after['forwards']-target)))
        row=dict(parent_label=rec['label'],original_parameters=x.tolist(),restored_parameters=xx.tolist(),original_option_loss=before['option_loss'],restored_option_loss=after['option_loss'],original_forward_gap=float(max(abs(before['forwards']-target))),restored_forward_gap=gap,eligible=gap<=1e-5,maximum_xi_unit_change=float(max(abs(xx[:7]-x[:7]))),restoration_nfev=opt.nfev,restoration_success=bool(opt.success),model=metadata(xx,p));rows.append(row)
        if rec['eligible']:candidates.append(dict(label=rec['label']+'_original',parameters=x.tolist(),option_loss=before['option_loss']))
        if row['eligible']:candidates.append(dict(label=rec['label']+'_restored',parameters=xx.tolist(),option_loss=after['option_loss']))
        dump(out/'restoration_ledger.json',rows);print('RESTORATION',json.dumps(row),flush=True);processed+=1
    assert candidates,'No forward-feasible candidate after restoration'
    best=min(candidates,key=lambda z:z['option_loss']);selected=np.array(best['parameters']);dump(out/'selection.json',dict(label=best['label'],parameters=selected.tolist(),model=metadata(selected,p),all_candidates=candidates,selected_before_new_extra_checking_draws=True,adopted=False));shutil.copy2(parent/'optimizer_ledger.json',out/'parent_optimizer_ledger.json')
    shutil.copy2(parent/'jacobian_checks.json',out/'jacobian_checks.json');shutil.copy2(parent/'reference_replay.json',out/'reference_replay.json')
    train=record(e,selected,out/'training','bank');del e
    ce=engine(dict(np.load(OLD_CHECK/'fresh_base_nodes.npz')),load_draws(OLD_CHECK/'fresh_base_cache'),PREVIOUS/'checking_fast',p,ends,market,forwards,controls,20261027379,out/'checking_bank');check=record(ce,selected,out/'fresh_check','bank')
    for prefix in ['prices_otm','prices_canonical','forwards','blocks','native_iv','score']:
        ext='npz' if prefix=='blocks' else ('json' if prefix=='score' else 'csv');shutil.copy2(PREVIOUS/'fresh_check'/f'{prefix}_two_speed.{ext}',out/'fresh_check'/f'{prefix}_two_speed.{ext}')
    comparison=paired(ce,out/'fresh_check','bank','two_speed');dump(out/'fresh_check/paired_losses.json',comparison)
    levels,floor,_=unpack(selected,p);pp=p.copy();pp[7:-1]=levels;pp[-1]=floor;dump(out/'experimental_model_bank.json',dict(parameters=pp.tolist(),ends_years=ends.tolist(),readout=metadata(selected,p),architecture=prior['architecture'],premium=prior['premium'],initial_state=prior['initial_state'],adopted=False))
    pd.DataFrame(dict(start_day=np.r_[0,ends[:-1]]*365,end_day=ends*365,active_in_fit=np.isin(np.arange(16),ACTIVE),bank=levels)).to_csv(out/'curve_levels.csv',index=False)
    for f in sources:assert hashlib.sha256(f.read_bytes()).hexdigest()==hashes[str(f)]
    dump(out/'run.json',dict(status='completed forward-restored pilot; audit follows',model_unchanged=True,training_score=train,checking_score=check,paired_losses=comparison,SPX_priced=False,baseline_replaced=False,wall_seconds=time.perf_counter()-begin));print('RESTORED BANK COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();run(args.parent.resolve(),args.output.resolve())
