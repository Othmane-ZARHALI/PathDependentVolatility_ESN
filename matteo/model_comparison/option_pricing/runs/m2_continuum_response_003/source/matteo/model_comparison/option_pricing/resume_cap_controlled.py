"""Resume ONLY fresh VIX pricing after the preserved cache-setup failure."""
import argparse,ctypes as ct,hashlib,json,time
from pathlib import Path
import numpy as np,pandas as pd
from .full_readout_fit import HERE,dump
from .m2_spx import ROOT,frozen_stage3
from .maturity_curve import CurveKernel
from .cap_tail_study import engine,tails
from .cap_controlled_study import configure
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from .audit_maturity_curve import replay
from esn_eq.pricing import DatedM2Pricer


def run(out):
    begin=time.perf_counter();d=json.loads((out/'design.json').read_text());prior=Path(d['parent']);pd0=json.loads((prior/'run.json').read_text());base=np.array(d['baseline']);selection=json.loads((out/'selection.json').read_text())['candidates'];cache=Path(d['cache']);fresh=out/'fresh_cache';check=out/'fresh_check';assert not list(check.glob('score_*.json'))
    for rel,sha in d['source_sha256'].items():
        snap=out/'source'/rel;assert hashlib.sha256(snap.read_bytes()).hexdigest()==sha
        if Path(rel).name!='cap_controlled_study.py':assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==sha
    fd=json.loads((fresh/'design.json').read_text())
    for row in fd['files']:assert hashlib.sha256((fresh/'production'/f"vix_d{row['days']}.npz").read_bytes()).hexdigest()==row['sha256']
    files=[Path(__file__),HERE/'cap_controlled_study.py',HERE/'cap_tail_study.py',out/'selection.json',out/'full/scores.json',fresh/'design.json',out/'fresh_setup_fix.diff'];captured=capture_sources(files);dest=out/'resume';dest.mkdir(exist_ok=False);save_sources(dest,captured,ROOT)
    dump(dest/'design.json',dict(status='resume only fresh VIX prices on existing frozen draws',source_sha256=source_hashes(captured,ROOT),original_failure_preserved=True,repeated_optimization=False,repeated_SPX=False,regenerated_random_streams=False,fresh_prices_before_failure=0))
    spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');q=pd.concat([spx,vix],ignore_index=True);f=pd.read_csv(cache/'forward_targets.csv');vc=json.loads((prior/'vix_controls.json').read_text())
    for row in selection:replay(out/'full',row['name'],np.array(row['parameters'])[7:-1],q)
    lib=out/'cap.dylib';k=CurveKernel(lib,lib,np.array(pd0['ends_years']));k.lib.configure.argtypes=[ct.c_int,ct.c_double];a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial=np.array(pd0['initial_state']);scores=[];ts=[];ms=[]
    for row in [dict(name='baseline',parameters=base.tolist(),continuation=0),*selection]:
        configure(k,row['continuation']);x=np.array(row['parameters']);e=engine(out,cache,k,dated,initial,spx,vix,f,vc,None,pd0);e.cache=fresh;e.folder=fresh/'production';print('FRESH CHECK',row['name'],flush=True)
        try:
            r=e.predict(x,True,False);scores.append(save(e,r,x,row['name'],check));t,m=tails(check,row['name'],vix,x,np.array(pd0['ends_days']));ts+=t;ms+=m
        finally:e.close()
        print('FRESH SCORE',row['name'],scores[-1]['original_otm_components'],flush=True)
    pd.DataFrame(ts).to_csv(check/'tails.csv',index=False);pd.DataFrame(ms).to_csv(check/'moments.csv',index=False);dump(check/'scores.json',scores);check_sources(captured)
    dump(out/'run.json',d|dict(status='completed controlled cap comparison; fresh setup failure recovered without repeating fits or SPX',selection=selection,full_scores=json.loads((out/'full/scores.json').read_text()),fresh_scores=scores,baseline_adopted_or_replaced=False,new_layer_implemented=False,resumed_only_fresh_VIX=True,resume_wall_seconds=time.perf_counter()-begin,original_failure='failure.json',original_sources='source',resume_sources='resume/source'))
    dump(dest/'run.json',dict(status='completed fresh-only continuation',seconds=time.perf_counter()-begin,scores=scores));print('COMPLETE FRESH RESUME',time.perf_counter()-begin,flush=True)

from .fit_parity_parameters_v3 import save_result as save
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);run(p.parse_args().run.resolve())
