"""Fixed-parameter comparison of common calls and larger wing path budgets."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2,frozen_stage3
from .full_readout_fit import HERE,NAMES,Kernel,dump
from .parity_adaptive_engine import ParityAdaptiveEngine,distances
from .parity_stock import StockValues
from .parity_adaptive import canonical_market
from .audit_continuous_strike_controls import shape
from .plot_joint_curve_fit import implied_vol
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.pricing import DatedM2Pricer


def ivs(q,prices):
    return np.array([implied_vol(float(x),float(t.forward),float(t.strike),float(t.T_years),'call')
                     if x>0 else np.nan for t,x in zip(q.itertuples(),prices)],float)


def record(engine,p,name,out):
    begin=time.perf_counter();r=engine.predict(p,True,derivatives=False)
    q,f=engine.loss.tables(r['prices'],r['errors'],r['forwards'],r['forward_errors'])
    q['actual_paths']=r['actual_paths'];q['effective_paths']=r['effective_paths']
    q['scaled_distance']=np.r_[*[distances(g) for _,g in engine.spx.groupby('expiry_days',sort=True)],*[distances(g) for _,g in engine.vix.groupby('expiry_days',sort=True)]]
    q['native_otm_model_price_U']=r['native_prices'];q['native_otm_mc_se_U']=r['native_errors']
    q['model_iv_percent']=ivs(q,r['prices']);q['iv_unavailable']=~np.isfinite(q.model_iv_percent)
    q['below_market_forward_intrinsic']=q.model_price_U<np.maximum(q.forward-q.strike,0)
    q.to_csv(out/f'prices_canonical_{name}.csv',index=False);f.to_csv(out/f'forwards_{name}.csv',index=False)
    otm,_=engine.otm_loss.tables(r['native_prices'],r['native_errors'],r['forwards'],r['forward_errors']);otm.to_csv(out/f'prices_otm_{name}.csv',index=False)
    np.savez_compressed(out/f'paths_{name}.npz',**r['raw'])
    tiers=[]
    for (instrument,days),g0 in engine.original_quotes.groupby(['instrument','expiry_days'],sort=False):
        g=g0.reset_index(drop=True);target=canonical_market(g)
        y=r['raw'][f'spx_{days}' if instrument=='SPX' else f'vix_values_{days}']
        for n in [512,1024,2048,4096]:
            if instrument=='SPX':sample=y[:n];price=sample.mean(0);se=sample.std(0,ddof=1)/np.sqrt(n)
            else:
                b=y.reshape(8,512,len(g))[:,:n//8].mean(1);price=b.mean(0);se=b.std(0,ddof=1)/np.sqrt(8)
            t=target.copy();t['candidate']=name;t['uniform_paths']=n;t['model_price_U']=price;t['mc_se_U']=se;t['model_iv_percent']=ivs(t,price);tiers.append(t)
    pd.concat(tiers,ignore_index=True).to_csv(out/f'uniform_budget_prices_{name}.csv',index=False)
    summary=dict(candidate=name,parameters=p.tolist(),canonical_loss=r['total'],canonical_components=r['components'],
                 original_otm_loss=r['original_otm_total'],original_otm_components=r['original_otm_components'],
                 inside_canonical_bidask=int(q.interval_residual.eq(0).sum()),inside_original_otm_bidask=int(otm.interval_residual.eq(0).sum()),
                 canonical_negative_prices=int(q.negative_estimate.sum()),native_otm_negative_prices=int(otm.negative_estimate.sum()),
                 iv_unavailable=int(q.iv_unavailable.sum()),normalizer=r['normalizer'],allocations=r['allocations'],moments=r['moments'],
                 shapes=[dict(instrument=ins,days=int(day),**shape(g,g.model_price_U.to_numpy())) for (ins,day),g in q.groupby(['instrument','expiry_days'])],
                 seconds=time.perf_counter()-begin)
    dump(out/f'score_{name}.json',summary);print('RECORD',json.dumps(summary),flush=True)
    return summary


def run(cache,controls,out):
    out.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    # Dependencies are ongoing authorized computations, not future scheduled work.
    while not (cache/'run.json').exists() or not (controls/'run.json').exists():
        if time.perf_counter()-begin>7200:raise TimeoutError('Cache/control preparation exceeded two hours')
        time.sleep(3)
    parent=HERE/'runs/m2_controlled_readout_fit_002';old=json.loads((parent/'run.json').read_text())
    selected=np.array(old['selected_parameters']);baseline=np.array(old['baseline'])
    spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');forwards=pd.read_csv(cache/'forward_targets.csv')
    sources=[Path(__file__),HERE/'parity_adaptive.py',HERE/'parity_adaptive_engine.py',HERE/'parity_stock.py',HERE/'parity_stock_kernel.cpp',HERE/'parity_hedge_controls.py',
             HERE/'controlled_readout_kernel.cpp',HERE/'controlled_readout_convolution_kernel.cpp',HERE/'controlled_readout_fast_fit.py',HERE/'all_maturity_core.py',
             HERE/'PARITY_ADAPTIVE_PLAN.md',HERE/'PARITY_ADAPTIVE_PLAN_V2.md',HERE/'PARITY_UNIFORM_PLAN.md',HERE/'plot_joint_curve_fit.py',HERE/'test_parity_adaptive.py',HERE/'test_parity_adaptive_engine.py',HERE/'test_parity_hedge_controls.py',HERE/'reusable_ou_initial.py',HERE/'test_reusable_ou_initial.py',parent/'run.json',cache/'run.json',cache/'design.json',
             controls/'run.json',controls/'coefficients.json',*sorted((M2/'src/esn_eq').glob('*.py'))]
    captured=capture_sources(sources);save_sources(out,captured,ROOT)
    lib=out/'convolution.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'controlled_readout_convolution_kernel.cpp'),'-o',str(lib)],check=True)
    stocklib=out/'stock.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'parity_stock_kernel.cpp'),'-o',str(stocklib)],check=True)
    a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial,_=dated.initial_state([1,0,0])
    design=dict(status='fixed-parameter estimator comparison',cache=str(cache),controls=str(controls),parameter_names=NAMES,selected=selected.tolist(),baseline=baseline.tolist(),
                source_sha256=source_hashes(captured,ROOT),architecture=asdict(a),premium=asdict(premium),initial_state=initial.tolist(),
                paths_per_strike=4096,allocation='uniform per user follow-up',diagnostic_prefixes=[512,1024,2048,4096],market_otm_inputs=len(spx)+len(vix),market_out_of_sample=False,
                vix_parity='Market puts converted using observed market F; native model puts reconstructed using model E[VIX]. Both objectives and forward errors retained.',
                changed_estimator='Canonical call payoff, independently piloted common hedge/Gaussian controls, 4096 paths at EVERY strike. Previous controlled losses are not direct comparators.',
                fitted_in_this_run=False)
    dump(out/'design.json',design)
    def digest_file(path):
        h=hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
        return h.hexdigest()
    immutable={str(p):digest_file(p) for root in [cache,controls] for p in sorted(root.rglob('*.np*'))}
    dump(out/'immutable_inputs_before.json',immutable)
    print('IMMUTABLE INPUTS FROZEN',len(immutable),flush=True)
    engine=ParityAdaptiveEngine(cache,Kernel(lib),dated,initial,spx,vix,forwards,frozen_controls=controls,stock_values=StockValues(stocklib))
    try:
        pilot=engine.fit_vix_controls(selected);dump(out/'vix_controls.json',pilot)
        scores=[record(engine,p,name,out) for p,name in [(selected,'selected'),(baseline,'baseline')]]
    finally:engine.close()
    check_sources(captured)
    for path,digest in immutable.items():
        if digest_file(Path(path))!=digest:raise RuntimeError('Immutable input changed: '+path)
    dump(out/'immutable_inputs_verified.json',dict(status='all unchanged',arrays=len(immutable)))
    dump(out/'run.json',design|dict(status='complete fixed-parameter estimator comparison',scores=scores,wall_seconds=time.perf_counter()-begin))
    print('COMPLETED',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',type=Path,required=True);a.add_argument('--controls',type=Path,required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args()
    run(args.cache.resolve(),args.controls.resolve(),args.output.resolve())
