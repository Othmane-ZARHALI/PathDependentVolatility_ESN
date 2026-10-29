"""Exploratory independent-pilot centering of exact weighted VIX estimators."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

from .m2_spx import ROOT,M2
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.outer_sampling import replication_statistics


def payoffs(values,quotes):
    sign=np.where(quotes.option_type.to_numpy()=='call',1.,-1.)
    return np.column_stack((values,values**2,np.maximum(sign*(values[:,None]-quotes.strike.to_numpy()),0.)))


def run(args):
    pilot=Path(__file__).with_name('runs')/'m2_vix_importance_sampling_001'
    sources=capture_sources([Path(__file__),Path(__file__).with_name('m2_spx.py'),
        *sorted((M2/'src/esn_eq').glob('*.py'))])
    inputs={str(p):p.read_bytes() for p in (args.sampling/'run.json',args.sampling/'outer_values.csv',pilot/'run.json',pilot/'outer_values.csv')}
    record=json.loads((args.sampling/'run.json').read_text())
    old=json.loads((pilot/'run.json').read_text())
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    if record['curve_levels']!=old['curve_levels'] or record['target_sha256']!=old['target_sha256'] or record['target_sha256']!=targets.artifact_sha256:
        raise ValueError('Control pilot and checking target differ')
    previous=pd.read_csv(pilot/'outer_values.csv'); fresh=pd.read_csv(args.sampling/'outer_values.csv')
    rows,reps,centers=[],[],[]
    for days in record['scope']['vix_fit_expiry_days']:
        q=targets.select('VIX',expiry_days=[days]).reset_index(drop=True)
        training=previous[(previous.days==days)&(previous.method=='mixture')]
        c=np.mean(training.weight.to_numpy()[:,None]*payoffs(training.vix.to_numpy(),q),axis=0)
        descriptors=[dict(metric='forward',strike=None,option_type=None,scale=.5),dict(metric='second_moment',strike=None,option_type=None,scale=None)]
        descriptors += [dict(metric='option',strike=float(x.strike),option_type=x.option_type,scale=float(x.price_scale_U)) for x in q.itertuples()]
        for j,desc in enumerate(descriptors):
            centers.append(dict(days=days,**desc,independent_pilot_center=float(c[j])))
        for method in fresh.method.unique():
            for power in record['arguments']['powers']:
                matrix=[]; masses=[]
                for replicate in range(record['arguments']['replications']):
                    frame=fresh[(fresh.days==days)&(fresh.method==method)&(fresh.replicate==replicate)&(fresh.point<2**power)].sort_values('point')
                    w=frame.weight.to_numpy()
                    if len(w)!=2**power:
                        raise ValueError('Incomplete checking block')
                    # c is fixed independently of these checking draws. No ratio or fitted normalization.
                    estimate=c+np.mean(w[:,None]*(payoffs(frame.vix.to_numpy(),q)-c),axis=0)
                    matrix.append(estimate); masses.append(w.mean())
                    for j,desc in enumerate(descriptors):
                        reps.append(dict(days=days,method=method,power=power,replicate=replicate,**desc,estimate=float(estimate[j]),weight_mean=float(w.mean())))
                mean,se=replication_statistics(matrix)
                for j,desc in enumerate(descriptors):
                    rows.append(dict(days=days,method=method,power=power,**desc,estimate=float(mean[j]),replication_se=float(se[j]),se_over_scale=float(se[j]/desc['scale']) if desc['scale'] else None))
    check_sources(sources)
    if any(Path(p).read_bytes()!=raw for p,raw in inputs.items()):
        raise RuntimeError('Input changed')
    args.output.mkdir(parents=True,exist_ok=False); save_sources(args.output,sources,ROOT)
    pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False)
    pd.DataFrame(reps).to_csv(args.output/'replications.csv',index=False)
    pd.DataFrame(centers).to_csv(args.output/'centers.csv',index=False)
    report=dict(status='exploratory mass-control diagnostic, preserves raw estimator',
        implementation_sha256=source_hashes(sources,ROOT),input_sha256={k:hashlib.sha256(v).hexdigest() for k,v in inputs.items()},
        estimator='c + mean(w*(H-c)); c estimated solely from independent old mixture pilot',
        identity='E_g[w]=1, so subtract c*(w-1); no self-normalization or market forward substitution',
        design_timing='Added after inspecting the completed 27-day fresh raw comparison; all control comparisons are exploratory, not pristine validation of an ex ante method.',
        limitations=['No negative controlled prices are clipped; negatives remain unresolved diagnostics.',
            'Errors conditional on one fitted pilot; no pilot-to-pilot assessment.',
            'Same approximation biases and rare-tail limitations as raw study remain.'])
    (args.output/'run.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sampling',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error('Output exists')
    run(args)
