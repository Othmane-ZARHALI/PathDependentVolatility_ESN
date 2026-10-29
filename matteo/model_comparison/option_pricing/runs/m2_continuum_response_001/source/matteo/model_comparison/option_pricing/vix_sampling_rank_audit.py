"""Paired fixed-curve ranking on fresh weighted VIX sampling observations."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .m2_spx import ROOT,M2
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.calibration_loss import JointPriceObjective,LossConfig


def transform_vix(values,old_level,new_level,floor=1e-4):
    """Only valid when the entire VIX window lies within one curve segment."""
    return np.sqrt(10000*floor+(new_level-floor)/(old_level-floor)*(np.asarray(values)**2-10000*floor))


def run(args):
    start_folder=Path(__file__).with_name('runs')/'m2_vix_curve_fit_001'
    sources=capture_sources([Path(__file__),Path(__file__).with_name('m2_spx.py'),
        *sorted((M2/'src/esn_eq').glob('*.py'))])
    names=('run.json','outer_values.csv')
    inputs={name:(args.sampling/name).read_bytes() for name in names}
    sampling=json.loads(inputs['run.json']); old=json.loads((start_folder/'run.json').read_text())
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    if sampling['target_sha256']!=targets.artifact_sha256 or sampling['curve_levels']!=old['selected_levels']:
        raise ValueError('Sampling target or anchor curve differs')
    family=NuisanceCurveFamily.from_targets(targets)
    quotes=targets.select('VIX').reset_index(drop=True)
    forwards=targets.vix_forwards.sort_values('expiry').reset_index(drop=True)
    objective=JointPriceObjective(quotes,forwards,LossConfig(**old['loss_config']))
    curves={name:np.asarray(old[name+'_levels']) for name in ('starting','selected')}
    if sampling['scope']['vix_fit_expiry_days']!=[27,55,118,153]:
        raise ValueError('Require full scope')
    # Independent cached-exposure check of the affine transformation, at every expiry.
    cache=np.load(start_folder/'exposures.npz'); transform_checks=[]
    for target in forwards.itertuples():
        manifest=old['exposure_manifest'][str(int(target.expiry))]
        segment=int(np.searchsorted(family.ends,target.T_years,side='right'))
        if target.T_years+manifest['delta']>=family.ends[segment]:
            raise ValueError('Single-level VIX transformation invalid across a curve knot')
        coefficients=cache[f'role0_expiry{int(target.expiry)}']
        selected=100*np.sqrt((manifest['floor_integral']+coefficients@(curves['selected']-1e-4))/manifest['delta'])
        starting=100*np.sqrt((manifest['floor_integral']+coefficients@(curves['starting']-1e-4))/manifest['delta'])
        mapped=transform_vix(selected,curves['selected'][segment],curves['starting'][segment])
        np.testing.assert_allclose(mapped,starting,rtol=2e-14,atol=1e-13)
        transform_checks.append(dict(days=int(target.expiry_days),segment=segment,max_error=float(np.max(np.abs(mapped-starting)))))
    raw=pd.read_csv(args.sampling/'outer_values.csv')
    R=sampling['arguments']['replications']; N=2**max(sampling['arguments']['powers'])
    scores,prices,forward_tables,summary,loo=[],[],[],[],[]
    for method in sorted(raw.method.unique()):
        data={}
        for candidate,levels in curves.items():
            p,f=np.zeros((R,len(quotes))),np.zeros((R,len(forwards)))
            for fi,target in enumerate(forwards.itertuples()):
                q=quotes[quotes.expiry==target.expiry]
                group=raw[(raw.days==target.expiry_days)&(raw.method==method)].sort_values(['replicate','point'])
                if len(group)!=R*N:
                    raise ValueError('Incomplete sampling block')
                v=group.vix.to_numpy().reshape(R,N); w=group.weight.to_numpy().reshape(R,N)
                segment=int(np.searchsorted(family.ends,target.T_years,side='right'))
                v=transform_vix(v,curves['selected'][segment],levels[segment])
                sign=np.where(q.option_type.to_numpy()=='call',1.,-1.)
                h=np.maximum(sign*(v[:,:,None]-q.strike.to_numpy()),0.)
                p[:,q.index]=(w[:,:,None]*h).mean(axis=1); f[:,fi]=(w*v).mean(axis=1)
            result=objective.evaluate(p.mean(axis=0),f.mean(axis=0),curve_penalty=family.penalty(levels))
            data[candidate]=(p,f,result)
            scores.append(dict(method=method,candidate=candidate,total=result['total'],**result['components']))
            table=result['quotes'].copy(); table['method']=method; table['candidate']=candidate
            table['sampling_se']=p.std(axis=0,ddof=1)/np.sqrt(R); prices.append(table)
            ft=forwards.copy(); ft['method']=method; ft['candidate']=candidate
            ft['model_forward']=f.mean(axis=0); ft['sampling_se']=f.std(axis=0,ddof=1)/np.sqrt(R); forward_tables.append(ft)
        change=data['selected'][2]['total']-data['starting'][2]['total']
        leave=[]
        for omitted in range(R):
            keep=np.arange(R)!=omitted
            values={name:objective.evaluate(d[0][keep].mean(axis=0),d[1][keep].mean(axis=0),curve_penalty=family.penalty(curves[name]))['total'] for name,d in data.items()}
            delta=values['selected']-values['starting']; leave.append(delta)
            loo.append(dict(method=method,omitted=omitted,selected_minus_starting=delta))
        jack=np.sqrt((R-1)/R*np.sum((leave-np.mean(leave))**2))
        summary.append(dict(method=method,selected_minus_starting=change,descriptive_jackknife_se=float(jack),
            leave_one_out_min=float(min(leave)),leave_one_out_max=float(max(leave))))
    check_sources(sources)
    if any((args.sampling/name).read_bytes()!=raw_bytes for name,raw_bytes in inputs.items()):
        raise RuntimeError('Sampling inputs changed')
    args.output.mkdir(parents=True,exist_ok=False); save_sources(args.output,sources,ROOT)
    record=dict(status='paired fixed-candidate loss audit, no fitting or ranking acceptance',
        sampling=str(args.sampling),input_sha256={k:hashlib.sha256(v).hexdigest() for k,v in inputs.items()},
        original_fit_sha256=hashlib.sha256((start_folder/'run.json').read_bytes()).hexdigest(),
        original_exposures_sha256=hashlib.sha256((start_folder/'exposures.npz').read_bytes()).hexdigest(),
        implementation_sha256=source_hashes(sources,ROOT),target_sha256=targets.artifact_sha256,
        summaries=summary,transform_checks=transform_checks,candidates={k:v.tolist() for k,v in curves.items()},
        loss_config=asdict(objective.config),limitations=['Eight-scramble jackknife is descriptive under nonsmooth interval losses and unresolved tails.',
            'Numerical errors and pilot variability remain; no final fit ranking or held-out claim.'])
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    pd.DataFrame(scores).to_csv(args.output/'scores.csv',index=False)
    pd.DataFrame(loo).to_csv(args.output/'leave_one_out.csv',index=False)
    pd.concat(prices,ignore_index=True).to_csv(args.output/'prices.csv',index=False)
    pd.concat(forward_tables,ignore_index=True).to_csv(args.output/'forwards.csv',index=False)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sampling',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error('Output exists')
    run(args)
