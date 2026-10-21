"""Record normalized-model VIX second-moment residuals without rescaling prices."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT,M2
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily


def run(args):
    sources=capture_sources([Path(__file__),Path(__file__).with_name('m2_spx.py'),*sorted((M2/'src/esn_eq').glob('*.py'))])
    inputs={name:(args.sampling/name).read_bytes() for name in ('run.json','outer_values.csv')}
    record=json.loads(inputs['run.json']);raw=pd.read_csv(args.sampling/'outer_values.csv')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19')
    if record['target_sha256']!=targets.artifact_sha256:
        raise ValueError('Target changed')
    curve=NuisanceCurveFamily.from_targets(targets).build(record['curve_levels']);rows=[];rep_rows=[]
    for days in record['scope']['vix_fit_expiry_days']:
        target_row=targets.vix_forwards[targets.vix_forwards.expiry_days==days].iloc[0]
        T=float(target_row.T_years);delta=30/365
        truth=10000*curve.integral(T,T+delta)/delta
        for method in raw.method.unique():
            frame=raw[(raw.days==days)&(raw.method==method)].copy()
            frame['raw_moment']=frame.weight*frame.vix**2
            # A diagnostic using the exact mass identity. Does not adjust prices.
            frame['centered_residual']=frame.weight*(frame.vix**2-truth)
            per=frame.groupby('replicate').agg(raw_moment=('raw_moment','mean'),centered_residual=('centered_residual','mean'),weight_mean=('weight','mean'))
            R=len(per)
            if R!=record['arguments']['replications']:
                raise ValueError('Incomplete replications')
            rows.append(dict(days=days,method=method,exact_model_target=truth,
                raw_second_moment=float(per.raw_moment.mean()),raw_residual=float(per.raw_moment.mean()-truth),raw_replication_se=float(per.raw_moment.std(ddof=1)/np.sqrt(R)),
                mass_centered_residual=float(per.centered_residual.mean()),centered_replication_se=float(per.centered_residual.std(ddof=1)/np.sqrt(R))))
            for replicate,row in per.iterrows():
                rep_rows.append(dict(days=days,method=method,replicate=int(replicate),exact_model_target=truth,**row.to_dict()))
    check_sources(sources)
    if any((args.sampling/name).read_bytes()!=value for name,value in inputs.items()):
        raise RuntimeError('Input changed')
    args.output.mkdir(parents=True,exist_ok=False);save_sources(args.output,sources,ROOT)
    pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False);pd.DataFrame(rep_rows).to_csv(args.output/'replications.csv',index=False)
    report=dict(status='numerical moment-identity diagnostic; no rescaling or normalization clearance',
        sampling=str(args.sampling),input_sha256={k:hashlib.sha256(v).hexdigest() for k,v in inputs.items()},
        implementation_sha256=source_hashes(sources,ROOT),target_sha256=targets.artifact_sha256,
        exact_model_identity='E[VIX_T^2] = (10000/delta) integral_T^(T+delta) xi(t) dt',
        centered_identity='E_g[w*(VIX^2-target)] = E_f[VIX^2]-target; E_g[w]=1',
        limitations=['Selected after inspecting the study second moments; exploratory, not a predeclared significance test.',
            'Reported SEs from eight scrambles may be unreliable for rare-tail integrands; no Gaussian confidence claim.',
            'A residual can reflect outer sampling, conditional quadrature, global normalizer or time integration errors.',
            'Gaussian likelihood-mass checks do not validate M_Q(t); neither identity is imposed by rescaling model prices.'])
    (args.output/'run.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--sampling',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Output exists')
    run(args)
