"""Freeze Gaussian path-score controls on design data, then replay fresh checks."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_volatility_sampling import ROOT, adjusted
from .spx_path_importance import weighted_controls
from .provenance import capture_sources, save_sources, source_hashes, check_sources


def dictionary(c, w, z):
    """13 zero-mean controls: original five + w*Z + w*(Z²-1)."""
    original=weighted_controls(np.zeros(len(w)),c,w)[1]
    return np.c_[original,w[:,None]*z,w[:,None]*(z*z-1.)]


def ridge_fit(y,c):
    centered=c-c.mean(axis=0); scale=centered.std(axis=0)
    scale=np.where(scale>1e-14,scale,1.)
    x=centered/scale
    return np.linalg.solve(x.T@x+np.eye(x.shape[1]),x.T@(y-y.mean()))/scale


def freeze(parent,output):
    # Enforce staging before checking archives exist; do not inspect partial checks.
    if (parent/'checking_draws.npz').exists() or (parent/'run.json').exists():
        raise ValueError('Freeze controls before completed checking artifacts exist')
    output.mkdir(parents=True,exist_ok=False)
    folder=Path(__file__).resolve().parent
    paths=[Path(__file__).resolve(),folder/'spx_volatility_sampling.py',folder/'spx_path_importance.py',
           parent.resolve()/'selection.json',parent.resolve()/'design_draws.npz']
    captured=capture_sources(paths)
    # Large input archives stay at parent; snapshot their hashes, not duplicate binaries.
    source_only={p:b for p,b in captured.items() if p.suffix!='.npz'}
    save_sources(output,source_only,ROOT)
    selection=json.loads((parent/'selection.json').read_text())
    # Forward/scale values are in the frozen market target loader; parent saves CSV at completion.
    from .m2_spx import ROOT as root
    from .spx_hedging_study import choose_quotes
    from esn_eq.market_targets import MarketTargets
    quotes=choose_quotes(MarketTargets.load(root/'data/processed/options/2026-02-19'))
    fits,design_rows={},[]
    with np.load(parent/'design_draws.npz') as draws:
        for name in ('equity','vol7','vol21'):
            prefix=f'{name}_r0_t0'; y,c,w,z=[draws[prefix+'_'+key] for key in ('y','c','w','z')]
            variants=[]
            for indices in (np.arange(len(w)),np.arange(len(w)//2),np.arange(len(w)//2,len(w))):
                variants.append(np.stack([ridge_fit(w[indices]*y[indices,j],dictionary(c[indices,j],w[indices],z[indices])) for j in range(y.shape[1])]))
            fits[name]=np.stack(variants)
            prefix=f'{name}_r1_t0';y,c,w,z=[draws[prefix+'_'+key] for key in ('y','c','w','z')]
            for variant,label in enumerate(('pooled','halfA','halfB')):
                x=np.stack([(w*y[:,j]-dictionary(c[:,j],w,z)@fits[name][variant,j])*row.forward for j,row in enumerate(quotes.itertuples())],axis=1)
                design_rows.append(dict(method=name,fit=label,variance=x.var(axis=0,ddof=1).tolist(),
                    central_score=float(np.log(np.maximum(x[:,[2,3,8,9]].var(axis=0,ddof=1)/quotes.price_scale_U.to_numpy()[[2,3,8,9]]**2,1e-24)).mean())))
    check_sources(captured)
    record=dict(status='frozen before completed checking archives',parent=str(parent.resolve()),
        frozen_at_utc=datetime.now(timezone.utc).isoformat(),selected=selection['selected'],
        dictionary='five original weighted controls plus four weighted Z and four weighted Z²-1 controls',
        regularization='ridge=1 on centered columns scaled by pilot standard deviation; constant columns use scale 1',
        coefficients={name:value.tolist() for name,value in fits.items()},design=design_rows,
        input_and_source_sha256=source_hashes(captured,ROOT),
        limitations=['No checking-based tuning, proposal reselection or clipping.','Primary proposal remains the original design-selected vol7/equity/vol21 winner.',
            'Ridge and dictionary change together relative to original OLS controls; their effects are not separately identified.'])
    (output/'frozen.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(f'Frozen score controls before completed checking data; primary proposal={record["selected"]}',flush=True)


def evaluate(parent,output):
    record=json.loads((output/'frozen.json').read_text());run=json.loads((parent/'run.json').read_text())
    quotes=pd.read_csv(parent/'contracts.csv');rows=[];archive={}
    for relative,digest in record['input_and_source_sha256'].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()!=digest:
            raise ValueError('Frozen score-control source/input changed')
    with np.load(parent/'checking_draws.npz') as draws:
        for name,fits in record['coefficients'].items():
            fits=np.asarray(fits)
            for trial in range(run['arguments']['trials']):
                prefix=f'{name}_r2_t{trial}';w=draws[prefix+'_w'];z=draws[prefix+'_z']
                for level in range(3):
                    y,c=draws[prefix+f'_g{level}_y'],draws[prefix+f'_g{level}_c']
                    values=[]
                    for j,row in enumerate(quotes.itertuples()):
                        control=dictionary(c[:,j],w,z)
                        x=np.stack([(w*y[:,j]-control@fits[k,j])*row.forward for k in range(3)],axis=1)
                        values.append(x)
                        for k,label in enumerate(('pooled','halfA','halfB')):
                            rows.append(dict(trial=trial,method=name,contract=j,steps=64*2**level,fit=label,
                                price=float(x[:,k].mean()),se=float(x[:,k].std(ddof=1)/np.sqrt(len(x))),scale=row.price_scale_U,
                                positive_payoffs=int((y[:,j]>0).sum())))
                    archive[prefix+f'_g{level}']=np.stack(values,axis=1)
    pd.DataFrame(rows).to_csv(output/'estimates.csv',index=False)
    np.savez_compressed(output/'checking_draws.npz',**archive)
    print(f'Evaluated {len(rows)} frozen score-control estimates; no reselection',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['freeze','evaluate'])
    parser.add_argument('parent',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args()
    (freeze if args.action=='freeze' else evaluate)(args.parent,args.output)
