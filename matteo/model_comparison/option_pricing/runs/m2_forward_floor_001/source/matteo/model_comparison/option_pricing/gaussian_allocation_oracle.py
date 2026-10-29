"""Known Gaussian expectations for the actual frozen allocation proposals."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import ndtr
from .m2_spx import ROOT,M2
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.stratified_sampling import StratifiedGaussianMixture
from esn_eq.outer_sampling import replication_statistics


def run(args):
    source=Path(__file__).with_name('runs')/'m2_vix_pilot_sampling_001/proposals.json'
    sources=capture_sources([Path(__file__),Path(__file__).with_name('m2_spx.py'),source,
        *sorted((M2/'src/esn_eq').glob('*.py'))])
    declarations=json.loads(source.read_text());rows=[];seeds=[]
    for days,record in declarations.items():
        mixture=DefensiveGaussianMixture(record['shifts'],record['probabilities'])
        stratifier=StratifiedGaussianMixture(mixture)
        for method_id,method in enumerate(('random_mixture','stratified')):
            estimates={p:[] for p in (7,14)}
            for replicate in range(8):
                seed=int(np.random.SeedSequence([20261029,int(days),method_id,replicate]).generate_state(1)[0])
                if method=='stratified':
                    draw=stratifier.draw(14,seed);z,w=draw.normals,draw.weights
                    masks={p:stratifier.prefix_mask(draw,p) for p in estimates}
                else:
                    z,w,_=mixture.draw(14,seed)
                    masks={p:np.arange(len(z))<2**p for p in estimates}
                h=np.column_stack((np.ones(len(z)),z[:,:8],z[:,:8]**2,np.maximum(z[:,0]-3,0)))
                for p,mask in masks.items():
                    estimates[p].append((w[mask,None]*h[mask]).mean(axis=0))
                seeds.append(dict(days=int(days),method=method,replicate=replicate,seed=seed))
            names=['mass']+[f'mean_{j}' for j in range(8)]+[f'second_{j}' for j in range(8)]+['call_z0_K3']
            truths=np.r_[1.,np.zeros(8),np.ones(8),np.exp(-4.5)/np.sqrt(2*np.pi)-3*ndtr(-3.)]
            for p,est in estimates.items():
                mean,se=replication_statistics(est)
                for j,name in enumerate(names):
                    rows.append(dict(days=int(days),method=method,power=p,metric=name,estimate=float(mean[j]),
                        truth=float(truths[j]),replication_se=float(se[j]),error=float(mean[j]-truths[j])))
    check_sources(sources);args.output.mkdir(parents=True,exist_ok=False);save_sources(args.output,sources,ROOT)
    pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False)
    record=dict(status='known-Gaussian oracle diagnostic for frozen actual proposals; no option pricing',
        powers=[7,14],replications=8,seed_records=seeds,implementation_sha256=source_hashes(sources,ROOT),
        limitations=['Finite randomization errors remain; correspondence on simple integrands does not certify nonlinear VIX payoffs.',
            'These mass/moment checks concern proposal likelihood correction, not the model date-zero M_Q normalizer.',
            'All oracle seeds are separate from pricing streams; component-prefix handling matches the pricing design.'])
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    s=pd.DataFrame(rows);print(s[(s.power==14)&(s.metric=='mass')].to_string(index=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Output exists')
    run(args)
