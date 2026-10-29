"""Measure concentration of the frozen model's Gaussian variance multiplier."""
from dataclasses import asdict
import argparse
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from scipy.special import ndtri
from scipy.stats import qmc
from .m2_spx import ROOT,M2,frozen_stage3
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.outer_sampling import replication_statistics
from esn_eq.pricing import DatedM2Pricer
from esn_eq.states import symmetric_square_root


def run(args):
    start=time.perf_counter();sources=capture_sources([Path(__file__),Path(__file__).with_name('m2_spx.py'),Path(__file__).with_name('active_scope.json'),*sorted((M2/'src/esn_eq').glob('*.py'))])
    architecture,parameters,premium,hashes=frozen_stage3();dated=DatedM2Pricer(architecture,premium);initial,_=dated.initial_state([0.,0.,0.])
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');scope=json.loads(Path(__file__).with_name('active_scope.json').read_text())
    if scope['vix_fit_expiry_days']!=[27,55,118,153] or scope['vix_deferred_expiry_days']:raise ValueError('Require full scope')
    rows,rep_rows,seeds,refs=[],[],[],[];thresholds=(100.,1000.,10000.,100000.)
    for days in scope['vix_fit_expiry_days']:
        expiry=float(targets.select('VIX',expiry_days=[days]).T_years.iloc[0]);t=expiry+15/365
        mean,cov=dated.gaussian_projection_law(t,initial);root=symmetric_square_root(cov)
        normalizers={order:float(dated.normalizer([t],initial,parameters,order)[0]) for order in (48,64,96,128)}
        refs.append(dict(days=days,time=t,**{f'M{q}':v for q,v in normalizers.items()}))
        estimates=[]
        for replicate in range(8):
            seed=int(np.random.SeedSequence([20261103,days,replicate]).generate_state(1)[0]);u=qmc.Sobol(4,scramble=True,bits=52,rng=seed).random_base2(18)
            z=ndtri(u[:,1:]);z[u[:,0]>.5]*=2.5
            log_ratio=-3*np.log(2.5)+.5*(1-2.5**-2)*np.sum(z*z,axis=1);w=np.exp(-np.logaddexp(-np.log(2.),-np.log(2.)+log_ratio))
            g=dated.multiplier(mean+z@root.T,parameters)
            vector=[w.mean(),np.mean(w*g)]
            for threshold in thresholds:vector.extend([np.mean(w*(g>threshold)),np.mean(w*g*(g>threshold))])
            estimates.append(vector);seeds.append(dict(days=days,replicate=replicate,seed=seed))
            for j,value in enumerate(vector):rep_rows.append(dict(days=days,replicate=replicate,metric=j,estimate=float(value)))
        means,se=replication_statistics(estimates)
        for j,threshold in enumerate(thresholds):
            rows.append(dict(days=days,time=t,threshold=threshold,probability=float(means[2+2*j]),probability_se=float(se[2+2*j]),
                multiplier_contribution=float(means[3+2*j]),contribution_se=float(se[3+2*j]),fraction_of_M128=float(means[3+2*j]/normalizers[128]),
                multiplier_mean=float(means[1]),mean_se=float(se[1]),weight_mean=float(means[0]),M128=normalizers[128]))
    check_sources(sources);args.output.mkdir(parents=True,exist_ok=False);save_sources(args.output,sources,ROOT)
    pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False);pd.DataFrame(rep_rows).to_csv(args.output/'replications.csv',index=False);pd.DataFrame(refs).to_csv(args.output/'normalizer_reference.csv',index=False)
    record=dict(status='fixed-model Gaussian multiplier tail diagnostic, no VIX price/fit',scope=scope,architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),
        initial_state=initial.tolist(),target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,implementation_sha256=source_hashes(sources,ROOT),seed_records=seeds,
        power=18,replications=8,proposal='half standard + half scale-2.5 3D Gaussian, exact density ratio bounded by 2',
        multiplier_upper_bound=float(np.exp(2*parameters.cap_level)),
        limitations=['One representative future date per maturity; not a full path/time law.', 'Tail fraction uses deterministic GH128 reference, not sample-weight normalization.',
            'Finite numerical estimates and sparse references do not certify all prices or parameter configurations.'],wall_seconds=time.perf_counter()-start)
    (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():parser.error('Output exists')
    run(args)
