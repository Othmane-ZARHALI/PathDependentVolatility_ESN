"""Low-dimensional population check of the implemented conditional GH rule.

Only E[VIX^2] is computed. Marginal Gaussian decompositions at each future
node bypass full-state outer sampling; they do not define VIX/stock paths.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
from scipy.special import ndtri,logsumexp
from scipy.stats import qmc
from .m2_spx import ROOT,M2,frozen_stage3
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.outer_sampling import GaussianOuterSampler,replication_statistics
from esn_eq.pricing import DatedM2Pricer
from esn_eq.states import symmetric_square_root


def discrete_gaussian(uniforms,order,scale=1.,broad=None):
    nodes,weights=hermgauss(order);nodes=np.sqrt(2)*nodes;weights=weights/np.sqrt(np.pi)
    mass=weights.sum();p=weights/mass;cdf=np.cumsum(p);cdf[-1]=1.
    index=np.searchsorted(cdf,uniforms,side='right')
    if broad is None:
        return nodes[index],float(mass**3)
    kappa=.5*(1-scale**-2)
    logZ=logsumexp(np.log(p)+kappa*nodes**2)
    tilted=np.exp(np.log(p)+kappa*nodes**2-logZ);cdf_tilt=np.cumsum(tilted);cdf_tilt[-1]=1.
    index[broad]=np.searchsorted(cdf_tilt,uniforms[broad],side='right')
    z=nodes[index]
    log_inner_ratio=kappa*np.sum(z*z,axis=1)-3*logZ
    return z,float(mass**3),log_inner_ratio


def composed_integrands(dated,parameters,initial,curve,expiry,times,time_weights,denominators,uniforms,orders=(32,48,64),floor=1e-4,proposal_scale=1.):
    """Average reproduces E_outer[GH_q E(G|X_T)] at each integration node.

    Future conditional means are themselves Gaussian with covariance B=P D
    Sigma_T D P'. Independent conditional innovations have covariance C.
    Replacing standard normals in the innovation with the GH discrete law
    reproduces its finite tensor rule, while keeping outer means continuous.
    """
    u=np.asarray(uniforms,float)
    dimension=6 if proposal_scale==1 else 7
    if u.ndim!=2 or u.shape[1]!=dimension or np.any((u<=0)|(u>=1)):raise ValueError('Invalid open-unit uniforms')
    sampler=GaussianOuterSampler(dated,expiry,initial)
    if proposal_scale==1:
        outer_z=ndtri(u[:,:3]);inner_z=ndtri(u[:,3:])
        discrete=[discrete_gaussian(u[:,3:],order) for order in orders]
        likelihood=[np.ones(len(u)) for _ in range(len(orders)+1)]
    else:
        broad=u[:,0]>.5;outer_z=ndtri(u[:,1:4]);inner_z=ndtri(u[:,4:])
        outer_z[broad]*=proposal_scale;inner_z[broad]*=proposal_scale
        kappa=.5*(1-proposal_scale**-2)
        log_outer_ratio=-3*np.log(proposal_scale)+kappa*np.sum(outer_z*outer_z,axis=1)
        log_inner_ratio=-3*np.log(proposal_scale)+kappa*np.sum(inner_z*inner_z,axis=1)
        likelihood=[np.exp(-np.logaddexp(-np.log(2.),-np.log(2.)+log_outer_ratio+log_inner_ratio))]
        discrete=[]
        for order in orders:
            z,mass,log_ratio=discrete_gaussian(u[:,4:],order,proposal_scale,broad)
            discrete.append((z,mass))
            likelihood.append(np.exp(-np.logaddexp(-np.log(2.),-np.log(2.)+log_outer_ratio+log_ratio)))
    out=np.zeros((len(u),len(denominators),len(orders)+1))
    delta=30/365
    for node,(t,w) in enumerate(zip(times,time_weights)):
        lag=float(t-expiry);decay=np.exp(-dated.layout.rates*lag)
        mean,_=dated.gaussian_projection_law(float(t),initial)
        _,C=dated.gaussian_projection_law(lag,np.zeros(len(initial)))
        B=dated.projection@(decay[:,None]*sampler.covariance*decay[None,:])@dated.projection.T
        center=mean+outer_z@symmetric_square_root(B).T;root=symmetric_square_root(C)
        values=[dated.multiplier(center+inner_z@root.T,parameters)]
        values += [mass*dated.multiplier(center+z@root.T,parameters) for z,mass in discrete]
        g=np.column_stack(values)*np.column_stack(likelihood)
        xi=float(curve(float(t)))
        out += (10000*w/delta)*(floor+(xi-floor)*g[:,None,:]/denominators[:,node,None])
    return out


def run(args):
    start=time.perf_counter();folder=Path(__file__).parent
    prior_path=folder/'runs/m2_vix_stratified_sampling_001/run.json';collapsed=args.collapsed
    source_paths=[Path(__file__),folder/'m2_spx.py',prior_path,collapsed/'run.json',collapsed/'nodes.npz',folder/'active_scope.json',*sorted((M2/'src/esn_eq').glob('*.py'))]
    sources=capture_sources(source_paths)
    previous=json.loads(prior_path.read_text());reference=json.loads((collapsed/'run.json').read_text())
    scope=json.loads((folder/'active_scope.json').read_text())
    if scope['vix_fit_expiry_days']!=[27,55,118,153] or scope['vix_deferred_expiry_days']:raise ValueError('Require full scope')
    targets=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');architecture,parameters,premium,hashes=frozen_stage3()
    dated=DatedM2Pricer(architecture,premium);initial,_=dated.initial_state([0.,0.,0.])
    expected=json.loads(json.dumps(dict(architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),target_sha256=targets.artifact_sha256)))
    for record in (previous,reference):
        if any(record[k]!=v for k,v in expected.items()):raise ValueError('Model or target changed')
    if reference['curve_levels']!=previous['curve_levels']:raise ValueError('Curve changed')
    curve=NuisanceCurveFamily.from_targets(targets).build(previous['curve_levels']);archive=np.load(collapsed/'nodes.npz')
    args.output.mkdir(parents=True,exist_ok=False);save_sources(args.output,sources,ROOT)
    (args.output/'status.json').write_text('{"status":"running"}\n')
    rows,rep_rows,differences,seeds=[],[],[],[]
    for days in scope['vix_fit_expiry_days']:
        expiry=float(targets.select('VIX',expiry_days=[days]).T_years.iloc[0]);truth=10000*curve.integral(expiry,expiry+30/365)/(30/365)
        for time_order in (8,16):
            prefix=f'days{days}_time{time_order}_';times=archive[prefix+'times'];weights=archive[prefix+'weights'];den=archive[prefix+'normalizers']
            matrices={p:[] for p in args.powers}
            for replicate in range(args.replications):
                seed=int(np.random.SeedSequence([args.seed,days,replicate]).generate_state(1)[0]);u=qmc.Sobol(6 if args.proposal_scale==1 else 7,scramble=True,bits=52,rng=seed).random_base2(max(args.powers))
                integrands=composed_integrands(dated,parameters,initial,curve,expiry,times,weights,den,u,proposal_scale=args.proposal_scale)
                for p in args.powers:
                    estimate=integrands[:2**p].mean(axis=0);matrices[p].append(estimate)
                    for j,norm_order in enumerate((48,64)):
                        for k,rule in enumerate(('continuous','GH32','GH48','GH64')):
                            rep_rows.append(dict(days=days,time_order=time_order,normalizer_order=norm_order,rule=rule,power=p,replicate=replicate,estimate=float(estimate[j,k]),target=truth))
                if time_order==8:seeds.append(dict(days=days,replicate=replicate,seed=seed))
            for p,matrix in matrices.items():
                x=np.array(matrix);mean,se=replication_statistics(x)
                diff=x[:,:,1:]-x[:,:,:1];dm,ds=replication_statistics(diff)
                for j,norm_order in enumerate((48,64)):
                    for k,rule in enumerate(('continuous','GH32','GH48','GH64')):
                        rows.append(dict(days=days,time_order=time_order,normalizer_order=norm_order,rule=rule,power=p,estimate=float(mean[j,k]),target=truth,residual=float(mean[j,k]-truth),replication_se=float(se[j,k])))
                        if k:differences.append(dict(days=days,time_order=time_order,normalizer_order=norm_order,rule=rule,power=p,GH_minus_continuous=float(dm[j,k-1]),paired_se=float(ds[j,k-1])))
            pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False)
            print(json.dumps(dict(days=days,time_order=time_order,elapsed=time.perf_counter()-start)),flush=True)
    check_sources(sources)
    report=dict(status='population conditional-GH moment diagnostic, not VIX option prices',arguments=vars(args)|dict(output=str(args.output),collapsed=str(args.collapsed)),scope=scope,
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),initial_state=initial.tolist(),curve_levels=previous['curve_levels'],target_sha256=targets.artifact_sha256,
        implementation_sha256=source_hashes(sources,ROOT),seed_records=seeds,conditional_orders=[32,48,64],normalizer_orders=[48,64],time_orders=[8,16],
        distinction='Continuous outer conditional-mean Gaussian, GH-discretized independent innovation; averaged integrands represent existing finite conditional GH rule.',
        limitations=['Does not construct the joint time law or price VIX futures/options.','Randomized six-dimensional numerical integral, not an exact population oracle.',
            'Correlated rule differences are paired; individual points within a scramble are not iid.','Inflated proposal uses exactly reweighted continuous outer and finite-GH inner laws; weights bounded by 2, never self normalized.',
            'High-order GH node probabilities below digital-uniform resolution are not individually resolved; readout is capped.',
            'No formal proof, fit or physical-feasibility claim.'],wall_seconds=time.perf_counter()-start)
    (args.output/'run.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    pd.DataFrame(rep_rows).to_csv(args.output/'replications.csv',index=False);pd.DataFrame(differences).to_csv(args.output/'conditional_changes.csv',index=False)
    (args.output/'status.json').write_text('{"status":"complete"}\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--collapsed',type=Path,required=True)
    parser.add_argument('--powers',type=int,nargs='+',default=[12,14,16]);parser.add_argument('--replications',type=int,default=8);parser.add_argument('--seed',type=int,default=20261031);parser.add_argument('--proposal-scale',type=float,default=1.)
    args=parser.parse_args()
    if args.output.exists() or args.replications<4 or min(args.powers)<4 or not np.isfinite(args.proposal_scale) or args.proposal_scale<1:parser.error('Require new output, >=4 replications and powers>=4')
    run(args)
