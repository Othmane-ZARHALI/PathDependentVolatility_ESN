"""Replay paired proposal/learning experiment and all descriptive warnings."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_call_tilt_study import PROPOSALS,METHODS,fit_controls,adjusted,summarize
from .spx_sampling_design import call_proposal


def pull(saved,prefix):
    return {key[len(prefix):]:saved[key] for key in saved.files if key.startswith(prefix) and not key.endswith('_adjusted')}


def audit(output):
    r=json.loads((output/'run.json').read_text());source=output/'source';e=pd.read_csv(output/'estimates.csv');costs=pd.read_csv(output/'costs.csv')
    for name,digest in r['implementation_sha256'].items():assert hashlib.sha256((source/name).read_bytes()).hexdigest()==digest
    seeds=[seed for expiry in r['expiry_records'] for row in expiry['seed_records'] for seed in row['seeds']];assert len(seeds)==len(set(seeds))
    old_run=json.loads((source/r['parent']/'run.json').read_text())
    old_seeds={seed for expiry in old_run['expiry_records'] for row in expiry['seed_records'] for seed in row['seeds']}
    assert not set(seeds)&old_seeds
    paired=[];learning=[];hits=[];moments=[];control_rows=[];sensitivity=[]
    for record in r['expiry_records']:
        day=record['expiry_days'];folder=output/f'expiry_{day}';quotes=pd.read_csv(folder/'contracts.csv');frozen=json.loads((folder/'frozen.json').read_text())
        old=json.loads((source/r['parent']/f'expiry_{day}/frozen.json').read_text())['coefficients']['base']
        def validate(data,label):
            w=data['w'];assert np.all((w>0)&(w<=2));assert np.isfinite(w).all() and np.all(data['y']>=0)
            np.testing.assert_allclose(w,call_proposal(2 if label=='original' else 3).weights(data['z']),rtol=1e-12,atol=1e-14)
            np.testing.assert_array_equal(data['b'][:,:,-1],np.broadcast_to(w[:,None]-1,data['y'].shape))
        with np.load(folder/'pilot_draws.npz') as saved:
            pilots={label:pull(saved,label+'_') for label in PROPOSALS}
            for label,data in pilots.items():
                validate(data,label);fits=fit_controls(data,old)
                for name,value in fits.items():np.testing.assert_allclose(value,frozen['coefficients'][label][name],rtol=1e-10,atol=1e-10)
                mid=len(data['y'])//2
                for j,q in enumerate(quotes.itertuples()):
                    hits.append(dict(expiry_days=day,proposal=label,contract=j,band=q.band,option_type=q.option_type,strike=q.strike,
                        full_hits=int((data['y'][:,j]>0).sum()),half_A_hits=int((data['y'][:mid,j]>0).sum()),half_B_hits=int((data['y'][mid:,j]>0).sum())))
        with np.load(folder/'checking_draws.npz') as saved:
            for trial in range(r['arguments']['trials']):
                arrays={};data_all={}
                for label in PROPOSALS:
                    prefix=f't{trial}_{label}_';data=pull(saved,prefix);validate(data,label);data_all[label]=data
                    x=adjusted(data,frozen['coefficients'][label],quotes.forward.to_numpy());arrays[label]=x
                    np.testing.assert_allclose(x,saved[prefix+'adjusted'],rtol=1e-11,atol=1e-10)
                    expected=pd.DataFrame(summarize(x,data,quotes,day,trial,label)).sort_values(['contract','method']).reset_index(drop=True)
                    actual=e[(e.expiry_days==day)&(e.trial==trial)&(e.proposal==label)].sort_values(['contract','method']).reset_index(drop=True)
                    pd.testing.assert_frame_equal(expected,actual,check_dtype=False,rtol=1e-10,atol=1e-10)
                    for j,q in enumerate(quotes.itertuples()):
                        delta=x[:,j,1]-x[:,j,4];se=float(delta.std(ddof=1)/np.sqrt(len(delta)))
                        learning.append(dict(expiry_days=day,trial=trial,proposal=label,contract=j,band=q.band,
                            difference=float(delta.mean()),se=se,z=float(delta.mean()/se) if se else None,
                            variance_gain_old128_to_new256=float(x[:,j,4].var(ddof=1)/x[:,j,1].var(ddof=1)) if x[:,j,1].var(ddof=1)>0 else None))
                        for method,k in [('ridge_A',2),('ridge_B',3)]:
                            sensitivity.append(dict(expiry_days=day,trial=trial,proposal=label,contract=j,band=q.band,method=method,
                                se_ratio=float(x[:,j,k].std(ddof=1)/x[:,j,1].std(ddof=1)) if x[:,j,1].std(ddof=1)>0 else None))
                        for k,name in enumerate(['stock','reference_payoff','deterministic_delta','spot_delta','mass']):
                            values=data['b'][:,j,k];mean=float(values.mean());se=float(values.std(ddof=1)/np.sqrt(len(values)))
                            control_rows.append(dict(expiry_days=day,trial=trial,proposal=label,contract=j,band=q.band,control=name,
                                mean=mean,se=se,z=mean/se if se else None))
                    monitor_values=np.c_[data['monitors'],data['w'][:,None]*data['z'],data['w'][:,None]*(data['z']**2-1)]
                    for k,name in enumerate(['mass','stock','reference',*[f'score_{i}' for i in range(4)],*[f'score_square_{i}' for i in range(4)]]):
                        values=monitor_values[:,k];target=1. if k<3 else 0.;mean=float(values.mean());se=float(values.std(ddof=1)/np.sqrt(len(values)))
                        moments.append(dict(expiry_days=day,trial=trial,proposal=label,statistic=name,mean=mean,se=se,expectation=target,z=(mean-target)/se if se else None))
                np.testing.assert_array_equal(data_all['original']['components'],data_all['call3']['components'])
                shift=np.zeros_like(data_all['original']['z']);shift[:,0]=(data_all['original']['components']==1)
                np.testing.assert_allclose(data_all['call3']['z']-data_all['original']['z'],shift,atol=1e-10,rtol=1e-10)
                for j,q in enumerate(quotes.itertuples()):
                    for k,method in enumerate(METHODS):
                        delta=arrays['call3'][:,j,k]-arrays['original'][:,j,k];se=float(delta.std(ddof=1)/np.sqrt(len(delta)))
                        paired.append(dict(expiry_days=day,trial=trial,contract=j,band=q.band,option_type=q.option_type,strike=q.strike,method=method,
                            difference=float(delta.mean()),se=se,z=float(delta.mean()/se) if se else None))
    summary=[]
    for (day,trial,method),g in e.groupby(['expiry_days','trial','method']):
        c=costs[(costs.expiry_days==day)&(costs.trial==trial)&(costs.role==1)].groupby('proposal').seconds.sum()
        multiplier=float(c['call3']/c['original'])
        for band in [None,0,1,2]:
            a=g if band is None else g[g.band==band];old=a[a.proposal=='original'];new=a[a.proposal=='call3']
            gain=float((old.variance/old.scale**2).mean()/(new.variance/new.scale**2).mean())
            summary.append(dict(expiry_days=day,trial=trial,method=method,band='all' if band is None else str(band),variance_gain=gain,
                pipeline_cost_multiplier=multiplier,variance_cost_gain=gain/multiplier))
    for name,values in [('comparisons',summary),('paired_differences',paired),('learning_differences',learning),('pilot_hits',hits),
        ('monitors',moments),('control_means',control_rows),('pilot_sensitivity',sensitivity)]:pd.DataFrame(values).to_csv(output/(name+'.csv'),index=False)
    primary=e[e.method=='ridge'];paired=pd.DataFrame(paired);monitors=pd.DataFrame(moments);control=pd.DataFrame(control_rows)
    stats=[]
    for (label,method),g in e[e.method.isin(['ridge','old128'])].groupby(['proposal','method']):
        stats.append(dict(proposal=label,method=method,negative=int((g.estimate<0).sum()),zero_raw_hits=int((g.nonzero_raw==0).sum()),
            sampling_passes=int(((g.estimate>=0)&(g.nonzero_raw>0)&(g.se<=g.scale/8)).sum()),median_se_over_scale=float((g.se/g.scale).median())))
    result=dict(auditor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),verification='Source hashes, fresh unique seeds, exact weights/mass, paired scores, frozen fits and all saved estimates replay.',
        estimates=len(e),stats=stats,primary_negative=primary[primary.estimate<0].to_dict(orient='records'),
        primary_paired_mean_warnings=paired[(paired.method=='ridge')&(paired.z.abs()>2)].to_dict(orient='records'),
        fixed_beta_paired_mean_warnings=paired[(paired.method=='old128')&(paired.z.abs()>2)].to_dict(orient='records'),
        bare_monitor_warnings=monitors[monitors.z.abs()>2].to_dict(orient='records'),control_mean_warning_count=int((control.z.abs()>2).sum()),
        wall_seconds=r['wall_seconds'],setup_seconds=r['setup_seconds'],
        limitations=['Observed ratios and correlated two-SE flags are descriptive, not independent hypothesis tests or convergence certificates.',
            'Primary compares new own pilots; identical old128 coefficients provide the sampling-law isolation.',
            'All band efficiencies share six-contract pipeline timings; no standalone band timing is measured.',
            'No checking-based selection, clipping, mass rescaling, model fit, timestep or normalizer refinement.'])
    (output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path);audit(parser.parse_args().output)
