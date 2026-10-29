"""Replay cutoff-3 panel fits, weights and statistics; retain all failures."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_central_panel import choose_panel,record_trial
from .spx_innovation_study import datasets,train,METHODS
from .spx_volatility_sampling import proposal


LABELS=['stock','reference_payoff','deterministic_delta','spot_delta','mass','ou_linear_C','ou_linear_Y','ou_linear_F']


def audit(output):
    r=json.loads((output/'run.json').read_text());e=pd.read_csv(output/'estimates.csv');costs=pd.read_csv(output/'costs.csv')
    saved_panel=pd.read_csv(output/'contracts.csv');source=output/'source'
    for name,digest in r['implementation_sha256'].items():
        assert hashlib.sha256((source/name).read_bytes()).hexdigest()==digest
    membership=source/'matteo/model_comparison/option_pricing/runs/spx_scaled_moneyness_003/training_membership.csv'
    regenerated=choose_panel(pd.read_csv(membership))
    pd.testing.assert_frame_equal(saved_panel,regenerated,check_dtype=False,rtol=1e-12,atol=1e-12)
    all_seeds=[seed for expiry in r['expiry_records'] for row in expiry['seed_records'] for seed in row['seeds']]
    assert len(all_seeds)==len(set(all_seeds))
    control_rows=[];monitors=[];paired_rows=[];half_rows=[];pilot_rows=[]
    for days,group in saved_panel.groupby('expiry_days'):
        folder=output/f'expiry_{int(days)}';quotes=pd.read_csv(folder/'contracts.csv')
        frozen=json.loads((folder/'frozen.json').read_text())
        def validate(data):
            w=data['w'];assert np.isfinite(w).all() and np.all((w>0)&(w<=2))
            np.testing.assert_allclose(proposal(1).weights(data['z']),w,rtol=1e-12,atol=1e-14)
            y,b,v=datasets(data)['price512'];assert np.all(y>=0)
            assert b.shape==(*y.shape,5) and v.shape==(*y.shape,3)
            np.testing.assert_array_equal(b[:,:,-1],np.broadcast_to(w[:,None]-1.,y.shape))
        with np.load(folder/'pilot_draws.npz') as saved:
            pilot={key:saved[key] for key in saved.files};validate(pilot)
            pilot_y=datasets(pilot)['price512'][0];middle=len(pilot_y)//2
            for j,quote in enumerate(quotes.itertuples()):
                pilot_rows.append(dict(expiry_days=int(days),contract=j,option_type=quote.option_type,strike=quote.strike,band=quote.band,
                    pilot_paths=len(pilot_y),full_hits=int((pilot_y[:,j]>0).sum()),
                    half_A_hits=int((pilot_y[:middle,j]>0).sum()),half_B_hits=int((pilot_y[middle:,j]>0).sum())))
            fitted=train(*datasets(pilot)['price512'])
            for name,value in fitted.items():np.testing.assert_allclose(value,frozen['coefficients'][name],rtol=1e-10,atol=1e-10)
        with np.load(folder/'checking_draws.npz') as saved:
            for trial in range(r['arguments']['trials']):
                prefix=f't{trial}_';data={key[len(prefix):]:saved[key] for key in saved.files if key.startswith(prefix) and not key.endswith('_adjusted')}
                validate(data)
                x,rows,paired=record_trial(data,frozen['coefficients'],quotes,trial)
                np.testing.assert_allclose(x,saved[prefix+'price512_adjusted'],rtol=1e-11,atol=1e-10)
                expected=pd.DataFrame(rows).sort_values(['trial','contract','method']).reset_index(drop=True)
                actual=e[(e.expiry_days==days)&(e.trial==trial)].sort_values(['trial','contract','method']).reset_index(drop=True)
                pd.testing.assert_frame_equal(expected,actual,check_dtype=False,rtol=1e-10,atol=1e-10)
                paired_rows +=[row|dict(expiry_days=int(days)) for row in paired]
                y,b,v=datasets(data)['price512'];controls=np.concatenate([b,v],axis=2)
                for contract in range(len(quotes)):
                    for k,label in enumerate(LABELS):
                        values=controls[:,contract,k];mean=float(values.mean());se=float(values.std(ddof=1)/np.sqrt(len(values)))
                        control_rows.append(dict(expiry_days=int(days),trial=trial,contract=contract,band=int(quotes.band.iloc[contract]),
                            control=label,mean=mean,se=se,z=mean/se if se else None,nonzero_mean_zero_se=bool(se==0 and mean!=0)))
                for k,label,target in [(0,'mass',1.),(1,'stock',1.),(3,'reference',1.)]:
                    values=data['monitors'][:,k];mean=float(values.mean());se=float(values.std(ddof=1)/np.sqrt(len(values)))
                    monitors.append(dict(expiry_days=int(days),trial=trial,statistic=label,mean=mean,se=se,expectation=target,z=(mean-target)/se if se else None))
                for k in range(3):
                    values=data['g3_innovation_moments'][:,k];mean=float(values.mean());se=float(values.std(ddof=1)/np.sqrt(len(values)))
                    monitors.append(dict(expiry_days=int(days),trial=trial,statistic=f'bare_linear_{k}',mean=mean,se=se,expectation=0.,z=mean/se if se else None))
    actual=pd.read_csv(output/'paired_differences.csv').sort_values(['expiry_days','trial','contract']).reset_index(drop=True)
    expected=pd.DataFrame(paired_rows).sort_values(['expiry_days','trial','contract']).reset_index(drop=True)
    pd.testing.assert_frame_equal(actual[expected.columns],expected,check_dtype=False,rtol=1e-10,atol=1e-10)
    summary=[]
    for keys,g in e[e.method.isin(['raw','legacy','base','augmented'])].groupby(['expiry_days','band','method']):
        days,band,method=keys;positive=g.estimate>0
        budget=r['arguments']['paths']*(g.se/(g.scale/8))**2
        summary.append(dict(expiry_days=int(days),band=int(band),method=method,observations=len(g),negative=int((g.estimate<0).sum()),
            sampling_passes=int(((g.estimate>=0)&(g.nonzero_raw>0)&(g.se<=g.scale/8)).sum()),
            min_raw_hits=int(g.nonzero_raw.min()),median_se_over_scale=float((g.se/g.scale).median()),
            median_se_over_price=float((g.loc[positive,'se']/g.loc[positive,'estimate']).median()) if positive.any() else None,
            largest_residual_share=float(g.largest_share_centered_squares.max()),
            median_plugin_paths=float(budget.median()),max_plugin_paths=float(budget.max())))
    for method in ['base_A','base_B','aug_A','aug_B']:
        primary='augmented' if method.startswith('aug') else 'base'
        merged=e[e.method==method].merge(e[e.method==primary],on=['expiry_days','trial','contract'],suffixes=('_half','_full'))
        for row in merged.itertuples():half_rows.append(dict(expiry_days=row.expiry_days,trial=row.trial,contract=row.contract,
            band=row.band_half,method=method,se_ratio=row.se_half/row.se_full if row.se_full>0 else None))
    comparisons=[]
    for (days,trial),g in e.groupby(['expiry_days','trial']):
        times=costs[(costs.expiry_days==days)&(costs.role==1)&(costs.trial==trial)]
        multiplier=float(times.augmented_seconds.sum()/times.ordinary_seconds.sum())
        for band in [None,0,1,2]:
            subset=g if band is None else g[g.band==band]
            b=subset[subset.method=='base'].sort_values('contract');a=subset[subset.method=='augmented'].sort_values('contract')
            variance_gain=float((b.variance/b.scale**2).mean()/(a.variance/a.scale**2).mean())
            comparisons.append(dict(expiry_days=int(days),trial=int(trial),band='all' if band is None else str(band),
                normalized_variance_gain=variance_gain,pipeline_runtime_multiplier=multiplier,
                variance_cost_gain=variance_gain/multiplier,
                cost_note='Expiry-wide six-contract pipeline; band ratios do not represent separately measured band-only costs.'))
    controls=pd.DataFrame(control_rows);monitors=pd.DataFrame(monitors);paired=pd.DataFrame(paired_rows)
    paired['z']=paired.difference/paired.se.replace(0,np.nan)
    independent=[]
    for (days,contract,method),g in e.groupby(['expiry_days','contract','method']):
        if len(g)==2 and method in ['base','augmented']:
            g=g.sort_values('trial');se=float(np.sqrt((g.se**2).sum()));difference=float(g.estimate.iloc[1]-g.estimate.iloc[0])
            independent.append(dict(expiry_days=int(days),contract=int(contract),method=method,difference=difference,se=se,z=difference/se if se else None))
    pd.DataFrame(pilot_rows).to_csv(output/'pilot_payoff_hits.csv',index=False)
    for name,frame in [('summary',pd.DataFrame(summary)),('pilot_sensitivity',pd.DataFrame(half_rows)),('comparisons',pd.DataFrame(comparisons)),
        ('control_means',controls),('monitors',monitors),('independent_trial_differences',pd.DataFrame(independent))]:frame.to_csv(output/(name+'.csv'),index=False)
    result=dict(auditor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),verification='Source snapshots, panel selection, unique seeds, exact weights/mass, independent pilot coefficients and all adjusted estimates/statistics replay.',
        estimates=len(e),paired_differences=len(paired),contract_control_means=len(controls),bare_monitors=len(monitors),
        primary_negative=e[(e.method.isin(['base','augmented']))&(e.estimate<0)].to_dict(orient='records'),
        paired_mean_warnings=paired[paired.z.abs()>2].to_dict(orient='records'),
        control_mean_warning_count=int((controls.z.abs()>2).sum()),bare_monitor_warnings=monitors[monitors.z.abs()>2].to_dict(orient='records'),
        pilot_max_se_ratio=float(pd.DataFrame(half_rows).se_ratio.max()),
        setup_seconds=r['setup_seconds'],wall_seconds=r['wall_seconds'],
        limitations=['All ratios/flags are descriptive, dependent across strikes/methods; small samples do not establish tail coverage.',
            'Required sample counts are plug-in fixed-variance forecasts, not prescribed allocations or accuracy guarantees.',
            'No stock/normalizer refinement, model fitting or clipping; negative estimates remain unusable as prices.'])
    (output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path);audit(parser.parse_args().output)
