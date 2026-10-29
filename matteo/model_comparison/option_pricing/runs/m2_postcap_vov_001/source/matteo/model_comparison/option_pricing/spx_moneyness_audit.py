"""Retrospective moneyness and positivity audit; no estimator selection/refitting."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_innovation_study import METHODS


def audit(source, output):
    output.mkdir(parents=True, exist_ok=False)
    files=['run.json','contracts.csv','estimates.csv','paired_differences.csv','checking_draws.npz']
    hashes={name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in files}
    record=json.loads((source/'run.json').read_text())
    quotes=pd.read_csv(source/'contracts.csv').reset_index(names='contract')
    quotes['forward_distance']=abs(quotes.K_over_F-1)
    quotes['group']=np.where(quotes.forward_distance<=.05,'near_ATM_within_5pct',
        np.where(quotes.forward_distance<=.25,'OTM_5_to_25pct','far_OTM_over_25pct'))
    quotes['OTM']=np.where(quotes.option_type=='call',quotes.K_over_F>1,quotes.K_over_F<1)
    estimates=pd.read_csv(source/'estimates.csv')
    estimates=estimates[estimates.target=='price512'].merge(
        quotes[['contract','K_over_F','forward_distance','group','OTM']],on='contract',validate='many_to_one')
    estimates['se_over_price']=estimates.se/estimates.estimate.abs().replace(0,np.nan)
    estimates['se_over_scale']=estimates.se/estimates.scale
    assert np.isfinite(estimates[['estimate','se']]).all().all()
    halves=[]
    for method in ['base_A','base_B','aug_A','aug_B']:
        primary='augmented' if method.startswith('aug') else 'base'
        merged=estimates[estimates.method==method].merge(estimates[estimates.method==primary],
            on=['trial','target','contract'],suffixes=('_half','_full'),validate='one_to_one')
        for row in merged.itertuples():
            halves.append(dict(trial=row.trial,contract=row.contract,group=row.group_half,method=method,
                se_ratio=row.se_half/row.se_full))
    tail=[]
    with np.load(source/'checking_draws.npz') as draws:
        for trial in range(record['arguments']['trials']):
            adjusted=draws[f't{trial}_price512_adjusted']
            raw=draws[f't{trial}_g3_y']*quotes.forward.to_numpy()[None,:]
            # Positive likelihood and positive-part payoff give pathwise positivity.
            assert np.isfinite(raw).all() and np.all(raw>=0)
            np.testing.assert_allclose(raw,adjusted[:,:,METHODS.index('raw')],rtol=1e-12,atol=1e-12)
            for contract in quotes.contract:
                for method in ['raw','base','augmented']:
                    values=adjusted[:,contract,METHODS.index(method)]
                    centered=values-values.mean();squares=centered**2;total=squares.sum()
                    tail.append(dict(trial=trial,contract=int(contract),group=quotes.loc[contract,'group'],method=method,
                        largest_share_centered_squares=float(squares.max()/total) if total else 0.,
                        largest_path=int(np.argmax(squares)),negative_path_contributions=int((values<0).sum())))
    estimates.to_csv(output/'estimates_by_moneyness.csv',index=False)
    pd.DataFrame(halves).to_csv(output/'pilot_sensitivity.csv',index=False)
    pd.DataFrame(tail).to_csv(output/'path_concentration.csv',index=False)
    primary=estimates[estimates.method.isin(['raw','legacy','base','augmented'])]
    summaries=[]
    for (group,method),g in primary.groupby(['group','method']):
        summaries.append(dict(group=group,method=method,count=len(g),negative=int((g.estimate<0).sum()),
            min_hits=int(g.nonzero_raw.min()),max_hits=int(g.nonzero_raw.max()),
            min_se_over_price=float(g.se_over_price.min()),max_se_over_price=float(g.se_over_price.max()),
            min_se_over_scale=float(g.se_over_scale.min()),max_se_over_scale=float(g.se_over_scale.max()),
            sampling_passes=int(((g.nonzero_raw>0)&(g.se<=g.scale/8)).sum())))
    result=dict(source=str(source),input_sha256=hashes,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        all_representatives_OTM=bool(quotes.OTM.all()),summaries=summaries,
        negative_estimates=estimates[estimates.estimate<0].to_dict(orient='records'),
        limitations=['Retrospective descriptive audit of existing independent checks; no new samples, fitting, clipping or estimator selection.',
            '5% and 25% forward-distance bands are reporting conventions, not adopted pricing/scope rules.',
            'Positive raw weighted payoffs guarantee nonnegative raw estimates; linear zero-mean controls do not preserve this pathwise.',
            'Positivity and pilot stability do not certify tail accuracy, timestep/normalizer convergence or market fit.',
            'All original SPX scope and all four VIX maturities retained.'])
    (output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(all_representatives_OTM=result['all_representatives_OTM'],summaries=summaries),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();audit(args.source,args.output)
