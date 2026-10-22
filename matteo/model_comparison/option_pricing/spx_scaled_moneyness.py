"""Training-only maturity-scaled SPX coverage; no pricing-performance selection."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from matteo.preprocessing.options import black76_implied_vol


CUTOFFS=(.5,1.,1.5,1.8,2.,2.5,3.,3.5)


def scaled_training_quotes(quotes):
    # Filter before choosing anchors or reading price fields for calculations.
    training=quotes[(quotes.instrument=='SPX')&(quotes.split=='train')].copy()
    if training.empty:raise ValueError('No SPX training quotes')
    if training.duplicated(['expiry','strike','option_type']).any():raise ValueError('Duplicate training contract')
    anchors=[];pieces=[]
    for expiry,group in training.groupby('expiry',sort=True):
        if group.T_years.nunique()!=1 or group.forward.nunique()!=1:raise ValueError('Inconsistent expiry inputs')
        if not np.isfinite(group[['strike','forward','T_years']]).all().all() or (group[['strike','forward','T_years']]<=0).any().any():
            raise ValueError('Invalid moneyness inputs')
        group=group.copy();group['log_moneyness']=np.log(group.strike/group.forward)
        anchor=group.assign(distance=abs(group.log_moneyness)).sort_values(['distance','strike','option_type']).iloc[0]
        vols={side:black76_implied_vol(float(anchor.forward),float(anchor.strike),float(anchor.T_years),
            float(anchor[f'price_{side}_U']),anchor.option_type) for side in ['bid','mid','ask']}
        if not all(np.isfinite(v) and v>0 for v in vols.values()):raise ValueError('Positive anchor IV required')
        if not vols['bid']<=vols['mid']<=vols['ask']:raise ValueError('Unordered anchor IVs')
        row=dict(expiry=int(expiry),expiry_days=int(anchor.expiry_days),T_years=float(anchor.T_years),
            forward=float(anchor.forward),anchor_strike=float(anchor.strike),anchor_option_type=anchor.option_type,
            anchor_log_moneyness=float(anchor.log_moneyness),anchor_source_row_ids=anchor.source_row_ids,
            **{f'iv_{side}':v for side,v in vols.items()})
        anchors.append(row)
        for side,v in vols.items():group[f'z_{side}']=abs(group.log_moneyness)/(v*np.sqrt(group.T_years))
        group['priority_c2']=group.z_mid<=2.
        pieces.append(group)
    return pd.concat(pieces).sort_values(['expiry','strike','option_type']),pd.DataFrame(anchors)


def run(input_path,output):
    output.mkdir(parents=True,exist_ok=False)
    quotes,anchors=scaled_training_quotes(pd.read_csv(input_path));rows=[]
    views=[('all_training',quotes),('short_training_29_57',quotes[quotes.expiry_days.isin([29,57])])]
    views += [(f'expiry_{int(days)}',g) for days,g in quotes.groupby('expiry_days')]
    for name,group in views:
        for cutoff in CUTOFFS:
            counts={side:int((group[f'z_{side}']<=cutoff).sum()) for side in ['bid','mid','ask']}
            rows.append(dict(view=name,cutoff=cutoff,total=len(group),included=counts['mid'],
                coverage=counts['mid']/len(group),included_bid_scale=counts['bid'],included_ask_scale=counts['ask']))
    coverage=pd.DataFrame(rows)
    for i,row in anchors.iterrows():
        width=2*row.iv_mid*np.sqrt(row.T_years)
        anchors.loc[i,'c2_K_over_F_low']=np.exp(-width);anchors.loc[i,'c2_K_over_F_high']=np.exp(width)
    quotes.to_csv(output/'training_membership.csv',index=False);anchors.to_csv(output/'atm_proxies.csv',index=False)
    coverage.to_csv(output/'coverage.csv',index=False)
    sources=[Path(__file__),Path('matteo/preprocessing/options.py'),input_path]
    record=dict(formula='z = abs(log(K/F(T))) / (sigma_ATM_proxy(T)*sqrt(T)); inclusive z <= c',
        anchor_policy='Nearest absolute log-moneyness SPX training quote per expiry; tie by strike then option_type. Invert undiscounted mid using existing Black76. Bid/ask proxy sensitivity retained. No call/put pooling.',
        selected_priority_cutoff=2.,reason='Round central-band cutoff meeting at least 50% training quote coverage at every training expiry; no model output used.',
        training_quotes=len(quotes),included_c2=int(quotes.priority_c2.sum()),
        sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        limitations=['Counts are training-only, quote-count coverage, not volume/open-interest weighting or held-out performance.',
            'ATM is approximated by a nearby training quote, not an exactly ATM quote; input IV/forward conventions remain preprocessing limitations.',
            'Two volatility units defines a central working band, broader than strict near ATM; no guaranteed Gaussian probability or payoff frequency.',
            'Long maturity bands are wider and include contracts not yet simulated; coverage is not numerical accuracy.',
            'All quotes/maturities retained in full scope; priority subset does not remove tail diagnostics or change fitting losses.'])
    (output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(coverage.to_string(index=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path('data/processed/options/2026-02-19/quotes_calibration.csv'))
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();run(args.input,args.output)
