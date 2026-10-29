"""Same-expiry forward conventions for Black IV; never changes option prices."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .plot_joint_curve_fit import implied_vol


def diagnose(run,output,names):
    output.mkdir(parents=True,exist_ok=False)
    rows=[];summary=[];sources=[]
    for name in names:
        qp=run/f'prices_fresh_pooled_{name}.csv';fp=run/f'forwards_fresh_pooled_{name}.csv'
        sources.extend([qp,fp]);q=pd.read_csv(qp);f=pd.read_csv(fp).set_index('expiry_days')
        q=q[q.instrument=='VIX']
        assert len(q)==175 and set(q.expiry_days)=={27,55,118,153}
        for days,group in q.groupby('expiry_days'):
            fm=float(f.loc[days,'model_forward']);market=float(f.loc[days,'forward_vix_points'])
            group_rows=[]
            for r in group.itertuples():
                inv=lambda price,forward: implied_vol(float(price),float(forward),float(r.strike),float(r.T_years),r.option_type)
                row=dict(candidate=name,expiry_days=int(days),strike=float(r.strike),option_type=r.option_type,
                    K_over_market_F=float(r.K_over_F),market_forward=market,model_forward=fm,
                    model_price_U=float(r.model_price_U),model_se_U=float(r.mc_se_U),
                    market_bid_iv=inv(r.price_bid_U,market),market_mid_iv=inv(r.price_mid_U,market),market_ask_iv=inv(r.price_ask_U,market),
                    model_iv_market_forward=inv(r.model_price_U,market),model_iv_own_forward=inv(r.model_price_U,fm))
                for convention,forward in [('market',market),('own',fm)]:
                    row[f'model_iv_{convention}_lower']=inv(r.model_price_U-2*r.mc_se_U,forward)
                    row[f'model_iv_{convention}_upper']=inv(r.model_price_U+2*r.mc_se_U,forward)
                rows.append(row);group_rows.append(row)
            frame=pd.DataFrame(group_rows);put=frame[frame.option_type=='put'].sort_values('strike').iloc[-1]
            call=frame[frame.option_type=='call'].sort_values('strike').iloc[0]
            summary.append(dict(candidate=name,expiry_days=int(days),model_forward=fm,market_forward=market,
                forward_error=fm-market,nearest_put_strike=put.strike,nearest_call_strike=call.strike,
                nearby_call_minus_put_market_forward=call.model_iv_market_forward-put.model_iv_market_forward,
                nearby_call_minus_put_own_forward=call.model_iv_own_forward-put.model_iv_own_forward,
                market_nearby_call_minus_put_iv=call.market_mid_iv-put.market_mid_iv,
                native_iv_unavailable=int(frame.model_iv_own_forward.isna().sum()),
                caveat='Nearby different strikes; not a same-strike parity test or exact discontinuity estimate.'))
    pd.DataFrame(rows).to_csv(output/'iv_conventions.csv',index=False)
    pd.DataFrame(summary).to_csv(output/'forward_and_nearby_iv.csv',index=False)
    # This is a conditional consistency check only: source premiums themselves
    # were reconstructed from IV using the supplied forward. It cannot discover
    # an independently observed alternative forward from actual cash premiums.
    rawpath=Path('data/processed/options/2026-02-19/quotes_all.csv');sources.append(rawpath)
    raw=pd.read_csv(rawpath);raw=raw[(raw.instrument=='VIX')&(raw.validity=='valid')]
    # Keep the market split: raw opposite-side or held-out quotes are not
    # needed to explain a forward-convention effect in training model prices.
    allowed=set(zip(q.expiry,q.strike,q.option_type))
    raw=raw[[key in allowed for key in zip(raw.expiry,raw.strike,raw.option_type)]]
    parity=[]
    for (expiry,strike),pair in raw.groupby(['expiry','strike']):
        if set(pair.option_type)!={'put','call'}:continue
        call=pair[pair.option_type=='call'].iloc[0];put=pair[pair.option_type=='put'].iloc[0]
        low=strike+call.price_bid_U-put.price_ask_U;high=strike+call.price_ask_U-put.price_bid_U
        parity.append(dict(expiry=int(expiry),strike=strike,supplied_forward=float(call.forward),
                           conditional_parity_forward_lower=low,conditional_parity_forward_upper=high,
                           supplied_forward_inside=bool(low<=call.forward<=high)))
    pd.DataFrame(parity,columns=['expiry','strike','supplied_forward',
        'conditional_parity_forward_lower','conditional_parity_forward_upper',
        'supplied_forward_inside']).to_csv(output/'conditional_market_parity.csv',index=False)
    manifest=dict(source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),*sources]},
        status='Display transformation only; prices, forwards, score and parameters unchanged.',
        native_forward_definition='Model E_Q[VIX_T] at the SAME option expiry T; not at T+30 days.',
        interval_caveat='CSV IV lower/upper endpoints transform price plus/minus 2 price SE with the forward held fixed. Native-forward endpoints omit estimated-forward uncertainty and price/forward covariance; they are not full native-IV uncertainty intervals and are not shaded in the plot.',
        parity_scope='Raw pairs restricted to the exact training option keys. OTM-only membership generally supplies no same-strike call/put pair. Historical diagnostic 001 used the broader raw valid VIX universe; it was never a model fit or selection score.',
        training_same_strike_pairs=len(parity),
        identification='Independent option-implied market forward requires independently quoted premiums and matching settlement/discount conventions. Current raw IVs do not supply those independent premiums.',
        summary=summary)
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(pd.DataFrame(summary)[['candidate','expiry_days','forward_error','nearby_call_minus_put_market_forward','nearby_call_minus_put_own_forward','native_iv_unavailable']].to_string(index=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--candidates',nargs='+',default=['baseline','selected'])
    args=parser.parse_args();diagnose(args.run,args.output,args.candidates)
