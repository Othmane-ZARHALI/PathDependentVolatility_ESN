"""Observed IV-error diagnostics by VIX strike region; no selection or smoothing."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from .plot_joint_curve_fit import implied_vol


def export(run):
    rows=[]
    for name in ['baseline','selected']:
        q=pd.read_csv(run/f'prices_canonical_{name}.csv');q=q[q.instrument=='VIX'].copy()
        q['market_iv_percent']=[implied_vol(float(t.price_mid_U),float(t.forward),float(t.strike),float(t.T_years),'call') for t in q.itertuples()]
        q['region']=np.where(q.K_over_F<=1,'K/F <= 1',np.where(q.K_over_F<=2,'1 < K/F <= 2','K/F > 2'))
        for (day,region),g in q.groupby(['expiry_days','region']):
            finite=np.isfinite(g.model_iv_percent)&np.isfinite(g.market_iv_percent)
            error=g.loc[finite,'model_iv_percent']-g.loc[finite,'market_iv_percent']
            rows.append(dict(candidate=name,days=int(day),region=region,quotes=len(g),defined_IVs=int(finite.sum()),
                             IV_MAE_percentage_points=float(abs(error).mean()) if len(error) else None,
                             IV_signed_error_percentage_points=float(error.mean()) if len(error) else None,
                             price_interval_MSE=float(np.mean(g.interval_residual**2))))
    pd.DataFrame(rows).to_csv(run/'vix_smile_diagnostics.csv',index=False);print('VIX SMILE DIAGNOSTICS SAVED')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--run',type=Path,required=True);args=a.parse_args();export(args.run.resolve())
