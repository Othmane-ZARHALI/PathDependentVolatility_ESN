"""Separate low-strike native-put error from forward-induced call-parity error.

This is a read-only diagnostic of saved prices, not new calibration or repricing.
Model-forward native-put IVs are displayed separately from market-forward call
IVs. Neither convention changes saved prices or repairs the forward mismatch.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .full_readout_fit import HERE,dump
from .plot_joint_curve_fit import implied_vol


def run(out):
    out.mkdir(parents=True,exist_ok=False)
    parent=HERE/'runs/m2_mixed_response_001';folder=parent/'fresh_check'
    paths=[folder/'prices_otm_additive_price.csv',folder/'prices_canonical_additive_price.csv',folder/'forwards_additive_price.csv',Path(__file__)]
    native=pd.read_csv(paths[0]);calls=pd.read_csv(paths[1]);forwards=pd.read_csv(paths[2]);rows=[]
    pd.testing.assert_series_equal(native.strike,calls.strike);pd.testing.assert_series_equal(native.expiry_days,calls.expiry_days)
    for i,r in native[native.option_type.eq('put')].iterrows():
        c=calls.loc[i];f=forwards[forwards.expiry_days.eq(r.expiry_days)].iloc[0]
        pe=r.model_price_U-r.price_mid_U;fe=f.model_forward-f.forward_vix_points;ce=c.model_price_U-c.price_mid_U
        assert abs(ce-pe-fe)<1e-10
        rows.append(dict(day=int(r.expiry_days),strike=r.strike,K_over_F=r.K_over_F,
            market_put=r.price_mid_U,model_put=r.model_price_U,put_SE=r.mc_se_U,
            native_put_error=pe,forward_error=fe,canonical_call_error=ce,
            parity_decomposition_error=ce-pe-fe,
            market_native_put_IV=implied_vol(r.price_mid_U,r.forward,r.strike,r.T_years,'put'),
            model_native_put_IV_model_forward=implied_vol(r.model_price_U,f.model_forward,r.strike,r.T_years,'put'),
            model_canonical_call_IV_market_forward=implied_vol(c.model_price_U,c.forward,c.strike,c.T_years,'call')))
    table=pd.DataFrame(rows);table.to_csv(out/'all_native_put_diagnostics.csv',index=False)
    lowest=table.sort_values(['day','strike']).groupby('day',sort=True).head(1);lowest.to_csv(out/'lowest_strikes.csv',index=False)
    inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (out/'source.py').write_bytes(Path(__file__).read_bytes())
    dump(out/'diagnostic.json',dict(status='PASS',parent=str(parent),input_sha256=inputs,
        identity='C_model-C_market = P_model-P_market + F_model-F_market, undiscounted prices',
        all_put_quotes=len(table),maximum_parity_decomposition_error=float(abs(table.parity_decomposition_error).max()),
        lowest_strikes=lowest.to_dict('records'),
        interpretation='Long-dated lowest-strike puts are underpriced despite a high plotted call IV; shorter lowest puts are overpriced. Enforce forward targets before diagnosing all displayed left ends as excessive lower tails.',
        limitations=['Same archived market targets and saved finite-sample prices.','Diagnostic identities and IV conventions do not produce new consistent model prices.','Lowest-strike conclusions do not characterize every low-strike quote.','No new calibration, floor change, model adoption or SPX computation.']))
    print(lowest.to_string(index=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();run(a.output)
