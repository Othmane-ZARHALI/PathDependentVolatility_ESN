"""Replay saved fits with direct nodal accumulation, fixed controls and exact budgets."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .vix_only_curve_fit import PRIOR,POST,ARCHIVE,ACTIVE,DAYS,DELTA,VIXFit,load_draws,loss_from_prices
from .postcap_vov_core import conditional_multiplier
from .parity_adaptive import gaussian_controls
from .full_readout_fit import dump
from .m2_spx import ROOT


def audit(run):
    design=json.loads((run/'design.json').read_text());checks=[]
    for group in ['source_sha256','input_sha256']:
        for path,h in design[group].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h,path
            if group=='source_sha256':assert hashlib.sha256((run/'source'/path).read_bytes()).hexdigest()==h,path
    assert design['SPX_priced'] is False and design['baseline_replaced'] is False
    d=json.loads((PRIOR/'run.json').read_text());p=np.array(d['selected']);ends=np.array(d['ends_years']);left=np.r_[0,ends[:-1]]
    selected=json.loads((run/'selection.json').read_text());controls=json.loads((PRIOR/'vix_controls.json').read_text())
    market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv')
    assert len(market)==217 and set(market.expiry_days)==set(DAYS) and p[4]>0
    for role,cache,nodepath in [('training',POST/'fresh_cache',POST/'fresh_check/base_conditional_nodes.npz'),('fresh_check',run/'fresh_cache',run/'fresh_conditional_nodes.npz')]:
        nodes=dict(np.load(nodepath));draws=load_draws(cache)
        for name,x in selected['parameters'].items():
            model=json.loads((run/f'experimental_model_{name}.json').read_text());pp=np.array(model['parameters'])
            assert len(pp[7:-1])==16;np.testing.assert_array_equal(pp[:7],p[:7]);assert pp[-1]==p[-1]
            off=np.setdiff1d(np.arange(16),ACTIVE);np.testing.assert_array_equal(pp[7+off],p[7+off])
            assert all(pp[7:-1]>=p[-1]) and max(pp[7:-1])<=.2
            stored=pd.read_csv(run/role/f'prices_otm_{name}.csv');fw=pd.read_csv(run/role/f'forwards_{name}.csv');errors=[]
            with np.load(run/role/f'paths_{name}.npz') as paths:
                for day,g in market.groupby('expiry_days',sort=True):
                    a=draws[day];rate=np.log(2)*365/model['half_life_days'];gamma=model['gamma'];T=day/365
                    z=np.sqrt(-np.expm1(-2*rate*T))*a['z_normals'];I=np.full(len(z),p[-1]*DELTA);weights=nodes[f'd{day}_weights'];times=nodes[f'd{day}_times']
                    for t,w,ratio in zip(times,weights,nodes[f'd{day}_ratios']):
                        seg=min(np.searchsorted(ends,t,side='right'),15)
                        I+=w*(pp[7+seg]-p[-1])*ratio*conditional_multiplier(z,T,t-T,gamma,rate)
                    v=100*np.sqrt(I/DELTA);np.testing.assert_allclose(v,paths[f'd{day}_V'],rtol=5e-14,atol=1e-10)
                    c=gaussian_controls(a['normals'],a['weights']);fixed=controls[str(day)]
                    b=np.array([np.interp(np.log(g.strike),np.log(fixed['anchors']),row) for row in fixed['beta']])
                    call=a['weights'][:,None]*np.maximum(v[:,None]-g.strike.to_numpy(),0)-c@b
                    forward=a['weights']*v-c@np.array(fixed['forward_beta'])
                    native=call-g.option_type.eq('put').to_numpy()[None,:]*(forward[:,None]-g.strike.to_numpy())
                    errors.append(float(np.max(abs(stored.loc[g.index,'model_price_U'].to_numpy()-native.mean(0)))))
                    assert abs(fw.loc[fw.expiry_days.eq(day),'model_forward'].iloc[0]-forward.mean())<1e-10
                    lengths=np.maximum(0,np.minimum(ends,T+DELTA)-np.maximum(left,T))
                    np.testing.assert_allclose(lengths.sum(),DELTA,rtol=1e-14)
                    exact=10000/DELTA*(lengths@pp[7:-1]);score=json.loads((run/role/f'score_{name}.json').read_text());row=next(r for r in score['moments'] if r['day']==day)
                    assert abs(exact-row['exact_second_moment'])<1e-9
            assert max(errors)<1e-9
            checks.append(dict(role=role,candidate=name,maximum_native_price_error=max(errors)))
    result=dict(status='PASS',price_path_replays=checks,source_files=len(design['source_sha256']),input_files=len(design['input_sha256']),all_16_pieces_preserved=True,off_window_pieces_unchanged=True,positive_beta_J=True,SPX_priced=False,baseline_replaced=False,limitations='Finite numerical replay, not population convergence or a global calibration certificate.')
    dump(run/'audit.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();audit(a.run)
