"""Local finite-sample feasibility diagnostic, not a population/model theorem."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from .full_readout_fit import SCALE,dump


def probe(run):
    d=json.loads((run/'design.json').read_text());q=pd.read_csv(run/'prices_canonical_baseline.csv');o=pd.read_csv(run/'prices_otm_baseline.csv');f=pd.read_csv(run/'forwards_baseline.csv')
    a=np.load(run/'jacobian_baseline.npz');P=a['price'];F=a['forward'];r=a['residuals'];J=a['residual'];n=int(q.instrument.eq('SPX').sum());v=q.instrument.eq('VIX').to_numpy();nf=len(f);nv=int(v.sum())
    cjoint=2*r@J;cV=2*r[n:n+nv]@J[n:n+nv];fg=2*r[n+nv:n+nv+nf]@J[n+nv:n+nv+nf]
    N=P.copy()
    for j,day in enumerate(f.expiry_days):
        mask=(o.instrument=='VIX')&(o.expiry_days==day)&(o.option_type=='put');N[mask]-=F[j]
    pn=o.model_price_U.to_numpy();outside=(pn<o.price_bid_U.to_numpy())|(pn>o.price_ask_U.to_numpy());rn=o.interval_residual.to_numpy()*np.sqrt(o.loss_weight.to_numpy())
    NJ=np.sqrt(o.loss_weight.to_numpy())[:,None]*outside[:,None]*N/o.price_scale_U.to_numpy()[:,None];nV=2*rn[v]@NJ[v]
    njoint=cjoint-cV+nV;G=np.stack([cjoint,njoint,cV+fg,nV+fg]);s=json.loads((run/'score_baseline.json').read_text())
    values=np.array([s['canonical_loss'],s['original_otm_loss'],s['canonical_components']['VIX']+s['canonical_components']['forwards'],s['original_otm_components']['VIX']+s['original_otm_components']['forwards']])
    margins=np.maximum(1e-4,np.array([.02,.02,.005,.005])*values)
    p=np.array(d['baseline']);lower=np.maximum((np.array(d['bounds_lower'])-p)/SCALE,-.1);upper=np.minimum((np.array(d['bounds_upper'])-p)/SCALE,.1)
    result=linprog(np.zeros(11),A_ub=G*SCALE,b_ub=-margins,bounds=list(zip(lower,upper)),method='highs')
    record=dict(status='local linear constraint feasibility only',feasible=bool(result.success),message=result.message,criteria_gradients=G.tolist(),margins=margins.tolist(),
                delta_parameters=(SCALE*result.x).tolist() if result.success else None,not_used_for_parameter_selection=True,
                limitation='Linearized finite4096-sample gradients at one starting vector; no model-family feasibility/infeasibility theorem.')
    dump(run/'validation/local_constraint_feasibility.json',record);print(json.dumps({k:v for k,v in record.items() if k!='criteria_gradients'},indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--run',type=Path,required=True);args=a.parse_args();probe(args.run.resolve())
