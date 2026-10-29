"""Finite-support VIX marginal feasibility, not a dynamic volatility model.

Fit every archived common-call bid/ask and observed VIX forward, optionally
holding the baseline's squared-VIX moment fixed. This bypasses any chosen SDE,
readout or Monte Carlo estimator. Grid refinements are diagnostics, not proofs
of continuum infeasibility. No SPX pricing or market data edits.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.optimize import linprog

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
source=ROOT/'runs/m2_postcap_vov_001/prices_canonical_baseline.csv'
q=pd.read_csv(source);d=json.loads((ROOT/'runs/m2_maturity_curve_fit_001/run.json').read_text());p=np.array(d['selected']);ends=np.array(d['ends_days']);records=[]
minimum=100*np.sqrt(p[-1])
for step in [.2,.05]:
 for top in [250.,1000.]:
  for day,g in q.groupby('expiry_days',sort=True):
   F=float(g.forward.iloc[0]);K=g.strike.to_numpy();scale=g.price_scale_U.to_numpy();bid=g.price_bid_U.to_numpy();ask=g.price_ask_U.to_numpy()
   y=np.unique(np.r_[minimum,np.arange(minimum,225,step),K,np.geomspace(225,top,100),top]);y=y[y>=minimum]
   P=np.maximum(y[None,:]-K[:,None],0);n=len(y)
   moments=np.maximum(0,np.minimum(ends,day+30)-np.maximum(np.r_[0,ends[:-1]],day));target=float(10000*moments@p[7:-1]/30)
   eq=np.stack([np.ones(n),y/F]);rhs=np.array([1.,1.])
   A=np.r_[P/scale[:,None],-P/scale[:,None]];b=np.r_[ask/scale,-bid/scale]
   solutions={}
   for fixed in [False,True]:
    eql=np.r_[eq,(y*y/target)[None,:]] if fixed else eq;br=np.r_[rhs,1.] if fixed else rhs
    res=linprog(np.r_[np.zeros(n),1.],A_ub=np.c_[A,-np.ones(len(A))],b_ub=b,A_eq=np.c_[eql,np.zeros(len(eql))],b_eq=br,bounds=[(0,None)]*(n+1),method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
    assert res.success,res.message
    probs=res.x[:-1];prices=P@probs;violation=np.maximum(bid-prices,np.maximum(prices-ask,0))/scale
    assert abs(probs.sum()-1)<1e-7 and abs(probs@y-F)<1e-6
    if fixed:assert abs(probs@(y*y)-target)<1e-5
    assert violation.max()<res.x[-1]+1e-6
    solutions['fixed_second' if fixed else 'free_second']=dict(minimum_max_quote_scale_relaxation=float(res.x[-1]),mean=float(probs@y),second=float(probs@(y*y)),minimum_probability=float(probs.min()),max_price_violation=float(violation.max()),active_atoms=int((probs>1e-10).sum()))
    if step==.05 and top==1000:
     role='fixed' if fixed else 'free';pd.DataFrame(dict(VIX=y,probability=probs)).to_csv(HERE/f'marginal_d{day}_{role}.csv',index=False)
     table=g.copy();table['marginal_call_price_U']=prices;table.to_csv(HERE/f'marginal_prices_d{day}_{role}.csv',index=False)
   minimum_second=linprog(y*y/target,A_ub=A,b_ub=b,A_eq=eq,b_eq=rhs,bounds=(0,None),method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
   low=float(minimum_second.fun*target) if minimum_second.success else None
   rec=dict(day=int(day),step=step,maximum_support=top,minimum_support=minimum,atoms=n,quotes=len(g),fixed_second_target=target,minimum_second_within_spreads_on_this_grid=low,**solutions);records.append(rec)
   print(day,'step',step,'top',top,'minimum E2',low,'fixed E2 slack',solutions['fixed_second']['minimum_max_quote_scale_relaxation'],flush=True)
(HERE/'distribution_feasibility.json').write_text(json.dumps(dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),model_parameters_source='m2_maturity_curve_fit_001/run.json selected',records=records,interpretation='Finite positive marginal laws only. Does not establish a jointly consistent time process, reservoir representation, SPX compatibility or unique tail law. No extrapolated second-moment upper bound is inferred.'),indent=2)+'\n')
