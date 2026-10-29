"""Exact rational payoff certificates for a VIX squared-moment obstruction.

LP proposes a static cash/forward/option portfolio. Fraction arithmetic then
checks its payoff is below VIX^2 on the ENTIRE nonnegative real axis: on each
strike interval the difference is a quadratic, whose minimum is at a boundary
or its stationary point. No grid-only infeasibility claim is made.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.optimize import linprog

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
SOURCE=ROOT/'runs/m2_postcap_vov_001/source/matteo/model_comparison/option_pricing/runs/m2_parity_adaptive_cache_001/vix_contracts.csv'
PARAM=ROOT/'runs/m2_maturity_curve_fit_001/run.json'

def frac(x):return Q(str(x))

def global_minimum(a,b,c,K):
    """Exact minimum of y^2-a-b*y-sum c_i*(y-K_i)+, y>=0."""
    knots=sorted(set([Q(0),*K]));best=None;arg=None
    for i,left in enumerate(knots):
        right=knots[i+1] if i+1<len(knots) else None
        slope=b+sum((cc for cc,kk in zip(c,K) if kk<=left),Q(0));critical=slope/2
        probes=[left]
        if critical>=left and (right is None or critical<=right):probes.append(critical)
        for y in probes:
            value=y*y-a-b*y-sum((cc*max(y-kk,Q(0)) for cc,kk in zip(c,K)),Q(0))
            if best is None or value<best:best,arg=value,y
    return best,arg


def certify(day,g,target,tol):
    F=float(g.forward.iloc[0]);K=g.strike.to_numpy();put=g.option_type.eq('put').to_numpy();sc=g.price_scale_U.to_numpy()
    bid=g.price_bid_U.to_numpy();ask=g.price_ask_U.to_numpy()
    # Read original native OTM prices directly from the archived quote input.
    # Store the declared decimal input values as exact rational numbers.
    FK=[frac(x) for x in K];FB=[frac(x) for x in bid];FA=[frac(x) for x in ask];FF=frac(F);FT=frac(tol)
    y=np.unique(np.r_[0.,np.arange(0,225,.1),K,np.geomspace(225,1000,100)]);history=[]
    for iteration in range(20):
        P=np.where(put[:,None],np.maximum(K[:,None]-y,0),np.maximum(y-K[:,None],0));n=len(y)
        A=np.r_[P/sc[:,None],-P/sc[:,None],(y/F)[None,:],(-y/F)[None,:]];rhs=np.r_[ask/sc,-bid/sc,(F+tol)/F,-(F-tol)/F]
        fit=linprog(y*y/1000,A_ub=A,b_ub=rhs,A_eq=np.ones((1,n)),b_eq=[1.],bounds=(0,None),method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
        assert fit.success,fit.message
        dual=fit.ineqlin.marginals;weights=[frac(x) for x in 1000*(dual[:len(K)]-dual[len(K):2*len(K)])/sc]
        a0=frac(1000*fit.eqlin.marginals[0]);b0=frac(1000*(dual[-2]-dual[-1])/F)
        # Rewrite native puts as calls plus cash/forward, using exact rationals.
        a=a0+sum((w*kk for w,kk,isput in zip(weights,FK,put) if isput),Q(0));b=b0-sum((w for w,isput in zip(weights,put) if isput),Q(0));c=weights
        minimum,arg=global_minimum(a,b,c,FK);history.append(dict(iteration=iteration,LP_second=float(fit.fun*1000),minimum_residual=float(minimum),argmin=float(arg),support_points=n))
        if minimum>=-Q(1,100):break
        y=np.unique(np.r_[y,float(arg)])
    else:raise RuntimeError('Cutting-plane did not converge')
    adjustment=min(minimum,Q(0))-Q(1,1000);a+=adjustment;a0+=adjustment
    checked,where=global_minimum(a,b,c,FK);assert checked>=Q(1,1000)
    cost=a0+b0*(FF-FT if b0>=0 else FF+FT)+sum((w*(bb if w>=0 else aa) for w,bb,aa in zip(weights,FB,FA)),Q(0))
    return dict(day=int(day),forward=float(F),forward_tolerance=float(tol),quotes=len(g),baseline_second_target=str(target),certified_lower_bound=str(cost),certified_lower_bound_float=float(cost),gap_over_baseline_float=float(cost-target),exact_global_minimum_residual=str(checked),payoff_cash=str(a),payoff_forward_coefficient=str(b),native_cash=str(a0),native_forward_coefficient=str(b0),native_option_coefficients=[str(x) for x in weights],call_coefficients=[str(x) for x in c],strikes=[str(x) for x in FK],native_sides=['put' if x else 'call' for x in put],native_bids=[str(x) for x in FB],native_asks=[str(x) for x in FA],history=history)


def verify(record):
    # Re-read ONLY serialized rational numbers, independently of LP primal/dual.
    a=Q(record['payoff_cash']);b=Q(record['payoff_forward_coefficient']);c=list(map(Q,record['call_coefficients']));K=list(map(Q,record['strikes']));w=list(map(Q,record['native_option_coefficients']));isput=[s=='put' for s in record['native_sides']]
    a0=Q(record['native_cash']);b0=Q(record['native_forward_coefficient'])
    assert a==a0+sum((cc*kk for cc,kk,p in zip(w,K,isput) if p),Q(0));assert b==b0-sum((cc for cc,p in zip(w,isput) if p),Q(0));assert w==c
    minimum,_=global_minimum(a,b,c,K);assert minimum>=0
    F=frac(record['forward']);tol=frac(record['forward_tolerance']);bid=list(map(Q,record['native_bids']));ask=list(map(Q,record['native_asks']))
    cost=a0+b0*(F-tol if b0>=0 else F+tol)+sum((cc*(bb if cc>=0 else aa) for cc,bb,aa in zip(w,bid,ask)),Q(0))
    assert cost==Q(record['certified_lower_bound']);assert minimum==Q(record['exact_global_minimum_residual'])
    return dict(day=record['day'],forward_tolerance=record['forward_tolerance'],global_payoff_dominance=True,exact_fraction_lower_bound_verified=True,baseline_excluded=Q(record['baseline_second_target'])<cost)


def main():
    q=pd.read_csv(SOURCE);d=json.loads(PARAM.read_text(),parse_float=str);p=list(map(frac,d['selected']));ends=list(map(frac,d['ends_days']));left=[Q(0),*ends[:-1]];records=[]
    for tol in [0.,.1,.5]:
        for day,g in q.groupby('expiry_days',sort=True):
            target=Q(10000,30)*sum((max(Q(0),min(r,Q(day+30))-max(l,Q(day)))*v for l,r,v in zip(left,ends,p[7:-1])),Q(0))
            rec=certify(day,g,target,tol);records.append(rec);print(day,'forward tolerance',tol,'certified min E[VIX2]',rec['certified_lower_bound_float'],'baseline',float(target),'excluded',rec['gap_over_baseline_float']>0,flush=True)
    result=dict(status='Exact rational static-payoff and worst-spread-price certificates',source=str(SOURCE),source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),parameters=str(PARAM),parameters_sha256=hashlib.sha256(PARAM.read_bytes()).hexdigest(),records=records,assumptions='VIX_T>=0, finite second moment, each original native option expectation within archived native bid/ask, mean VIX in indicated interval around observed forward. Every payoff inequality checked globally, not just on a finite grid. These bounds need not be sharp.',scope='No dynamic model, time consistency, SPX consistency or full market feasibility theorem follows. No Lean certificate.')
    path=HERE/'moment_certificates.json';path.write_text(json.dumps(result,indent=2)+'\n');saved=json.loads(path.read_text());checks=[verify(r) for r in saved['records']];(HERE/'moment_certificate_audit.json').write_text(json.dumps(dict(status='PASS',arithmetic='fractions.Fraction exact rational arithmetic',records=checks),indent=2)+'\n')

if __name__=='__main__':main()
