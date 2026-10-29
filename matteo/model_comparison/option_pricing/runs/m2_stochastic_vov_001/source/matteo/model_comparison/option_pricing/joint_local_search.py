"""Bounded, normalized local search; no full pricer calls inside optimization."""
import time
import numpy as np
from scipy.optimize import minimize
from .calibration_policy import LOWER, UPPER
from .full_readout_fit import SCALE


def propose(model, radius=.05, maxiter=200):
    anchor = model.evaluate(model.parameters)
    normalization = np.maximum(1., anchor['criteria'])
    # A small second-order reserve: unlike the old percentage gate it shrinks
    # quadratically with the radius. The anchor itself need not satisfy it.
    margins = normalization * min(.01, .003*(radius/.05)**2)
    limits = anchor['criteria'] - margins
    lo = np.maximum((LOWER-model.parameters)/SCALE, -radius)
    hi = np.minimum((UPPER-model.parameters)/SCALE, radius)
    history = []
    memo = {}
    best = [np.inf, None, None]
    begin = time.perf_counter()

    def current(z):
        p = model.parameters + SCALE*np.asarray(z)
        key = tuple(p)
        if key not in memo:
            r = model.evaluate(p)
            slack = (limits-r['criteria'])/normalization
            feasible = np.min(slack) >= -1e-9
            if feasible and r['merit'] < best[0]:
                best[:] = [r['merit'], p.copy(), r]
            history.append(dict(parameters=p.tolist(),criteria=r['criteria'].tolist(),
                                minimum_constraint_slack=float(slack.min()),merit=float(r['merit'])))
            memo.clear();memo[key]=r
        return memo[key]

    fit = minimize(lambda z:current(z)['merit']/normalization[3], np.clip(np.zeros(11), lo, hi),
                   jac=lambda z:current(z)['merit_gradient']*SCALE/normalization[3],
                   method='SLSQP', bounds=list(zip(lo,hi)),
                   constraints=[dict(type='ineq',fun=lambda z:(limits-current(z)['criteria'])/normalization,
                                     jac=lambda z:-current(z)['criteria_gradients']*SCALE[None,:]/normalization[:,None])],
                   options=dict(maxiter=maxiter,ftol=1e-11))
    if best[1] is not None:
        # Rounding at an active SLSQP bound must not make a new model infeasible.
        candidate=np.clip(best[1],LOWER,UPPER)
        value=model.evaluate(candidate)
        if np.min((limits-value['criteria'])/normalization)<-1e-8:
            candidate=None;value=None
    else:candidate=None;value=None
    return dict(parameters=candidate, prediction=value, history=history, seconds=time.perf_counter()-begin,
                success=bool(fit.success), message=str(fit.message), nfev=fit.nfev,
                full_pricer_calls=0, radius=radius, reserve_margins=margins.tolist())
