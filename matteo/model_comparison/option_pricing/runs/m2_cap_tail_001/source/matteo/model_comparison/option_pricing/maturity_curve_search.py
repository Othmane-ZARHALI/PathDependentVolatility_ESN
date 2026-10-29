"""Fast local search with arbitrary curve size and hard xi_i>=floor constraints."""
import time
import numpy as np
from scipy.optimize import minimize
from .maturity_curve import bounds,validate


def propose(model,radius=.1,curve_only=False,maxiter=200):
    n=len(model.parameters)-8;lower,upper,scale=bounds(n);base=model.evaluate(model.parameters)
    norm=np.maximum(1.,base['criteria']);limits=base['criteria']-norm*min(.002,.0005*(radius/.1)**2)
    lo=np.maximum((lower-model.parameters)/scale,-radius);hi=np.minimum((upper-model.parameters)/scale,radius)
    if curve_only:lo[:7]=0;hi[:7]=0
    history=[];memo={};best=[np.inf,None,None];begin=time.perf_counter()
    def current(z):
        p=model.parameters+scale*np.asarray(z);key=tuple(p)
        if key not in memo:
            r=model.evaluate(p);slack=(limits-r['criteria'])/norm
            feasible=np.min(slack)>=-1e-8 and np.min(p[7:-1]-p[-1])>=-1e-12
            if feasible and r['merit']<best[0]:best[:]=[r['merit'],p.copy(),r]
            history.append(dict(parameters=p.tolist(),criteria=r['criteria'].tolist(),minimum_slack=float(slack.min())))
            memo.clear();memo[key]=r
        return memo[key]
    relation=np.zeros((n,n+8));relation[:,7:-1]=np.eye(n);relation[:,-1]=-1
    fit=minimize(lambda z:current(z)['merit']/norm[3],np.zeros(len(lo)),jac=lambda z:current(z)['merit_gradient']*scale/norm[3],
                 method='SLSQP',bounds=list(zip(lo,hi)),constraints=[
                    dict(type='ineq',fun=lambda z:(limits-current(z)['criteria'])/norm,jac=lambda z:-current(z)['criteria_gradients']*scale[None,:]/norm[:,None]),
                    dict(type='ineq',fun=lambda z:(relation@(model.parameters+scale*z))/.05,jac=lambda z:relation*scale[None,:]/.05)],
                 options=dict(maxiter=maxiter,ftol=1e-10))
    candidate=best[1];prediction=None
    if candidate is not None:
        candidate=np.clip(candidate,lower,upper);candidate[7:-1]=np.maximum(candidate[7:-1],candidate[-1]);validate(candidate)
        prediction=model.evaluate(candidate)
        if np.min((limits-prediction['criteria'])/norm)<-1e-7:candidate=None;prediction=None
    return dict(parameters=candidate,prediction=prediction,history=history,seconds=time.perf_counter()-begin,nfev=int(fit.nfev),success=bool(fit.success),message=str(fit.message),radius=radius,curve_only=curve_only)
