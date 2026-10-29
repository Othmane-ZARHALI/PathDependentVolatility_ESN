"""Independent mean-one OU exponential outside the unchanged capped G.

This module is value-only. It does not implement readout-parameter derivatives.
E[L_u | Z_T] is analytic, so no extra conditional Gaussian quadrature is needed.
"""
import numpy as np
from .stochastic_vov_core import LayerKernel, LayerEngine, outer_z, arr
from .maturity_curve import window_nodes

SPECS=[('baseline',0.,np.log(2)*365/30,30,1.)]+[
    (f'g{gamma:g}_h{half}',gamma,np.log(2)*365/half,half,float(np.exp(gamma)))
    for gamma in [.25,.5,1.] for half in [7,30,90]]


def conditional_multiplier(z, expiry, lag, gamma, rate):
    """Z0=0, dZ=-aZdt+sqrt(2a)dB; q(T)=1-exp(-2aT)."""
    decay=np.exp(-rate*lag);q=-np.expm1(-2*rate*expiry)
    return np.exp(gamma*decay*np.asarray(z)-.5*gamma**2*decay**2*q)


class PostKernel(LayerKernel):
    def __init__(self,library,ends):
        super().__init__(library,ends)
        self.lib.poststock_curve.argtypes=self.lib.zstock_curve.argtypes
        self.lib.zstock_curve=self.lib.poststock_curve
        self.spec=SPECS[0]
        self.configure(0.)


class PostEngine(LayerEngine):
    pass


def exposures(e,k,p,norm,cache,specs=SPECS,time_order=8,order=32,days=(27,55,118,153)):
    result={s[0]:{} for s in specs};draws={};seeds={};node_records={}
    for day in days:
        with np.load(cache/'production'/f'vix_d{day}.npz') as saved:a={key:saved[key] for key in saved.files}
        if 'z_normals' in a:eps=a['z_normals'];records=[]
        else:eps,records=outer_z(cache,day)
        a['z_normals']=eps;draws[day]=a;seeds[str(day)]=records
        n=len(eps);times,weights=window_nodes(day/365,e.ends,e.tau,time_order)
        ratios=[]
        for t in times:
            ratios.append(k.expectation(*e.law(float(t-day/365),a['states'],order),p)[:,0]/np.exp(norm(t)[0]))
        ratios=np.array(ratios);node_records[str(day)]=dict(times=times,weights=weights,ratios=ratios)
        for name,gamma,rate,_,_ in specs:
            z=np.sqrt(-np.expm1(-2*rate*day/365))*eps
            expo=np.zeros((n,len(e.ends)))
            for t,w,ratio in zip(times,weights,ratios):
                j=min(np.searchsorted(e.ends,t,side='right'),len(e.ends)-1)
                expo[:,j]+=w*ratio*conditional_multiplier(z,day/365,t-day/365,gamma,rate)
            result[name][day]=expo
        print('POSTCAP EXPOSURES',day,'nodes',len(times),'GH',order,'time',time_order,flush=True)
    return result,draws,seeds,node_records
