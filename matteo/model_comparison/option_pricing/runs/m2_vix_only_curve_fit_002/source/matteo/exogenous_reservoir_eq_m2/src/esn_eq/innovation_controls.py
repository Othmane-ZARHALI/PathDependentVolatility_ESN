"""Full-state OU innovation recording and predictable volatility martingales."""
from dataclasses import dataclass
import numpy as np
from scipy.special import expit,ndtr
from .states import joint_transition_covariance


@dataclass(frozen=True)
class ProjectedInnovations:
    times: np.ndarray
    values: np.ndarray  # steps x paths x (C,Y,F), target-centered even under proposal
    covariance: np.ndarray  # steps x 3 x 3
    equity_covariance: np.ndarray  # steps x 3


class InnovationRecorder:
    """Observe full states, retain only innovations on declared nested grids.

    eps=P*(X_next-D*X-b) uses the full multi-rate conditional mean, never a
    Markov approximation to the three projections. Center against the TARGET
    Q drift, not the proposal component shift. Observer inputs are not mutated.
    """
    def __init__(self,pricer,fine_times,grids,paths):
        self.pricer=pricer;self.fine_times=np.array(fine_times,float,copy=True)
        if (self.fine_times.ndim!=1 or len(self.fine_times)<2 or self.fine_times[0]!=0 or
                not np.isfinite(self.fine_times).all() or np.any(np.diff(self.fine_times)<=0) or paths<2):
            raise ValueError('Invalid fine event grid or path count')
        self.entries={};self.next_step=0;self.transitions={}
        for name,values in grids.items():
            times=np.array(values,float,copy=True);indices=np.searchsorted(self.fine_times,times)
            if (times.ndim!=1 or len(times)<2 or times[0]!=0 or times[-1]!=self.fine_times[-1] or
                    np.any(np.diff(times)<=0) or np.any(indices>=len(self.fine_times)) or
                    not np.array_equal(self.fine_times[indices],times)):
                raise ValueError('Innovation grids must be exact nested event subsets')
            covariance=[];cross=[]
            for dt in np.diff(times):
                if dt not in self.transitions:
                    q=joint_transition_covariance(pricer.layout,float(dt));decay=np.exp(-pricer.layout.rates*dt)
                    shift=pricer.q_drift*(-np.expm1(-pricer.layout.rates*dt))/pricer.layout.rates
                    self.transitions[dt]=(decay,shift,pricer.projection@q[1:,1:]@pricer.projection.T,pricer.projection@q[1:,0])
                _,_,r,a=self.transitions[dt];covariance.append(r);cross.append(a)
            self.entries[name]=dict(times=times,indices=indices,index=0,start=None,
                values=np.empty((len(times)-1,paths,3)),covariance=np.asarray(covariance),cross=np.asarray(cross))

    def observe(self,step,before,after):
        if step!=self.next_step:raise ValueError('Observer must be consumed once in order')
        for e in self.entries.values():
            k=e['index']
            if k>=len(e['times'])-1:continue
            if step==e['indices'][k]:e['start']=before
            if step+1==e['indices'][k+1]:
                dt=e['times'][k+1]-e['times'][k];decay,shift,_,_=self.transitions[dt]
                mean=(e['start']*decay+shift)@self.pricer.projection.T
                e['values'][k]=after@self.pricer.projection.T-mean;e['index']+=1;e['start']=None
        self.next_step+=1

    def result(self,name):
        e=self.entries[name]
        if self.next_step!=len(self.fine_times)-1 or e['index']!=len(e['times'])-1:
            raise ValueError('Incomplete innovation recording')
        return ProjectedInnovations(e['times'],e['values'],e['covariance'],e['cross'])


def readout_gradient(projections,parameters):
    """Gradient of log G, without constructing quadratic predictors."""
    c,y,f=np.asarray(projections).T
    score=parameters.cluster_loading*c+parameters.feedback_curvature*(y-parameters.feedback_shift)**2+parameters.spike_curvature*(f-parameters.spike_shift)**2
    slope=expit(parameters.cap_sharpness*(parameters.cap_level-score))
    first=np.stack([np.full_like(c,parameters.cluster_loading),2*parameters.feedback_curvature*(y-parameters.feedback_shift),2*parameters.spike_curvature*(f-parameters.spike_shift)],axis=1)
    return 2*slope[:,None]*first


def readout_derivatives(projections,parameters):
    """Smooth-readout gradient log G and half diagonal Hessian G/G, before exponential guards."""
    c,y,f=np.asarray(projections).T
    score=parameters.cluster_loading*c+parameters.feedback_curvature*(y-parameters.feedback_shift)**2+parameters.spike_curvature*(f-parameters.spike_shift)**2
    slope=expit(parameters.cap_sharpness*(parameters.cap_level-score))
    first=np.stack([np.full_like(c,parameters.cluster_loading),2*parameters.feedback_curvature*(y-parameters.feedback_shift),2*parameters.spike_curvature*(f-parameters.spike_shift)],axis=1)
    second=np.array([0.,2*parameters.feedback_curvature,2*parameters.spike_curvature])
    gradient=2*slope[:,None]*first
    curvature=2*slope[:,None]*second-2*parameters.cap_sharpness*(slope*(1-slope))[:,None]*first**2
    return gradient,.5*(curvature+gradient**2)


@dataclass(frozen=True)
class InnovationHedges:
    terminal_ratios: dict
    stock_hedges: np.ndarray
    volatility_hedges: np.ndarray  # paths x contracts x 3 or 9 (linear/square/mixed)


def spx_innovation_controls(pricer,cache,innovations,parameters,curve,maturities,strikes,kinds,order=48,variance_floor=1e-4,*,linear_only=False):
    """Existing adapted stock hedges plus nine predictable OU-innovation gains.

    Each volatility gain multiplies a start-of-step Black log-variance vega by
    a readout derivative and one centered noise: eps_i, eps_i²-R_ii, or
    eps_i*dW-Cov(eps_i,dW). Their target conditional means are zero, including
    when OU/equity shocks are correlated. Under proposal sampling, multiply the
    complete gains by the full-path likelihood. No stock update is changed. With linear_only, compute only three linear
    gains: no Hessians, squares, mixed noises or nine-column accumulation.
    """
    maturities,strikes=np.asarray(maturities,float),np.asarray(strikes,float);kinds=np.asarray(kinds)
    if (maturities.ndim!=1 or not len(maturities) or strikes.shape!=maturities.shape or kinds.shape!=maturities.shape or
            not np.isfinite(strikes).all() or np.any(strikes<=0) or not np.isin(kinds,['call','put']).all() or
            not set(maturities).issubset(set(cache.times[1:])) or not np.array_equal(cache.times,innovations.times)):
        raise ValueError('Invalid contracts or innovation grid')
    if cache.architecture!=pricer.architecture or cache.premium!=pricer.premium:raise ValueError('Incompatible cache')
    if not np.isfinite(variance_floor) or variance_floor<0:raise ValueError('Invalid variance floor')
    if any(0<t<cache.times[-1] and t not in cache.times for t in curve.ends):raise ValueError('Curve knots must be exact events')
    xi=curve(cache.times[:-1]);normalizer=pricer.normalizer(cache.times[:-1],cache.initial_state,parameters,order)
    if np.any(xi<=variance_floor):raise ValueError('Curve below variance floor')
    variance=variance_floor+(xi[:,None]-variance_floor)*pricer.multiplier(cache.projections,parameters)/normalizer[:,None]
    if not np.isfinite(variance).all() or np.any(variance<=0):raise FloatingPointError('Invalid variance paths')
    dt=np.diff(cache.times);accumulated=np.r_[0.,np.cumsum(xi*dt)];total=accumulated[np.searchsorted(cache.times,maturities)]
    paths=cache.equity_increments.shape[1]
    if innovations.values.shape!=(len(dt),paths,3):raise ValueError('Invalid innovation shape')
    log_stock=np.zeros(paths);hedges=np.zeros((paths,len(strikes),2));vol=np.zeros((paths,len(strikes),3 if linear_only else 9));terminal={}
    for step,t in enumerate(cache.times[:-1]):
        active=np.flatnonzero(maturities>t);current=np.exp(log_stock);dw=cache.equity_increments[step]
        log_next=log_stock-.5*variance[step]*dt[step]+np.sqrt(variance[step])*dw
        following=np.exp(log_next)
        remaining=total[active]-accumulated[step]
        d1=(log_stock[:,None]-np.log(strikes[active])+.5*remaining)/np.sqrt(remaining)
        hedges[:,active,0]+=(ndtr(d1)-(kinds[active]=='put'))*(following-current)[:,None]
        spot=remaining*(variance[step,:,None]/xi[step]);spot_d1=(log_stock[:,None]-np.log(strikes[active])+.5*spot)/np.sqrt(spot)
        hedges[:,active,1]+=(ndtr(spot_d1)-(kinds[active]=='put'))*(following-current)[:,None]
        vega=.5*current[:,None]*np.sqrt(remaining)*np.exp(-.5*d1*d1)/np.sqrt(2*np.pi)
        eps=innovations.values[step]
        if linear_only:
            noise=readout_gradient(cache.projections[step],parameters)*eps
        else:
            gradient,hessian=readout_derivatives(cache.projections[step],parameters)
            squares=eps**2-np.diag(innovations.covariance[step])
            mixed=eps*dw[:,None]-innovations.equity_covariance[step]
            noise=np.c_[gradient*eps,hessian*squares,gradient*mixed]
        vol[:,active,:]+=vega[:,:,None]*noise[:,None,:]
        log_stock=log_next
        if cache.times[step+1] in maturities:terminal[float(cache.times[step+1])]=following.copy()
    if not np.isfinite(vol).all() or not np.isfinite(hedges).all() or not all(np.isfinite(x).all() for x in terminal.values()):
        raise FloatingPointError('Nonfinite innovation control or stock path')
    return InnovationHedges(terminal,hedges,vol)
