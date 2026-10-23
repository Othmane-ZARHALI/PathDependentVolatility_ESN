"""Positive partial conditional SPX integration under the existing discrete law.

Condition on jointly sampled OU-driving shocks, not merely (C,Y,F), and integrate
an independent primitive equity direction. No importance sampler is accepted.
Covariance reconstruction against the existing numerical sampler is audited.
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import log_ndtr
from .states import joint_transition_covariance,symmetric_square_root
from .pricing import PreparedQPaths


def independent_equity_fraction(layout):
    """Squared equity projection onto the primitive-driver nullspace.

    Rank choice is on the constant Brownian-loading matrix, never on a fine-dt
    OU covariance. No OU state mode is dropped from the simulation.
    """
    loadings=layout.driver_loadings
    # Preserve structural endpoints exactly: sqrt(epsilon) leakage otherwise
    # creates spurious retained equity noise in the independent Black oracle.
    if np.all(loadings[:,0]==0):return 1.
    if np.any((np.abs(loadings[:,0])==1)&np.all(loadings[:,1:]==0,axis=1)):return 0.
    _,singular,vh=np.linalg.svd(loadings,full_matrices=True)
    rank=int(np.sum(singular>1e-12*singular[0])) if len(singular) else 0
    residual=vh[rank:].T@vh[rank:,0]
    fraction=float(residual@residual)
    if np.max(np.abs(loadings@residual),initial=0)>1e-12:raise FloatingPointError('Primitive nullspace mismatch')
    return float(np.clip(fraction,0.,1.))


def partial_transition(pricer,dt):
    """Joint correlated-equity/OU root plus independent equity variance.

    The starting covariance is the one reconstructed from the existing
    prepare() factor, including its PSD roundoff. The new factorization changes
    neither Q drift nor intended covariance. Tiny PSD reconstruction errors are
    retained and checked in standardized covariance units, not silently ignored.
    """
    if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite dt required')
    intended=joint_transition_covariance(pricer.layout,dt);cross=intended[1:,0]
    old_root=symmetric_square_root(intended[1:,1:]-np.outer(cross,cross)/dt)
    implemented=intended.copy();implemented[1:,1:]=old_root@old_root.T+np.outer(cross,cross)/dt
    fraction=independent_equity_fraction(pricer.layout);independent=dt*fraction
    correlated=implemented.copy();correlated[0,0]-=independent
    eigenvalues=np.linalg.eigvalsh(correlated);root=symmetric_square_root(correlated)
    reconstructed=root@root.T;reconstructed[0,0]+=independent
    scales=np.sqrt(np.maximum(np.diag(implemented),np.finfo(float).tiny))
    standardized=float(np.max(abs(reconstructed-implemented)/np.outer(scales,scales)))
    if standardized>2e-9:raise FloatingPointError('Conditional transition covariance does not match legacy sampler')
    audit=dict(dt=float(dt),independent_fraction=fraction,correlated_min_eigenvalue=float(eigenvalues.min()),
        max_abs_legacy_covariance_error=float(np.max(abs(reconstructed-implemented))),max_standardized_covariance_error=standardized,
        max_abs_intended_covariance_error=float(np.max(abs(reconstructed-intended))))
    return root,independent,audit


def _log_difference(a,b):
    """Nonnegative exp(a)-exp(b) for analytically ordered Black terms."""
    if np.any(b>a+1e-9):raise FloatingPointError('Unordered Black tail terms')
    delta=np.minimum(b-a,0.)
    return np.exp(a)*(-np.expm1(delta))


def lognormal_payoffs(mean,variance,strikes,kinds):
    """E[(exp(N(mean,var))-K)+] or put; positive log-domain Black tails.

    Zero variance uses intrinsic payoff. Underflow to zero can still occur for
    numerically unrepresentable tails and is not an accuracy certificate.
    """
    mean,variance,strikes=np.broadcast_arrays(np.asarray(mean,float),np.asarray(variance,float),np.asarray(strikes,float))
    kinds=np.broadcast_to(np.asarray(kinds),mean.shape)
    if (not np.isfinite(mean).all() or not np.isfinite(variance).all() or not np.isfinite(strikes).all() or
            np.any(variance<0) or np.any(strikes<=0) or not np.isin(kinds,['call','put']).all()):raise ValueError('Invalid conditional Black input')
    forward=np.exp(mean+.5*variance)
    if not np.isfinite(forward).all():raise FloatingPointError('Nonfinite conditional forward')
    result=np.array(np.maximum(np.where(kinds=='call',forward-strikes,strikes-forward),0.),copy=True)
    positive=variance>0
    if positive.any():
        m,v,k,kind=mean[positive],variance[positive],strikes[positive],kinds[positive]
        logf=m+.5*v;logk=np.log(k);sd=np.sqrt(v);d2=(m-logk)/sd;d1=d2+sd
        call_otm=logk>=logf
        a=np.where(call_otm,logf+log_ndtr(d1),logk+log_ndtr(-d2))
        b=np.where(call_otm,logk+log_ndtr(d2),logf+log_ndtr(-d1))
        # Finite input can produce -inf log probabilities at extreme standardized tails.
        otm=np.zeros_like(a);finite=np.isfinite(a);otm[finite]=_log_difference(a[finite],b[finite])
        intrinsic=np.where(kind=='call',np.maximum(np.exp(logf)-k,0.),np.maximum(k-np.exp(logf),0.))
        result[positive]=otm+intrinsic
    if not np.isfinite(result).all() or np.any(result<0):raise FloatingPointError('Invalid conditional payoff')
    return result


@dataclass(frozen=True)
class ConditionalSPXPaths:
    raw_payoffs: np.ndarray
    conditional_payoffs: np.ndarray
    raw_stock: np.ndarray
    conditional_stock: np.ndarray
    log_means: np.ndarray
    log_variances: np.ndarray
    transition_audit: tuple
    cache: PreparedQPaths | None


def simulate_conditional_spx(pricer,times,initial,parameters,curve,maturities,strikes,kinds,paths,seed,
        order=48,variance_floor=1e-4,*,keep_cache=False):
    """Streaming target-Q coupled raw and partial-conditional payoffs.

    Given retained Gaussian shocks U, log S_T = m_T + sum sqrt(V_t)*nu*sqrt(dt)*Z_t.
    All V_t are known from the U-driven full state. Thus conditional log variance
    is sum V_t*nu²*dt. This conditions on the retained drivers (possibly more
    information than the state path), not a Markov approximation in projections.
    """
    times=np.asarray(times,float);initial=pricer._state(initial)
    maturities,strikes=np.asarray(maturities,float),np.asarray(strikes,float);kinds=np.asarray(kinds)
    if (times.ndim!=1 or len(times)<2 or times[0]!=0 or not np.isfinite(times).all() or np.any(np.diff(times)<=0) or paths<2 or
            maturities.ndim!=1 or not len(maturities) or strikes.shape!=maturities.shape or kinds.shape!=maturities.shape or
            not set(maturities).issubset(set(times[1:])) or np.any(strikes<=0) or not np.isfinite(strikes).all() or
            not np.isin(kinds,['call','put']).all()):raise ValueError('Invalid conditional simulation inputs')
    if not np.isfinite(variance_floor) or variance_floor<0:raise ValueError('Invalid variance floor')
    if any(0<t<times[-1] and t not in times for t in curve.ends):raise ValueError('Insert exact curve knots')
    xi=curve(times[:-1]);normalizer=pricer.normalizer(times[:-1],initial,parameters,order)
    if np.any(xi<=variance_floor):raise ValueError('Curve below variance floor')
    nstate=len(initial);state=np.broadcast_to(initial,(paths,nstate)).copy();rng=np.random.default_rng(seed)
    raw_log=np.zeros(paths);mean=np.zeros(paths);conditional_var=np.zeros(paths)
    raw_stock=np.empty((paths,len(strikes)));means=np.empty_like(raw_stock);variances=np.empty_like(raw_stock)
    transitions={};audits=[]
    features=np.empty((len(times)-1,paths,3)) if keep_cache else None
    equity=np.empty((len(times)-1,paths)) if keep_cache else None
    for step,dt in enumerate(np.diff(times)):
        if dt not in transitions:
            root,independent,audit=partial_transition(pricer,float(dt));audits.append(audit)
            decay=np.exp(-pricer.layout.rates*dt);shift=pricer.q_drift*(-np.expm1(-pricer.layout.rates*dt))/pricer.layout.rates
            transitions[dt]=(root,independent,decay,shift)
        root,independent,decay,shift=transitions[dt]
        projection=state@pricer.projection.T
        variance=variance_floor+(xi[step]-variance_floor)*pricer.multiplier(projection,parameters)/normalizer[step]
        if not np.isfinite(variance).all() or np.any(variance<=0):raise FloatingPointError('Invalid conditional variance path')
        z=rng.standard_normal((paths,nstate+2));joint=z[:,:nstate+1]@root.T
        dw=joint[:,0]+np.sqrt(independent)*z[:,-1]
        increment=-.5*variance*dt
        raw_log+=increment+np.sqrt(variance)*dw
        mean+=increment+np.sqrt(variance)*joint[:,0]
        conditional_var+=variance*independent
        if keep_cache:features[step]=projection;equity[step]=dw
        state=state*decay+shift+joint[:,1:]
        expiry=times[step+1];columns=np.flatnonzero(maturities==expiry)
        if len(columns):
            raw_stock[:,columns]=np.exp(raw_log)[:,None];means[:,columns]=mean[:,None];variances[:,columns]=conditional_var[:,None]
    sign=np.where(kinds=='call',1.,-1.)
    raw=np.maximum(sign[None,:]*(raw_stock-strikes[None,:]),0.)
    conditional=lognormal_payoffs(means,variances,strikes[None,:],kinds[None,:])
    forwards=np.exp(means+.5*variances)
    if not np.isfinite(raw).all() or not np.isfinite(forwards).all():raise FloatingPointError('Nonfinite conditional simulation output')
    cache=PreparedQPaths(times.copy(),features,equity,initial.copy(),{},seed,pricer.architecture,pricer.premium) if keep_cache else None
    return ConditionalSPXPaths(raw,conditional,raw_stock,forwards,means,variances,tuple(audits),cache)
