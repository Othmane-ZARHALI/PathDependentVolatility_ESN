"""Endpoint-preserving orthogonal Gaussian reflections for coupled SPX paths."""
from dataclasses import dataclass
import time
import numpy as np
from .pricing import PreparedQPaths
from .path_importance import ImportancePaths
from .volatility_path_importance import GaussianPathTilt
from .states import joint_transition_covariance,symmetric_square_root


def joint_root(pricer,dt):
    covariance=joint_transition_covariance(pricer.layout,float(dt));cross=covariance[1:,0]
    root=np.zeros_like(covariance);root[0,0]=np.sqrt(dt);root[1:,0]=cross/np.sqrt(dt)
    root[1:,1:]=symmetric_square_root(covariance[1:,1:]-np.outer(cross,cross)/dt)
    return root


def endpoint_reflector(root1,root2,decay2):
    """T orthogonal, T²=I, M*T=M for M=[diag(1,decay2)*L1,L2].

    SVD needs no covariance inverse. All returned right-singular rows are
    retained, including a null-row completion if rank deficient. This always
    contains the endpoint row space, so endpoints and Gaussian law are fixed;
    at numerical rank deficiency it may reflect only part of the null space.
    No claim of optimal variance reduction or exact conditional-mean reflection
    in that rank-deficient case is made.
    """
    a=np.r_[1.,np.asarray(decay2,float)]
    if root1.shape!=root2.shape or root1.shape!=(len(a),len(a)):
        raise ValueError('Incompatible joint roots and decay')
    endpoint=np.c_[a[:,None]*root1,root2]
    _,singular,right=np.linalg.svd(endpoint,full_matrices=False)
    reflection=2*right.T@right-np.eye(endpoint.shape[1])
    diagnostics=dict(endpoint_residual=float(np.max(np.abs(endpoint@reflection-endpoint))),
        orthogonal_residual=float(np.max(np.abs(reflection@reflection.T-np.eye(len(reflection))))),
        involution_residual=float(np.max(np.abs(reflection@reflection-np.eye(len(reflection))))),
        smallest_singular=float(singular[-1]),largest_singular=float(singular[0]))
    if max(diagnostics[k] for k in ('endpoint_residual','orthogonal_residual','involution_residual'))>1e-10:
        raise FloatingPointError('Invalid endpoint reflection')
    return reflection,diagnostics


@dataclass(frozen=True)
class BridgePaths:
    original: ImportancePaths
    reflected: PreparedQPaths
    reflected_scores: np.ndarray
    endpoint_error: float
    original_seconds: float
    reflection_seconds: float


class EndpointBridge:
    """Reusable kernels; each draw returns iid pairs, NOT two independent samples.

    The original path is sampled under q, with w=f/q. Its reflected path is
    T*u; target f is invariant under T. Therefore the ORIGINAL weight is used
    for BOTH payoffs: E_q[w(u)*P(T*u)]=E_f[P]. Evaluating f/q at reflected
    scores would be incorrect unless the proposal were also transformed.
    """
    def __init__(self,pricer,times,basis):
        self.pricer=pricer;self.times=np.array(times,float,copy=True);self.basis=np.array(basis,float,copy=True)
        d=len(pricer.layout.rates)+1
        if (self.times.ndim!=1 or len(self.times)<3 or (len(self.times)-1)%2 or
                self.times[0]!=0 or not np.isfinite(self.times).all() or np.any(np.diff(self.times)<=0) or
                self.basis.shape[:2]!=(len(self.times)-1,d)):
            raise ValueError('Require paired fine steps and a matching Gaussian basis')
        self.times.setflags(write=False);self.basis.setflags(write=False)
        self.transitions={};self.reflections={};self.diagnostics=[]
        for dt in np.diff(self.times):
            if dt not in self.transitions:
                decay=np.exp(-pricer.layout.rates*dt)
                shift=pricer.q_drift*(-np.expm1(-pricer.layout.rates*dt))/pricer.layout.rates
                self.transitions[dt]=(joint_root(pricer,dt),decay,shift)
        for first in range(0,len(self.times)-1,2):
            key=tuple(np.diff(self.times[first:first+3]))
            if key not in self.reflections:
                l1,_,_=self.transitions[key[0]];l2,decay,_=self.transitions[key[1]]
                reflection,diagnostic=endpoint_reflector(l1,l2,decay)
                self.reflections[key]=reflection;self.diagnostics.append(diagnostic)

    def draw(self,initial,paths,seed,mixture,checkpoints=()):
        if paths<2:raise ValueError('Need at least two independent pairs')
        p=self.pricer;initial=p._state(initial);times=self.times
        checkpoints=set(map(float,checkpoints))
        if not checkpoints.issubset(set(times)):raise ValueError('Checkpoints must be exact events')
        component_seed,path_seed=np.random.SeedSequence(seed).spawn(2)
        tilt=GaussianPathTilt(self.basis,mixture,paths,component_seed)
        actual_seed=int(path_seed.generate_state(1)[0]);rng=np.random.default_rng(actual_seed)
        state=np.broadcast_to(initial,(paths,len(initial))).copy();anti=state.copy()
        features=np.empty((len(times)-1,paths,3));anti_features=np.empty_like(features)
        equity=np.empty((len(times)-1,paths));anti_equity=np.empty_like(equity)
        saved={0.:state.copy()} if 0. in checkpoints else {};anti_saved={0.:anti.copy()} if 0. in checkpoints else {}
        anti_scores=np.zeros_like(tilt.scores);error=0.;original_seconds=reflection_seconds=0.;d=len(initial)+1
        for first in range(0,len(times)-1,2):
            begin=time.perf_counter()
            u1=tilt.transform(first,rng.standard_normal((paths,d)))
            u2=tilt.transform(first+1,rng.standard_normal((paths,d)))
            original_seconds+=time.perf_counter()-begin
            key=tuple(np.diff(times[first:first+3]))
            begin=time.perf_counter();reflected=np.c_[u1,u2]@self.reflections[key].T
            reflection_seconds+=time.perf_counter()-begin
            for step,u,v in ((first,u1,reflected[:,:d]),(first+1,u2,reflected[:,d:])):
                root,decay,shift=self.transitions[times[step+1]-times[step]]
                begin=time.perf_counter();features[step]=state@p.projection.T
                innovation=u@root.T;equity[step]=innovation[:,0]
                state=state*decay+shift+innovation[:,1:]
                if times[step+1] in checkpoints:saved[float(times[step+1])]=state.copy()
                original_seconds+=time.perf_counter()-begin
                begin=time.perf_counter();anti_features[step]=anti@p.projection.T
                innovation=v@root.T;anti_equity[step]=innovation[:,0]
                anti=anti*decay+shift+innovation[:,1:];anti_scores+=v@self.basis[step]
                if times[step+1] in checkpoints:anti_saved[float(times[step+1])]=anti.copy()
                reflection_seconds+=time.perf_counter()-begin
            error=max(error,float(np.max(np.abs(state-anti))))
        if error>5e-9:raise FloatingPointError(f'Coarse OU endpoints drifted: {error}')
        cache=PreparedQPaths(times,features,equity,initial,saved,actual_seed,p.architecture,p.premium)
        anti_cache=PreparedQPaths(times,anti_features,anti_equity,initial,anti_saved,actual_seed,p.architecture,p.premium)
        original=ImportancePaths(cache,mixture.weights(tilt.scores),tilt.scores,tilt.components)
        return BridgePaths(original,anti_cache,anti_scores,error,original_seconds,reflection_seconds)
