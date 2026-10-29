"""Defensive covariance inflation in deterministic readout-score directions."""
from dataclasses import dataclass
import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc


def score_subspace(dated,parameters,sampler,window,rank=4):
    """Leading eigenspace of a window-averaged quadratic-score sensitivity.

    For A_s=P exp(-Lambda*s) L, use A_s' diag(0,b,c) A_s plus
    (A_s' grad(score at mean_t))(...)'. No option payoff or checking draw
    selects directions. This heuristic is not an optimal-variance assertion.
    """
    if not isinstance(rank,int) or not 1<=rank<=len(sampler.mean):
        raise ValueError('Invalid score rank')
    if parameters.echo_loading!=0 or parameters.orthogonal_curvature!=0:
        raise ValueError('Score geometry supports Gaussian three-projection M2 only')
    curvature=np.array([0.,parameters.feedback_curvature,parameters.spike_curvature])
    gram=np.zeros((len(sampler.mean),len(sampler.mean)))
    for t,w in zip(window.times,window.weights):
        A=(dated.projection*np.exp(-dated.layout.rates*(t-window.expiry)))@sampler.root
        # Mean of eta_t implied by the terminal distribution, not a fitted state.
        lag=float(t-window.expiry);decay=np.exp(-dated.layout.rates*lag)
        shift=dated.q_drift*(-np.expm1(-dated.layout.rates*lag))/dated.layout.rates
        mu=dated.projection@(decay*sampler.mean+shift)
        gradient=np.array([parameters.cluster_loading,2*parameters.feedback_curvature*(mu[1]-parameters.feedback_shift),2*parameters.spike_curvature*(mu[2]-parameters.spike_shift)])
        g=A.T@gradient
        gram+=(w/window.delta)*((A.T*curvature)@A+np.outer(g,g))
    values,vectors=np.linalg.eigh((gram+gram.T)/2)
    ordering=np.argsort(-values,kind='stable');values=values[ordering];basis=vectors[:,ordering[:rank]]
    for j in range(rank):basis[:,j]*=np.sign(basis[np.argmax(np.abs(basis[:,j])),j])
    basis.setflags(write=False);values.setflags(write=False)
    return basis,values


@dataclass(frozen=True)
class InflatedGaussianDraw:
    normals: np.ndarray
    weights: np.ndarray
    broad: np.ndarray


class LowRankInflatedGaussian:
    """g=.5 phi_I+.5 phi_{I+(scale^2-1)UU'}, with full dimension retained.

    f/g = 1/[.5+.5 scale^(-rank) exp((1-scale^-2)||zU||^2/2)].
    Rotate Gaussian Sobol coordinates to put the subspace first; orthogonal
    complement coordinates remain standard Gaussian. No self normalization.
    """
    def __init__(self,basis,scale=2.5):
        U=np.array(basis,float,copy=True)
        if (U.ndim!=2 or not 1<=U.shape[1]<=U.shape[0] or not np.isfinite(U).all()
                or not np.allclose(U.T@U,np.eye(U.shape[1]),rtol=0,atol=2e-12)
                or not np.isfinite(scale) or scale<=1):
            raise ValueError('Require orthonormal subspace and finite scale > 1')
        complement=np.linalg.svd(U.T,full_matrices=True)[2][U.shape[1]:].T
        rotation=np.column_stack((U,complement))
        for j in range(rotation.shape[1]):rotation[:,j]*=np.sign(rotation[np.argmax(np.abs(rotation[:,j])),j])
        U.setflags(write=False);rotation.setflags(write=False)
        self.basis,self.rotation,self.scale=U,rotation,float(scale)

    def weights(self,z):
        z=np.asarray(z,float)
        if z.ndim!=2 or z.shape[1]!=len(self.basis) or not np.isfinite(z).all():
            raise ValueError('Invalid full whitened states')
        log_ratio=-self.basis.shape[1]*np.log(self.scale)+.5*(1-self.scale**-2)*np.sum((z@self.basis)**2,axis=1)
        return np.exp(-np.logaddexp(-np.log(2.),-np.log(2.)+log_ratio))

    def draw(self,power,seed):
        if not isinstance(power,int) or power<1:raise ValueError('Require integer power>=1')
        u=qmc.Sobol(len(self.basis)+1,scramble=True,bits=52,rng=seed).random_base2(power)
        if np.any((u<=0)|(u>=1)):raise FloatingPointError('Sobol endpoint')
        coordinates=ndtri(u[:,1:]);broad=u[:,0]>.5
        coordinates[:,:self.basis.shape[1]]*=np.where(broad,self.scale,1.)[:,None]
        z=coordinates@self.rotation.T
        return InflatedGaussianDraw(z,self.weights(z),broad)
