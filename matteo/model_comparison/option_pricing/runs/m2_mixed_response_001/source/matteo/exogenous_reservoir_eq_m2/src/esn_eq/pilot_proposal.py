"""Frozen pilot-only Gaussian-mixture design; never called on checking draws."""
import numpy as np
from .importance_sampling import DefensiveGaussianMixture


def fit_payoff_mixture(normals, likelihood_weights, targets, coordinates=8, max_shift=3.):
    """Approximate payoff-weighted pilot laws using two shifts per target.

    Each target is nonnegative. In a fixed leading-coordinate subspace, estimate
    its weighted mean and covariance. The two centers are mean +/- spread along
    the largest covariance eigenvector, with spread sqrt(max(eigenvalue-1,.25)).
    Norms are capped at max_shift; other coordinates remain zero. This is a
    heuristic proposal, not a fitted pricing model or an optimal-variance claim.
    Half the probability remains unshifted. Pilot mass normalization is solely
    for estimating proposal moments, never for evaluation likelihood weights.
    """
    z,w,h = np.asarray(normals,float),np.asarray(likelihood_weights,float),np.asarray(targets,float)
    if (z.ndim != 2 or w.shape != (len(z),) or h.ndim != 2 or h.shape[0] != len(z)
            or h.shape[1]<1 or not np.isfinite(z).all() or not np.isfinite(w).all()
            or not np.isfinite(h).all() or np.any(w<=0) or np.any(h<0)
            or not isinstance(coordinates,int) or coordinates<1
            or not np.isfinite(max_shift) or max_shift<=0):
        raise ValueError('Invalid independent pilot arrays or proposal bounds')
    d = min(coordinates,z.shape[1])
    x = z[:,:d]
    shifts,records = [np.zeros(z.shape[1])],[]
    for column in range(h.shape[1]):
        mass = w*h[:,column]
        total = mass.sum()
        if not np.isfinite(total) or total<=0:
            raise ValueError('Pilot target has no finite positive mass; cannot learn its tail')
        p = mass/total
        mean = p @ x
        centered = x-mean
        cov = (centered.T*p) @ centered
        values,vectors = np.linalg.eigh(cov)
        direction = vectors[:,-1]
        # Fix the eigenvector sign for stable serialized proposals.
        direction *= np.sign(direction[np.argmax(np.abs(direction))])
        spread = np.sqrt(max(float(values[-1])-1.,.25))
        unclipped = []
        for sign in (-1.,1.):
            center = mean+sign*spread*direction
            norm = float(np.linalg.norm(center))
            unclipped.append(norm)
            center *= min(1.,max_shift/norm) if norm else 1.
            shift = np.zeros(z.shape[1]); shift[:d] = center
            shifts.append(shift)
        records.append(dict(target=column,positive_pilot_count=int(np.sum(mass>0)),
            contribution_ess=float(total**2/np.sum(mass**2)),max_pilot_share=float(p.max()),
            mean=mean.tolist(),largest_covariance_eigenvalue=float(values[-1]),
            spread=float(spread),unclipped_shift_norms=unclipped))
    probabilities = np.r_[.5,np.full(len(shifts)-1,.5/(len(shifts)-1))]
    return DefensiveGaussianMixture(shifts,probabilities),records
