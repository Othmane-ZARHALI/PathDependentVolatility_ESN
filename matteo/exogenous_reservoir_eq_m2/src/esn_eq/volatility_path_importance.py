"""Orthonormal full-path Gaussian scores targeting equity and OU excursions."""
import numpy as np

from .path_importance import ImportancePaths
from .states import joint_transition_covariance, symmetric_square_root


def volatility_path_basis(pricer, times, anchor):
    """Scores: W_H/sqrt(H), residual C_H, residual Y_anchor, residual F_anchor.

    OU scores retain only innovations conditional on the equity increments.
    The exact OU semigroup maps those innovations to their projected endpoint.
    Gram-Schmidt gives unit independent scores without changing the full path
    dimension. The anchor must be an exact event. Each score defines a
    deterministic Cameron-Martin direction in transition-normal coordinates.
    """
    times = np.asarray(times, float)
    if (times.ndim != 1 or len(times) < 2 or times[0] != 0 or
            not np.isfinite(times).all() or np.any(np.diff(times) <= 0) or
            anchor <= 0 or anchor not in times):
        raise ValueError("Require a finite event grid with a positive exact anchor")
    basis = np.zeros((len(times)-1, len(pricer.layout.rates)+1, 4))
    basis[:, 0, 0] = np.sqrt(np.diff(times)/times[-1])
    roots = {}
    for step, dt in enumerate(np.diff(times)):
        if dt not in roots:
            cov = joint_transition_covariance(pricer.layout, float(dt))
            cross = cov[1:, 0]
            roots[dt] = symmetric_square_root(cov[1:, 1:]-np.outer(cross,cross)/dt)
        for axis, (projection, endpoint) in enumerate(((0,times[-1]),(1,anchor),(2,anchor)), 1):
            if times[step+1] <= endpoint:
                a = pricer.projection[projection]*np.exp(-pricer.layout.rates*(endpoint-times[step+1]))
                basis[step, 1:, axis] = a@roots[dt]
    for axis in range(1, 4):
        for previous in range(axis):
            basis[:,:,axis] -= np.sum(basis[:,:,axis]*basis[:,:,previous])*basis[:,:,previous]
        length = np.linalg.norm(basis[:,:,axis])
        if length < 1e-10:
            raise ValueError("Degenerate conditional OU score")
        basis[:,:,axis] /= length
    basis.setflags(write=False)
    return basis


class GaussianPathTilt:
    """One-use mixture draw, accumulating scores during exact state generation."""
    def __init__(self, basis, mixture, paths, seed):
        self.basis = np.array(basis, float, copy=True)
        if (self.basis.ndim != 3 or not np.isfinite(self.basis).all() or
                self.basis.shape[2] != mixture.shifts.shape[1]):
            raise ValueError("Invalid path basis")
        flat = self.basis.reshape(-1, self.basis.shape[2])
        if not np.allclose(flat.T@flat, np.eye(flat.shape[1]), atol=1e-10, rtol=0):
            raise ValueError("Path score directions must be orthonormal")
        self.basis.setflags(write=False)
        self.components = np.random.default_rng(seed).choice(len(mixture.shifts), size=paths, p=mixture.probabilities)
        self.shifts = mixture.shifts[self.components]
        self.scores = np.zeros((paths, flat.shape[1]))
        self.next_step = 0

    def transform(self, step, normals):
        if step != self.next_step or normals.shape != (len(self.scores), self.basis.shape[1]):
            raise ValueError("Path tilt must be consumed once in order")
        shifted = normals+self.shifts@self.basis[step].T
        self.scores += shifted@self.basis[step]
        self.next_step += 1
        return shifted


def prepare_volatility_paths(pricer, times, initial, paths, seed, basis, mixture, checkpoints=(), *, state_observer=None):
    """Exact f/q weights from ALL orthonormal path scores; no mass rescaling."""
    if np.asarray(basis).shape[:2] != (len(times)-1, len(pricer.layout.rates)+1):
        raise ValueError("Path basis does not match grid/state dimensions")
    component_seed, path_seed = np.random.SeedSequence(seed).spawn(2)
    tilt = GaussianPathTilt(basis, mixture, paths, component_seed)
    options = dict(gaussian_tilt=tilt)
    if state_observer is not None:
        options['state_observer'] = state_observer
    cache = pricer.prepare(times, initial, paths, int(path_seed.generate_state(1)[0]), checkpoints, **options)
    if tilt.next_step != len(times)-1:
        raise ValueError("Incomplete path tilt")
    return ImportancePaths(cache, mixture.weights(tilt.scores), tilt.scores, tilt.components)
