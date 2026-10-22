"""Gaussian terminal OU sampling for conditional VIX outer expectations.

Randomized Sobol uncertainty is measured across independent scrambles, never
as though points within one scramble were independent Monte Carlo draws.
"""
import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc

from .states import joint_transition_covariance


def pca_root(covariance, directions):
    """Full-rank Gaussian factor, ordered by projected variance contribution.

    All eigenvectors are retained, including unused modes. Reordering columns
    changes Sobol's coordinate allocation but preserves the Gaussian covariance.
    Directions are deterministic decayed C/Y/F projections, not fitted payoffs.
    """
    covariance, directions = np.asarray(covariance, float), np.asarray(directions, float)
    if (covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]
            or not np.isfinite(covariance).all() or directions.ndim != 2
            or directions.shape[1] != len(covariance) or not np.isfinite(directions).all()
            or not np.allclose(covariance, covariance.T, atol=1e-14, rtol=1e-12)):
        raise ValueError("Require symmetric finite covariance and compatible directions")
    values, vectors = np.linalg.eigh(covariance)
    if values.min() < -1e-11*max(1., float(values.max())):
        raise ValueError("Covariance is not positive semidefinite")
    # Fix eigenvector signs to make the coordinate convention reproducible.
    signs = np.sign(vectors[np.argmax(np.abs(vectors), axis=0), np.arange(len(values))])
    root = vectors*(signs*np.sqrt(np.maximum(values, 0.)))
    scores = np.sum((directions @ root)**2, axis=0)
    order = np.argsort(-scores, kind="stable")
    return root[:, order], scores[order]


class GaussianOuterSampler:
    """Full conditional X_T distribution from dated X_0, with no time stepping."""
    def __init__(self, dated, expiry, initial_state, delta=30/365):
        if not np.isfinite([expiry, delta]).all() or expiry < 0 or delta <= 0:
            raise ValueError("Require nonnegative expiry and positive finite window")
        initial = dated._state(initial_state)
        rates = dated.layout.rates
        self.mean = initial*np.exp(-rates*expiry)+dated.q_drift*(-np.expm1(-rates*expiry))/rates
        self.covariance = joint_transition_covariance(dated.layout, expiry)[1:, 1:]
        lags = np.r_[0., np.geomspace(min(1/float(rates.max()), delta), delta, 8)]
        directions = np.concatenate([dated.projection*np.exp(-rates*lag) for lag in lags])
        self.root, self.coordinate_scores = pca_root(self.covariance, directions)
        self.direction_lags = lags
        for array in (self.mean, self.covariance, self.root, self.coordinate_scores, self.direction_lags):
            array.setflags(write=False)

    def draw(self, power, seed, method="sobol"):
        """Return 2**power states; fixed seed preserves prefixes on refinement.

        Sobol uses 52-bit digital points and inverse-normal transformation.
        Endpoint draws fail explicitly rather than clipping extreme normals.
        """
        if not isinstance(power, int) or power < 1:
            raise ValueError("Require integer power >= 1")
        if method == "sobol":
            uniforms = qmc.Sobol(len(self.mean), scramble=True, bits=52, rng=seed).random_base2(power)
            if np.any((uniforms <= 0) | (uniforms >= 1)):
                raise FloatingPointError("Inverse-normal Sobol transform encountered an endpoint")
            normals = ndtri(uniforms)
        elif method == "iid":
            normals = np.random.default_rng(seed).standard_normal((2**power, len(self.mean)))
        else:
            raise ValueError("Unknown outer sampling method")
        return self.mean+normals @ self.root.T


def replication_statistics(estimates):
    """Mean and SE of the average across independent complete replications.

    Axis zero is the replication unit. Replication count, not points inside a
    Sobol net, sets the uncertainty degrees of freedom. A zero observed SE
    supplies no guarantee for an unobserved rare-event payoff.
    """
    estimates = np.asarray(estimates, float)
    if estimates.ndim < 1 or len(estimates) < 2 or not np.isfinite(estimates).all():
        raise ValueError("Require at least two finite independent replication estimates")
    return estimates.mean(axis=0), estimates.std(axis=0, ddof=1)/np.sqrt(len(estimates))
