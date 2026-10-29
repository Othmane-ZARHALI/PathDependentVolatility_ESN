"""Defensive Gaussian mixtures in the full terminal OU whitened coordinates."""
import numpy as np
from scipy.special import logsumexp, ndtri
from scipy.stats import qmc


class DefensiveGaussianMixture:
    """Sample g(z)=sum pi_k phi(z-m_k), returning exact f/g weights.

    Component zero is unshifted, hence weights are bounded by 1/pi_0.
    Weights are never divided by their sample average. Uncertainty must be
    measured across independent complete scrambles, including mixture draws.
    """
    def __init__(self, shifts, probabilities):
        shifts = np.array(shifts, dtype=float, copy=True)
        probabilities = np.array(probabilities, dtype=float, copy=True)
        if (shifts.ndim != 2 or shifts.shape[1] < 1 or
                probabilities.shape != (len(shifts),) or
                not np.isfinite(shifts).all() or not np.isfinite(probabilities).all() or
                np.any(probabilities <= 0) or not np.isclose(probabilities.sum(), 1., atol=1e-14, rtol=0) or
                np.any(shifts[0] != 0)):
            raise ValueError("Require finite shifts, unshifted component zero and positive probabilities summing to one")
        shifts.setflags(write=False)
        probabilities.setflags(write=False)
        self.shifts, self.probabilities = shifts, probabilities

    def weights(self, normals):
        normals = np.asarray(normals, float)
        if normals.ndim != 2 or normals.shape[1] != self.shifts.shape[1] or not np.isfinite(normals).all():
            raise ValueError("Invalid whitened states")
        log_ratios = normals @ self.shifts.T - .5*np.sum(self.shifts**2, axis=1)
        return np.exp(-logsumexp(log_ratios+np.log(self.probabilities),axis=1))

    def draw(self, power, seed):
        if not isinstance(power, int) or power < 1:
            raise ValueError("Require integer power >= 1")
        u = qmc.Sobol(self.shifts.shape[1]+1,scramble=True,bits=52,rng=seed).random_base2(power)
        if np.any((u <= 0) | (u >= 1)):
            raise FloatingPointError("Sobol endpoint")
        component = np.searchsorted(np.cumsum(self.probabilities),u[:,0],side="right")
        z = ndtri(u[:,1:])+self.shifts[component]
        return z,self.weights(z),component

    @classmethod
    def symmetric_leading(cls, dimension, magnitude=2.):
        """Predeclared shifts +/-2 in each of the first two PCA coordinates.

        No payoff fitting or validation-sample tuning; full dimension retained.
        Half the proposal remains the original Gaussian distribution.
        """
        if not isinstance(dimension,int) or dimension < 1 or not np.isfinite(magnitude) or magnitude <= 0:
            raise ValueError("Invalid dimension or shift")
        modes = min(2,dimension)
        shifts = np.zeros((1+2*modes,dimension))
        for j in range(modes):
            shifts[1+2*j,j],shifts[2+2*j,j] = magnitude,-magnitude
        return cls(shifts,np.r_[.5,np.full(2*modes,.5/(2*modes))])
