"""Reusable conditional VIX exposures for piecewise-constant curve profiling.

At fixed readout, dated state and Q dynamics, the date-zero normalizer does
not depend on xi. Conditional integrated variance is affine in its segment
levels. This cache preserves the Gaussian/time approximations used to build
it; it supplies no additional numerical accuracy guarantee.
"""
from dataclasses import dataclass

import numpy as np

from .vix_pricing import ConditionalVIXResult


@dataclass(frozen=True)
class VIXCurveExposures:
    coefficients: np.ndarray  # (outer states, curve segments), units of years
    floor_integral: float
    ends: tuple[float, ...]
    delta: float
    order: int
    owner: object
    configuration: tuple

    def evaluate(self, levels):
        if self.owner._configuration() != self.configuration:
            raise ValueError("Curve exposures belong to a changed dated model")
        levels = np.asarray(levels,float)
        floor = self.owner.variance_floor
        if (levels.shape != (len(self.ends),) or not np.isfinite(levels).all()
                or np.any(levels <= floor)):
            raise ValueError("Require one finite level above the floor per curve segment")
        integrated = self.floor_integral+self.coefficients @ (levels-floor)
        if not np.isfinite(integrated).all() or np.any(integrated <= 0):
            raise FloatingPointError("Invalid cached conditional variance integral")
        return ConditionalVIXResult(100*np.sqrt(integrated/self.delta),integrated,
                                    "cached_gaussian_quadrature",None,None)


def prepare_curve_exposures(pricer, states, window, order=32):
    """Integrate E[G|X_T]/M0 separately on each declared curve segment.

    The window must have been constructed by this pricer using these same
    segment boundaries. Candidate levels change through evaluate(), without
    mutating the template pricer or resetting its date-zero normalizer.
    """
    if window.owner is not pricer or window.configuration != pricer._configuration():
        raise ValueError("Window belongs to a different date-zero model")
    states = np.asarray(states,float)
    if (states.ndim != 2 or states.shape[1] != len(pricer.initial_state)
            or len(states) < 1 or not np.isfinite(states).all()
            or not isinstance(order,int) or order < 2):
        raise ValueError("Require finite full states and integer quadrature order >=2")
    ends = pricer.curve.ends
    segments = np.minimum(np.searchsorted(ends,window.times,side="right"),len(ends)-1)
    coefficients = np.zeros((len(states),len(ends)))
    for time,weight,normalizer,segment in zip(window.times,window.weights,window.normalizer,segments):
        mean,cov = pricer.dated.gaussian_projection_law(float(time-window.expiry),states)
        expected = pricer.dated.gaussian_expectation(mean,cov,pricer.parameters,order)
        coefficients[:,segment] += weight*expected/normalizer
    if not np.isfinite(coefficients).all() or np.any(coefficients < 0):
        raise FloatingPointError("Invalid curve exposure coefficients")
    coefficients.setflags(write=False)
    return VIXCurveExposures(coefficients,float(pricer.variance_floor*window.weights.sum()),
                             ends,window.delta,order,pricer,pricer._configuration())
