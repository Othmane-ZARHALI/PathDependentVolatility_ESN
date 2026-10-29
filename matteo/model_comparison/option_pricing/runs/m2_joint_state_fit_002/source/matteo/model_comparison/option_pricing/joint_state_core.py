"""Restricted-state repricing on fixed primitive Gaussian draws.

The candidate OU path is the zero-initial-state path plus exp(-lambda*t)*X0.
G, the dated normalizer and conditional VIX exposures MUST be recomputed.
Frozen zero-state anchor controls remain functions of the original primitives.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

from esn_eq.pricing import DatedM2Pricer


def projected_initial_shift(pricer, times, initial):
    times = np.asarray(times, float)
    initial = pricer._state(initial)
    return (np.exp(-times[:, None]*pricer.layout.rates)*initial) @ pricer.projection.T


class AdaptiveDatedNormalizer:
    """Positive log-GH48 interpolation, with explicit sampled-error diagnostics.

    This is a NEW deterministic approximation, not a GH convergence certificate.
    Adaptive midpoints plus independent quarter-grid checks compare to direct GH
    at the SAME order. Zero-state full stock-grid values are separately checked.
    """
    def __init__(self, pricer, parameters, initial, horizon, *, order=48, tolerance=2e-7):
        self.pricer, self.parameters = pricer, parameters
        self.initial = pricer._state(initial)
        self.horizon, self.order, self.tolerance = float(horizon), order, tolerance
        self.tau = 1/float(max(pricer.layout.rates))
        high = np.log1p(horizon/self.tau)
        cache = {}

        def direct(s):
            todo = [float(x) for x in s if float(x) not in cache]
            if todo:
                times = self.tau*np.expm1(np.asarray(todo))
                values = DatedM2Pricer.normalizer(pricer, times, self.initial, parameters, order)
                cache.update(zip(todo, np.log(values)))
            return np.array([cache[float(x)] for x in s])

        knots = np.linspace(0, high, 97)
        for iteration in range(12):
            self.interpolator = CubicSpline(knots, direct(knots), extrapolate=False)
            probes = (knots[:-1]+knots[1:])/2
            errors = np.abs(np.expm1(self.interpolator(probes)-direct(probes)))
            bad = errors > tolerance/2
            if not bad.any():
                break
            knots = np.unique(np.r_[knots, probes[bad]])
            if len(knots) > 4096:
                raise RuntimeError('Declared normalizer knot budget exhausted')
        else:
            raise RuntimeError('Normalizer interpolation did not meet sampled tolerance')
        probes = np.linspace(0, high, 514)[1:-1]
        errors = np.abs(np.expm1(self.interpolator(probes)-direct(probes)))
        if errors.max() > tolerance:
            raise RuntimeError(f'Independent normalizer interpolation screen failed: {errors.max()}')
        self.audit = dict(order=order, tolerance=tolerance, knots=len(knots),
                          direct_evaluations=len(cache), midpoint_max_relative_error=float(np.max(np.abs(np.expm1(
                              self.interpolator((knots[:-1]+knots[1:])/2)-direct((knots[:-1]+knots[1:])/2))))),
                          independent_max_relative_error=float(errors.max()),
                          note='Sampled same-GH-order interpolation error; not GH bias or a uniform bound.')
        self.knots, self.logvalues = knots, direct(knots)

    def __call__(self, times):
        times = np.asarray(times, float)
        if not np.isfinite(times).all() or np.any(times < 0) or np.any(times > self.horizon+1e-12):
            raise ValueError('Normalizer times outside declared support')
        value = np.exp(self.interpolator(np.log1p(np.minimum(times,self.horizon)/self.tau)))
        if not np.isfinite(value).all() or np.any(value <= 0):
            raise FloatingPointError('Invalid interpolated dated normalizer')
        return value


class StateDatedPricer(DatedM2Pricer):
    """Candidate-specific normalizer; full-state Gaussian dynamics are unchanged."""
    def normalizer(self, times, initial, parameters, order):
        if (not np.array_equal(initial,self.state_normalizer.initial)
                or parameters != self.state_normalizer.parameters or order != self.state_normalizer.order):
            raise ValueError('Attempted reuse of another candidate normalizer')
        return self.state_normalizer(times)
