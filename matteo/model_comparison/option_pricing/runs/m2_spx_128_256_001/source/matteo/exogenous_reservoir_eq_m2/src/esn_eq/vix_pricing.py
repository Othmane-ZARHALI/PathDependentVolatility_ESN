"""Conditional continuous-window VIX valuation for Gaussian baseline M2.

Both engines average future variance before taking its square root. The
date-zero normalizer is fixed across all conditional states/continuations.
Composite midpoint integration respects every deterministic curve jump.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .pricing import DatedM2Pricer, ForwardVarianceCurve
from .states import joint_transition_covariance, symmetric_square_root


@dataclass(frozen=True)
class VIXWindow:
    expiry: float
    delta: float
    times: np.ndarray
    weights: np.ndarray
    normalizer: np.ndarray
    xi: np.ndarray
    owner: object
    configuration: tuple


@dataclass(frozen=True)
class ConditionalVIXResult:
    vix: np.ndarray
    integrated_variance: np.ndarray
    engine: str
    inner_paths: int | None
    integrated_variance_se: np.ndarray | None


class ConditionalVIXPricer:
    def __init__(self, dated: DatedM2Pricer, parameters, curve: ForwardVarianceCurve,
                 initial_state, normalizer_order: int = 32, variance_floor: float = 1e-4):
        self.dated, self.parameters, self.curve = dated, parameters, curve
        self.initial_state = dated._state(initial_state)
        if not np.isfinite(variance_floor) or variance_floor < 0:
            raise ValueError("Variance floor must be finite and nonnegative")
        self.normalizer_order, self.variance_floor = normalizer_order, variance_floor

    def _configuration(self):
        return (self.dated.architecture, self.dated.premium, self.parameters, self.curve,
                self.initial_state.tobytes(), self.normalizer_order, self.variance_floor)

    def window(self, expiry: float, max_step: float = 1/504, delta: float = 30/365) -> VIXWindow:
        if (not np.isfinite([expiry, max_step, delta]).all()
                or expiry < 0 or max_step <= 0 or delta <= 0):
            raise ValueError("Require nonnegative expiry and positive finite window/step")
        end = expiry + delta
        self.curve([expiry, end])
        events = np.unique([expiry, end, *[k for k in self.curve.ends if expiry < k < end]])
        lefts, rights = [], []
        for a, b in zip(events[:-1], events[1:]):
            count = max(1, int(np.ceil((b-a)/max_step)))
            grid = np.linspace(a, b, count+1)
            lefts.extend(grid[:-1]); rights.extend(grid[1:])
        lefts, rights = np.asarray(lefts), np.asarray(rights)
        times, weights = (lefts+rights)/2, rights-lefts
        xi = self.curve(times)
        if np.any(xi <= self.variance_floor):
            raise ValueError("Forward variance must exceed the floor throughout VIX window")
        normalizer = self.dated.normalizer(times, self.initial_state, self.parameters, self.normalizer_order)
        for array in (times, weights, normalizer, xi):
            array.setflags(write=False)
        return VIXWindow(expiry, delta, times, weights, normalizer, xi, self, self._configuration())

    def evaluate(self, states, window: VIXWindow, engine="gaussian_quadrature",
                 order: int = 20, inner_paths: int = 128, seed: int = 0,
                 batch_size: int = 16) -> ConditionalVIXResult:
        if window.owner is not self or window.configuration != self._configuration():
            raise ValueError("Window belongs to a different date-zero model/normalizer")
        states = np.asarray(states, dtype=float)
        if (states.ndim != 2 or states.shape[1] != len(self.initial_state)
                or len(states) < 1 or not np.isfinite(states).all() or batch_size < 1):
            raise ValueError("Require a finite batch of full conditional OU states")
        if engine == "gaussian_quadrature":
            integrated = np.zeros(len(states))
            for time, weight, xi, normalizer in zip(window.times, window.weights, window.xi, window.normalizer):
                mean, cov = self.dated.gaussian_projection_law(float(time-window.expiry), states)
                expected = self.dated.gaussian_expectation(mean, cov, self.parameters, order)
                integrated += weight * (self.variance_floor + (xi-self.variance_floor)*expected/normalizer)
            inner = None
            inner_se = None
        elif engine == "nested_mc":
            if not isinstance(inner_paths, int) or inner_paths < 2:
                raise ValueError("At least two inner continuations required")
            integrated, inner_se = self._nested(states, window, inner_paths, seed, batch_size)
            inner = inner_paths
        else:
            raise ValueError("Unknown VIX engine")
        if not np.isfinite(integrated).all() or np.any(integrated <= 0):
            raise FloatingPointError("Nonpositive or nonfinite conditional variance integral")
        return ConditionalVIXResult(100*np.sqrt(integrated/window.delta), integrated, engine, inner, inner_se)

    def _nested(self, states, window, inner, seed, batch_size):
        """Full-state conditional OU paths; each outer state has independent inner streams.

        Streams keyed by outer index and integration node preserve inner prefixes
        across count refinements, with the same window and outer-state ordering.
        """
        output = np.empty(len(states))
        errors = np.empty(len(states))
        rates = self.dated.layout.rates
        previous = window.expiry
        transitions = []
        for time in window.times:
            dt = float(time-previous)
            root = symmetric_square_root(joint_transition_covariance(self.dated.layout, dt)[1:, 1:])
            decay = np.exp(-rates*dt)
            shift = self.dated.q_drift*(-np.expm1(-rates*dt))/rates
            transitions.append((root, decay, shift))
            previous = time
        for start in range(0, len(states), batch_size):
            stop = min(start+batch_size, len(states))
            current = np.repeat(states[start:stop, None, :], inner, axis=1)
            path_integral = np.zeros((stop-start, inner))
            for node, ((root, decay, shift), weight, xi, normalizer) in enumerate(zip(
                    transitions, window.weights, window.xi, window.normalizer)):
                normals = np.stack([np.random.default_rng(np.random.SeedSequence(
                    [seed, outer, node])).standard_normal((inner, len(rates))) for outer in range(start, stop)])
                current = current*decay+shift+normals @ root.T
                multiplier = self.dated.multiplier(current @ self.dated.projection.T, self.parameters)
                path_integral += weight*(self.variance_floor+(xi-self.variance_floor)*multiplier/normalizer)
            output[start:stop] = path_integral.mean(axis=1)
            errors[start:stop] = path_integral.std(axis=1, ddof=1)/np.sqrt(inner)
        return output, errors

    def moment_target(self, window: VIXWindow) -> float:
        if window.owner is not self or window.configuration != self._configuration():
            raise ValueError("Window belongs to another model")
        return 1e4*self.curve.integral(window.expiry, window.expiry+window.delta)/window.delta
