"""Dated baseline-M2 pricing in forward units with deterministic Q normalization.

OU transitions are exact; stock integration uses adapted exponential Euler.
Gaussian quadrature is a numerical approximation, not an exact normalizer.
Only the Gaussian M2 cell is supported (zero echo/orthogonal loadings).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product

import numpy as np
from numpy.polynomial.hermite import hermgauss

from .config import ArchitectureConfig, FloatArray, ReadoutParameters, RiskPremium
from .features import FeatureEngine, FeaturePaths
from .readout import TemperedExponentialQuadratic
from .states import (build_mode_layout, joint_transition_covariance,
                     primitive_risk_prices, stationary_covariance, symmetric_square_root)


@lru_cache(maxsize=8)
def _normal_rule(order: int):
    nodes, weights = hermgauss(order)
    index = np.array(list(product(range(order), repeat=3)))
    return np.sqrt(2.) * nodes[index], np.prod(weights[index] / np.sqrt(np.pi), axis=1)


@dataclass(frozen=True)
class ForwardVarianceCurve:
    """Piecewise-constant xi with explicit finite support; no silent extrapolation."""

    ends: tuple[float, ...]
    levels: tuple[float, ...]

    def __post_init__(self):
        ends, levels = np.asarray(self.ends), np.asarray(self.levels)
        if (ends.ndim != 1 or not len(ends) or len(ends) != len(levels)
                or not np.isfinite(ends).all() or not np.isfinite(levels).all()
                or np.any(np.diff(np.r_[0., ends]) <= 0) or np.any(levels <= 0)):
            raise ValueError("Require positive increasing endpoints and positive finite levels")

    def __call__(self, times) -> FloatArray:
        times = np.asarray(times, dtype=float)
        if not np.isfinite(times).all() or np.any(times < 0) or np.any(times > self.ends[-1]):
            raise ValueError("Curve does not cover requested times")
        indices = np.minimum(np.searchsorted(self.ends, times, side="right"), len(self.ends) - 1)
        return np.asarray(self.levels)[indices]

    def integral(self, start: float, end: float) -> float:
        self(np.asarray([start, end]))
        if end < start:
            raise ValueError("Integration endpoints are reversed")
        left = np.r_[0., self.ends[:-1]]
        lengths = np.maximum(0., np.minimum(end, self.ends) - np.maximum(start, left))
        return float(lengths @ np.asarray(self.levels))


def dated_grid(events, max_step: float) -> FloatArray:
    """Include every event exactly, subdividing intervals without rounded expiries."""
    events = np.unique(np.r_[0., np.asarray(events, dtype=float)])
    if not np.isfinite(events).all() or events[0] < 0 or events[-1] <= 0:
        raise ValueError("Require finite nonnegative events and a positive horizon")
    if not np.isfinite(max_step) or max_step <= 0:
        raise ValueError("max_step must be finite and positive")
    pieces = [np.array([0.])]
    for left, right in zip(events[:-1], events[1:]):
        count = max(1, int(np.ceil((right - left) / max_step)))
        pieces.append(np.linspace(left, right, count + 1)[1:])
    return np.concatenate(pieces)


@dataclass(frozen=True)
class PreparedQPaths:
    """Readout-independent projections and equity increments; dated initial state."""

    times: FloatArray
    projections: FloatArray  # start-of-step, (steps, paths, C/Y/F)
    equity_increments: FloatArray
    initial_state: FloatArray
    checkpoint_states: dict[float, FloatArray]
    seed: int
    architecture: ArchitectureConfig
    premium: RiskPremium


@dataclass(frozen=True)
class SPXPathResult:
    """All maturity ratios use the same stock paths and deterministic carry removal."""

    terminal_ratios: dict[float, FloatArray]
    normalizer: FloatArray
    variance_mean: FloatArray
    variance_mean_se: FloatArray
    cap_fraction: float
    quadrature_order: int


def coarsen_paths(cache: PreparedQPaths, times) -> PreparedQPaths:
    """Couple time-grid diagnostics by summing the same Brownian increments.

    Exact OU semigroup transitions make stored fine-grid states valid on any
    subset grid. Only the adapted stock-integration rule changes.
    """
    times = np.asarray(times, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or times[0] != cache.times[0]
            or times[-1] != cache.times[-1] or np.any(np.diff(times) <= 0)):
        raise ValueError("Coarse grid must span the same interval with increasing events")
    indices = np.searchsorted(cache.times, times)
    if np.any(indices >= len(cache.times)) or not np.array_equal(cache.times[indices], times):
        raise ValueError("Coarse events must be exact members of the cached grid")
    equity = np.add.reduceat(cache.equity_increments, indices[:-1], axis=0)
    return PreparedQPaths(times.copy(), cache.projections[indices[:-1]].copy(), equity,
        cache.initial_state.copy(), {t: x.copy() for t, x in cache.checkpoint_states.items() if t in times},
        cache.seed, cache.architecture, cache.premium)


def payoff_statistics(draws, strikes, option_types, forward_scale: float = 1.) -> tuple[FloatArray, FloatArray]:
    """Undiscounted prices and pathwise MC SE; no IV inversion or failed-point mask."""
    draws, strikes = np.asarray(draws, float), np.asarray(strikes, float)
    kinds = np.asarray(option_types)
    if (draws.ndim != 1 or len(draws) < 2 or strikes.ndim != 1
            or kinds.shape != strikes.shape or not np.isfinite(draws).all()
            or not np.isfinite(strikes).all() or np.any(draws < 0) or np.any(strikes <= 0)
            or not np.isin(kinds, ["call", "put"]).all()
            or not np.isfinite(forward_scale) or forward_scale <= 0):
        raise ValueError("Invalid payoff inputs")
    means, errors = [], []
    for strike, kind in zip(strikes, kinds):
        payoff = np.maximum(draws - strike, 0.) if kind == "call" else np.maximum(strike - draws, 0.)
        means.append(forward_scale * payoff.mean())
        errors.append(forward_scale * payoff.std(ddof=1) / np.sqrt(len(draws)))
    return np.asarray(means), np.asarray(errors)


class DatedM2Pricer:
    """Exact Gaussian state law from fixed full X0 with existing P-to-Q drift convention."""

    def __init__(self, architecture: ArchitectureConfig, premium: RiskPremium):
        if architecture.echo.enabled:
            raise ValueError("Dated pricing currently supports Gaussian M2 only")
        self.architecture, self.premium = architecture, premium
        self.layout = build_mode_layout(architecture)
        engine = FeatureEngine(architecture, self.layout)
        self.projection = np.zeros((3, len(self.layout.rates)))
        for i, name in enumerate(("cluster", "feedback", "spike")):
            self.projection[i, self.layout.slices[name]] = engine.weights[name]
        self.stationary = stationary_covariance(self.layout)
        effective = self.layout.driver_loadings @ primitive_risk_prices(self.layout, premium)
        self.q_drift = -np.sqrt(2 * self.layout.rates) * effective
        self.readout = TemperedExponentialQuadratic()

    def initial_state(self, projections) -> tuple[FloatArray, float]:
        """Minimum stationary-Mahalanobis-norm state for declared (C0,Y0,F0)."""
        u = np.asarray(projections, dtype=float)
        covariance = self.projection @ self.stationary @ self.projection.T
        if u.shape != (3,) or not np.isfinite(u).all():
            raise ValueError("Require three finite projected state coordinates (C,Y,F)")
        if np.linalg.cond(covariance) > 1e10:
            raise ValueError("Projected state prior is singular or ill-conditioned")
        solved = np.linalg.solve(covariance, u)
        return self.stationary @ self.projection.T @ solved, float(u @ solved)

    def _state(self, initial) -> FloatArray:
        initial = np.array(initial, dtype=float, copy=True)
        if initial.shape != self.layout.rates.shape or not np.isfinite(initial).all():
            raise ValueError("Initial full state has wrong dimension or nonfinite values")
        return initial

    def gaussian_projection_law(self, elapsed: float, initial) -> tuple[FloatArray, FloatArray]:
        """Conditional C/Y/F mean and covariance; full state, not projections, is Markov."""
        initial = np.asarray(initial, float)
        if (initial.shape[-1:] != self.layout.rates.shape or not np.isfinite(initial).all()
                or not np.isfinite(elapsed) or elapsed < 0):
            raise ValueError("Invalid conditional state or elapsed time")
        rates = self.layout.rates
        decay = np.exp(-rates * elapsed)
        shift = self.q_drift * (-np.expm1(-rates * elapsed)) / rates
        mean = (decay * initial + shift) @ self.projection.T
        covariance = self.stationary * (-np.expm1(-(rates[:, None] + rates[None, :]) * elapsed))
        return mean, self.projection @ covariance @ self.projection.T

    def multiplier(self, projections, parameters: ReadoutParameters) -> FloatArray:
        if parameters.echo_loading != 0 or parameters.orthogonal_curvature != 0:
            raise ValueError("Three-projection quadrature requires zero echo/orthogonal loadings")
        projections = np.asarray(projections)
        zeros = np.zeros(projections.shape[:-1])
        features = FeaturePaths(projections[..., 1], projections[..., 2], projections[..., 0], zeros, zeros)
        return self.readout.multiplier(features, parameters)

    def gaussian_expectation(self, mean, covariance, parameters, order: int) -> FloatArray:
        """Tensor Gauss-Hermite expectation of G; handle degenerate endpoint directly."""
        if not isinstance(order, int) or order < 2:
            raise ValueError("Quadrature order must be an integer >=2")
        mean = np.asarray(mean, dtype=float)
        if np.max(np.abs(covariance)) == 0:
            return self.multiplier(mean, parameters)
        points, weight = _normal_rule(order)
        offset = points @ symmetric_square_root(covariance).T
        # Chunking avoids allocating all outer-states x quadrature nodes at once.
        flat = mean.reshape(-1, 3)
        result = np.empty(len(flat))
        for start in range(0, len(flat), 32):
            values = self.multiplier(flat[start:start + 32, None, :] + offset[None, :, :], parameters)
            result[start:start + 32] = values @ weight
        return result.reshape(mean.shape[:-1])

    def normalizer(self, times, initial, parameters, order: int) -> FloatArray:
        """Date-zero conditional normalizer; independent of stock simulation draws."""
        initial = self._state(initial)
        values = []
        for time in times:
            mean, covariance = self.gaussian_projection_law(float(time), initial)
            values.append(float(self.gaussian_expectation(mean, covariance, parameters, order)))
        values = np.asarray(values)
        if np.any(~np.isfinite(values)) or np.any(values <= 0):
            raise FloatingPointError("Invalid Gaussian normalizer")
        return values

    def prepare(self, times, initial, paths: int, seed: int, checkpoints=(), *, equity_drift=None, gaussian_tilt=None) -> PreparedQPaths:
        """Cache exact OU/stock shocks, optionally under a constant pathwise shift.

        ``equity_drift`` is one deterministic Brownian drift per path. It changes
        the sampling law, including correlated OU innovations; callers MUST
        supply likelihood weights for target-Q expectations. The Q drift itself
        is unchanged. None preserves the original target sampler.
        Alternatively ``gaussian_tilt.transform(step, normals)`` may shift the
        independent joint-transition normals and accumulate path scores. It
        must return a finite array of the same shape. Its caller owns the
        corresponding full-path likelihood. Do not combine the two interfaces.
        """
        times = np.array(times, dtype=float, copy=True)
        if (times.ndim != 1 or len(times) < 2 or times[0] != 0
                or not np.isfinite(times).all() or np.any(np.diff(times) <= 0) or paths < 2):
            raise ValueError("Require increasing grid starting at zero and at least two paths")
        checkpoints = set(map(float, checkpoints))
        if not checkpoints.issubset(set(times)):
            raise ValueError("Every checkpoint must be an exact grid event")
        initial = self._state(initial)
        if gaussian_tilt is not None and equity_drift is not None:
            raise ValueError("Use only one Gaussian shift interface")
        if equity_drift is not None:
            equity_drift = np.asarray(equity_drift, float)
            if equity_drift.shape != (paths,) or not np.isfinite(equity_drift).all():
                raise ValueError("Require one finite equity drift per path")
        state = np.broadcast_to(initial, (paths, len(initial))).copy()
        features = np.empty((len(times) - 1, paths, 3))
        equity = np.empty((len(times) - 1, paths))
        saved = {0.: state.copy()} if 0. in checkpoints else {}
        rng = np.random.default_rng(seed)
        transitions = {}
        for step, dt in enumerate(np.diff(times)):
            features[step] = state @ self.projection.T
            if dt not in transitions:
                covariance = joint_transition_covariance(self.layout, float(dt))
                cross = covariance[1:, 0]
                conditional = covariance[1:, 1:] - np.outer(cross, cross) / dt
                transitions[dt] = symmetric_square_root(conditional), cross
            root, cross = transitions[dt]
            normals = rng.standard_normal((paths, len(initial) + 1))
            if gaussian_tilt is not None:
                transformed = np.asarray(gaussian_tilt.transform(step, normals), float)
                if transformed.shape != normals.shape or not np.isfinite(transformed).all():
                    raise ValueError("Invalid transformed transition normals")
                normals = transformed
            equity[step] = np.sqrt(dt) * normals[:, 0]
            if equity_drift is not None:
                equity[step] += equity_drift * dt
            innovation = normals[:, 1:] @ root.T + equity[step, :, None] * cross / dt
            decay = np.exp(-self.layout.rates * dt)
            shift = self.q_drift * (-np.expm1(-self.layout.rates * dt)) / self.layout.rates
            state = state * decay + shift + innovation
            if times[step + 1] in checkpoints:
                saved[float(times[step + 1])] = state.copy()
        return PreparedQPaths(times, features, equity, initial, saved, seed,
                              self.architecture, self.premium)

    def spx(self, cache: PreparedQPaths, parameters: ReadoutParameters,
            curve: ForwardVarianceCurve, maturities, order: int = 20,
            variance_floor: float = 1e-4) -> SPXPathResult:
        """Price normalized stock paths, applying no sample-dependent normalizations."""
        times = cache.times
        if cache.architecture != self.architecture or cache.premium != self.premium:
            raise ValueError("Path cache uses a different architecture or Q risk premium")
        maturities = set(map(float, maturities))
        if not maturities or not maturities.issubset(set(times[1:])):
            raise ValueError("Maturities must be positive exact grid events")
        if not np.isfinite(variance_floor) or variance_floor < 0:
            raise ValueError("Variance floor must be finite and nonnegative")
        # No interval may straddle a deterministic curve jump.
        if any(0 < t < times[-1] and t not in times for t in curve.ends):
            raise ValueError("Insert every relevant curve knot into the time grid")
        curve(np.asarray([0., times[-1]]))
        xi = curve(times[:-1])
        if np.any(xi <= variance_floor):
            raise ValueError("Forward variance must exceed the variance floor")
        normalizer = self.normalizer(times[:-1], cache.initial_state, parameters, order)
        multiplier = self.multiplier(cache.projections, parameters)
        variance = variance_floor + (xi[:, None] - variance_floor) * multiplier / normalizer[:, None]
        if not np.isfinite(variance).all():
            raise FloatingPointError("Nonfinite variance")
        log_ratio = np.zeros(cache.equity_increments.shape[1])
        ratios = {}
        for step, dt in enumerate(np.diff(times)):
            log_ratio += -0.5 * variance[step] * dt + np.sqrt(variance[step]) * cache.equity_increments[step]
            if times[step + 1] in maturities:
                ratio = np.exp(log_ratio)
                if not np.isfinite(ratio).all():
                    raise FloatingPointError("Nonfinite terminal stock ratio")
                ratios[float(times[step + 1])] = ratio.copy()
        raw_score = (parameters.cluster_loading * cache.projections[..., 0]
                     + parameters.feedback_curvature * (cache.projections[..., 1] - parameters.feedback_shift)**2
                     + parameters.spike_curvature * (cache.projections[..., 2] - parameters.spike_shift)**2)
        return SPXPathResult(ratios, normalizer, variance.mean(axis=1),
                             variance.std(axis=1, ddof=1) / np.sqrt(variance.shape[1]),
                             float(np.mean(raw_score > parameters.cap_level)), order)
