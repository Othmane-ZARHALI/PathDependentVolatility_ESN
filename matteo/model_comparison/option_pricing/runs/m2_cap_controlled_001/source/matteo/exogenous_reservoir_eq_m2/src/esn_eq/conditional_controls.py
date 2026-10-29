"""Conditional expectations of the existing discrete SPX controls.

All means below condition on retained shocks U. Their target-Q expectations,
not their conditional values, are the known control means. The proposal tilts
only retained standard normals; independent equity normals are integrated out.
"""
from dataclasses import dataclass

import numpy as np
from scipy.special import ndtr

from .conditional_spx import partial_transition, lognormal_payoffs
from .volatility_path_importance import GaussianPathTilt


def stock_delta_mean(mean, variance, strikes, remaining, kinds):
    """E[e^L delta_Black(e^L; K,W)], L~N(mean,variance), W>0.

    Exponential tilting shifts the Gaussian mean by its variance. Integrating
    Phi then adds variance to W. This is not delta evaluated at E[e^L].
    """
    mean, variance, strikes, remaining = np.broadcast_arrays(
        np.asarray(mean, float), np.asarray(variance, float),
        np.asarray(strikes, float), np.asarray(remaining, float))
    kinds = np.broadcast_to(np.asarray(kinds), mean.shape)
    if (any(not np.isfinite(x).all() for x in (mean, variance, strikes, remaining))
            or np.any(variance < 0) or np.any(strikes <= 0) or np.any(remaining <= 0)
            or not np.isin(kinds, ['call', 'put']).all()):
        raise ValueError('Invalid Gaussian delta inputs')
    d = (mean+variance-np.log(strikes)+.5*remaining)/np.sqrt(remaining+variance)
    # Evaluate put directly to avoid cancellation in Phi(d)-1.
    delta = np.where(kinds == 'call', ndtr(d), -ndtr(-d))
    return np.exp(mean+.5*variance)*delta


def retained_path_basis(pricer, times, anchor):
    """Four orthonormal scores in the NEW retained-normal coordinates.

    Retained W_H, then residual projected C_H, Y_anchor and F_anchor. This is a
    new defensive proposal, analogous to vol7; it is not the marginal density
    of the old full-path vol7 proposal and cannot use that proposal's weights.
    """
    times = np.asarray(times, float)
    if (times.ndim != 1 or len(times) < 2 or times[0] != 0
            or not np.isfinite(times).all() or np.any(np.diff(times) <= 0)
            or anchor <= 0 or anchor not in times):
        raise ValueError('Require a finite increasing grid with exact anchor')
    basis = np.zeros((len(times)-1, len(pricer.layout.rates)+1, 4))
    roots = {}
    for step, dt in enumerate(np.diff(times)):
        if dt not in roots:
            roots[dt] = partial_transition(pricer, float(dt))[0]
        root = roots[dt]
        basis[step, :, 0] = root[0]
        for axis, (projection, endpoint) in enumerate(
                ((0, times[-1]), (1, anchor), (2, anchor)), 1):
            if times[step+1] <= endpoint:
                loading = pricer.projection[projection]*np.exp(
                    -pricer.layout.rates*(endpoint-times[step+1]))
                basis[step, :, axis] = loading@root[1:]
    for axis in range(4):
        for previous in range(axis):
            basis[:, :, axis] -= np.sum(basis[:, :, axis]*basis[:, :, previous])*basis[:, :, previous]
        length = np.linalg.norm(basis[:, :, axis])
        if length < 1e-10:
            raise ValueError('Degenerate retained score direction')
        basis[:, :, axis] /= length
    basis.setflags(write=False)
    return basis


@dataclass(frozen=True)
class ConditionalControlPaths:
    payoff: np.ndarray
    controls: np.ndarray  # stock-1, reference payoff-Black, two adapted gains
    stock: np.ndarray
    reference_stock: np.ndarray
    weights: np.ndarray
    scores: np.ndarray
    components: np.ndarray
    transition_audit: tuple
    raw_payoff: np.ndarray | None
    raw_controls: np.ndarray | None


def simulate_conditional_controls(pricer, times, initial, parameters, curve,
        maturities, strikes, kinds, paths, seed, order=48, variance_floor=1e-4,
        *, basis=None, mixture=None, include_raw=False):
    """Stream exact conditional controls under a fixed discrete approximation.

    At step k let d=-V dt/2+sqrt(V)*a, q=nu² V dt, L~N(m,v).
    Independence of the eliminated next equity shock gives
    E[delta(Y_k)(Y_next-Y_k)|U] = expm1(d+q/2)*E[e^L delta(e^L)|U].
    Future U determines V, but never restricts eliminated independent shocks.
    Original zero means hold after averaging U under target Q.
    """
    times = np.asarray(times, float)
    maturities, strikes = np.asarray(maturities, float), np.asarray(strikes, float)
    kinds = np.asarray(kinds)
    initial = pricer._state(initial)
    if (times.ndim != 1 or len(times) < 2 or times[0] != 0
            or not np.isfinite(times).all() or np.any(np.diff(times) <= 0)
            or paths < 2 or maturities.ndim != 1 or not len(maturities)
            or strikes.shape != maturities.shape or kinds.shape != maturities.shape
            or np.any(strikes <= 0) or not np.isfinite(strikes).all()
            or not np.isin(kinds, ['call', 'put']).all()
            or not set(maturities).issubset(set(times[1:]))):
        raise ValueError('Invalid conditional-control inputs')
    if not np.isfinite(variance_floor) or variance_floor < 0:
        raise ValueError('Invalid variance floor')
    if any(0 < t < times[-1] and t not in times for t in curve.ends):
        raise ValueError('Insert exact curve knots')
    if (basis is None) != (mixture is None):
        raise ValueError('Supply both retained proposal basis and mixture')
    if basis is not None and np.asarray(basis).shape[:2] != (len(times)-1, len(initial)+1):
        raise ValueError('Retained basis/grid mismatch')
    component_seed, path_seed, raw_seed = np.random.SeedSequence(seed).spawn(3)
    tilt = GaussianPathTilt(basis, mixture, paths, component_seed) if basis is not None else None
    rng = np.random.default_rng(path_seed)
    raw_rng = np.random.default_rng(raw_seed)
    dt_all = np.diff(times)
    xi = curve(times[:-1])
    if np.any(xi <= variance_floor):
        raise ValueError('Curve below variance floor')
    norm = pricer.normalizer(times[:-1], initial, parameters, order)
    accumulated = np.r_[0., np.cumsum(xi*dt_all)]
    total = accumulated[np.searchsorted(times, maturities)]
    known = lognormal_payoffs(-.5*total, total, strikes, kinds)
    state = np.broadcast_to(initial, (paths, len(initial))).copy()
    m, v, mr, vr = (np.zeros(paths) for _ in range(4))
    raw_log, raw_ref = np.zeros(paths), np.zeros(paths)
    gains = np.zeros((paths, len(strikes), 2))
    raw_gains = np.zeros_like(gains) if include_raw else None
    output = [np.empty((paths, len(strikes))) for _ in range(6)]
    payoff, stock, reference_payoff, reference_stock, raw_payoff, raw_stock = output
    raw_reference_payoff = np.empty_like(payoff) if include_raw else None
    transitions, audits = {}, []
    for step, dt in enumerate(dt_all):
        if dt not in transitions:
            root, independent, audit = partial_transition(pricer, float(dt))
            decay = np.exp(-pricer.layout.rates*dt)
            shift = pricer.q_drift*(-np.expm1(-pricer.layout.rates*dt))/pricer.layout.rates
            transitions[dt] = root, independent, decay, shift
            audits.append(audit)
        root, independent, decay, shift = transitions[dt]
        projection = state@pricer.projection.T
        variance = variance_floor+(xi[step]-variance_floor)*pricer.multiplier(projection, parameters)/norm[step]
        if not np.isfinite(variance).all() or np.any(variance <= 0):
            raise FloatingPointError('Invalid conditional variance')
        z = rng.standard_normal((paths, len(initial)+1))
        if tilt is not None:
            z = tilt.transform(step, z)
        joint = z@root.T
        d = -.5*variance*dt+np.sqrt(variance)*joint[:, 0]
        q = variance*independent
        active = np.flatnonzero(maturities > times[step])
        remaining = total[active]-accumulated[step]
        if include_raw:
            dw = joint[:, 0]+np.sqrt(independent)*raw_rng.standard_normal(paths)
            next_log = raw_log-.5*variance*dt+np.sqrt(variance)*dw
            raw_increment = np.exp(next_log)-np.exp(raw_log)
        for control in (0, 1):
            w = np.broadcast_to(remaining, (paths, len(active)))
            if control == 1:
                w = w*(variance[:, None]/xi[step])
            sd = stock_delta_mean(m[:, None], v[:, None], strikes[active], w, kinds[active])
            gains[:, active, control] += np.expm1(d+.5*q)[:, None]*sd
            if include_raw:
                d1 = (raw_log[:, None]-np.log(strikes[active])+.5*w)/np.sqrt(w)
                delta = np.where(kinds[active] == 'call', ndtr(d1), -ndtr(-d1))
                raw_gains[:, active, control] += delta*raw_increment[:, None]
        m += d
        v += q
        mr += -.5*xi[step]*dt+np.sqrt(xi[step])*joint[:, 0]
        vr += xi[step]*independent
        if include_raw:
            raw_log = next_log
            raw_ref += -.5*xi[step]*dt+np.sqrt(xi[step])*dw
        state = state*decay+shift+joint[:, 1:]
        columns = np.flatnonzero(maturities == times[step+1])
        if len(columns):
            payoff[:, columns] = lognormal_payoffs(m[:, None], v[:, None], strikes[columns], kinds[columns])
            stock[:, columns] = np.exp(m+.5*v)[:, None]
            reference_payoff[:, columns] = lognormal_payoffs(mr[:, None], vr[:, None], strikes[columns], kinds[columns])
            reference_stock[:, columns] = np.exp(mr+.5*vr)[:, None]
            if include_raw:
                raw_stock[:, columns] = np.exp(raw_log)[:, None]
                sign = np.where(kinds[columns] == 'call', 1., -1.)
                raw_payoff[:, columns] = np.maximum(sign*(np.exp(raw_log)[:, None]-strikes[columns]), 0.)
                raw_reference_payoff[:, columns] = np.maximum(sign*(np.exp(raw_ref)[:, None]-strikes[columns]), 0.)
    controls = np.concatenate([(stock-1.)[:, :, None], (reference_payoff-known)[:, :, None], gains], axis=2)
    if include_raw:
        raw_controls = np.concatenate([(raw_stock-1.)[:, :, None], (raw_reference_payoff-known)[:, :, None], raw_gains], axis=2)
    else:
        raw_controls = None
    weights = mixture.weights(tilt.scores) if tilt is not None else np.ones(paths)
    scores = tilt.scores if tilt is not None else np.empty((paths, 0))
    components = tilt.components if tilt is not None else np.zeros(paths, int)
    if not all(np.isfinite(a).all() for a in (payoff, controls, stock, reference_stock, weights)):
        raise FloatingPointError('Nonfinite conditional controls')
    return ConditionalControlPaths(payoff, controls, stock, reference_stock,
        weights, scores, components, tuple(audits), raw_payoff if include_raw else None, raw_controls)
