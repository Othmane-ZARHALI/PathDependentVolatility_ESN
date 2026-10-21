"""SPX controls with known expectations and independently trained coefficients."""
from __future__ import annotations

import numpy as np
from scipy.special import ndtr


def deterministic_stock_controls(cache, curve, maturities):
    """Constant deterministic-volatility stocks on the SAME equity Brownian path."""
    curve([0., cache.times[-1]])
    xi = curve(cache.times[:-1])
    if any(0 < t < cache.times[-1] and t not in cache.times for t in curve.ends):
        raise ValueError("Curve jumps must be exact cache events")
    requested = set(map(float, maturities))
    if not requested.issubset(set(cache.times[1:])):
        raise ValueError("Control maturities must be exact cache events")
    log_stock = np.zeros(cache.equity_increments.shape[1])
    result = {}
    for i, dt in enumerate(np.diff(cache.times)):
        log_stock += -.5*xi[i]*dt+np.sqrt(xi[i])*cache.equity_increments[i]
        if cache.times[i+1] in requested:
            result[float(cache.times[i+1])] = np.exp(log_stock)
    return result


def _control_data(model, reference, strikes, kinds, total_variance):
    model, reference, strikes = (np.asarray(x, dtype=float) for x in (model, reference, strikes))
    kinds = np.asarray(kinds)
    if (model.ndim != 1 or reference.shape != model.shape or len(model) < 2
            or strikes.ndim != 1 or kinds.shape != strikes.shape
            or any(not np.isfinite(x).all() for x in (model, reference, strikes))
            or np.any(model < 0) or np.any(reference <= 0) or np.any(strikes <= 0)
            or not np.isin(kinds, ["call", "put"]).all()
            or not np.isfinite(total_variance) or total_variance <= 0):
        raise ValueError("Invalid stock/control/payoff inputs")
    std = np.sqrt(total_variance)
    for strike, kind in zip(strikes, kinds):
        sign = 1 if kind == "call" else -1
        target = np.maximum(sign*(model-strike), 0)
        reference_payoff = np.maximum(sign*(reference-strike), 0)
        d1 = (-np.log(strike)+.5*total_variance)/std
        d2 = d1-std
        known = (ndtr(d1)-strike*ndtr(d2) if kind == "call"
                 else strike*ndtr(-d2)-ndtr(-d1))
        yield target, np.stack((model-1., reference_payoff-known), axis=1)


def fit_payoff_controls(pilot_model, pilot_reference, strikes, kinds, total_variance):
    """Fit on independent pilot paths; caller must not reuse evaluation paths."""
    coefficients = []
    for target, controls in _control_data(pilot_model, pilot_reference, strikes, kinds, total_variance):
        centered = controls-controls.mean(axis=0)
        beta = np.linalg.lstsq(centered, target-target.mean(), rcond=1e-10)[0]
        coefficients.append(beta)
    return np.asarray(coefficients)


def controlled_payoff_statistics(model, reference, strikes, kinds, total_variance,
                                 coefficients, forward_scale=1.):
    """Unclipped unbiased estimates conditional on an independent frozen pilot.

    SE is conditional on the pilot, normalizer and time grid. Negative estimates
    remain visible and require numerical resolution before interval fitting.
    """
    coefficients = np.asarray(coefficients, dtype=float)
    if (coefficients.shape != (len(strikes), 2) or not np.isfinite(coefficients).all()
            or not np.isfinite(forward_scale) or forward_scale <= 0):
        raise ValueError("Invalid frozen coefficients or forward scale")
    means, errors = [], []
    for beta, (target, controls) in zip(coefficients, _control_data(model, reference, strikes, kinds, total_variance)):
        adjusted = target-controls @ beta
        means.append(forward_scale*adjusted.mean())
        errors.append(forward_scale*adjusted.std(ddof=1)/np.sqrt(len(adjusted)))
    return np.asarray(means), np.asarray(errors)
