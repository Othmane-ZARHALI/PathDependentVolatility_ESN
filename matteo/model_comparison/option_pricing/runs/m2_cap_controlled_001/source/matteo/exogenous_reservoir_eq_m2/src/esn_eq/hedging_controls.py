"""Predictable delta-like martingale controls for dated SPX simulation.

For the exponential-Euler stock, each adapted stock increment has conditional
mean zero. Therefore sum(delta_t * (Y_next-Y_t)) has expectation zero for any
bounded, adapted delta_t. The Black-inspired functions here are controls, not
claims about the true stochastic-volatility option delta.
"""
from dataclasses import dataclass

import numpy as np
from scipy.special import ndtr


@dataclass(frozen=True)
class HedgedSPXPaths:
    terminal_ratios: dict
    hedges: np.ndarray  # paths x contracts x two predictable controls


def spx_hedging_controls(pricer, cache, parameters, curve, maturities, strikes, kinds,
                         order=48, variance_floor=1e-4):
    """Two controls: deterministic-xi Black delta, and a spot-variance proxy.

    The second remaining-variance proxy is W(t,T)*V_t/xi(t), evaluated at the
    start of the step. Both multiply the model stock's actual discrete
    martingale increment. Neither control uses future states or shocks.
    """
    maturities,strikes = np.asarray(maturities,float),np.asarray(strikes,float)
    kinds = np.asarray(kinds)
    if (maturities.ndim != 1 or not len(maturities) or strikes.shape != maturities.shape
            or kinds.shape != maturities.shape or not np.isfinite(maturities).all()
            or not np.isfinite(strikes).all() or np.any(strikes <= 0)
            or not np.isin(kinds,["call","put"]).all()
            or not set(maturities).issubset(set(cache.times[1:]))):
        raise ValueError("Require valid contracts with exact cached maturities")
    if cache.architecture != pricer.architecture or cache.premium != pricer.premium:
        raise ValueError("Incompatible path cache")
    if not np.isfinite(variance_floor) or variance_floor < 0:
        raise ValueError("Invalid variance floor")
    if any(0 < t < cache.times[-1] and t not in cache.times for t in curve.ends):
        raise ValueError("Curve knots must be exact integration events")
    xi = curve(cache.times[:-1])
    if np.any(xi <= variance_floor):
        raise ValueError("Curve levels must exceed variance floor")
    normalizer = pricer.normalizer(cache.times[:-1],cache.initial_state,parameters,order)
    variance = variance_floor+(xi[:,None]-variance_floor)*pricer.multiplier(cache.projections,parameters)/normalizer[:,None]
    if not np.isfinite(variance).all() or np.any(variance <= 0):
        raise FloatingPointError("Invalid variance paths")
    dt = np.diff(cache.times)
    accumulated = np.r_[0.,np.cumsum(xi*dt)]
    end_indices = np.searchsorted(cache.times,maturities)
    total_variance = accumulated[end_indices]
    log_stock = np.zeros(cache.equity_increments.shape[1])
    hedges = np.zeros((len(log_stock),len(strikes),2))
    terminal = {}
    for step,time in enumerate(cache.times[:-1]):
        active = np.flatnonzero(maturities > time)
        current = np.exp(log_stock)
        log_next = log_stock-.5*variance[step]*dt[step]+np.sqrt(variance[step])*cache.equity_increments[step]
        following = np.exp(log_next)
        stock_increment = following-current
        remaining = total_variance[active]-accumulated[step]
        if np.any(remaining <= 0):
            raise FloatingPointError("Nonpositive remaining deterministic variance")
        for control in (0,1):
            w = np.broadcast_to(remaining,(len(current),len(active)))
            if control == 1:
                w = w*(variance[step,:,None]/xi[step])
            d1 = (log_stock[:,None]-np.log(strikes[active])+.5*w)/np.sqrt(w)
            delta = ndtr(d1)-(kinds[active] == "put")
            hedges[:,active,control] += delta*stock_increment[:,None]
        log_stock = log_next
        if cache.times[step+1] in maturities:
            terminal[float(cache.times[step+1])] = following.copy()
    if not np.isfinite(hedges).all() or any(not np.isfinite(y).all() for y in terminal.values()):
        raise FloatingPointError("Nonfinite stock or hedging controls")
    return HedgedSPXPaths(terminal,hedges)
