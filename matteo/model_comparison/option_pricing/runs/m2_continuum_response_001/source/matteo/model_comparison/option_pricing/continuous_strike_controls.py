"""Continuous strike interpolation of FROZEN zero-mean random controls.

This interpolates control variates path by path, never market/model prices.
Deterministic interpolation of known-zero controls remains known-zero.
Continuity does not imply finite-sample positivity, monotonicity or convexity.
"""
import numpy as np


def interpolation_weights(strikes, anchors):
    strikes, anchors = np.asarray(strikes, float), np.asarray(anchors, float)
    if (strikes.ndim != 1 or anchors.ndim != 1 or not len(anchors)
            or not np.isfinite(strikes).all() or not np.isfinite(anchors).all()
            or np.any(strikes <= 0) or np.any(anchors <= 0)
            or np.any(np.diff(anchors) <= 0)):
        raise ValueError('Require positive strikes and strictly increasing anchors')
    out=np.zeros((len(strikes),len(anchors)))
    if len(anchors)==1:out[:,0]=1;return out
    right=np.clip(np.searchsorted(anchors,strikes,side='right'),1,len(anchors)-1)
    left=right-1;fraction=np.clip((strikes-anchors[left])/(anchors[right]-anchors[left]),0,1)
    out[np.arange(len(strikes)),left]=1-fraction
    out[np.arange(len(strikes)),right]=fraction
    return out


def correction_draws(strikes, anchors, frozen_zero_mean_draws):
    draws=np.asarray(frozen_zero_mean_draws,float)
    if draws.ndim!=2 or draws.shape[1]!=len(anchors) or not np.isfinite(draws).all():
        raise ValueError('Require path-by-anchor frozen control draws')
    return draws@interpolation_weights(strikes,anchors).T


def affine_correction_draws(strikes, anchors, frozen_zero_mean_draws):
    """Global affine control: preserves convexity, not positivity/monotonicity.

    Least-squares weights depend only on declared anchor strikes. No price,
    payoff or evaluation-path information selects the weights.
    """
    strikes,anchors=np.asarray(strikes,float),np.asarray(anchors,float)
    draws=np.asarray(frozen_zero_mean_draws,float)
    interpolation_weights(strikes,anchors)  # shared input validation
    if draws.ndim!=2 or draws.shape[1]!=len(anchors) or len(anchors)<2:
        raise ValueError('Require at least two anchor control columns')
    center=float(np.mean(anchors));scale=float(np.ptp(anchors))
    design=np.column_stack([np.ones(len(anchors)),(anchors-center)/scale])
    target=np.column_stack([np.ones(len(strikes)),(strikes-center)/scale])
    return draws@(target@np.linalg.pinv(design)).T


def warm_start(previous, lower, upper, names):
    """Keep the old parameters, changing only beta_J to its positive bound."""
    previous=np.asarray(previous,float);lower=np.asarray(lower,float);upper=np.asarray(upper,float)
    if previous.shape!=lower.shape or previous.shape!=upper.shape or len(names)!=len(previous):
        raise ValueError('Incompatible parameter vectors')
    result=previous.copy();j=names.index('spike_shift')
    if lower[j]<=0:raise ValueError('beta_J lower bound must be strictly positive')
    result[j]=max(result[j],lower[j])
    if not np.isfinite(result).all() or np.any(result<lower) or np.any(result>upper):
        raise ValueError('Previous fit is outside another bound; do not silently change it')
    return result
