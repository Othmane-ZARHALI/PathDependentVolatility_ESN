"""Defensive importance sampling of the complete correlated SPX/OU path."""
from dataclasses import dataclass

import numpy as np

from .importance_sampling import DefensiveGaussianMixture
from .pricing import PreparedQPaths


@dataclass(frozen=True)
class ImportancePaths:
    cache: PreparedQPaths
    weights: np.ndarray
    score: np.ndarray
    components: np.ndarray


def prepare_importance_paths(pricer, times, initial, paths, seed, magnitude=2., checkpoints=()):
    """Half target, quarter each Brownian shift +/- magnitude/sqrt(horizon).

    For Z=W_H/sqrt(H), q/f=sum pi exp(m Z-m²/2). This is the
    FULL discrete joint path likelihood: conditional OU residuals are unchanged
    and the stock/OU cross innovations shift together. Weights are f/q, never
    self-normalized, and bounded by two. The same full-horizon weight also
    applies to earlier expiries and any coupled coarsening of this cache.
    Paths are iid; no Sobol or antithetic uncertainty claim is made here.
    """
    times = np.asarray(times, float)
    if (times.ndim != 1 or len(times) < 2 or times[0] != 0 or
            not np.isfinite(times).all() or np.any(np.diff(times) <= 0) or
            not np.isfinite(magnitude) or magnitude <= 0):
        raise ValueError("Require increasing finite grid and positive shift magnitude")
    mixture = DefensiveGaussianMixture([[0.], [magnitude], [-magnitude]], [.5, .25, .25])
    component_seed, path_seed = np.random.SeedSequence(seed).spawn(2)
    component = np.random.default_rng(component_seed).choice(3, size=paths, p=mixture.probabilities)
    drift = mixture.shifts[component, 0]/np.sqrt(times[-1])
    cache = pricer.prepare(times, initial, paths, int(path_seed.generate_state(1)[0]),
                           checkpoints, equity_drift=drift)
    score = cache.equity_increments.sum(axis=0)/np.sqrt(times[-1])
    return ImportancePaths(cache, mixture.weights(score[:, None]), score, component)
