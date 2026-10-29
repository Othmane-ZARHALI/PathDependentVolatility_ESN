"""Hard admissibility condition for any number of forward-variance pieces.

No monotonicity condition across maturities is imposed. The calibration policy
sets the separate 5% minimum-volatility bound on the variance floor.
"""
import numpy as np


def validate_forward_variance_levels(levels, variance_floor):
    values = np.asarray(levels, dtype=float)
    floor = np.asarray(variance_floor, dtype=float)
    if floor.ndim != 0 or not np.isfinite(floor) or floor <= 0:
        raise ValueError('Require a finite positive variance floor')
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError('Require a nonempty finite vector of forward-variance levels')
    if np.any(values < floor):
        raise ValueError('Every forward-variance level xi_i must be >= variance floor')
    return values
