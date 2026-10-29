"""Current calibration policy; historical run bounds are preserved unchanged.

User decision, 24 October 2026: minimum annualized volatility about 5%,
with a calibratable floor in the 5--6% volatility range.
This is a calibration constraint, not a ban on repricing historical models.
"""
import numpy as np
from .controlled_readout_fast_fit import LOWER as HISTORICAL_LOWER, UPPER as HISTORICAL_UPPER

MINIMUM_VOLATILITY = .05
MAXIMUM_FLOOR_VOLATILITY = .06
MINIMUM_VARIANCE_FLOOR = MINIMUM_VOLATILITY ** 2
MAXIMUM_VARIANCE_FLOOR = MAXIMUM_FLOOR_VOLATILITY ** 2
LOWER = HISTORICAL_LOWER.copy()
UPPER = HISTORICAL_UPPER.copy()
LOWER[10] = MINIMUM_VARIANCE_FLOOR
UPPER[10] = MAXIMUM_VARIANCE_FLOOR


def validate(parameters):
    p = np.asarray(parameters, dtype=float)
    if p.shape != (11,) or not np.isfinite(p).all():
        raise ValueError('Require eleven finite calibration parameters')
    if np.any(p < LOWER) or np.any(p > UPPER):
        raise ValueError('Parameters outside current calibration bounds (volatility floor 5--6%)')
    if np.any(p[7:10] < p[10]) or p[4] <= 0:
        raise ValueError('Require xi >= variance floor and beta_J > 0')
    return p
