"""Declared nuisance inputs for the one-date M2 calibration pilot.

These construct candidate inputs and separate penalties; they do not fit
parameters, clear input compatibility flags or establish physical feasibility.
"""
from dataclasses import dataclass

import numpy as np

from .config import RiskPremium
from .pricing import ForwardVarianceCurve


@dataclass(frozen=True)
class NuisanceCurveFamily:
    ends: tuple[float, ...]
    lower: tuple[float, ...]
    upper: tuple[float, ...]
    smoothness_weight: float = .1
    curve_id: str = "nuisance_piecewise_constant_v1"

    def __post_init__(self):
        lower, upper = np.asarray(self.lower), np.asarray(self.upper)
        ForwardVarianceCurve(self.ends, self.lower)
        if (lower.shape != upper.shape or not np.isfinite(upper).all()
                or np.any(upper <= lower) or not np.isfinite(self.smoothness_weight)
                or self.smoothness_weight < 0):
            raise ValueError("Invalid nuisance-curve bounds or smoothness weight")

    @classmethod
    def from_targets(cls, targets):
        rows = targets.curves[targets.curves.curve_id == "nuisance_piecewise_constant_v1"].sort_values("start_years")
        expected_penalty = "0.1*sum((xi_next-xi_current)^2)"
        if (rows.empty or rows.start_years.iloc[0] != 0
                or not np.array_equal(rows.start_years.to_numpy()[1:], rows.end_years.to_numpy()[:-1])
                or not rows.basis.eq("constant_xi_on_segment").all()
                or not rows.source_split.eq("train_only").all()
                or not rows.penalty.eq(expected_penalty).all()):
            raise ValueError("Unexpected nuisance-curve contract; do not evaluate metadata expressions")
        return cls(tuple(rows.end_years),tuple(rows.coefficient_lower),tuple(rows.coefficient_upper))

    def build(self, levels):
        levels = np.asarray(levels, float)
        if (levels.shape != (len(self.ends),) or not np.isfinite(levels).all()
                or np.any(levels < self.lower) or np.any(levels > self.upper)):
            raise ValueError("Candidate variance levels outside declared bounds")
        return ForwardVarianceCurve(self.ends,tuple(levels))

    def penalty(self, levels):
        self.build(levels)
        return float(self.smoothness_weight*np.sum(np.diff(levels)**2))


@dataclass(frozen=True)
class RestrictedInitialState:
    projections_C_Y_F: tuple[float, float, float]
    full_state: np.ndarray
    prior_squared_norm: float
    penalty: float
    convention: str = "minimum stationary-P prior norm; latent, not historically filtered"


def restricted_initial_state(dated, projections_C_Y_F, weight=.1):
    if not np.isfinite(weight) or weight < 0:
        raise ValueError("State penalty weight must be finite and nonnegative")
    state, norm = dated.initial_state(projections_C_Y_F)
    state.setflags(write=False)
    return RestrictedInitialState(tuple(map(float,projections_C_Y_F)),state,norm,weight*norm/3)


def risk_price_scenarios(equity=.35):
    """The five predeclared primitive unspanned scenarios, retaining equity risk.

    State drift still includes rho*equity in scenario zero; zero unspanned
    prices do not mean a zero total OU drift shift.
    """
    if not np.isfinite(equity):
        raise ValueError("Equity risk price must be finite")
    return tuple(RiskPremium(equity=equity,feedback_idiosyncratic=f,spike_idiosyncratic=j)
                 for f,j in ((0.,0.),(.25,0.),(-.25,0.),(0.,.25),(0.,-.25)))
