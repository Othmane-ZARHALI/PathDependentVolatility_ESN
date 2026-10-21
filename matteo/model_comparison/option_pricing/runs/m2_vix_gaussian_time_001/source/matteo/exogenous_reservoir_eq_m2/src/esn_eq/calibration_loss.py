"""Frozen common price loss; independent of either model's pricing engine.

The loss is a descriptive pilot objective, not a likelihood or an accuracy
certificate. Prices are undiscounted values at actual quoted strikes. Market
scales remain fixed when model parameters or numerical budgets change.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .market_targets import interval_residuals


@dataclass(frozen=True)
class LossConfig:
    spx_weight: float = 1.
    vix_weight: float = 1.
    forward_weight: float = 1.
    forward_scale: float = .5
    midpoint_weight: float = 0.

    def __post_init__(self):
        weights = [self.spx_weight, self.vix_weight, self.forward_weight, self.midpoint_weight]
        if (not np.isfinite(weights).all() or min(weights) < 0
                or not np.isfinite(self.forward_scale) or self.forward_scale <= 0):
            raise ValueError("Require finite nonnegative weights and a positive forward scale")


class JointPriceObjective:
    """Freeze selected rows and return a weighted least-squares residual vector.

    Prediction arrays must follow the supplied row order. Default construction
    rejects held-out rows; evaluation on a different split requires an explicit
    new objective with that split. Forwards are supplied once per expiry.
    """
    def __init__(self, quotes, forwards, config=LossConfig(), split="train"):
        self._quotes = quotes.copy(deep=True).reset_index(drop=True)
        self._forwards = forwards.copy(deep=True).reset_index(drop=True)
        self.config = config
        q, f = self._quotes, self._forwards
        if (q.empty or not q.split.eq(split).all() or not q.instrument.isin(["SPX", "VIX"]).all()
                or not q.option_type.isin(["put", "call"]).all()
                or q.duplicated(["instrument", "expiry", "strike", "option_type"]).any()):
            raise ValueError("Require unique contracts from one explicitly selected split")
        interval_residuals(q.price_bid_U, q.price_bid_U, q.price_ask_U, q.price_scale_U)
        if (not np.isfinite(q[["strike", "forward", "T_years"]]).all().all()
                or (q[["strike", "forward", "T_years"]] <= 0).any().any()):
            raise ValueError("Require positive finite strike, market forward and maturity")
        for _, group in q.groupby(["instrument", "expiry"]):
            if group.forward.nunique() != 1 or group.T_years.nunique() != 1:
                raise ValueError("Inconsistent market forward/maturity within expiry")
        if (f.expiry.duplicated().any()
                or not np.isfinite(f.forward_vix_points).all() or (f.forward_vix_points <= 0).any()):
            raise ValueError("Require one positive finite VIX forward per expiry")
        self._weights = np.zeros(len(q))
        for instrument, weight in (("SPX", config.spx_weight), ("VIX", config.vix_weight)):
            rows = q[q.instrument == instrument]
            if rows.empty and weight > 0:
                raise ValueError("Positive family weight requires corresponding quotes")
            for _, group in rows.groupby("expiry"):
                self._weights[group.index] = weight/(rows.expiry.nunique()*len(group))
        if f.empty and config.forward_weight > 0:
            raise ValueError("Positive forward weight requires forward targets")

    def evaluate(self, prices, model_forwards, *, state_penalty=0., curve_penalty=0.):
        q, f, cfg = self._quotes.copy(), self._forwards.copy(), self.config
        prices, model_forwards = np.asarray(prices, float), np.asarray(model_forwards, float)
        if prices.shape != (len(q),) or model_forwards.shape != (len(f),):
            raise ValueError("Prediction shape must match frozen quote/forward order")
        if (not np.isfinite(model_forwards).all() or np.any(model_forwards <= 0)
                or not np.isfinite([state_penalty, curve_penalty]).all()
                or min(state_penalty, curve_penalty) < 0):
            raise ValueError("Invalid forward prediction or regularisation penalty")
        q["model_price_U"] = prices
        q["interval_residual"] = interval_residuals(prices, q.price_bid_U, q.price_ask_U, q.price_scale_U)
        q["midpoint_residual"] = (prices-(q.price_bid_U+q.price_ask_U)/2)/q.price_scale_U
        q["loss_weight"] = self._weights
        f["model_forward"] = model_forwards
        f["error_points"] = model_forwards-f.forward_vix_points
        f["residual"] = f.error_points/cfg.forward_scale
        components = {}
        for instrument in ("SPX", "VIX"):
            rows = q[q.instrument == instrument]
            components[instrument] = float(np.sum(rows.loss_weight*rows.interval_residual**2))
        components["midpoint"] = float(cfg.midpoint_weight*np.sum(self._weights*q.midpoint_residual**2))
        components["forwards"] = float(cfg.forward_weight*np.mean(f.residual**2)) if len(f) else 0.
        components["state_penalty"], components["curve_penalty"] = float(state_penalty), float(curve_penalty)
        residuals = np.r_[np.sqrt(self._weights)*q.interval_residual,
                          np.sqrt(cfg.midpoint_weight*self._weights)*q.midpoint_residual,
                          np.sqrt(cfg.forward_weight/max(len(f), 1))*f.residual,
                          np.sqrt([state_penalty, curve_penalty])]
        blocks = []
        for keys, rows in q.groupby(["instrument", "expiry", "option_type"]):
            blocks.append(dict(instrument=keys[0], expiry=keys[1], option_type=keys[2],
                quotes=len(rows), interval_mse=float(np.mean(rows.interval_residual**2)),
                inside_fraction=float(np.mean(rows.interval_residual == 0)),
                max_abs_scaled_distance=float(np.max(np.abs(rows.interval_residual)))))
        return dict(total=sum(components.values()), components=components, residuals=residuals,
                    quotes=q, forwards=f, blocks=pd.DataFrame(blocks),
                    numerical_status="not assessed by loss; supplied by pricing diagnostics")
