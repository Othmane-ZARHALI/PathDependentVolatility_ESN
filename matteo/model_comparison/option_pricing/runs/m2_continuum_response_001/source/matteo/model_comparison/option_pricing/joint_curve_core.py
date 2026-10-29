"""Fast fixed-shock curve repricing and explicitly signed exploratory loss.

The OU path law and G/M0 do not depend on deterministic xi at fixed state,
readout and Q dynamics. Recompute the stock, never reuse its terminal value
after xi changes. Anchor controls are fixed functions of these same primitive
paths with known target means, hence remain valid for every candidate xi.
"""
from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class CurvePathBatch:
    multiplier_ratio: np.ndarray
    equity: np.ndarray
    dt: np.ndarray
    segment: np.ndarray
    floor: float = 1e-4

    def terminal(self, levels, derivative=True, chunk=2048):
        levels = np.asarray(levels, float)
        if (levels.shape != (3,) or not np.isfinite(levels).all() or np.any(levels <= self.floor)
                or self.multiplier_ratio.shape != self.equity.shape
                or self.equity.shape[0] != len(self.dt)
                or self.segment.shape != self.dt.shape):
            raise ValueError('Invalid frozen path/curve inputs')
        logstock = np.zeros(self.equity.shape[1])
        gradient = np.zeros((len(logstock), 3))
        for start in range(0, len(self.dt), chunk):
            end = min(start+chunk, len(self.dt))
            a = self.multiplier_ratio[start:end]
            dw = self.equity[start:end]
            dt = self.dt[start:end, None]
            segment = self.segment[start:end]
            v = self.floor+(levels[segment, None]-self.floor)*a
            if not np.isfinite(v).all() or np.any(v <= 0):
                raise FloatingPointError('Invalid repriced variance')
            sv = np.sqrt(v)
            logstock += (-.5*v*dt+sv*dw).sum(axis=0)
            if derivative:
                terms = a*(-.5*dt+.5*dw/sv)
                for j in np.unique(segment):
                    gradient[:, j] += terms[segment == j].sum(axis=0)
        stock = np.exp(logstock)
        if not np.isfinite(stock).all():
            raise FloatingPointError('Nonfinite repriced stock')
        return stock, stock[:, None]*gradient


def path_payoffs(stock, derivative, quotes):
    sign = np.where(quotes.option_type.to_numpy() == 'call', 1., -1.)
    distance = sign*(stock[:, None]-quotes.K_over_F.to_numpy())
    h = np.maximum(distance, 0.)*quotes.forward.to_numpy()
    dh = (distance > 0)[:, :, None]*sign[None, :, None]*derivative[:, None, :]*quotes.forward.to_numpy()[None, :, None]
    return h, dh


def donor_indices(quotes, representatives):
    """Market-only nearest same-side anchor hedge; all controls retain mean zero."""
    result = []
    for q in quotes.itertuples():
        ix = np.flatnonzero(representatives.option_type.to_numpy() == q.option_type)
        if not len(ix):
            raise ValueError('Missing same-side hedge representative')
        distance = abs(np.log(representatives.K_over_F.to_numpy()[ix]/q.K_over_F))
        result.append(int(ix[np.argmin(distance)]))
    return np.array(result)


class ExploratoryJointLoss:
    """Same interval formula and family weights, extended to SIGNED MC estimates.

    Production JointPriceObjective rightly rejects negative model prices. This
    extension scores finite signed estimates without clipping, while reporting
    their invalid-price status. It supplies no qualified-price/fit clearance.
    """
    def __init__(self, quotes, forwards):
        from esn_eq.calibration_loss import JointPriceObjective
        # Reuse all production validation, including split and forward matching.
        JointPriceObjective(quotes, forwards)
        self.quotes = quotes.copy().reset_index(drop=True)
        self.forwards = forwards.copy().reset_index(drop=True)
        self.weights = np.zeros(len(quotes))
        for instrument in ['SPX', 'VIX']:
            q = self.quotes[self.quotes.instrument == instrument]
            for _, group in q.groupby('expiry'):
                self.weights[group.index] = 1/(q.expiry.nunique()*len(group))
        self.sqrtweights = np.sqrt(self.weights)

    def evaluate(self, prices, forwards, levels, price_jac=None, forward_jac=None):
        prices, forwards, levels = np.asarray(prices), np.asarray(forwards), np.asarray(levels)
        q = self.quotes
        if (prices.shape != (len(q),) or forwards.shape != (len(self.forwards),)
                or not all(np.isfinite(x).all() for x in (prices, forwards, levels))):
            raise ValueError('Invalid signed numerical predictions')
        distance = np.where(prices < q.price_bid_U, prices-q.price_bid_U,
            np.where(prices > q.price_ask_U, prices-q.price_ask_U, 0.))
        r = distance/q.price_scale_U.to_numpy()
        fr = (forwards-self.forwards.forward_vix_points.to_numpy())/.5
        # Vector representation of 0.1*sum(diff(levels)^2), same scalar penalty.
        penalty_residual = np.sqrt(.1)*np.diff(levels)
        residual = np.r_[self.sqrtweights*r, fr/np.sqrt(len(fr)), penalty_residual]
        components = {instrument: float(np.sum(self.weights[q.instrument == instrument]*r[q.instrument == instrument]**2))
            for instrument in ['SPX', 'VIX']}
        components |= dict(midpoint=0., forwards=float(np.mean(fr**2)), state_penalty=0., curve_penalty=float(np.sum(penalty_residual**2)))
        result = dict(total=float(residual@residual), components=components, residuals=residual,
            negative_count=int(np.sum(prices < 0)), invalid_forward_count=int(np.sum(forwards <= 0)))
        if price_jac is not None:
            outside = (prices < q.price_bid_U.to_numpy())|(prices > q.price_ask_U.to_numpy())
            pj = self.sqrtweights[:, None]*outside[:, None]*price_jac/q.price_scale_U.to_numpy()[:, None]
            fj = forward_jac/(.5*np.sqrt(len(fr)))
            cj = np.sqrt(.1)*np.array([[-1., 1., 0.], [0., -1., 1.]])
            result['jacobian'] = np.concatenate([pj, fj, cj])
        return result

    def tables(self, prices, errors, forwards, forward_errors):
        q, f = self.quotes.copy(), self.forwards.copy()
        q['model_price_U'], q['mc_se_U'], q['loss_weight'] = prices, errors, self.weights
        q['interval_residual'] = np.where(prices < q.price_bid_U, prices-q.price_bid_U,
            np.where(prices > q.price_ask_U, prices-q.price_ask_U, 0.))/q.price_scale_U
        q['negative_estimate'] = prices < 0
        q['nominal_se_screen'] = (prices > 0)&(errors <= q.price_scale_U.to_numpy()/8)
        f['model_forward'], f['mc_se'] = forwards, forward_errors
        f['error_points'] = forwards-f.forward_vix_points
        return q, f
