"""Finite-grid interval feasibility and loss benchmarks for OTM option views."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, linprog, minimize
from scipy.sparse import coo_matrix, diags, vstack


def _call_intervals(group: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    strikes = group.strike.to_numpy(dtype=float)
    forwards = group.forward.to_numpy(dtype=float)
    shift = np.where(group.option_type.eq("put"), forwards - strikes, 0.0)
    return (strikes, group.price_bid_U.to_numpy(dtype=float) + shift,
            group.price_ask_U.to_numpy(dtype=float) + shift,
            group.price_scale_U.to_numpy(dtype=float))


def _curve_constraints(strikes: np.ndarray) -> tuple[coo_matrix, np.ndarray, np.ndarray]:
    n = len(strikes)
    if n < 2:
        return coo_matrix((0, n)).tocsr(), np.empty(0), np.empty(0)
    gaps = np.diff(strikes)
    rows, cols, data, lower, upper = [], [], [], [], []
    for i, gap in enumerate(gaps):
        rows.extend((len(lower), len(lower)))
        cols.extend((i, i + 1))
        data.extend((-1 / gap, 1 / gap))
        lower.append(-1.0)
        upper.append(0.0)
    for i in range(n - 2):
        left, right = gaps[i], gaps[i + 1]
        row = len(lower)
        rows.extend((row, row, row))
        cols.extend((i, i + 1, i + 2))
        data.extend((1 / left, -1 / left - 1 / right, 1 / right))
        lower.append(0.0)
        upper.append(np.inf)
    matrix = coo_matrix((data, (rows, cols)), shape=(len(lower), n)).tocsr()
    return matrix, np.asarray(lower), np.asarray(upper)


def _minimum_violation(strikes: np.ndarray, lower: np.ndarray, upper: np.ndarray,
                       scales: np.ndarray, forward: float) -> tuple[float, np.ndarray, str]:
    """Solve the fixed-scale squared interval-distance problem on a convex call grid."""
    matrix, lb, ub = _curve_constraints(strikes)
    intrinsic = np.maximum(forward - strikes, 0.0)
    start = intrinsic.copy()  # Feasible under all discrete no-arbitrage constraints.
    n = len(strikes)

    def value_and_grad(curve: np.ndarray) -> tuple[float, np.ndarray]:
        residual = np.minimum(curve - lower, 0.0) + np.maximum(curve - upper, 0.0)
        scaled = residual / scales
        return float(np.mean(scaled**2)), 2 * residual / (n * scales**2)

    result = minimize(lambda x: value_and_grad(x), start, jac=True, method="SLSQP",
                      bounds=Bounds(intrinsic, np.full(n, forward)),
                      constraints=[LinearConstraint(matrix, lb, ub)],
                      options={"ftol": 1e-10, "maxiter": 500})
    if not result.success:
        def hessian(curve: np.ndarray):
            outside = (curve < lower) | (curve > upper)
            return diags(2 * outside / (n * scales**2), format="csr")

        result = minimize(lambda x: value_and_grad(x), start, jac=True,
                          hess=hessian, method="trust-constr",
                          bounds=Bounds(intrinsic, np.full(n, forward)),
                          constraints=[LinearConstraint(matrix, lb, ub)],
                          options={"gtol": 1e-10, "xtol": 1e-12,
                                   "maxiter": 1000, "sparse_jacobian": True})
    if not result.success:
        return np.nan, result.x, f"optimizer_failed:{result.message}"
    return float(result.fun), result.x, "optimal"


def surface_audit(calibration: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Check full diagnostic and fit-training grids without changing any quotes."""
    records, violations = [], []
    if calibration.empty:
        return pd.DataFrame(), pd.DataFrame()
    for (instrument, expiry), all_group in calibration.groupby(["instrument", "expiry"], sort=True):
        views = [("all", all_group)]
        train = all_group[(all_group.split == "train") & (all_group.maturity_role == "fit")]
        if not train.empty:
            views.append(("training", train))
        for view, group in views:
            group = group.sort_values("strike")
            strikes, lower, upper, scales = _call_intervals(group)
            forward = float(group.forward.iloc[0])
            matrix, slope_lb, slope_ub = _curve_constraints(strikes)
            intrinsic = np.maximum(forward - strikes, 0.0)
            finite_upper = np.isfinite(slope_ub)
            finite_lower = np.isfinite(slope_lb)
            inequality = vstack((matrix[finite_upper], -matrix[finite_lower])).tocsr()
            rhs = np.concatenate((slope_ub[finite_upper], -slope_lb[finite_lower]))
            price_lb = np.maximum(lower, intrinsic)
            price_ub = np.minimum(upper, forward)
            lp = (linprog(np.zeros(len(strikes)), A_ub=inequality, b_ub=rhs,
                          bounds=list(zip(price_lb, price_ub)), method="highs")
                  if np.all(price_lb <= price_ub) else None)
            feasible = bool(lp is not None and lp.success)
            if feasible:
                loss, curve, optimizer = 0.0, lp.x, "feasible_lp"
            else:
                loss, curve, optimizer = _minimum_violation(strikes, lower, upper, scales, forward)
            records.append(dict(instrument=instrument, expiry=int(expiry), view=view,
                                quote_count=len(group), feasible=feasible,
                                minimum_squared_interval_loss=loss,
                                loss_status=optimizer))
            if not feasible:
                residual = np.maximum(lower - curve, 0) + np.maximum(curve - upper, 0)
                order = np.argsort(-np.abs(residual / scales), kind="stable")
                for rank, idx in enumerate(order[:20], 1):
                    if residual[idx] <= 1e-10 * forward:
                        break
                    row = group.iloc[int(idx)]
                    violations.append(dict(instrument=instrument, expiry=int(expiry), view=view,
                                           rank=rank, strike=float(strikes[idx]), K_over_F=float(row.K_over_F),
                                           source_row_ids=row.source_row_ids,
                                           gap_U=float(residual[idx]), scaled_gap=float(residual[idx] / scales[idx])))
    return pd.DataFrame(records), pd.DataFrame(violations)
