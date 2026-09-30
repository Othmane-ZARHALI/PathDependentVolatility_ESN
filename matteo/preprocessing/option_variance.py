"""Observed-strike SPX variance diagnostics and conservative VIX moment gates."""

from __future__ import annotations

import numpy as np
import pandas as pd


def _strike_spacing(strikes: np.ndarray) -> np.ndarray:
    if len(strikes) < 2:
        return np.zeros(len(strikes))
    result = np.empty(len(strikes))
    result[0] = strikes[1] - strikes[0]
    result[-1] = strikes[-1] - strikes[-2]
    result[1:-1] = (strikes[2:] - strikes[:-2]) / 2
    return result


def cboe_style_term(quotes: pd.DataFrame, allowed_strikes: set[float] | None = None) -> dict:
    """Date-level, observed-grid diagnostic; no zero-bid cutoff is inferred."""
    forward = float(quotes.forward.iloc[0])
    listed = np.sort(quotes.strike.astype(float).unique())
    below = listed[listed <= forward]
    if not len(below):
        return dict(status="no_K0", variance_mid=np.nan, variance_bid=np.nan,
                    variance_ask=np.nan, total_variance_mid=np.nan, quote_count=0)
    k0 = float(below[-1])
    valid = quotes[quotes.validity == "valid"]
    if allowed_strikes is not None:
        valid = valid[valid.strike.astype(float).isin(allowed_strikes)]
    if valid.empty:
        return dict(status="no_valid_quotes", variance_mid=np.nan, variance_bid=np.nan,
                    variance_ask=np.nan, total_variance_mid=np.nan, quote_count=0)
    maturity = float(valid.T_years.iloc[0])
    strike_values = np.sort(valid.strike.astype(float).unique())
    included = []
    for strike in strike_values:
        kind = "put" if strike < k0 else "call" if strike > k0 else "at_k0"
        rows = valid[valid.strike.astype(float) == strike]
        if kind == "at_k0":
            if set(rows.option_type) != {"call", "put"}:
                continue
            price = rows[["price_bid_U", "price_mid_U", "price_ask_U"]].astype(float).mean().to_dict()
        else:
            selected = rows[rows.option_type == kind]
            if selected.empty:
                continue
            price = selected.iloc[0][["price_bid_U", "price_mid_U", "price_ask_U"]].astype(float).to_dict()
        included.append((float(strike), price))
    strikes = np.asarray([x[0] for x in included])
    has_k0 = k0 in strikes
    if len(strikes) < 2 or not has_k0:
        return dict(status="incomplete_k0_pair_or_grid", K0=k0, variance_mid=np.nan,
                    variance_bid=np.nan, variance_ask=np.nan, total_variance_mid=np.nan,
                    quote_count=len(strikes), strike_min=float(strikes.min()) if len(strikes) else np.nan,
                    strike_max=float(strikes.max()) if len(strikes) else np.nan)
    spacings = _strike_spacing(strikes)
    correction = ((forward / k0) - 1) ** 2 / maturity
    variances = {}
    for side, field in (("bid", "price_bid_U"), ("mid", "price_mid_U"), ("ask", "price_ask_U")):
        prices = np.asarray([x[1][field] for x in included], dtype=float)
        variances[f"variance_{side}"] = float(2 / maturity * np.sum(spacings / strikes**2 * prices) - correction)
    gaps = np.diff(strikes)
    return dict(status="diagnostic_complete", K0=k0, **variances,
                total_variance_mid=maturity * variances["variance_mid"],
                quote_count=len(strikes), strike_min=float(strikes[0]), strike_max=float(strikes[-1]),
                max_strike_gap=float(gaps.max()), wing_truncation="observed_grid_endpoints",
                interval_extremes_surface_consistent="unverified")


def vix_option_moment_bound(group: pd.DataFrame, pivot: float, model_mean: float) -> dict:
    """Monotonicity lower bound, using covered finite-strike segments only."""
    puts = group[(group.option_type == "put") & (group.strike < pivot)].sort_values("strike")
    calls = group[(group.option_type == "call") & (group.strike > pivot)].sort_values("strike")
    put_integral = 0.0
    call_integral = 0.0
    if not puts.empty:
        strikes = puts.strike.to_numpy(dtype=float)
        bids = puts.price_bid_U.to_numpy(dtype=float)
        endpoints = np.r_[strikes[1:], pivot]
        put_integral = float(np.sum((endpoints - strikes) * bids))
    if not calls.empty:
        strikes = calls.strike.to_numpy(dtype=float)
        bids = calls.price_bid_U.to_numpy(dtype=float)
        starts = np.r_[pivot, strikes[:-1]]
        call_integral = float(np.sum((strikes - starts) * bids))
    bound = 2 * pivot * model_mean - pivot**2 + 2 * (put_integral + call_integral)
    return dict(moment_lower_bound=max(bound, 0.0), put_integral_lower=put_integral,
                call_integral_lower=call_integral, put_segments=len(puts), call_segments=len(calls),
                omitted_put_tail=f"[0,{float(puts.strike.min())}]" if not puts.empty else "all",
                omitted_call_tail=f"[{float(calls.strike.max())},infinity)" if not calls.empty else "all")


def _piecewise_curve(terms: pd.DataFrame, side: str, horizon: float, config) -> list[dict]:
    """Exact W interpolation and constant-last-forward-variance extrapolation."""
    terms = terms.sort_values("T_years")
    knots = np.r_[0.0, terms.T_years.to_numpy(dtype=float)]
    total = np.r_[0.0, terms.T_years.to_numpy(dtype=float) *
                  terms[f"variance_{side}"].to_numpy(dtype=float)]
    slopes = np.diff(total) / np.diff(knots)
    if not np.isfinite(slopes).all() or (slopes < config.variance_floor + config.variance_margin).any():
        return []
    curve = []
    for i, (left, right, w_left, slope) in enumerate(zip(knots[:-1], knots[1:], total[:-1], slopes)):
        curve.append(dict(curve_id=f"training_{side}_piecewise_v1", segment=i,
                          start_years=float(left), end_years=float(right),
                          W_start=float(w_left), xi_annual_variance=float(slope),
                          basis="constant_xi_on_segment", status="fitted_training_only",
                          floor=config.variance_floor, margin=config.variance_margin,
                          knot_derivative="right_continuous_except_H_left",
                          extrapolation="last_fitted_slope_constant_to_H",
                          smoothness_penalty=0.0, integration_error=0.0,
                          horizon_years=horizon, source_expiry_days="|".join(map(str, terms.expiry_days)),
                          source_split="train_only"))
    if horizon > knots[-1]:
        curve.append(dict(curve_id=f"training_{side}_piecewise_v1", segment=len(slopes),
                          start_years=float(knots[-1]), end_years=horizon,
                          W_start=float(total[-1]), xi_annual_variance=float(slopes[-1]),
                          basis="constant_xi_on_segment", status="extrapolation_sensitivity_required",
                          floor=config.variance_floor, margin=config.variance_margin,
                          knot_derivative="right_continuous_except_H_left",
                          extrapolation="last_fitted_slope_constant_to_H",
                          smoothness_penalty=0.0, integration_error=0.0,
                          horizon_years=horizon, source_expiry_days="|".join(map(str, terms.expiry_days)),
                          source_split="train_only"))
    frame = pd.DataFrame(curve)
    if any(abs(total_variance_at(frame, float(t)) - float(w)) > 1e-12
           for t, w in zip(knots, total)):
        raise ValueError("Forward-variance integration failed to recover training terms")
    return curve


def forward_variance_at(curve: pd.DataFrame, time: float) -> float:
    """Evaluate the right-continuous derivative of a fitted CSV curve."""
    if "curve_id" in curve and curve.curve_id.nunique() != 1:
        raise ValueError("Select exactly one curve_id before evaluation")
    segment = curve[(curve.start_years <= time) & (time <= curve.end_years)]
    if segment.empty or "xi_annual_variance" not in curve:
        raise ValueError("Fitted curve does not cover requested maturity")
    return float(segment.iloc[-1].xi_annual_variance)


def total_variance_at(curve: pd.DataFrame, time: float) -> float:
    """Integrate the fitted forward-variance curve exactly from zero."""
    if "curve_id" in curve and curve.curve_id.nunique() != 1:
        raise ValueError("Select exactly one curve_id before evaluation")
    segment = curve[(curve.start_years <= time) & (time <= curve.end_years)]
    if segment.empty or "xi_annual_variance" not in curve:
        raise ValueError("Fitted curve does not cover requested maturity")
    row = segment.iloc[-1]
    return float(row.W_start + row.xi_annual_variance * (time - row.start_years))


def _nuisance_family(horizon: float, config) -> pd.DataFrame:
    """Low-dimensional profile family; coefficients are deliberately unset."""
    breakpoints = sorted(set(min(x, horizon) for x in (0.0, 0.25, 0.75, horizon)))
    return pd.DataFrame([
        dict(curve_id="nuisance_piecewise_constant_v1", segment=i,
             start_years=left, end_years=right, basis="constant_xi_on_segment",
             coefficient=f"xi_{i}", coefficient_lower=config.variance_floor + config.variance_margin,
             coefficient_upper=config.nuisance_xi_upper,
             knot_derivative="right_continuous_except_H_left",
             extrapolation="last_segment_constant_to_H", W_start=0.0 if i == 0 else np.nan,
             penalty=f"{config.nuisance_penalty}*sum((xi_next-xi_current)^2)",
             status="unfitted_requires_calibration_profile", horizon_years=horizon,
             source_split="train_only")
        for i, (left, right) in enumerate(zip(breakpoints[:-1], breakpoints[1:]))])


def variance_outputs(quotes: pd.DataFrame, calibration: pd.DataFrame, config) -> dict:
    valid = quotes[quotes.validity == "valid"]
    vix = quotes[quotes.instrument == "VIX"]
    forward_rows = []
    for expiry, group in vix.groupby("expiry_serial", sort=True):
        forward_rows.append(dict(expiry=int(float(expiry)), expiry_date=group.expiry_date.iloc[0],
                                 expiry_days=int(float(expiry) - float(group.reference.iloc[0])),
                                 T_years=float(group.T_years.iloc[0]),
                                 forward_vix_points=float(group.forward.iloc[0]),
                                 source_provenance="likely_traded_future_user_reported",
                                 quote_interval_available=False))
    forwards = pd.DataFrame(forward_rows)
    spx = quotes[quotes.instrument == "SPX"]
    terms = []
    for expiry, group in spx.groupby("expiry_serial", sort=True):
        days = int(float(expiry) - float(group.reference.iloc[0]))
        full = cboe_style_term(group)
        full.update(expiry=int(float(expiry)), expiry_date=group.expiry_date.iloc[0],
                    expiry_days=days, T_years=float(group.T_years.iloc[0]),
                    view="all_quote_diagnostic", training_source="all_valid_two_sided")
        terms.append(full)
        selected = calibration[(calibration.instrument == "SPX") &
                               (calibration.expiry == int(float(expiry))) &
                               (calibration.split == "train") & (calibration.maturity_role == "fit")]
        if selected.empty:
            continue
        allowed = set(selected.strike.astype(float))
        training = cboe_style_term(group, allowed)
        training.update(expiry=int(float(expiry)), expiry_date=group.expiry_date.iloc[0],
                        expiry_days=days, T_years=float(group.T_years.iloc[0]),
                        view="training_only", training_source="train_strikes_both_option_types_only")
        terms.append(training)
    terms_frame = pd.DataFrame(terms)
    horizon = max(float(spx.T_years.astype(float).max()),
                  float(vix.T_years.astype(float).max()) + config.vix_window_days / config.day_count)
    training_terms = terms_frame[(terms_frame.view == "training_only") &
                                 (terms_frame.status == "diagnostic_complete")].copy()
    training_expiries = set(training_terms.expiry_days)
    expected_expiries = set(config.spx_fit_maturities_days)
    training_surfaces_ok = bool(calibration[(calibration.instrument == "SPX") &
                                            (calibration.maturity_role == "fit")]
                                .training_surface_feasible.eq(True).all())
    curve_rows = []
    if training_expiries == expected_expiries and training_surfaces_ok:
        for side in ("bid", "mid", "ask"):
            curve_rows.extend(_piecewise_curve(training_terms, side, horizon, config))
    if curve_rows:
        forward_variance = pd.DataFrame(curve_rows)
        curve_status = "fitted_training_only_with_tail_sensitivity_and_moment_gates"
    else:
        forward_variance = _nuisance_family(horizon, config)
        curve_status = "unfitted_nuisance_family_requires_profile"
    checks = []
    for _, target in forwards.iterrows():
        expiry = int(target.expiry)
        quoted = float(target.forward_vix_points)
        training = calibration[(calibration.instrument == "VIX") &
                               (calibration.expiry == expiry) & (calibration.split == "train")]
        for curve_id, curve in forward_variance.groupby("curve_id", sort=True):
            fitted = curve_id.startswith("training_")
            window = (10000 / (config.vix_window_days / config.day_count) *
                      (total_variance_at(curve, float(target.T_years) + config.vix_window_days / config.day_count)
                       - total_variance_at(curve, float(target.T_years)))) if fitted else np.nan
            for error in (-config.vix_forward_error_band, 0.0, config.vix_forward_error_band):
                mean = quoted + error
                bound = vix_option_moment_bound(training, quoted, mean)
                jensen_gap = window - mean**2 if fitted else np.nan
                status = ("unresolved_curve_not_fitted" if not fitted else
                          "moment_lower_bound_exceeds_curve" if window + 1e-8 < bound["moment_lower_bound"] else
                          "forward_jensen_bound_exceeds_curve" if jensen_gap < -1e-8 else
                          "necessary_moment_gate_passes")
                checks.append(dict(expiry=expiry, expiry_date=target.expiry_date,
                                   curve_id=curve_id, quote_split="training_only", quoted_forward=quoted,
                                   assumed_model_mean=mean, forward_error=error, pivot=quoted,
                                   **bound, window_second_moment=window,
                                   moment_gap=window - bound["moment_lower_bound"] if fitted else np.nan,
                                   forward_jensen_gap=jensen_gap,
                                   compatibility_status=status))
    checks_frame = pd.DataFrame(checks)
    conflicts_at_quoted_forward = int(((checks_frame.forward_error == 0) &
                                       checks_frame.compatibility_status.str.endswith("exceeds_curve")).sum())
    conflict_expiries = sorted(checks_frame.loc[
        (checks_frame.forward_error == 0) &
        checks_frame.compatibility_status.str.endswith("exceeds_curve"), "expiry"].astype(int).unique().tolist())
    conflict_days = forwards.set_index("expiry").loc[conflict_expiries, "expiry_days"].astype(int).tolist()
    if conflicts_at_quoted_forward and curve_rows:
        nuisance = _nuisance_family(horizon, config)
        forward_variance = pd.concat((forward_variance, nuisance), ignore_index=True)
        extra = []
        for _, target in forwards.iterrows():
            quoted = float(target.forward_vix_points)
            training = calibration[(calibration.instrument == "VIX") &
                                   (calibration.expiry == int(target.expiry)) &
                                   (calibration.split == "train")]
            for error in (-config.vix_forward_error_band, 0.0, config.vix_forward_error_band):
                bound = vix_option_moment_bound(training, quoted, quoted + error)
                extra.append(dict(expiry=int(target.expiry), expiry_date=target.expiry_date,
                                  curve_id="nuisance_piecewise_constant_v1", quote_split="training_only",
                                  quoted_forward=quoted, assumed_model_mean=quoted + error,
                                  forward_error=error, pivot=quoted, **bound,
                                  window_second_moment=np.nan, moment_gap=np.nan,
                                  forward_jensen_gap=np.nan,
                                  compatibility_status="unresolved_curve_not_fitted"))
        checks_frame = pd.concat((checks_frame, pd.DataFrame(extra)), ignore_index=True)
        curve_status = "fitted_curves_conflict_nuisance_profile_required"
    diagnostic = terms_frame[(terms_frame.view == "all_quote_diagnostic") &
                             (terms_frame.status == "diagnostic_complete")].sort_values("expiry_days")
    stored_vix_spot = float(vix.spot.iloc[0])
    spot_comparison = {"status": "unavailable", "stored_vix_spot": stored_vix_spot}
    if (diagnostic.expiry_days == 29).any() and (diagnostic.expiry_days == 57).any():
        near = diagnostic[diagnostic.expiry_days.isin((29, 57))].set_index("expiry_days")
        w30 = float(near.loc[29, "total_variance_mid"] +
                    (near.loc[57, "total_variance_mid"] - near.loc[29, "total_variance_mid"]) / 28)
        if w30 > 0:
            spot_comparison = dict(status="approximate_date_level_interpolation_only",
                                   approximate_vix_30d=100 * np.sqrt(w30 / (30 / config.day_count)),
                                   stored_vix_spot=stored_vix_spot)
    return dict(vix_forwards=forwards, spx_variance_terms=terms_frame,
                forward_variance=forward_variance, cross_market_checks=checks_frame,
                audit=dict(spot_comparison=spot_comparison,
                           curve_status=curve_status,
                           cross_market_conflicts_at_quoted_forward=conflicts_at_quoted_forward,
                           cross_market_conflict_expiries=conflict_expiries,
                           cross_market_conflict_expiry_days=conflict_days,
                           overall_compatibility_status=(
                               "unresolved_input_conflict_requires_nuisance_curve_profile"
                               if conflicts_at_quoted_forward else "necessary_gates_passed_only"),
                           variance_strip_note="actual observed strikes; all valid two-sided quotes; no inferred zero-bid cutoff",
                           moment_bound_note="finite-strike monotonicity lower bound; omitted tails nonnegative"))
