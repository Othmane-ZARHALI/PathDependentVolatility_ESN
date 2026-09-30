"""Audited SPX/VIX snapshot preprocessing in undiscounted forward-price units.

The functions do no I/O except ``load_quotes`` and ``run``.  Importing this
module never processes the market snapshot.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import shlex
import sys

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import ndtr


VERSION = "1.0.0"
EXCEL_EPOCH = datetime(1899, 12, 30)
REQUIRED = (
    "reference", "expiry", "callput", "strike", "exercise", "bid", "ask",
    "volume_bid", "volume_ask", "spot", "forward", "moneyness", "expiry_yf",
)
NUMERIC = tuple(c for c in REQUIRED if c != "exercise")


@dataclass(frozen=True)
class OptionConfig:
    quote_unit: str = "iv_decimal"
    day_count: int = 365
    agreement_tolerance: float = 1e-10
    short_expiry_days: int = 2
    variance_floor: float = 1e-4
    variance_margin: float = 1e-8
    price_scale_fraction: float = 1e-5
    vix_window_days: int = 30
    vix_forward_error_band: float = 0.5
    nuisance_xi_upper: float = 0.25
    nuisance_penalty: float = 0.1
    spx_final_maturity_days: int = 119
    spx_fit_maturities_days: tuple[int, ...] = (29, 57, 211, 302)
    strike_cutoff: str = "all_valid_observed_strikes_no_zero_bid_inference"

    def __post_init__(self) -> None:
        if self.quote_unit != "iv_decimal":
            raise ValueError("Only explicitly confirmed iv_decimal is supported")
        if (self.day_count <= 0 or self.variance_floor <= 0 or self.variance_margin <= 0 or
                self.price_scale_fraction <= 0 or self.vix_forward_error_band < 0 or
                self.nuisance_xi_upper <= self.variance_floor + self.variance_margin or
                self.nuisance_penalty < 0):
            raise ValueError("Invalid preprocessing configuration")


def excel_date(serial: int) -> str:
    return (EXCEL_EPOCH + timedelta(days=int(serial))).date().isoformat()


def black76_forward(forward: float, strike: float, maturity: float,
                    volatility: float, option_type: str) -> float:
    """Black-76 forward value with an intrinsic limit and stable OTM evaluation."""
    if not all(np.isfinite(v) for v in (forward, strike, maturity, volatility)):
        raise ValueError("Nonfinite Black input")
    if min(forward, strike, maturity) <= 0 or volatility < 0:
        raise ValueError("Nonpositive forward/strike/maturity or negative volatility")
    if option_type not in ("call", "put"):
        raise ValueError("Unknown option type")
    intrinsic = max(forward - strike, 0.0) if option_type == "call" else max(strike - forward, 0.0)
    total_vol = volatility * np.sqrt(maturity)
    if total_vol < 1e-12:
        return intrinsic
    d1 = (np.log(forward / strike) + 0.5 * total_vol**2) / total_vol
    d2 = d1 - total_vol
    if strike >= forward:
        otm_call = forward * ndtr(d1) - strike * ndtr(d2)
        price = otm_call if option_type == "call" else otm_call + strike - forward
    else:
        otm_put = strike * ndtr(-d2) - forward * ndtr(-d1)
        price = otm_put if option_type == "put" else otm_put + forward - strike
    return max(float(price), intrinsic)


def black76_implied_vol(forward: float, strike: float, maturity: float,
                        price: float, option_type: str, upper_vol: float = 10.0) -> float:
    """Numerical inversion for a quote-conversion check, not a data-unit guess."""
    intrinsic = max(forward - strike, 0) if option_type == "call" else max(strike - forward, 0)
    if price < intrinsic or price > (forward if option_type == "call" else strike):
        raise ValueError("Price outside elementary option bounds")
    if price == intrinsic:
        return 0.0
    if black76_forward(forward, strike, maturity, upper_vol, option_type) < price:
        raise ValueError("IV inversion bracket too small")
    return float(brentq(lambda vol: black76_forward(forward, strike, maturity, vol, option_type) - price,
                        0.0, upper_vol, xtol=1e-12))


def load_quotes(path: Path, instrument: str, config: OptionConfig) -> pd.DataFrame:
    """Read immutable CSV bytes and reject structural errors before conversion."""
    if instrument not in ("SPX", "VIX"):
        raise ValueError("instrument must be SPX or VIX")
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = set(REQUIRED) - set(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if frame.empty:
        raise ValueError(f"{path}: empty input")
    if frame.columns.duplicated().any():
        raise ValueError(f"{path}: duplicate columns")
    frame = frame.loc[:, REQUIRED].copy()
    raw = frame.copy()
    for column in NUMERIC:
        frame[column] = pd.to_numeric(frame[column].replace("", np.nan), errors="coerce")
    structural = ("reference", "expiry", "callput", "strike", "spot", "forward", "moneyness", "expiry_yf")
    for column in structural:
        if not np.isfinite(frame[column].to_numpy(dtype=float)).all():
            raise ValueError(f"{path}: invalid {column}")
    for column in ("volume_bid", "volume_ask"):
        populated = raw[column].ne("")
        if not np.isfinite(frame.loc[populated, column].to_numpy(dtype=float)).all():
            raise ValueError(f"{path}: invalid {column} type")
    for column in ("reference", "expiry"):
        if not (frame[column] == frame[column].round()).all():
            raise ValueError(f"{path}: noninteger {column}")
    if frame.reference.nunique() != 1:
        raise ValueError("A run requires one reference date per source")
    if not (frame.reference < frame.expiry).all():
        raise ValueError(f"{path}: reference must precede expiry")
    if not frame.callput.isin((-1, 1)).all():
        raise ValueError(f"{path}: callput must be +1 or -1")
    if not raw.exercise.str.lower().eq("european").all():
        raise ValueError(f"{path}: only European exercise is supported")
    for column in ("spot", "forward", "strike", "expiry_yf"):
        if not (frame[column] > 0).all():
            raise ValueError(f"{path}: {column} must be positive")
    if frame.duplicated(("reference", "expiry", "callput", "strike")).any():
        raise ValueError(f"{path}: duplicate contract key")
    for expiry, group in frame.groupby("expiry", sort=True):
        for column in ("spot", "forward"):
            values = group[column].to_numpy(dtype=float)
            if not np.allclose(values, values[0], atol=config.agreement_tolerance, rtol=0):
                raise ValueError(f"{path}: inconsistent {column} at expiry {expiry}")
    expected_t = (frame.expiry - frame.reference) / config.day_count
    if not np.allclose(frame.expiry_yf, expected_t, atol=config.agreement_tolerance, rtol=0):
        raise ValueError(f"{path}: expiry_yf disagrees with ACT/{config.day_count}")
    if not np.allclose(frame.moneyness, frame.strike / frame.forward,
                       atol=config.agreement_tolerance, rtol=0):
        raise ValueError(f"{path}: moneyness disagrees with strike/forward")
    raw.insert(0, "raw_row_id", [f"{instrument}:{i + 1}" for i in range(len(raw))])
    raw.insert(1, "instrument", instrument)
    return raw


def convert_quotes(raw: pd.DataFrame, config: OptionConfig) -> pd.DataFrame:
    """Preserve source strings and attach price intervals and exclusion reasons."""
    rows = []
    for source in raw.to_dict("records"):
        reference = int(float(source["reference"]))
        expiry = int(float(source["expiry"]))
        forward = float(source["forward"])
        strike = float(source["strike"])
        maturity = float(source["expiry_yf"])
        option_type = "call" if float(source["callput"]) == 1 else "put"
        row = dict(source)
        row.update(reference_serial=reference, expiry_serial=expiry,
                   reference_date=excel_date(reference), expiry_date=excel_date(expiry),
                   T_years=maturity, option_type=option_type, K_over_F=strike / forward,
                   source_bid=source["bid"], source_ask=source["ask"],
                   quote_unit=config.quote_unit, price_bid_U=np.nan, price_ask_U=np.nan,
                   price_mid_U=np.nan, source_width=np.nan, price_width_U=np.nan,
                   short_expiry=(expiry - reference <= config.short_expiry_days))
        reasons = []
        values = {}
        for side in ("bid", "ask"):
            raw_value = source[side]
            if raw_value == "":
                reasons.append(f"missing_{side}")
            else:
                try:
                    value = float(raw_value)
                except ValueError:
                    reasons.append(f"non_numeric_{side}")
                else:
                    if not np.isfinite(value):
                        reasons.append(f"nonfinite_{side}")
                    elif value < 0:
                        reasons.append(f"negative_{side}")
                    elif value == 0:
                        reasons.append(f"zero_{side}")
                    else:
                        values[side] = value
        if len(values) == 2:
            row["source_width"] = values["ask"] - values["bid"]
            if values["ask"] < values["bid"]:
                reasons.append("crossed_source_quote")
        tol = 1e-10 * max(forward, strike, 1)
        intrinsic = max((forward - strike) if option_type == "call" else (strike - forward), 0)
        upper = forward if option_type == "call" else strike
        for side, sigma in values.items():
            try:
                price = black76_forward(forward, strike, maturity, sigma, option_type)
            except (ValueError, OverflowError):
                reasons.append(f"black_conversion_failed_{side}")
                continue
            row[f"price_{side}_U"] = price
            if price < intrinsic - tol or price > upper + tol:
                reasons.append(f"price_bound_violation_{side}")
        if len(values) == 2 and np.isfinite(row["price_bid_U"]) and np.isfinite(row["price_ask_U"]):
            bid, ask = row["price_bid_U"], row["price_ask_U"]
            row["price_width_U"] = ask - bid
            if bid > ask + tol:
                reasons.append("crossed_price_quote")
            if not reasons:
                row["price_mid_U"] = (bid + ask) / 2
        row["validity"] = "valid" if not reasons else "excluded"
        row["reason_codes"] = "|".join(reasons)
        rows.append(row)
    return pd.DataFrame(rows)


def parity_checks(quotes: pd.DataFrame) -> pd.DataFrame:
    valid = quotes[quotes.validity == "valid"].copy()
    valid["_expiry"] = valid.expiry.astype(float).astype(int)
    valid["_strike"] = valid.strike.astype(float)
    records = []
    keys = ["instrument", "_expiry", "_strike"]
    for key, group in valid.groupby(keys, sort=True):
        if set(group.option_type) != {"call", "put"}:
            continue
        call = group[group.option_type == "call"].iloc[0]
        put = group[group.option_type == "put"].iloc[0]
        shift = float(call.forward) - float(call.strike)
        lower = max(call.price_bid_U, put.price_bid_U + shift)
        upper = min(call.price_ask_U, put.price_ask_U + shift)
        records.append(dict(instrument=key[0], expiry=key[1], strike=key[2], K_over_F=call.K_over_F,
                            call_row_id=call.raw_row_id, put_row_id=put.raw_row_id,
                            intersects=lower <= upper + 1e-10 * float(call.forward),
                            gap_U=max(lower - upper, 0.0), intersection_bid_U=lower,
                            intersection_ask_U=upper))
    return pd.DataFrame(records)


def calibration_view(quotes: pd.DataFrame, parity: pd.DataFrame,
                     config: OptionConfig) -> pd.DataFrame:
    """Select OTM intervals, intersect ATM put/call, and assign strike roles."""
    valid = quotes[quotes.validity == "valid"].copy()
    valid["_expiry"] = valid.expiry.astype(float).astype(int)
    valid["_strike"] = valid.strike.astype(float)
    records = []
    for (instrument, expiry), group in valid.groupby(["instrument", "_expiry"], sort=True):
        forward = float(group.forward.iloc[0])
        expiry_days = int(float(expiry) - float(group.reference.iloc[0]))
        maturity_role = ("fit" if instrument == "VIX" or expiry_days in config.spx_fit_maturities_days
                         else "final_maturity_test" if expiry_days == config.spx_final_maturity_days
                         else "short_expiry_diagnostic" if expiry_days <= config.short_expiry_days
                         else "maturity_diagnostic" if expiry_days < 730 else "extrapolation_diagnostic")
        for strike, pair in group.groupby("_strike", sort=True):
            k = float(strike)
            chosen = pair[(pair.option_type == ("put" if k < forward else "call"))] if k != forward else pair
            if chosen.empty:
                continue
            if k == forward and len(chosen) == 2:
                match = parity[(parity.instrument == instrument) &
                               (parity.expiry.astype(float) == float(expiry)) &
                               (parity.strike.astype(float) == k)]
                if match.empty or not bool(match.iloc[0].intersects):
                    continue
                p = match.iloc[0]
                source_ids = f"{p.put_row_id}|{p.call_row_id}"
                bid, ask = float(p.intersection_bid_U), float(p.intersection_ask_U)
                kind = "atm_intersection"
            else:
                q = chosen.iloc[0]
                source_ids = q.raw_row_id
                bid, ask = float(q.price_bid_U), float(q.price_ask_U)
                kind = q.option_type
            records.append(dict(instrument=instrument, expiry=int(float(expiry)),
                                expiry_date=excel_date(int(float(expiry))), expiry_days=expiry_days,
                                T_years=float(group.T_years.iloc[0]), strike=k, forward=forward,
                                K_over_F=k / forward, option_type=kind, source_row_ids=source_ids,
                                price_bid_U=bid, price_ask_U=ask, price_mid_U=(bid + ask) / 2,
                                maturity_role=maturity_role))
    result = pd.DataFrame(records)
    if result.empty:
        return result
    result = result.sort_values(["instrument", "expiry", "strike"], kind="stable").reset_index(drop=True)
    result["strike_rank"] = result.groupby(["instrument", "expiry"]).cumcount()
    result["split"] = np.where(result.strike_rank % 10 == 4, "validation",
                               np.where(result.strike_rank % 10 == 9, "final_test", "train"))
    result.loc[result.maturity_role == "final_maturity_test", "split"] = "final_maturity_test"
    result.loc[~result.maturity_role.isin(("fit", "final_maturity_test")), "split"] = "diagnostic"
    result["price_scale_U"] = np.maximum((result.price_ask_U - result.price_bid_U) / 2,
                                           config.price_scale_fraction * result.forward)
    return result


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, float_format="%.17g", lineterminator="\n")


def _json(path: Path, payload: dict) -> None:
    def clean(value):
        if isinstance(value, dict):
            return {key: clean(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [clean(item) for item in value]
        if isinstance(value, (float, np.floating)) and not np.isfinite(value):
            return None
        if isinstance(value, np.generic):
            return value.item()
        return value
    path.write_text(json.dumps(clean(payload), indent=2, sort_keys=True, allow_nan=False) + "\n")


def _roundtrip_audit(quotes: pd.DataFrame) -> dict:
    errors = []
    unresolved = 0
    for row in quotes[quotes.validity == "valid"].itertuples(index=False):
        for side in ("bid", "ask"):
            source = float(getattr(row, f"source_{side}"))
            price = float(getattr(row, f"price_{side}_U"))
            try:
                recovered = black76_implied_vol(float(row.forward), float(row.strike),
                                                 float(row.T_years), price, row.option_type)
            except ValueError:
                unresolved += 1
                continue
            if recovered == 0 and source > 0:
                unresolved += 1  # Intrinsic numerical limit: IV not identifiable.
            else:
                errors.append(abs(recovered - source))
    return dict(converted_sides=len(quotes[quotes.validity == "valid"]) * 2,
                identifiable_sides=len(errors), numerically_unidentifiable_sides=unresolved,
                max_abs_iv_error=max(errors, default=None), declared_tolerance=1e-7)


def _spread_audit(quotes: pd.DataFrame) -> list[dict]:
    result = []
    for (instrument, expiry), group in quotes.groupby(["instrument", "expiry_date"], sort=True):
        valid = group[group.validity == "valid"]
        entry = dict(instrument=instrument, expiry_date=expiry, source_rows=len(group),
                     retained=len(valid), excluded=len(group) - len(valid),
                     min_K_over_F=float(group.K_over_F.min()), max_K_over_F=float(group.K_over_F.max()))
        for column in ("source_width", "price_width_U"):
            values = valid[column].dropna().astype(float)
            for percentile in (0.1, 0.5, 0.9, 0.99):
                entry[f"{column}_q{int(percentile * 100):02d}"] = (
                    float(values.quantile(percentile)) if len(values) else None)
        result.append(entry)
    return result


def run(spx_path: Path, vix_path: Path, output: Path, config: OptionConfig) -> dict:
    from matteo.preprocessing.option_surface import surface_audit
    from matteo.preprocessing.option_variance import variance_outputs

    spx_path, vix_path, output = Path(spx_path), Path(vix_path), Path(output)
    source = {"SPX": spx_path, "VIX": vix_path}
    raw = {name: load_quotes(path, name, config) for name, path in source.items()}
    if raw["SPX"].reference.iloc[0] != raw["VIX"].reference.iloc[0]:
        raise ValueError("SPX and VIX reference dates disagree")
    reference_date = excel_date(int(float(raw["SPX"].reference.iloc[0])))
    quotes = pd.concat([convert_quotes(raw[name], config) for name in ("SPX", "VIX")], ignore_index=True)
    parity = parity_checks(quotes)
    calibration = calibration_view(quotes, parity, config)
    surfaces, violations = surface_audit(calibration)
    calibration = calibration.merge(surfaces[["instrument", "expiry", "view", "feasible"]]
                                    .query("view == 'all'")[["instrument", "expiry", "feasible"]],
                                    on=["instrument", "expiry"], how="left")
    calibration["surface_status"] = np.where(calibration.feasible, "interval_surface_feasible",
                                               "interval_surface_infeasible")
    calibration = calibration.drop(columns="feasible")
    training_status = surfaces[surfaces.view == "training"][["instrument", "expiry", "feasible"]].rename(
        columns={"feasible": "training_surface_feasible"})
    calibration = calibration.merge(training_status, on=["instrument", "expiry"], how="left")
    variance = variance_outputs(quotes, calibration, config)
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(quotes, output / "quotes_all.csv")
    _write_csv(calibration, output / "quotes_calibration.csv")
    _write_csv(variance["vix_forwards"], output / "vix_forwards.csv")
    for key in ("spx_variance_terms", "forward_variance", "cross_market_checks"):
        _write_csv(variance[key], output / f"{key}.csv")
    audit = dict(reference_date=reference_date, source_rows={n: len(f) for n, f in raw.items()},
                 valid_rows=quotes.groupby("instrument").validity.apply(lambda x: int((x == "valid").sum())).to_dict(),
                 reason_counts={name: int(sum(name in codes.split("|") for codes in quotes.reason_codes))
                                for name in sorted({r for codes in quotes.reason_codes for r in codes.split("|") if r})},
                 by_expiry=quotes.groupby(["instrument", "expiry_date", "validity"]).size().rename("count").reset_index().to_dict("records"),
                 parity=parity.to_dict("records"), surface_checks=surfaces.to_dict("records"),
                 surface_violations=violations.to_dict("records"),
                 spread_and_wing_summary=_spread_audit(quotes),
                 iv_roundtrip=_roundtrip_audit(quotes), variance_audit=variance["audit"])
    audit["gates"] = dict(
        row_reconciliation=all(len(raw[name]) == int((quotes.instrument == name).sum())
                               for name in raw),
        excluded_rows_have_reason=bool(quotes.loc[quotes.validity == "excluded", "reason_codes"].ne("").all()),
        iv_roundtrip_within_tolerance=(audit["iv_roundtrip"]["numerically_unidentifiable_sides"] == 0 and
                                       audit["iv_roundtrip"]["max_abs_iv_error"] <=
                                       audit["iv_roundtrip"]["declared_tolerance"]),
        parity_conflict_pairs=int((~parity.intersects).sum()) if not parity.empty else 0,
        infeasible_surfaces=int((~surfaces.feasible).sum()) if not surfaces.empty else 0,
        cross_market_status=variance["audit"]["overall_compatibility_status"])
    _json(output / "audit.json", audit)
    inputs = {}
    for name, path in source.items():
        data = path.read_bytes()
        inputs[name] = {"path": str(path.resolve()), "sha256": hashlib.sha256(data).hexdigest(),
                        "byte_count": len(data)}
    modules = (Path(__file__), Path(__file__).with_name("option_surface.py"),
               Path(__file__).with_name("option_variance.py"))
    code_hashes = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in modules}
    command = " ".join(map(shlex.quote, ("python", "-m", "matteo.preprocessing.options",
                                           "--spx", str(spx_path), "--vix", str(vix_path),
                                           "--output", str(output), "--quote-unit", config.quote_unit)))
    manifest = dict(preprocessing_version=VERSION, code_sha256=code_hashes,
                    run_time_utc=datetime.now(timezone.utc).isoformat(),
                    reference_date=reference_date, input_files=inputs, configuration=asdict(config),
                    date_conversion="Excel 1900 serial, epoch 1899-12-30; date-level only",
                    price_unit="undiscounted_forward_value", quote_unit=config.quote_unit,
                    vendor_iv_forward_and_day_count_match="unverified",
                    vix_forward_provenance="likely_traded_future_user_reported; venue unverified",
                    vix_spot_replication="approximate diagnostic only",
                    command=command)
    _json(output / "manifest.json", manifest)
    return audit


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spx", type=Path, required=True)
    parser.add_argument("--vix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quote-unit", choices=["iv_decimal"], default="iv_decimal")
    args = parser.parse_args(argv)
    run(args.spx, args.vix, args.output, OptionConfig(quote_unit=args.quote_unit))


if __name__ == "__main__":
    main(sys.argv[1:])
