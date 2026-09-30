from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from CrossSectionalData.market_vix import black76_undisc as historical_black76

from matteo.preprocessing.options import (
    OptionConfig, black76_forward, black76_implied_vol, calibration_view, convert_quotes, load_quotes,
    parity_checks, run,
)
from matteo.preprocessing.option_surface import surface_audit
from matteo.preprocessing.option_variance import (
    cboe_style_term, forward_variance_at, total_variance_at, vix_option_moment_bound,
)


def source_row(instrument="SPX", strike="100", callput="1", bid="0.2", ask="0.21",
               expiry="46101", forward="100"):
    return dict(raw_row_id=f"{instrument}:1", instrument=instrument, reference="46072",
                expiry=expiry, callput=callput, strike=strike, exercise="european",
                bid=bid, ask=ask, volume_bid="", volume_ask="", spot="100" if instrument == "SPX" else "20",
                forward=forward, moneyness=str(float(strike) / float(forward)),
                expiry_yf=str((int(expiry) - 46072) / 365))


def test_black_parity_limits_and_extreme_strikes():
    for strike in (1e-5, 80, 100, 120, 1e7):
        call = black76_forward(100, strike, 0.25, 0.4, "call")
        put = black76_forward(100, strike, 0.25, 0.4, "put")
        assert call - put == pytest.approx(100 - strike, abs=1e-7)
    assert black76_forward(100, 90, 1e-12, 1e-10, "call") == 10
    assert black76_forward(100, 110, 1e-12, 1e-10, "put") == 10
    with pytest.raises(ValueError):
        black76_forward(0, 100, 1, 0.2, "call")


def test_black_matches_historical_reference_away_from_limits():
    for maturity in (1 / 365, 30 / 365, 1.0):
        for strike in (70, 100, 140):
            for is_call in (True, False):
                actual = black76_forward(100, strike, maturity, 0.35,
                                         "call" if is_call else "put")
                assert actual == pytest.approx(
                    historical_black76(100, strike, maturity, 0.35, is_call), abs=1e-11)
                intrinsic = max((100 - strike) if is_call else (strike - 100), 0)
                if actual > intrinsic + 1e-8:
                    assert black76_implied_vol(100, strike, maturity, actual,
                                                "call" if is_call else "put") == pytest.approx(0.35, abs=1e-7)


def test_quote_reasons_and_separate_prices():
    rows = [source_row(), source_row(bid="", ask="0.3"),
            source_row(bid="0", ask="0.3"), source_row(bid="0.3", ask="0.2"),
            source_row(bid="-0.1", ask="0.2"), source_row(bid="nan", ask="0.2")]
    converted = convert_quotes(pd.DataFrame(rows), OptionConfig())
    assert converted.iloc[0].validity == "valid"
    assert converted.iloc[0].price_mid_U == pytest.approx(
        (converted.iloc[0].price_bid_U + converted.iloc[0].price_ask_U) / 2)
    assert converted.iloc[1].reason_codes == "missing_bid"
    assert np.isnan(converted.iloc[1].price_mid_U)
    assert np.isfinite(converted.iloc[1].price_ask_U)
    assert "zero_bid" in converted.iloc[2].reason_codes
    assert "crossed_source_quote" in converted.iloc[3].reason_codes
    assert "negative_bid" in converted.iloc[4].reason_codes
    assert "nonfinite_bid" in converted.iloc[5].reason_codes
    equal = convert_quotes(pd.DataFrame([source_row(bid="0.2", ask="0.2")]), OptionConfig())
    assert equal.iloc[0].validity == "valid"
    assert equal.iloc[0].price_width_U == 0


def test_schema_key_and_moneyness_rejected(tmp_path: Path):
    path = tmp_path / "quotes.csv"
    row = source_row()
    pd.DataFrame([row, row]).drop(columns=["raw_row_id", "instrument"]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="duplicate contract key"):
        load_quotes(path, "SPX", OptionConfig())
    bad = dict(row, moneyness="0.99")
    pd.DataFrame([bad]).drop(columns=["raw_row_id", "instrument"]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="moneyness"):
        load_quotes(path, "SPX", OptionConfig())


def test_parity_and_infeasible_surface():
    rows = [source_row(strike="100", callput="1"), source_row(strike="100.0", callput="-1")]
    rows[1]["raw_row_id"] = "SPX:2"
    quotes = convert_quotes(pd.DataFrame(rows), OptionConfig())
    parity = parity_checks(quotes)
    assert len(parity) == 1 and bool(parity.iloc[0].intersects)
    view = calibration_view(quotes, parity, OptionConfig())
    assert len(view) == 1 and view.iloc[0].option_type == "atm_intersection"
    surface = pd.DataFrame([
        dict(instrument="SPX", expiry=46101, strike=90.0, forward=100.0, option_type="call",
             price_bid_U=19.0, price_ask_U=20.0, price_scale_U=0.5, K_over_F=0.9,
             source_row_ids="a", split="train", maturity_role="fit"),
        dict(instrument="SPX", expiry=46101, strike=100.0, forward=100.0, option_type="call",
             price_bid_U=1.0, price_ask_U=2.0, price_scale_U=0.5, K_over_F=1.0,
             source_row_ids="b", split="train", maturity_role="fit"),
    ])
    checks, violations = surface_audit(surface)
    assert not checks.feasible.any()
    assert (checks.minimum_squared_interval_loss > 0).all()
    assert not violations.empty


def test_surface_convexity_on_uneven_strikes():
    surface = pd.DataFrame([
        dict(instrument="VIX", expiry=46099, strike=k, forward=50.0, option_type="call",
             price_bid_U=p, price_ask_U=p + 0.01, price_scale_U=0.01,
             K_over_F=k / 50, source_row_ids=str(k), split="train", maturity_role="fit")
        for k, p in ((50.0, 12.0), (60.0, 11.0), (75.0, 1.0))
    ])
    checks, _ = surface_audit(surface)
    assert not checks.feasible.any()
    assert (checks.minimum_squared_interval_loss > 0).all()


def test_vix_moment_bound_uses_monotonic_endpoint_bids():
    group = pd.DataFrame([dict(option_type="put", strike=10, price_bid_U=1),
                          dict(option_type="put", strike=15, price_bid_U=2),
                          dict(option_type="call", strike=25, price_bid_U=3),
                          dict(option_type="call", strike=30, price_bid_U=2)])
    result = vix_option_moment_bound(group, pivot=20, model_mean=20)
    assert result["put_integral_lower"] == 15
    assert result["call_integral_lower"] == 25
    assert result["moment_lower_bound"] == 480


def test_piecewise_forward_variance_evaluation():
    curve = pd.DataFrame([dict(start_years=0.0, end_years=1.0, W_start=0.0,
                               xi_annual_variance=0.04),
                          dict(start_years=1.0, end_years=2.0, W_start=0.04,
                               xi_annual_variance=0.06)])
    assert forward_variance_at(curve, 1.0) == 0.06
    assert total_variance_at(curve, 1.5) == pytest.approx(0.07)
    with pytest.raises(ValueError):
        forward_variance_at(curve, 2.1)


def test_held_out_k0_makes_training_strip_incomplete():
    rows = []
    for strike in (90, 100, 110):
        for kind in ("call", "put"):
            rows.append(dict(validity="valid", forward=105.0, T_years=0.1,
                             strike=float(strike), option_type=kind,
                             price_bid_U=1.0, price_mid_U=1.1, price_ask_U=1.2))
    result = cboe_style_term(pd.DataFrame(rows), {90.0, 110.0})
    assert result["K0"] == 100
    assert result["status"] == "incomplete_k0_pair_or_grid"


def test_miniature_snapshot_outputs(tmp_path: Path):
    config = OptionConfig()
    sources = {}
    for instrument, forward, expiry, strikes in (("SPX", "100", "46101", (90, 100, 110)),
                                                 ("VIX", "20", "46099", (18, 20, 22))):
        rows = []
        for strike in strikes:
            for callput in (-1, 1):
                rows.append(source_row(instrument, str(strike), str(callput), expiry=expiry,
                                       forward=forward))
        path = tmp_path / f"{instrument}.csv"
        pd.DataFrame(rows).drop(columns=["raw_row_id", "instrument"]).to_csv(path, index=False)
        sources[instrument] = path
    output = tmp_path / "processed"
    audit = run(sources["SPX"], sources["VIX"], output, config)
    assert audit["source_rows"] == {"SPX": 6, "VIX": 6}
    assert audit["gates"]["row_reconciliation"]
    assert audit["gates"]["iv_roundtrip_within_tolerance"]
    quotes = pd.read_csv(output / "quotes_all.csv")
    assert len(quotes) == 12
    assert set(quotes.validity) == {"valid"}
    assert len(pd.read_csv(output / "vix_forwards.csv")) == 1
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["quote_unit"] == "iv_decimal"
    assert len(manifest["input_files"]["SPX"]["sha256"]) == 64
    assert len(manifest["code_sha256"]["options.py"]) == 64
