"""Common quote interface for model-independent bid/ask price comparisons.

Loads existing preprocessing artifacts; never re-pools put/call IVs or modifies
training/test roles. Input compatibility flags remain attached to the snapshot.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .pricing import ForwardVarianceCurve


def interval_residuals(prices, bid, ask, scales) -> np.ndarray:
    """Signed interval distance / frozen scale. Failed prices are errors, never zero."""
    prices, bid, ask, scales = (np.asarray(v, dtype=float) for v in (prices, bid, ask, scales))
    if (len({v.shape for v in (prices, bid, ask, scales)}) != 1
            or any(not np.isfinite(v).all() for v in (prices, bid, ask, scales))
            or np.any(scales <= 0) or np.any(bid < 0) or np.any(ask < bid)
            or np.any(prices < 0)):
        raise ValueError("Finite nonnegative prices, ordered intervals and positive scales required")
    return np.where(prices < bid, prices - bid, np.where(prices > ask, prices - ask, 0.)) / scales


@dataclass(frozen=True)
class MarketTargets:
    quotes: pd.DataFrame
    vix_forwards: pd.DataFrame
    curves: pd.DataFrame
    audit: dict
    manifest: dict
    artifact_sha256: dict[str, str]

    @classmethod
    def load(cls, directory: Path):
        directory = Path(directory)
        names = ("quotes_calibration.csv", "vix_forwards.csv", "forward_variance.csv",
                 "audit.json", "manifest.json")
        hashes = {name: hashlib.sha256((directory / name).read_bytes()).hexdigest() for name in names}
        quotes, forwards, curves = (pd.read_csv(directory / name) for name in names[:3])
        audit, manifest = (json.loads((directory / name).read_text()) for name in names[3:])
        required = ["T_years", "strike", "forward", "price_bid_U", "price_ask_U", "price_scale_U"]
        if (not np.isfinite(quotes[required].to_numpy(dtype=float)).all()
                or (quotes[["T_years", "strike", "forward", "price_scale_U"]] <= 0).any().any()
                or (quotes.price_bid_U < 0).any() or (quotes.price_ask_U < quotes.price_bid_U).any()
                or not quotes.option_type.isin(["put", "call"]).all()
                or not quotes.instrument.isin(["SPX", "VIX"]).all()
                or not quotes.split.isin(["train", "validation", "final_test", "diagnostic", "final_maturity_test"]).all()):
            raise ValueError("Invalid common option targets")
        if quotes.duplicated(["instrument", "expiry", "strike", "option_type"]).any():
            raise ValueError("Duplicate option contract in common targets")
        if not np.allclose(quotes.K_over_F, quotes.strike / quotes.forward, atol=1e-10, rtol=0):
            raise ValueError("Stored moneyness differs from actual strike/forward")
        for _, group in quotes.groupby(["instrument", "expiry"]):
            if group.T_years.nunique() != 1 or group.forward.nunique() != 1:
                raise ValueError("Inconsistent maturity/forward within expiry")
        if (forwards.expiry.duplicated().any()
                or not np.isfinite(forwards[["T_years", "forward_vix_points"]]).all().all()
                or (forwards[["T_years", "forward_vix_points"]] <= 0).any().any()):
            raise ValueError("Invalid VIX forward targets")
        for expiry, group in quotes[quotes.instrument == "VIX"].groupby("expiry"):
            match = forwards[forwards.expiry == expiry]
            if (len(match) != 1 or not np.isclose(match.T_years.iloc[0], group.T_years.iloc[0], rtol=0, atol=1e-12)
                    or not np.isclose(match.forward_vix_points.iloc[0], group.forward.iloc[0], rtol=0, atol=1e-10)):
                raise ValueError("VIX option and forward targets disagree")
        if (manifest.get("price_unit") != "undiscounted_forward_value"
                or audit.get("reference_date") != manifest.get("reference_date")):
            raise ValueError("Unexpected price units or inconsistent reference dates")
        return cls(quotes, forwards, curves, audit, manifest, hashes)

    def select(self, instrument: str, split: str = "train", expiry_days=None) -> pd.DataFrame:
        """Explicit split selection; never return held-out quotes by default."""
        if instrument not in {"SPX", "VIX"} or split not in set(self.quotes.split):
            raise ValueError("Unknown instrument or split")
        rows = self.quotes[(self.quotes.instrument == instrument) & (self.quotes.split == split)]
        if expiry_days is not None:
            rows = rows[rows.expiry_days.isin(expiry_days)]
        return rows.sort_values(["expiry", "strike", "option_type"]).copy()

    def curve(self, curve_id: str) -> ForwardVarianceCurve:
        """Load a fitted diagnostic curve; compatibility flags are not cleared here."""
        rows = self.curves[self.curves.curve_id == curve_id].sort_values("start_years")
        if rows.empty or not np.isfinite(rows.xi_annual_variance).all():
            raise ValueError("Curve is absent or requires nuisance-parameter fitting")
        starts, ends = rows.start_years.to_numpy(), rows.end_years.to_numpy()
        if starts[0] != 0 or not np.array_equal(starts[1:], ends[:-1]):
            raise ValueError("Curve segments must be contiguous from date zero")
        return ForwardVarianceCurve(tuple(ends), tuple(rows.xi_annual_variance))
