"""Load the historical OHLC files and form aligned returns and variance.

The asset selection and Garman--Klass proxy follow Othmane's
``bench9_real_leverage_compute.py`` and ``bench9_real_acf_compute.py``.
Inputs live in the repository's ignored ``data/raw`` directory.
"""

from pathlib import Path
import pickle

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
STOCK_FILE = RAW_DIR / "list_dics_SP500_ohlcvol_19940101_filtered.pkl"
INDEX_FILE = RAW_DIR / "indices_data_all_quotations1990-2023.xlsx"

STOCK_NAMES = ("LHX", "BK", "HUM", "TSN", "IEX", "KMB")
INDEX_TICKERS = {"SP500": "^GSPC", "RUSSELL2000": "^RUT", "NASDAQ100": "^NDX"}
ASSET_NAMES = STOCK_NAMES + tuple(INDEX_TICKERS)
EXCLUDED_STOCKS = {
    "GPC", "QCOM", "INTC", "LUV", "NI",
    "LH", "SPGI", "DOV", "F", "UDR", "SO", "CTAS", "TGT", "BDX",
}


def _check_ohlc(frame, name):
    """Reject bad dates or prices before lagged calculations can hide them."""
    frame = frame.copy()
    frame.index = pd.to_datetime(frame.index)
    if frame.index.hasnans or frame.index.has_duplicates:
        raise ValueError(f"{name}: missing or duplicate dates")
    frame = frame.sort_index()
    frame = frame.astype(float)
    values = frame[["O", "H", "L", "C"]].to_numpy()
    if not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError(f"{name}: nonfinite or nonpositive OHLC price")
    return frame


def _stock_frame(entry):
    frame = pd.DataFrame(
        {key: entry[column].iloc[:, 0].to_numpy() for key, column in
         (("O", "Open"), ("H", "High"), ("L", "Low"), ("C", "Close"))},
        index=entry["Date"],
    )
    return _check_ohlc(frame, entry["Name"])


def load_stocks(path=STOCK_FILE):
    """Repeat the seed-2026 stock sample used for Othmane's benchmark."""
    with Path(path).open("rb") as handle:
        panel = pickle.load(handle)  # Only load a trusted local pickle.
    order = np.random.default_rng(2026).permutation(len(panel))
    selected = {}
    for idx in order:
        entry = panel[idx]
        name = entry["Name"]
        if name in EXCLUDED_STOCKS:
            continue
        frame = _stock_frame(entry)
        if len(frame) < 6 * 250:
            continue
        selected[name] = frame
        if len(selected) == len(STOCK_NAMES):
            break
    if tuple(selected) != STOCK_NAMES:
        raise ValueError(f"Stock sample changed: expected {STOCK_NAMES}, got {tuple(selected)}")
    return selected


def load_indices(path=INDEX_FILE):
    """Read the workbook's two header rows and the three index OHLC sets."""
    raw = pd.read_excel(path, sheet_name="Stock_Data", header=None)
    category = raw.iloc[0].ffill()
    ticker = raw.iloc[1]
    dates = pd.to_datetime(raw.iloc[3:, 0].to_numpy())
    out = {}
    for name, symbol in INDEX_TICKERS.items():
        columns = {}
        for label, short in (("Open", "O"), ("High", "H"),
                             ("Low", "L"), ("Close", "C")):
            matches = [col for col in range(1, raw.shape[1])
                       if category.iloc[col] == label and ticker.iloc[col] == symbol]
            if len(matches) != 1:
                raise ValueError(f"{name}: expected one {label} column, got {len(matches)}")
            columns[short] = pd.to_numeric(raw.iloc[3:, matches[0]], errors="coerce").to_numpy()
        frame = pd.DataFrame(columns, index=dates).dropna()
        frame = frame.loc[frame.index >= "1994-01-01"]
        out[name] = _check_ohlc(frame, name)
    return out


def load_bench9(stock_path=STOCK_FILE, index_path=INDEX_FILE):
    """Return the six selected stocks and three indices in benchmark order."""
    stocks = load_stocks(stock_path)
    indices = load_indices(index_path)
    return {name: (stocks if name in stocks else indices)[name] for name in ASSET_NAMES}


def garman_klass_variance(frame):
    """Daily OHLC variance proxy, including the first available day."""
    open_, high, low, close = (frame[key].to_numpy(dtype=float)
                               for key in ("O", "H", "L", "C"))
    variance = (0.5 * np.log(high / low) ** 2
                - (2 * np.log(2) - 1) * np.log(close / open_) ** 2)
    variance = np.maximum(variance, 1e-12)
    if not np.isfinite(variance).all():
        raise ValueError("Nonfinite variance")
    return variance


def returns_and_variance(frame):
    """Close-to-close log returns and same-day Garman--Klass variance.

    The first variance observation is removed to align day t's variance
    with the log return from day t-1 to day t, as in Othmane's scripts.
    """
    close = frame["C"].to_numpy(dtype=float)
    returns = np.diff(np.log(close))
    if not np.isfinite(returns).all():
        raise ValueError("Nonfinite return")
    return returns, garman_klass_variance(frame)[1:]
