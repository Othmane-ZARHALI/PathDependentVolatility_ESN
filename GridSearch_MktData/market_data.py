"""
market_data.py (reconstructed)
======================================
Loads the S&P 500 OHLC panel and builds per-asset stylised-fact
references from 6 non-overlapping time windows, per esn_market_protocol
Section 1.3: "Five randomly selected tickers from the S&P 500 panel
(list_dics_SP500_ohlcvol_19940101_filtered.pkl, 243 tickers, daily
OHLCV, 1994-2025) ... Each asset's reference is built by splitting its
full history into 6 non-overlapping time windows, computing the 11
stylised facts once per window ... Variance proxy: Garman-Klass."
"""
import pickle
import numpy as np
import pandas as pd

RAW_PKL_PATH = "/mnt/user-data/uploads/list_dics_SP500_ohlcvol_19940101_filtered.pkl"
ASSET_SAMPLE_SEED = 11  # reproduces GPC, QCOM, INTC, LUV, NI as the first 5

_raw_cache = None


def _load_raw():
    global _raw_cache
    if _raw_cache is None:
        with open(RAW_PKL_PATH, "rb") as f:
            _raw_cache = pickle.load(f)
    return _raw_cache


def _ticker_frame(entry):
    """Convert one raw dict entry (single-column Open/High/Low/Close
    DataFrames + Date index) into a plain O/H/L/C DataFrame."""
    idx = entry["Date"]
    df = pd.DataFrame({
        "O": entry["Open"].iloc[:, 0].values,
        "H": entry["High"].iloc[:, 0].values,
        "L": entry["Low"].iloc[:, 0].values,
        "C": entry["Close"].iloc[:, 0].values,
    }, index=idx)
    return df


def _stats_for_window(win, stat_func):
    """Build daily log-returns x and Garman-Klass variance v for one
    window, then apply stat_func(x, v). Returns None if the window is
    too short or degenerate."""
    if len(win) < 60:
        return None
    O, H, L, C = win["O"].values, win["H"].values, win["L"].values, win["C"].values
    v = 0.5 * (np.log(H / L)) ** 2 - (2 * np.log(2) - 1) * (np.log(C / O)) ** 2
    v = np.maximum(v, 1e-12)
    x = np.diff(np.log(np.maximum(C, 1e-12)))
    v = v[1:]  # align with x (x has len n-1)
    if len(x) < 60 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(v)):
        return None
    try:
        return stat_func(x, v)
    except Exception:
        return None


def build_market_regime_cache(stat_func, verbose=False, n_assets=5):
    """Select the first n_assets admissible tickers from the seed-11
    permutation and build each one's list of 6-window stat dicts."""
    raw = _load_raw()
    rng = np.random.default_rng(ASSET_SAMPLE_SEED)
    order = rng.permutation(len(raw))
    cache = {}
    for idx in order:
        entry = raw[idx]
        name = entry["Name"]
        df = _ticker_frame(entry)
        if len(df) < 6 * 250:
            continue
        n = len(df)
        window_len = n // 6
        stats_list = []
        for w in range(6):
            win = df.iloc[w * window_len:(w + 1) * window_len]
            st = _stats_for_window(win, stat_func)
            if st is not None:
                stats_list.append(st)
        if len(stats_list) >= 5:
            cache[name] = stats_list
            if verbose:
                print(f"{name}: {len(stats_list)} windows ok")
        if len(cache) >= n_assets:
            break
    return cache
