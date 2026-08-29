"""
bench9_real_leverage_compute.py
======================================
Computes the real leverage correlation profile, corr(x_t, v_{t+L}) for
lags L=0,...,40, directly from real market data for the 9-asset bench
used throughout the leverage-effect, Taylor-effect, Zumbach-effect, and
kurtosis-effect notes: 6 stocks (LHX, BK, HUM, TSN, IEX, KMB) and 3
indices (SP500, RUSSELL2000, NASDAQ100).

Produces bench9_real_leverage_wide.json: a flat dict, asset name ->
list of 41 leverage-correlation values (lag 0 through lag 40).

STOCK SELECTION: 6 tickers from the S&P 500 panel
(list_dics_SP500_ohlcvol_19940101_filtered.pkl), drawn via a fixed
seed (2026), excluding every ticker already used elsewhere in this
note family -- the companion note's own 5 assets (GPC, QCOM, INTC,
LUV, NI, its own seed 11) and the roughness-exercise note's own 9
assets (LH, SPGI, DOV, F, UDR, SO, CTAS, TGT, BDX) -- so the 9-asset
bench is disjoint from every earlier calibration sample in this whole
project. Verified to reproduce the exact 6 tickers used throughout
this note family: LHX, BK, HUM, TSN, IEX, KMB.

INDEX SELECTION: the only 3 real indices with available OHLC data
(indices_data_all_quotations1990-2023.xlsx): SP500 (^GSPC), RUSSELL2000
(^RUT), NASDAQ100 (^NDX).

VARIANCE PROXY: Garman-Klass, v_t = 0.5*(ln(H_t/L_t))^2 -
(2*ln(2)-1)*(ln(C_t/O_t))^2, paired with daily log-return
x_t = ln(C_t/C_{t-1}) -- the same proxy used throughout this whole
note family (esn_market_protocol.tex, §1.3).

Requires: list_dics_SP500_ohlcvol_19940101_filtered.pkl and
indices_data_all_quotations1990-2023.xlsx (both under
/mnt/user-data/uploads/ in this project's own environment).
"""
import pickle
import json
import numpy as np
import pandas as pd

RAW_PKL_PATH = "/mnt/user-data/uploads/list_dics_SP500_ohlcvol_19940101_filtered.pkl"
INDICES_XLSX_PATH = "/mnt/user-data/uploads/indices_data_all_quotations1990-2023.xlsx"

STOCK_SAMPLE_SEED = 2026
EXCLUDED_TICKERS = {
    "GPC", "QCOM", "INTC", "LUV", "NI",              # main protocol (seed 11)
    "LH", "SPGI", "DOV", "F", "UDR", "SO", "CTAS", "TGT", "BDX",  # roughness exercise
}
N_STOCKS = 6
MAX_LAG = 40
NAMES9 = ["LHX", "BK", "HUM", "TSN", "IEX", "KMB",
          "SP500", "RUSSELL2000", "NASDAQ100"]


def _ticker_frame(entry):
    """Convert one raw dict entry (single-column Open/High/Low/Close
    DataFrames + Date index) into a plain O/H/L/C DataFrame."""
    idx = entry["Date"]
    return pd.DataFrame({
        "O": entry["Open"].iloc[:, 0].values,
        "H": entry["High"].iloc[:, 0].values,
        "L": entry["Low"].iloc[:, 0].values,
        "C": entry["Close"].iloc[:, 0].values,
    }, index=idx)


def select_six_stocks():
    with open(RAW_PKL_PATH, "rb") as f:
        raw = pickle.load(f)
    rng = np.random.default_rng(STOCK_SAMPLE_SEED)
    order = rng.permutation(len(raw))
    selected = {}
    for idx in order:
        entry = raw[idx]
        name = entry["Name"]
        if name in EXCLUDED_TICKERS:
            continue
        df = _ticker_frame(entry)
        if len(df) < 6 * 250:
            continue
        selected[name] = df
        if len(selected) >= N_STOCKS:
            break
    return selected


def load_three_indices():
    """SP500 (^GSPC), RUSSELL2000 (^RUT), NASDAQ100 (^NDX) -- the only
    3 real indices with OHLC data available. The workbook's own header
    spans rows 0-1 (category e.g. 'Open'/'High'/'Low'/'Close' merged
    across 3 ticker columns each, row 0; ticker, row 1), with column 0
    the date (label 'Date' in row 2, data from row 3)."""
    xls = pd.ExcelFile(INDICES_XLSX_PATH)
    raw = pd.read_excel(xls, sheet_name=0, header=None)
    category_row = raw.iloc[0].ffill()
    ticker_row = raw.iloc[1]
    dates = pd.to_datetime(raw.iloc[3:, 0].values)

    ticker_map = {"SP500": "^GSPC", "RUSSELL2000": "^RUT", "NASDAQ100": "^NDX"}
    out = {}
    for name, ticker in ticker_map.items():
        cols = {}
        for cat in ["Open", "High", "Low", "Close"]:
            for col_idx in range(1, raw.shape[1]):
                if category_row[col_idx] == cat and ticker_row[col_idx] == ticker:
                    cols[cat[0]] = raw.iloc[3:, col_idx].astype(float).values
                    break
        df = pd.DataFrame(dict(O=cols["O"], H=cols["H"], L=cols["L"], C=cols["C"]),
                           index=dates).dropna()
        # Restrict to 1994 onward, matching the S&P 500 stock panel's own
        # start date (list_dics_SP500_ohlcvol_19940101_filtered.pkl) --
        # the indices file itself starts in 1990, but every other note in
        # this family uses the 1994-2023 overlap so index and stock
        # history lengths are comparable.
        df = df[df.index >= "1994-01-01"]
        out[name] = df
    return out


def real_x_v(df):
    """Daily log-return x_t and Garman-Klass variance v_t, aligned
    (v_t uses day t's own O/H/L/C; x_t = log-return into day t)."""
    O, H, L, C = df["O"].values, df["H"].values, df["L"].values, df["C"].values
    v = 0.5 * (np.log(H / L)) ** 2 - (2 * np.log(2) - 1) * (np.log(C / O)) ** 2
    v = np.maximum(v, 1e-12)
    x = np.diff(np.log(np.maximum(C, 1e-12)))
    v = v[1:]  # align with x (length n-1)
    return x, v


def leverage_profile(x, v, lags):
    """corr(x_t, v_{t+L}) for each lag L (L=0 is the contemporaneous
    correlation)."""
    out = []
    for lag in lags:
        if lag == 0:
            out.append(float(np.corrcoef(x, v)[0, 1]))
        else:
            out.append(float(np.corrcoef(x[:-lag], v[lag:])[0, 1]))
    return out


if __name__ == "__main__":
    stocks = select_six_stocks()
    print("Selected stocks:", list(stocks.keys()))
    indices = load_three_indices()
    print("Loaded indices:", list(indices.keys()))

    all_dfs = dict(stocks)
    all_dfs.update(indices)

    lags = list(range(0, MAX_LAG + 1))
    results = {}
    for name in NAMES9:
        x, v = real_x_v(all_dfs[name])
        results[name] = leverage_profile(x, v, lags)
        print(f"{name}: n={len(x)}  leverage(lag0)={results[name][0]:.4f}")

    json.dump(results, open("bench9_real_leverage_wide.json", "w"))
    print("\nSaved bench9_real_leverage_wide.json")
