"""
bench9_real_acf_compute.py
======================================
Computes the real autocorrelation of |x_t| and x_t^2 (the Taylor
effect's own 2 underlying curves) at lags L=1,...,40, directly from
real market data for the same 9-asset bench used throughout the
leverage-effect, Taylor-effect, Zumbach-effect, and kurtosis-effect
notes: 6 stocks (LHX, BK, HUM, TSN, IEX, KMB) and 3 indices (SP500,
RUSSELL2000, NASDAQ100).

Produces 2 files, each a flat dict (asset name -> list of 40 values,
lag 1 through lag 40; lag 0, the trivial acf(x,0)=1, is excluded, as
throughout the Taylor-effect note):
  - bench9_real_abs_acf_40.json : acf(|x_t|, L)
  - bench9_real_sq_acf_40.json  : acf(x_t^2, L)

STOCK/INDEX SELECTION: identical to
bench9_real_leverage_compute.py -- 6 stocks via a fixed seed (2026)
from the S&P 500 panel, excluding every ticker already used elsewhere
in this note family, plus the 3 real indices with available OHLC data
(SP500=^GSPC, RUSSELL2000=^RUT, NASDAQ100=^NDX), restricted to 1994
onward to match the stock panel's own start date. Reproduces the
exact 9 assets used throughout: LHX, BK, HUM, TSN, IEX, KMB, SP500,
RUSSELL2000, NASDAQ100.

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
        # start date -- the indices file itself starts in 1990, but every
        # other note in this family uses the 1994-2023 overlap so index
        # and stock history lengths are comparable.
        df = df[df.index >= "1994-01-01"]
        out[name] = df
    return out


def real_x(df):
    """Daily log-return x_t."""
    C = df["C"].values
    return np.diff(np.log(np.maximum(C, 1e-12)))


def acf(series, lags):
    """Sample autocorrelation of `series` at each lag in `lags`."""
    return [float(np.corrcoef(series[:-lag], series[lag:])[0, 1]) for lag in lags]


if __name__ == "__main__":
    stocks = select_six_stocks()
    print("Selected stocks:", list(stocks.keys()))
    indices = load_three_indices()
    print("Loaded indices:", list(indices.keys()))

    all_dfs = dict(stocks)
    all_dfs.update(indices)

    lags = list(range(1, MAX_LAG + 1))
    abs_results = {}
    sq_results = {}
    for name in NAMES9:
        x = real_x(all_dfs[name])
        abs_results[name] = acf(np.abs(x), lags)
        sq_results[name] = acf(x ** 2, lags)
        print(f"{name}: n={len(x)}  acf(|x|,1)={abs_results[name][0]:.4f}  "
              f"acf(x^2,1)={sq_results[name][0]:.4f}")

    json.dump(abs_results, open("bench9_real_abs_acf_40.json", "w"))
    json.dump(sq_results, open("bench9_real_sq_acf_40.json", "w"))
    print("\nSaved bench9_real_abs_acf_40.json and bench9_real_sq_acf_40.json")
