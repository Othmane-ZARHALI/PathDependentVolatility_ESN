"""Reconstruct the real moment-scaling Hurst targets for the nine assets."""

import json

import numpy as np

from .market_data import ASSET_NAMES, PROCESSED_DIR, garman_klass_variance, load_bench9


MOMENT_ORDERS = (0.5, 1.0, 1.5, 2.0)
REAL_LAGS = tuple(range(2, 41))  # Exclude lag 1 for the noisy real volatility proxy.


def estimate_hurst(log_volatility, lags=REAL_LAGS, orders=MOMENT_ORDERS):
    """Slope in q of log-log structure-function slopes zeta(q)."""
    series = np.asarray(log_volatility, dtype=float)
    lags = np.asarray(lags, dtype=int)
    if len(series) <= int(lags.max()) or not np.isfinite(series).all():
        raise ValueError("Log-volatility series is too short or nonfinite")
    log_lags = np.log(lags)
    zeta = []
    for order in orders:
        moments = [np.mean(np.abs(series[lag:] - series[:-lag]) ** order)
                   for lag in lags]
        if min(moments) <= 0:
            raise ValueError("A structure-function moment is zero")
        zeta.append(np.polyfit(log_lags, np.log(moments), 1)[0])
    return float(np.polyfit(orders, zeta, 1)[0])


def compute_targets(frames=None):
    """Estimate H from the full daily log-Garman--Klass-volatility series."""
    frames = load_bench9() if frames is None else frames
    return {name: estimate_hurst(0.5 * np.log(garman_klass_variance(frames[name])))
            for name in ASSET_NAMES}


if __name__ == "__main__":
    destination = PROCESSED_DIR / "bench9_real_hurst.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(compute_targets(), indent=2) + "\n")
    print(destination)
