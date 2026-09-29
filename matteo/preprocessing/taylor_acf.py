"""Generate absolute-return and squared-return ACF targets at lags 1-40."""

import json

import numpy as np

from .market_data import ASSET_NAMES, PROCESSED_DIR, load_bench9, returns_and_variance


def autocorrelation_profile(values, max_lag=40):
    """Pearson correlation of a series with itself at each positive lag."""
    series = np.asarray(values, dtype=float)
    if len(series) <= max_lag + 1:
        raise ValueError("Series must be longer than max_lag + 1")
    return [float(np.corrcoef(series[:-lag], series[lag:])[0, 1])
            for lag in range(1, max_lag + 1)]


def compute_targets(frames=None, max_lag=40):
    """Return the two JSON-ready Taylor-effect curves for nine assets."""
    frames = load_bench9() if frames is None else frames
    absolute, squared = {}, {}
    for name in ASSET_NAMES:
        returns, _ = returns_and_variance(frames[name])
        absolute[name] = autocorrelation_profile(np.abs(returns), max_lag)
        squared[name] = autocorrelation_profile(returns ** 2, max_lag)
    return absolute, squared


if __name__ == "__main__":
    absolute, squared = compute_targets()
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    for name, values in (("bench9_real_abs_acf_40.json", absolute),
                         ("bench9_real_sq_acf_40.json", squared)):
        destination = PROCESSED_DIR / name
        destination.write_text(json.dumps(values, indent=2) + "\n")
        print(destination)
