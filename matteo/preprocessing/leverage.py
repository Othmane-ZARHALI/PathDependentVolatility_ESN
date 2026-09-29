"""Generate the nine-asset real leverage profiles at lags 0 through 40."""

import json

import numpy as np

from .market_data import ASSET_NAMES, PROCESSED_DIR, load_bench9, returns_and_variance


def leverage_profile(returns, variance, max_lag=40):
    """Corr(return[t], variance[t + lag]), including contemporaneous lag 0."""
    x = np.asarray(returns, dtype=float)
    v = np.asarray(variance, dtype=float)
    if len(x) != len(v) or len(x) <= max_lag + 1:
        raise ValueError("Aligned series must be longer than max_lag + 1")
    return [float(np.corrcoef(x if lag == 0 else x[:-lag],
                              v if lag == 0 else v[lag:])[0, 1])
            for lag in range(max_lag + 1)]


def compute_targets(frames=None, max_lag=40):
    """Build the JSON-ready leverage target for every benchmark asset."""
    frames = load_bench9() if frames is None else frames
    return {name: leverage_profile(*returns_and_variance(frames[name]), max_lag)
            for name in ASSET_NAMES}


if __name__ == "__main__":
    destination = PROCESSED_DIR / "bench9_real_leverage_wide.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(compute_targets(), indent=2) + "\n")
    print(destination)
