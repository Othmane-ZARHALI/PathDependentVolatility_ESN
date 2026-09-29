"""Reconstruct finite-sample corrected excess kurtosis of daily returns."""

import json

import numpy as np

from .market_data import ASSET_NAMES, PROCESSED_DIR, load_bench9, returns_and_variance


def corrected_excess_kurtosis(returns):
    """Adjusted Fisher-Pearson kurtosis (SciPy fisher=True, bias=False)."""
    x = np.asarray(returns, dtype=float)
    n = len(x)
    if n < 4 or not np.isfinite(x).all():
        raise ValueError("Need at least four finite returns")
    centered = x - x.mean()
    second = np.mean(centered ** 2)
    if second <= 0:
        raise ValueError("Return variance must be positive")
    raw_excess = np.mean(centered ** 4) / second ** 2 - 3.0
    return float((n - 1) * ((n + 1) * raw_excess + 6.0) / ((n - 2) * (n - 3)))


def compute_targets(frames=None):
    """Compute one real corrected excess-kurtosis target per asset."""
    frames = load_bench9() if frames is None else frames
    return {name: corrected_excess_kurtosis(returns_and_variance(frames[name])[0])
            for name in ASSET_NAMES}


if __name__ == "__main__":
    destination = PROCESSED_DIR / "bench9_real_kurtosis.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(compute_targets(), indent=2) + "\n")
    print(destination)
