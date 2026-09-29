"""Reconstruct the two real Zumbach correlation curves at even lags 2-40."""

import json

import numpy as np

from .market_data import ASSET_NAMES, PROCESSED_DIR, load_bench9, returns_and_variance


LAGS = tuple(range(2, 41, 2))


def zumbach_components(returns, variance, lags=LAGS):
    """Corr(past return squared, future variance) and its reverse-time peer."""
    x = np.asarray(returns, dtype=float)
    v = np.asarray(variance, dtype=float)
    if len(x) != len(v) or len(x) <= 2 * max(lags) + 30:
        raise ValueError("Aligned series is too short for Zumbach windows")
    cumulative_x = np.concatenate(([0.0], np.cumsum(x)))
    cumulative_v = np.concatenate(([0.0], np.cumsum(v)))
    forward_variance, forward_return = [], []
    for lag in lags:
        starts = lag + np.arange(len(x) - 2 * lag)
        past_return = cumulative_x[starts] - cumulative_x[starts - lag]
        future_return = cumulative_x[starts + lag] - cumulative_x[starts]
        past_variance = (cumulative_v[starts] - cumulative_v[starts - lag]) / lag
        future_variance = (cumulative_v[starts + lag] - cumulative_v[starts]) / lag
        forward_variance.append(float(np.corrcoef(past_return ** 2, future_variance)[0, 1]))
        forward_return.append(float(np.corrcoef(past_variance, future_return ** 2)[0, 1]))
    return forward_variance, forward_return


def compute_targets(frames=None, lags=LAGS):
    """Return the JSON structure consumed by joint calibration."""
    frames = load_bench9() if frames is None else frames
    output = {"lags": list(lags), "pR2_fV": {}, "pV_fR2": {}}
    for name in ASSET_NAMES:
        returns, variance = returns_and_variance(frames[name])
        output["pR2_fV"][name], output["pV_fR2"][name] = zumbach_components(
            returns, variance, lags)
    return output


if __name__ == "__main__":
    destination = PROCESSED_DIR / "bench9_real_zumbach_profile.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(compute_targets(), indent=2) + "\n")
    print(destination)
