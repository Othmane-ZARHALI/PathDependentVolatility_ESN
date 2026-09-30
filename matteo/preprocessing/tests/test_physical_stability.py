"""Checks for the reusable OHLC stability estimators and block boundaries."""

import numpy as np

from matteo.preprocessing.hurst import estimate_hurst
from matteo.preprocessing.kurtosis import corrected_excess_kurtosis
from matteo.preprocessing.leverage import leverage_profile
from matteo.preprocessing.physical_stability import (
    _acf, _corr, hurst, sample_blocks, statistics, zumbach,
)
from matteo.preprocessing.taylor_acf import autocorrelation_profile
from matteo.preprocessing.zumbach import zumbach_components


def test_statistics_match_existing_extractors():
    rng = np.random.default_rng(421)
    returns = rng.normal(size=600) * 0.01
    variance = np.exp(rng.normal(size=600) - 10)
    values = statistics(returns, variance)
    forward, backward = zumbach_components(returns, variance, (5, 10, 20))
    abs_acf = np.array(autocorrelation_profile(abs(returns), 20))
    sq_acf = np.array(autocorrelation_profile(returns**2, 20))
    assert np.isclose(values[0], estimate_hurst(0.5 * np.log(variance)))
    assert np.isclose(values[1], corrected_excess_kurtosis(returns))
    assert np.isclose(values[2], np.mean(leverage_profile(returns, variance, 5)[1:]))
    assert np.allclose(values[3:6], np.array(forward) - backward)
    assert np.isclose(values[6], np.mean(abs_acf - sq_acf))
    assert np.isclose(values[4], zumbach(returns, variance, 10))


def test_block_labels_exclude_artificial_level_jump():
    rng = np.random.default_rng(229)
    values = np.r_[rng.normal(size=126), rng.normal(size=126) + 100]
    segments = np.repeat((0, 1), 126)
    assert _acf(values, 1) > 0.9
    expected = _corr(np.r_[values[:125], values[126:251]],
                     np.r_[values[1:126], values[127:252]])
    assert np.isclose(_acf(values, 1, segments), expected)

    indices, labels = sample_blocks(500, 126, rng)
    assert len(indices) == len(labels) == 500
    assert np.all(np.diff(indices)[np.diff(labels) == 0] == 1)


def test_hurst_uses_only_within_block_pairs():
    rng = np.random.default_rng(17)
    log_vol = rng.normal(size=252)
    labels = np.repeat((0, 1), 126)
    shifted = log_vol.copy()
    shifted[126:] += 100
    assert np.isclose(hurst(log_vol, labels), hurst(shifted, labels))
    assert not np.isclose(hurst(log_vol), hurst(shifted))
