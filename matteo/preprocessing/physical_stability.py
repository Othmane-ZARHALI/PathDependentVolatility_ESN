"""Physical-statistic uncertainty for Othmane's nine daily OHLC series.

Run from the repository root with ``python -m matteo.preprocessing.physical_stability``.
The bootstrap samples contiguous daily tuples, then excludes every lagged
comparison that crosses a sampled-block join. The primary variance measure is
Othmane's daily Garman--Klass proxy; an overnight-gap sensitivity is separate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from .hurst import MOMENT_ORDERS, REAL_LAGS
from .market_data import (
    ASSET_NAMES, INDEX_FILE, INDEX_TICKERS, ROOT, STOCK_FILE, STOCK_NAMES,
    garman_klass_variance,
    load_indices, load_stocks,
)


OUT_DIR = ROOT / "matteo" / "exogenous_reservoir_eq_m2" / "experiments" / "physical_stability_001"
METRICS = (
    "hurst", "kurtosis", "leverage_1_5", "zumbach_5", "zumbach_10",
    "zumbach_20", "taylor_gap_1_20", "clustering_1_20", "return_acf_1",
)


def observations(frame):
    """Return aligned close-to-close returns and same-day OHLC log ratios."""
    open_, high, low, close = (frame[c].to_numpy(float) for c in ("O", "H", "L", "C"))
    returns = np.diff(np.log(close))
    high_low = np.log(high[1:] / low[1:])
    close_open = np.log(close[1:] / open_[1:])
    overnight = np.log(open_[1:] / close[:-1])
    assert np.allclose(returns, overnight + close_open, rtol=0, atol=1e-12)
    return returns, high_low, close_open, overnight


def variance(high_low, close_open, overnight=None):
    raw = 0.5 * high_low**2 - (2 * np.log(2) - 1) * close_open**2
    if overnight is not None:
        raw = raw + overnight**2
    return np.maximum(raw, 1e-12)


def parkinson_variance(high_low):
    """High/low-only sensitivity, independent of the recorded open."""
    return np.maximum(high_low**2 / (4 * np.log(2)), 1e-12)


def _mask(segments, lag, n):
    return slice(None) if segments is None else segments[:-lag] == segments[lag:]


def _corr(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    xx = x - x.mean()
    yy = y - y.mean()
    scale = np.sqrt(np.dot(xx, xx) * np.dot(yy, yy))
    return float(np.dot(xx, yy) / scale) if scale > 0 else np.nan


def _acf(x, lag, segments=None):
    valid = _mask(segments, lag, len(x))
    return _corr(x[:-lag][valid], x[lag:][valid])


def hurst(log_vol, segments=None):
    """Othmane's four-order structure-function slope, lags 2--40."""
    log_lags = np.log(REAL_LAGS)
    slopes = []
    for order in MOMENT_ORDERS:
        moments = []
        for lag in REAL_LAGS:
            valid = _mask(segments, lag, len(log_vol))
            diff = np.abs((log_vol[lag:] - log_vol[:-lag])[valid])
            moments.append(np.mean(diff**order))
        slopes.append(np.polyfit(log_lags, np.log(moments), 1)[0])
    return float(np.polyfit(MOMENT_ORDERS, slopes, 1)[0])


def zumbach(returns, var, horizon, segments=None):
    """Othmane-style squared cumulative return versus mean daily GK variance."""
    n = len(returns)
    t = np.arange(horizon, n - horizon)
    if segments is not None:
        good = segments[t - horizon] == segments[t + horizon - 1]
        t = t[good]
    if len(t) < 30:
        return np.nan
    cr = np.r_[0.0, np.cumsum(returns)]
    cv = np.r_[0.0, np.cumsum(var)]
    past_r = cr[t] - cr[t - horizon]
    future_r = cr[t + horizon] - cr[t]
    past_v = (cv[t] - cv[t - horizon]) / horizon
    future_v = (cv[t + horizon] - cv[t]) / horizon
    return _corr(past_r**2, future_v) - _corr(past_v, future_r**2)


def statistics(returns, var, segments=None):
    n = len(returns)
    log_vol = 0.5 * np.log(var)
    centered = returns - returns.mean()
    m2 = np.mean(centered**2)
    raw_kurt = np.mean(centered**4) / m2**2 - 3
    kurt = (n - 1) * ((n + 1) * raw_kurt + 6) / ((n - 2) * (n - 3))
    leverage = np.mean([
        _corr(returns[:-lag][_mask(segments, lag, n)],
              var[lag:][_mask(segments, lag, n)])
        for lag in range(1, 6)
    ])
    gap = np.mean([
        _acf(np.abs(returns), lag, segments) - _acf(returns**2, lag, segments)
        for lag in range(1, 21)
    ])
    clustering = np.mean([_acf(log_vol, lag, segments) for lag in range(1, 21)])
    return np.array((
        hurst(log_vol, segments), float(kurt), leverage,
        *(zumbach(returns, var, lag, segments) for lag in (5, 10, 20)),
        gap, clustering, _acf(returns, 1, segments),
    ), dtype=float)


def sample_blocks(n, length, rng):
    """Fixed-length moving blocks; labels expose artificial joins."""
    count = int(np.ceil(n / length))
    starts = rng.integers(0, n - length + 1, size=count)
    indices = (starts[:, None] + np.arange(length)).ravel()[:n]
    segments = np.repeat(np.arange(count), length)[:n]
    return indices, segments


def bootstrap(data, length, reps, rng):
    returns, high_low, close_open, _ = data
    n = len(returns)
    draws = np.empty((reps, len(METRICS)))
    for b in range(reps):
        indices, segments = sample_blocks(n, length, rng)
        var = variance(high_low[indices], close_open[indices])
        draws[b] = statistics(returns[indices], var, segments)
    if not np.isfinite(draws).all():
        raise ValueError("A bootstrap statistic was nonfinite")
    return draws


def heterogeneity(points, draws):
    """Centered null bootstrap of inverse-variance between-period dispersion."""
    se = draws.std(axis=1, ddof=1)
    weights = 1 / se**2
    common = np.sum(weights * points, axis=0) / weights.sum(axis=0)
    observed = np.sum(weights * (points - common)**2, axis=0)
    centered = draws - draws.mean(axis=1, keepdims=True)
    null_points = common[None, None, :] + centered
    null_common = np.sum(weights[:, None, :] * null_points, axis=0) / weights.sum(axis=0)
    null_q = np.sum(weights[:, None, :] * (null_points - null_common[None])**2, axis=0)
    p_value = (1 + np.sum(null_q >= observed, axis=0)) / (draws.shape[1] + 1)
    return common, observed, p_value


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def analyze_rolling(asset, frame, args, output):
    """Rolling estimates and pointwise uncertainty, never independent tests."""
    rows = []
    proxy_rows = []
    quality_rows = []
    for window_years in args.rolling_window_years:
        final_start = args.last_year - window_years + 1
        if final_start < args.first_year:
            raise ValueError(f"{asset}: {window_years}-year rolling window exceeds the date range")
        for start in range(args.first_year, final_start + 1, args.rolling_step_years):
            end = start + window_years - 1
            selected = frame[(frame.index.year >= start) & (frame.index.year <= end)]
            returns, high_low, close_open, overnight = observations(selected)
            if len(returns) < max(args.min_observations, args.rolling_block_length):
                raise ValueError(f"{asset} rolling {start}-{end}: only {len(returns)} returns")
            gk = variance(high_low, close_open)
            with_overnight = variance(high_low, close_open, overnight)
            parkinson = parkinson_variance(high_low)
            point = statistics(returns, gk)
            overnight_point = statistics(returns, with_overnight)
            parkinson_point = statistics(returns, parkinson)
            # The stream for each window is independent of disjoint-block settings.
            seed = np.random.SeedSequence((args.seed, ASSET_NAMES.index(asset),
                                           1, window_years, start))
            draws = bootstrap((returns, high_low, close_open, overnight),
                              args.rolling_block_length, args.rolling_replicates,
                              np.random.default_rng(seed))
            common = {
                "asset": asset, "window_years": window_years,
                "window_start_year": start, "window_end_year": end,
                "first_date": str(selected.index.min().date()),
                "last_date": str(selected.index.max().date()),
                "n_returns": len(returns),
            }
            quality_rows.append({
                **common,
                "open_equal_previous_close_share": float(np.mean(np.abs(overnight) < 1e-10)),
                "overnight_variance_share": float(np.sum(overnight**2) / np.sum(with_overnight)),
            })
            for j, metric in enumerate(METRICS):
                low, high = np.quantile(draws[:, j], (0.025, 0.975))
                rows.append({
                    **common, "metric": metric, "estimate": point[j],
                    "bootstrap_mean": draws[:, j].mean(),
                    "bootstrap_se": draws[:, j].std(ddof=1),
                    "ci_low": low, "ci_high": high,
                    "block_length": args.rolling_block_length,
                    "replicates": args.rolling_replicates,
                })
                proxy_rows.append({
                    **common, "metric": metric, "gk": point[j],
                    "gk_plus_overnight_sq": overnight_point[j],
                    "parkinson": parkinson_point[j],
                })
        print(asset, f"{window_years}-year rolling windows done", flush=True)
    write_csv(output / "rolling_estimates.csv", rows,
              ("asset", "window_years", "window_start_year", "window_end_year",
               "first_date", "last_date", "n_returns", "metric", "estimate",
               "bootstrap_mean", "bootstrap_se", "ci_low", "ci_high",
               "block_length", "replicates"))
    write_csv(output / "rolling_proxy_sensitivity.csv", proxy_rows,
              ("asset", "window_years", "window_start_year", "window_end_year",
               "first_date", "last_date", "n_returns", "metric", "gk",
               "gk_plus_overnight_sq", "parkinson"))
    write_csv(output / "rolling_source_checks.csv", quality_rows,
              ("asset", "window_years", "window_start_year", "window_end_year",
               "first_date", "last_date", "n_returns",
               "open_equal_previous_close_share", "overnight_variance_share"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assets", nargs="+", default=["SP500"],
                        choices=(*ASSET_NAMES, "all"))
    parser.add_argument("--first-year", type=int, default=1994)
    parser.add_argument("--last-year", type=int, default=2023)
    parser.add_argument("--period-years", type=int, default=5)
    parser.add_argument("--include-partial", action="store_true",
                        help="Include a final period shorter than --period-years")
    parser.add_argument("--min-observations", type=int, default=500)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--block-lengths", type=int, nargs="+", default=(126, 63))
    parser.add_argument("--seed", type=int, default=20260930)
    parser.add_argument("--rolling", action="store_true",
                        help="Also compute rolling windows with pointwise bootstrap intervals")
    parser.add_argument("--rolling-window-years", type=int, nargs="+", default=(5, 10))
    parser.add_argument("--rolling-step-years", type=int, default=1)
    parser.add_argument("--rolling-replicates", type=int, default=500)
    parser.add_argument("--rolling-block-length", type=int, default=126)
    parser.add_argument("--output", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    if args.replicates < 100 or any(x < 60 for x in args.block_lengths):
        parser.error("Use at least 100 replicates and block lengths of at least 60 days")
    if len(set(args.block_lengths)) != len(args.block_lengths):
        parser.error("Block lengths must be distinct")
    if args.period_years < 1 or args.last_year < args.first_year:
        parser.error("Invalid year range or period length")
    if args.rolling and (args.rolling_replicates < 100 or args.rolling_block_length < 60
                         or args.rolling_step_years < 1
                         or any(x < 1 for x in args.rolling_window_years)
                         or len(set(args.rolling_window_years)) != len(args.rolling_window_years)):
        parser.error("Invalid rolling windows, step, bootstrap length or replicate count")
    if "all" in args.assets and len(args.assets) != 1:
        parser.error("Use --assets all alone")

    assets = ASSET_NAMES if args.assets == ["all"] else tuple(dict.fromkeys(args.assets))
    frames = {}
    if any(asset in STOCK_NAMES for asset in assets):
        frames.update(load_stocks())
    if any(asset in INDEX_TICKERS for asset in assets):
        frames.update(load_indices())
    for asset in assets:
        analyze_asset(asset, frames[asset], args)


def analyze_asset(asset, frame, args):
    output = args.output / asset
    output.mkdir(parents=True, exist_ok=True)
    # Asset-specific streams give identical results in single- and multi-asset runs.
    rng = np.random.default_rng(np.random.SeedSequence((args.seed, ASSET_NAMES.index(asset))))

    periods = []
    estimate_rows = []
    proxy_rows = []
    quality_rows = []
    yearly_rows = []
    covariance_rows = []
    draws_by_length = {length: [] for length in args.block_lengths}
    points = []
    _, full_hl, full_co, full_overnight = observations(frame)
    full_gk = variance(full_hl, full_co)
    for year in range(args.first_year, args.last_year + 1):
        use = frame.index[1:].year == year
        if not use.any():
            continue
        yearly_rows.append({
            "asset": asset, "year": year, "n_returns": int(use.sum()),
            "open_equal_previous_close_share": float(np.mean(np.abs(full_overnight[use]) < 1e-10)),
            "overnight_variance_share": float(np.sum(full_overnight[use]**2) /
                                              np.sum(full_gk[use] + full_overnight[use]**2)),
        })
    for start in range(args.first_year, args.last_year + 1, args.period_years):
        end = min(start + args.period_years - 1, args.last_year)
        if end - start + 1 < args.period_years and not args.include_partial:
            continue
        selected = frame[(frame.index.year >= start) & (frame.index.year <= end)]
        if len(selected) - 1 < args.min_observations:
            raise ValueError(f"{asset} {start}-{end}: only {len(selected)-1} returns")
        if len(selected) - 1 < max(args.block_lengths):
            raise ValueError(f"{asset} {start}-{end}: bootstrap block exceeds the period")
        data = observations(selected)
        returns, high_low, close_open, overnight = data
        primary = variance(high_low, close_open)
        adjusted = variance(high_low, close_open, overnight)
        parkinson = parkinson_variance(high_low)
        raw_gk = 0.5 * high_low**2 - (2 * np.log(2) - 1) * close_open**2
        # Check our direct GK construction against the established code.
        assert np.allclose(primary, garman_klass_variance(selected)[1:])
        point = statistics(returns, primary)
        points.append(point)
        period = f"{start}-{end}"
        periods.append(period)
        quality_rows.append({
            "asset": asset, "period": period, "first_date": str(selected.index.min().date()),
            "last_date": str(selected.index.max().date()), "n_returns": len(returns),
            "gk_floored_days": int(np.sum(raw_gk <= 1e-12)),
            "open_equal_previous_close_share": float(np.mean(np.abs(overnight) < 1e-10)),
            "overnight_variance_share": float(np.sum(overnight**2) / np.sum(adjusted)),
        })
        adjusted_point = statistics(returns, adjusted)
        parkinson_point = statistics(returns, parkinson)
        for j, metric in enumerate(METRICS):
            proxy_rows.append({"asset": asset, "period": period, "metric": metric,
                               "gk": point[j], "gk_plus_overnight_sq": adjusted_point[j],
                               "parkinson": parkinson_point[j],
                               "overnight_difference": adjusted_point[j] - point[j],
                               "parkinson_difference": parkinson_point[j] - point[j]})
        for length in args.block_lengths:
            draws = bootstrap(data, length, args.replicates, rng)
            draws_by_length[length].append(draws)
            cov = np.cov(draws, rowvar=False)
            for i, metric_i in enumerate(METRICS):
                for j, metric_j in enumerate(METRICS):
                    covariance_rows.append({
                        "asset": asset, "period": period, "block_length": length,
                        "metric_i": metric_i, "metric_j": metric_j, "covariance": cov[i, j],
                    })
            for j, metric in enumerate(METRICS):
                low, high = np.quantile(draws[:, j], (0.025, 0.975))
                estimate_rows.append({
                    "asset": asset, "period": period, "n_returns": len(returns),
                    "metric": metric, "block_length": length, "replicates": args.replicates,
                    "estimate": point[j], "bootstrap_mean": draws[:, j].mean(),
                    "bootstrap_se": draws[:, j].std(ddof=1), "ci_low": low, "ci_high": high,
                })
        print(asset, period, "n=", len(returns), "done", flush=True)

    hetero_rows = []
    if not points:
        raise ValueError(f"{asset}: no eligible periods")
    points = np.stack(points)
    if len(points) > 1:
        for length, blocks in draws_by_length.items():
            draws = np.stack(blocks)
            common, q, p = heterogeneity(points, draws)
            for j, metric in enumerate(METRICS):
                hetero_rows.append({
                    "asset": asset, "metric": metric, "block_length": length,
                    "common_inverse_variance": common[j], "observed_range": np.ptp(points[:, j]),
                    "q_statistic": q[j], "p_common_statistic": p[j],
                })

    write_csv(output / "period_estimates.csv", estimate_rows,
              ("asset", "period", "n_returns", "metric", "block_length", "replicates",
               "estimate", "bootstrap_mean", "bootstrap_se", "ci_low", "ci_high"))
    write_csv(output / "heterogeneity.csv", hetero_rows,
              ("asset", "metric", "block_length", "common_inverse_variance",
               "observed_range", "q_statistic", "p_common_statistic"))
    write_csv(output / "proxy_sensitivity.csv", proxy_rows,
              ("asset", "period", "metric", "gk", "gk_plus_overnight_sq",
               "parkinson", "overnight_difference", "parkinson_difference"))
    write_csv(output / "source_checks.csv", quality_rows,
              ("asset", "period", "first_date", "last_date", "n_returns",
               "gk_floored_days", "open_equal_previous_close_share",
               "overnight_variance_share"))
    write_csv(output / "source_years.csv", yearly_rows,
              ("asset", "year", "n_returns", "open_equal_previous_close_share",
               "overnight_variance_share"))
    write_csv(output / "bootstrap_covariance.csv", covariance_rows,
              ("asset", "period", "block_length", "metric_i", "metric_j", "covariance"))
    if args.rolling:
        analyze_rolling(asset, frame, args, output)
    source = STOCK_FILE if asset in STOCK_NAMES else INDEX_FILE
    (output / "run.json").write_text(json.dumps({
        "asset": asset, "periods": periods, "seed": args.seed,
        "first_year": args.first_year, "last_year": args.last_year,
        "period_years": args.period_years, "include_partial": args.include_partial,
        "min_observations": args.min_observations,
        "bootstrap_replicates": args.replicates, "block_lengths": args.block_lengths,
        "rolling": args.rolling,
        "rolling_window_years": args.rolling_window_years if args.rolling else [],
        "rolling_step_years": args.rolling_step_years if args.rolling else None,
        "rolling_replicates": args.rolling_replicates if args.rolling else None,
        "rolling_block_length": args.rolling_block_length if args.rolling else None,
        "primary_proxy": "same-day Garman-Klass variance, floor 1e-12",
        "sensitivity_proxy": "Garman-Klass plus squared overnight open/previous-close log gap",
        "second_sensitivity_proxy": "Parkinson high-low variance, floor 1e-12",
        "bootstrap": "fixed moving blocks of aligned daily tuples; exclude cross-block lag comparisons",
        "interval": "95% percentile bootstrap; point and bootstrap mean both reported",
        "heterogeneity": "inverse-variance Q with centered within-period bootstrap null",
        "source_file": str(source), "source_sha256": file_hash(source),
        "source_first": str(frame.index.min().date()), "source_last": str(frame.index.max().date()),
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
