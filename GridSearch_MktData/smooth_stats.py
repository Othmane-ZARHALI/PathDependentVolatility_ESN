"""
smooth_stats.py (reconstructed)
======================================
Implements compute_statistics_smooth exactly per the formulas documented
in esn_market_protocol.tex Section 1.4 ("Smooth surrogate statistics"):

  - H_hat: log-spaced (15 points, log-uniform in [1, n/4]) structure-
    function regression, replacing the hard 7-fixed-dyadic-lag version.
  - q995_vol_ann: Harrell-Davis estimator (Beta-weighted average of
    every order statistic), replacing the hard single order statistic.
  - max_vol_ann, max_ret_acf: log-sum-exp soft-max (beta = 8/std of the
    values maximised over), replacing the hard max.
  - taylor_frac: sigmoid relaxation (kappa=200) of the indicator
    1[abs_acf(L) > sq_acf(L)].
  - kurtosis: tanh-winsorised 4th moment (clip scale cs = 12*sqrt(Var(x))).

mean_vol_ann, mean_vol_acf, taylor_gap, zumbach, and leverage are
already smooth (means/mean-correlations) and are computed identically
to the hard version.
"""
import math
import numpy as np
from scipy.special import betainc

TRADING_DAYS = 252.0


def _corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() < 1e-14 or b.std() < 1e-14:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def _harrell_davis(sorted_vals, q):
    n = len(sorted_vals)
    a = (n + 1) * q
    b = (n + 1) * (1 - q)
    edges = np.arange(0, n + 1) / n
    cdf = betainc(a, b, edges)
    w = cdf[1:] - cdf[:-1]
    return float(np.sum(w * sorted_vals))


def _softmax(vals, beta):
    vals = np.asarray(vals, float)
    m = np.max(vals)
    return float(m + (1.0 / beta) * math.log(np.mean(np.exp(beta * (vals - m)))))


def compute_statistics_smooth(daily_x, daily_var):
    x = np.asarray(daily_x, float)
    v = np.asarray(daily_var, float)
    sd = np.sqrt(np.maximum(v, 1e-30))
    lv = np.log(v + 1e-30)
    n = len(x)
    lg = np.arange(1, 21)

    ret_acf = np.array([_corr(x[:-L], x[L:]) for L in lg])
    abs_acf = np.array([_corr(np.abs(x[:-L]), np.abs(x[L:])) for L in lg])
    sq_acf = np.array([_corr(x[:-L] ** 2, x[L:] ** 2) for L in lg])
    vol_acf = np.array([_corr(v[:-L], v[L:]) for L in lg])
    lev = np.array([_corr(x[:-L], v[L:]) for L in lg])

    ann = math.sqrt(TRADING_DAYS)

    # ---- H_hat: 15 log-spaced lags in [1, n/4] ----
    lag_max = max(2, n // 4)
    lags = np.unique(np.round(np.geomspace(1, lag_max, 15)).astype(int))
    lags = lags[lags >= 1]
    xs, ys = [], []
    for L in lags:
        if L < n:
            d = lv[L:] - lv[:-L]
            vv = float(np.mean(d * d))
            if vv > 1e-30 and np.isfinite(vv):
                xs.append(math.log(L))
                ys.append(math.log(vv))
    H_smooth = 0.5 * float(np.polyfit(xs, ys, 1)[0]) if len(xs) >= 3 else np.nan

    # ---- q995_vol_ann: Harrell-Davis ----
    sd_sorted = np.sort(sd)
    q995_smooth = ann * _harrell_davis(sd_sorted, 0.995)

    # ---- max_vol_ann, max_ret_acf: log-sum-exp soft-max ----
    beta_vol = 8.0 / max(float(np.std(sd)), 1e-12)
    max_vol_smooth = ann * _softmax(sd, beta_vol)

    abs_ret_acf = np.abs(ret_acf)
    beta_racf = 8.0 / max(float(np.nanstd(abs_ret_acf)), 1e-12)
    max_ret_acf_smooth = _softmax(abs_ret_acf[np.isfinite(abs_ret_acf)], beta_racf)

    # ---- taylor_frac: sigmoid relaxation ----
    kappa = 200.0
    diffs = abs_acf - sq_acf
    taylor_frac_smooth = float(np.nanmean(1.0 / (1.0 + np.exp(-kappa * diffs))))

    # ---- kurtosis: tanh-winsorised 4th moment ----
    c = x - x.mean()
    var = float(np.var(c))
    cs = 12.0 * math.sqrt(max(var, 1e-30))
    c_soft = cs * np.tanh(c / cs)
    kurt_smooth = float(np.mean(c_soft ** 4) / (var ** 2 + 1e-30) - 3.0)

    # ---- already-smooth terms (identical to hard version) ----
    zvals = []
    for L in [5, 10, 20]:
        m = n - 2 * L
        if m > 30:
            csx = np.concatenate([[0.0], np.cumsum(x)])
            csv = np.concatenate([[0.0], np.cumsum(v)])
            idx = L + np.arange(m)
            pR = csx[idx] - csx[idx - L]
            fR = csx[idx + L] - csx[idx]
            pV = (csv[idx] - csv[idx - L]) / L
            fV = (csv[idx + L] - csv[idx]) / L
            zvals.append(_corr(pR ** 2, fV) - _corr(pV, fR ** 2))

    return dict(
        H_hat=float(H_smooth),
        mean_vol_ann=float(sd.mean()) * ann,
        q995_vol_ann=float(q995_smooth),
        max_vol_ann=float(max_vol_smooth),
        mean_vol_acf=float(np.nanmean(vol_acf)),
        taylor_gap=float(np.nanmean(abs_acf) - np.nanmean(sq_acf)),
        taylor_frac=float(taylor_frac_smooth),
        zumbach=float(np.nanmean(zvals)),
        leverage=float(np.nanmean(lev)),
        kurtosis=float(kurt_smooth),
        max_ret_acf=float(max_ret_acf_smooth),
    )
