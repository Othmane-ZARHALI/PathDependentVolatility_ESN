"""
market_vix.py
======================================
Model-free VIX / VVIX replication from raw SPX / VIX option-chain
snapshots (CBOE VIX white-paper methodology, Black-76 undiscounted
prices built from mid-IV quotes since the input files quote implied
vols, not dollar premia).

Inputs
------
SPX_CSV, VIX_CSV : option-chain CSVs with columns
    reference, expiry, callput (+1/-1), strike, bid, ask (mid-IV,
    decimal), spot, forward, moneyness, expiry_yf
Both files are single-snapshot (one `reference` date), multi-expiry.

Public functions
-----------------
compute_cboe_vix(csv_path, ref, near_expiry, next_expiry)
    -> dict with sigma1^2, sigma2^2, VIX_30d (CBOE 2-step interpolation)

excel_to_date(serial)
    -> datetime, Excel-serial-date helper (1899-12-30 epoch)
"""
import numpy as np
import pandas as pd
from scipy.stats import norm
from datetime import datetime, timedelta


def excel_to_date(serial):
    return datetime(1899, 12, 30) + timedelta(days=int(serial))


def build_moneyness_grid(lo=0.60, hi=1.60, atm_lo=0.85, atm_hi=1.20,
                          atm_step=0.005, wing_step=0.02):
    """ATM-weighted moneyness grid: dense (atm_step) inside
    [atm_lo,atm_hi], coarse (wing_step) in the wings. Weighting the
    inversion toward the ATM region matters for 2 reasons: (i) it's
    where real strikes cluster and where the market itself is most
    liquid/most economically meaningful, and (ii) the CBOE discretized
    sum (cboe_variance_from_grid) uses local strike spacing DeltaK_i as
    its own quadrature weight, so a denser ATM grid directly sharpens
    the variance integral exactly where curvature (and MC precision)
    matters most, without paying for that density in the flat wings."""
    left = np.arange(lo, atm_lo, wing_step)
    mid = np.arange(atm_lo, atm_hi + 1e-9, atm_step)
    right = np.arange(atm_hi + wing_step, hi + 1e-9, wing_step)
    grid = np.unique(np.round(np.concatenate([left, mid, right]), 6))
    return grid


def black76_undisc(F, K, T, sigma, is_call):
    """Undiscounted (forward-measure) Black-76 price, i.e. e^{RT}*Q(K)."""
    if sigma <= 0 or T <= 0:
        return max(0.0, (F - K) if is_call else (K - F))
    d1 = (np.log(F / K) + 0.5 * sigma**2 * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    if is_call:
        return F * norm.cdf(d1) - K * norm.cdf(d2)
    return K * norm.cdf(-d2) - F * norm.cdf(-d1)


def cboe_variance_from_grid(F, T, strikes, Q_forward):
    """CORE CBOE discretized-strike replication (Eq. in the turn
    derivation), shared by BOTH the market path (Q_forward built from
    mid-IV via Black-76) and the model path (Q_forward built directly
    from Monte-Carlo undiscounted option prices) -- same summation
    logic, same strike-grid truncation rule, applied to whichever
    quote source is passed in. `strikes` must be sorted ascending;
    `Q_forward[i]` is the forward-value (undiscounted, i.e. already
    == e^{RT}Q(K)) OTM price at strikes[i], or NaN if unavailable."""
    strikes = np.asarray(strikes, float)
    Q_forward = np.asarray(Q_forward, float)
    below = strikes[strikes < F]
    K0 = below.max() if len(below) else strikes.min()

    idx0 = int(np.where(strikes == K0)[0][0])
    keep = {idx0}
    for i in range(idx0 - 1, -1, -1):
        if np.isnan(Q_forward[i]):
            break
        keep.add(i)
    for i in range(idx0 + 1, len(strikes)):
        if np.isnan(Q_forward[i]):
            break
        keep.add(i)
    keep = sorted(keep)
    Ks = strikes[keep]
    Qs = Q_forward[keep]
    n = len(Ks)
    dK = np.zeros(n)
    dK[0] = Ks[1] - Ks[0]
    dK[-1] = Ks[-1] - Ks[-2]
    for i in range(1, n - 1):
        dK[i] = (Ks[i + 1] - Ks[i - 1]) / 2

    term = np.sum(dK / Ks**2 * Qs)
    sigma2 = (2 / T) * term - (1 / T) * (F / K0 - 1) ** 2
    return dict(sigma2=sigma2, K0=float(K0), n_strikes=n, kmin=float(Ks.min()), kmax=float(Ks.max()))


def _cboe_sigma2_fromIV(df_e, F, T):
    """One expiry's CBOE-style sigma^2, built from mid-IV quotes via
    Black-76 (e^{RT}Q(K) = Black76_undiscounted(K) exactly under this
    convention, so R never needs to be estimated separately), then
    handed to the shared cboe_variance_from_grid() core."""
    calls = df_e[df_e.callput == 1].groupby('strike').agg(
        cbid=('bid', 'first'), cask=('ask', 'first'))
    puts = df_e[df_e.callput == -1].groupby('strike').agg(
        pbid=('bid', 'first'), pask=('ask', 'first'))
    strikes = sorted(set(calls.index) | set(puts.index))
    below = [k for k in strikes if k < F]
    K0 = max(below) if below else min(strikes)

    Qs = []
    for K in strikes:
        c_iv = p_iv = np.nan
        if K in calls.index and not (pd.isna(calls.loc[K, 'cbid']) or pd.isna(calls.loc[K, 'cask'])):
            c_iv = (calls.loc[K, 'cbid'] + calls.loc[K, 'cask']) / 2
        if K in puts.index and not (pd.isna(puts.loc[K, 'pbid']) or pd.isna(puts.loc[K, 'pask'])):
            p_iv = (puts.loc[K, 'pbid'] + puts.loc[K, 'pask']) / 2
        if K < K0:
            Q = black76_undisc(F, K, T, p_iv, False) if not pd.isna(p_iv) else np.nan
        elif K > K0:
            Q = black76_undisc(F, K, T, c_iv, True) if not pd.isna(c_iv) else np.nan
        else:
            vals = []
            if not pd.isna(c_iv):
                vals.append(black76_undisc(F, K, T, c_iv, True))
            if not pd.isna(p_iv):
                vals.append(black76_undisc(F, K, T, p_iv, False))
            Q = np.mean(vals) if vals else np.nan
        Qs.append(Q)

    res = cboe_variance_from_grid(F, T, strikes, Qs)
    return res


def compute_cboe_vix(csv_path, near_expiry, next_expiry):
    """Full 2-step CBOE VIX (or VVIX, if fed a VIX-option chain)
    replication: near-term + next-term sigma^2, interpolated to a
    constant 30-day maturity."""
    df = pd.read_csv(csv_path)
    ref = int(df.reference.iloc[0])
    out = {}
    for e in (near_expiry, next_expiry):
        df_e = df[df.expiry == e]
        F = float(df_e.forward.iloc[0])
        T = (e - ref) / 365
        res = _cboe_sigma2_fromIV(df_e, F, T)
        res.update(T=T, F=F, expiry=int(e), expiry_date=str(excel_to_date(e).date()))
        out[str(e)] = res

    T1, T2 = out[str(near_expiry)]['T'], out[str(next_expiry)]['T']
    s1, s2 = out[str(near_expiry)]['sigma2'], out[str(next_expiry)]['sigma2']
    N_T1, N_T2, N_30, N_365 = T1 * 365, T2 * 365, 30, 365
    vix2 = (T1 * s1 * (N_T2 - N_30) / (N_T2 - N_T1)
            + T2 * s2 * (N_30 - N_T1) / (N_T2 - N_T1)) * (N_365 / N_30)
    return dict(reference_date=str(excel_to_date(ref).date()), near=out[str(near_expiry)],
                next=out[str(next_expiry)], index_30d=100 * np.sqrt(vix2))


if __name__ == "__main__":
    import json
    spx = compute_cboe_vix("/mnt/user-data/uploads/1788726546415_SPX_data_19_02_2026.csv", 46101, 46129)
    vix = compute_cboe_vix("/mnt/user-data/uploads/1788726546415_VIX_data_19_02_2026.csv", 46099, 46127)
    print("VIX  (from SPX chain) 30d:", round(spx['index_30d'], 3))
    print("VVIX (from VIX chain) 30d:", round(vix['index_30d'], 3))
