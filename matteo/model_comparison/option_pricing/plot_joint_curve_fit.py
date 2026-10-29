"""Rebuild an inspectable chart from saved fresh checks; no new model pricing.

Run with the option-pricing Python environment (numpy, pandas, scipy).
All four expiries of both instruments remain visible. IV is a display transform
using the frozen MARKET forward for both candidates, including VIX. Prices that
do not admit a finite Black IV are recorded, never clipped into admissibility.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import ndtr

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "runs/m2_joint_curve_fit_002"
OUT = ROOT / "visualizations/joint-curve-fit-002"


def black_price(vol: float, f: float, k: float, t: float, side: str) -> float:
    sign = 1 if side == "call" else -1
    if vol == 0:
        return max(sign * (f - k), 0.0)
    s = vol * np.sqrt(t)
    d1 = np.log(f / k) / s + s / 2
    return float(sign * (f * ndtr(sign * d1) - k * ndtr(sign * (d1 - s))))


def implied_vol(price: float, f: float, k: float, t: float, side: str) -> float | None:
    intrinsic = max(f - k, 0.0) if side == "call" else max(k - f, 0.0)
    upper = f if side == "call" else k
    if not intrinsic <= price < upper:
        return None
    if price == intrinsic:
        return 0.0
    try:
        return 100 * brentq(lambda v: black_price(v, f, k, t, side) - price,
                           0, 64, xtol=1e-12)
    except ValueError:
        return None


def static_figures(data: dict) -> list[Path]:
    """Matplotlib exports provide an offline, reviewable view of the same data."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter, MaxNLocator

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    colors = {"market": "#30343b", "anchor": "#c17621", "selected": "#246cb1"}
    outputs = []

    def save(fig, stem):
        for suffix in ("png", "svg"):
            p = OUT / f"{stem}.{suffix}"
            fig.savefig(p, dpi=160, bbox_inches="tight")
            outputs.append(p)
        plt.close(fig)

    def curve_plot(ax, facet, mode):
        rows = facet["points"]
        x = np.array([p["x"] for p in rows])
        if mode == "iv":
            vals = {k: np.array([p["iv_" + k] if p["iv_" + k] is not None else np.nan
                                 for p in rows])
                    for k in ("bid", "ask", "mid", "anchor", "selected", "anchorlo", "anchorhi", "selectedlo", "selectedhi")}
            ylabel = "Annualized Black IV (%)"
        else:
            mid = np.array([p["mid"] for p in rows])
            scale = np.array([p["scale"] for p in rows])
            vals = {k: (np.array([p[k] for p in rows]) - mid) / scale
                    for k in ("bid", "ask", "mid", "anchor", "selected")}
            for key, se in (("anchor", "ase"), ("selected", "sse")):
                error = np.array([p[se] for p in rows]) * 2 / scale
                vals[key + "lo"], vals[key + "hi"] = vals[key] - error, vals[key] + error
            ylabel = "Midpoint error (quote scales)"
        # Keep call/put sides separate, especially when the VIX model and market
        # forwards differ. Joining their IVs would invent an intervening curve.
        for j, side in enumerate(("put", "call")):
            mask = np.array([p["side"] == side for p in rows])
            xx = x[mask]
            ax.fill_between(xx, vals["bid"][mask], vals["ask"][mask], color=colors["market"],
                            alpha=.16, label="Market bid–ask" if j == 0 else None)
            ax.plot(xx, vals["mid"][mask], color=colors["market"], lw=1)
            for key, style, label in (("anchor", "--", "Anchor ±2 sampling SE"), ("selected", "-", "Joint fit ±2 sampling SE")):
                ax.fill_between(xx, vals[key + "lo"][mask], vals[key + "hi"][mask], color=colors[key], alpha=.12)
                ax.plot(xx, vals[key][mask], style, color=colors[key], lw=1.7, label=label if j == 0 else None)
        ax.set_xscale("log")
        ax.set_xlim(x.min() * .995, x.max() * 1.005)
        ticks = np.linspace(x.min(), x.max(), 5)
        ax.set_xticks(ticks, [f"{v:.2f}" for v in ticks])
        ax.set_xticks([], minor=True)
        ax.set_xlabel("Strike / market forward (log scale)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{facet['instrument']} · {facet['days']} days · {len(rows)} options", loc="left")
        ax.grid(alpha=.2)
        ax.axvline(1, color="0.65", lw=.7)
        ax.yaxis.set_major_locator(MaxNLocator(5))

    for inst in ("SPX", "VIX"):
        facets = [f for f in data["facets"] if f["instrument"] == inst]
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharey=True, layout="constrained")
        for ax, facet in zip(axes.flat, facets):
            curve_plot(ax, facet, "iv")
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncols=3)
        fig.suptitle(f"{inst}: market vs anchor and joint fit — pooled fresh checks\n"
                     "Fixed market forwards; calls/puts separated. Bands: sampling error only; numerical bias unresolved.", fontsize=13)
        save(fig, inst.lower() + "-implied-volatility")

    fig, axes = plt.subplots(4, 2, figsize=(12, 14), layout="constrained")
    for ax, facet in zip(axes.flat, data["facets"]):
        curve_plot(ax, facet, "residual")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncols=3)
    fig.suptitle("Price errors relative to market midpoint — all objective options\n"
                 "Inside gray band = inside bid/ask. Frozen quote scales; separate y scales. Bands: ±2 sampling SE only.", fontsize=12)
    save(fig, "option-price-errors")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), layout="constrained")
    f = data["forwards"]
    days = [r["x"] for r in f]
    axes[0].plot(days, [r["market"] for r in f], "D-", color=colors["market"], label="Market")
    for key, style, se in (("anchor", "o--", "ase"), ("selected", "o-", "sse")):
        axes[0].errorbar(days, [r[key] for r in f], yerr=[2 * r[se] for r in f],
                         fmt=style, color=colors[key], capsize=3, label=key.title() if key == "anchor" else "Joint fit")
    axes[0].set(xlabel="Maturity (days)", ylabel="VIX forward (points)", title="All four VIX forwards · ±2 sampling SE")
    axes[0].legend()
    ends = data["curves"]["ends_days"]
    x = [0, ends[0], ends[1], 365]
    for key, style in (("anchor", "--"), ("selected", "-")):
        v = 100 * np.sqrt(data["curves"][key])
        axes[1].step(x, [*v, v[-1]], where="post", ls=style, color=colors[key],
                     label=key.title() if key == "anchor" else "Joint fit")
    axes[1].set(xlabel="Time (days)", ylabel="Square root of annual variance (%)", title="Fitted variance input · first year")
    axes[1].legend()
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle("What the curve optimization changed", fontsize=15)
    fig.supxlabel("The variance input is not a VIX forward. Stock-grid, normalizer and VIX quadrature errors are not included in the bands.", fontsize=10)
    save(fig, "forwards-and-variance-input")
    return outputs


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data: dict = {"facets": [], "forwards": [], "curves": {}, "source": str(RUN)}
    files = [RUN / f"{kind}_fresh_pooled_{candidate}.csv"
             for kind in ("prices", "forwards") for candidate in ("anchor", "selected")]
    files += [RUN / "scores.json", Path(__file__), OUT / "joint-fit-template.html"]
    a = pd.read_csv(files[0])
    b = pd.read_csv(files[1])
    keys = ["instrument", "expiry_days", "strike", "option_type"]
    pd.testing.assert_frame_equal(a[keys], b[keys])
    for col in ["forward", "price_bid_U", "price_ask_U", "price_mid_U", "price_scale_U"]:
        np.testing.assert_array_equal(a[col], b[col])
    assert len(a) == 882 and (a.instrument == "SPX").sum() == 707
    assert set(a[a.instrument == "VIX"].expiry_days) == {27, 55, 118, 153}
    invalid = []
    max_repricing_error = 0.0
    for (instrument, days), rows in a.groupby(["instrument", "expiry_days"], sort=True):
        points = []
        for i, r in rows.sort_values("K_over_F").iterrows():
            s = b.loc[i]
            p = dict(x=float(r.K_over_F), k=float(r.strike), side=r.option_type,
                     f=float(r.forward), t=float(r.T_years), scale=float(r.price_scale_U),
                     mid=float(r.price_mid_U), bid=float(r.price_bid_U), ask=float(r.price_ask_U),
                     anchor=float(r.model_price_U), selected=float(s.model_price_U),
                     ase=float(r.mc_se_U), sse=float(s.mc_se_U))
            # IV and band endpoints are each inverted independently. Missing bounds
            # remain null; neither signed prices nor bounds are silently truncated.
            prices = {key: p[key] for key in ("bid", "ask", "mid", "anchor", "selected")}
            for candidate, se in (("anchor", p["ase"]), ("selected", p["sse"])):
                prices[candidate + "lo"] = p[candidate] - 2 * se
                prices[candidate + "hi"] = p[candidate] + 2 * se
            for key, price in prices.items():
                iv = implied_vol(price, p["f"], p["k"], p["t"], p["side"])
                p["iv_" + key] = iv
                if iv is None:
                    invalid.append(dict(instrument=instrument, days=int(days), strike=p["k"],
                                        field=key, price=price))
                else:
                    error = abs(black_price(iv / 100, p["f"], p["k"], p["t"], p["side"]) - price)
                    max_repricing_error = max(max_repricing_error, error)
            points.append(p)
        data["facets"].append(dict(instrument=instrument, days=int(days), points=points))
    af = pd.read_csv(files[2])
    bf = pd.read_csv(files[3])
    pd.testing.assert_frame_equal(af[["expiry_days", "forward_vix_points"]],
                                  bf[["expiry_days", "forward_vix_points"]])
    for i, r in af.iterrows():
        data["forwards"].append(dict(x=int(r.expiry_days), market=float(r.forward_vix_points),
                                     anchor=float(r.model_forward), selected=float(bf.loc[i, "model_forward"]),
                                     ase=float(r.mc_se), sse=float(bf.loc[i, "mc_se"])))
    scores = json.loads((RUN / "scores.json").read_text())
    for candidate in ("anchor", "selected"):
        data["curves"][candidate] = next(r["levels"] for r in scores if r["candidate"] == candidate)
    data["curves"]["ends_days"] = [0.25 * 365, 0.75 * 365, 5.832876712328767 * 365]
    data["invalid_iv"] = invalid
    data["validation"] = dict(rows=len(a), max_iv_repricing_error_U=max_repricing_error,
                              missing_iv_values=len(invalid))
    assert max_repricing_error < 1e-7
    raw = json.dumps(data, separators=(",", ":"), allow_nan=False)
    (OUT / "plot-data.json").write_text(raw + "\n")
    template = (OUT / "joint-fit-template.html").read_text()
    assert template.count("__PLOT_DATA__") == 1
    chart = template.replace("__PLOT_DATA__", raw)
    (OUT / "joint-fit.html").write_text(chart)
    outputs = static_figures(data)
    manifest = dict(source_hashes={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in files},
                    output_sha256=hashlib.sha256(chart.encode()).hexdigest(),
                    validation=data["validation"],
                    interpretation="Pooled fresh checks; fixed model/grid/quadrature. Bands are +/-2 sampling SE only.",
                    iv_convention="Undiscounted Black, fixed market forward, annualized percent; missing values retained in data.")
    manifest["static_output_hashes"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(data["validation"]))


if __name__ == "__main__":
    main()
