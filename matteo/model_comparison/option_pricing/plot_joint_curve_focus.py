"""Focused SPX display: retain cutoff 3 and omit four user-deferred expiries.

Does not price, optimize, score held-out contracts, or change the original plots.
The 119-day reserved market curve is now visually exposed at user request.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from plot_joint_curve_fit import implied_vol, black_price

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[2]
PREVIOUS = ROOT / "visualizations/joint-curve-fit-002"
OUT = ROOT / "visualizations/joint-curve-focus-004"
MARKET = PROJECT / "data/processed/options/2026-02-19/quotes_calibration.csv"
SCOPE = ROOT / "active_scope.json"
ROLES = {
    "fit": "used for fitting",
    "final_maturity_test": "reserved maturity · market only",
    "short_expiry_diagnostic": "1-day excluded · market only",
    "maturity_diagnostic": "maturity diagnostic · market only",
    "extrapolation_diagnostic": "extrapolation · market only",
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = json.loads((PREVIOUS / "plot-data.json").read_text())
    quotes = pd.read_csv(MARKET)
    scope = json.loads(SCOPE.read_text())
    maturity_scope = scope["spx_maturity_scope"]
    kept = maturity_scope["analysis_expiry_days"]
    cutoff = scope["spx_research_priority"]["criterion"]["cutoff"]
    assert cutoff == 3 and kept == [29,57,119,211,302,666,1030]
    spx = quotes[quotes.instrument == "SPX"]
    coverage = []
    roundtrip = 0.0
    for days, rows in spx.groupby("expiry_days", sort=True):
        if int(days) not in kept:
            continue
        existing = next((f for f in data["facets"] if f["instrument"] == "SPX" and f["days"] == days), None)
        role = rows.maturity_role.iloc[0]
        assert rows.maturity_role.nunique() == 1
        info = dict(days=int(days), role=role, eligible_otm_quotes=len(rows),
                    training_quotes=int((rows.split == "train").sum()),
                    objective_quotes=len(existing["points"]) if existing else 0,
                    model_overlay=existing is not None)
        coverage.append(info)
        if existing:
            existing["role_label"] = "fit · cutoff 3"
            continue
        # Additional slices are visual diagnostics, never new fitted members.
        # Their proxy uses the nearest-ATM market quote from that slice because
        # those maturities have no training split. It does not alter the frozen
        # training-only proxies or membership at the four fitted maturities.
        rows = rows.assign(distance=np.abs(np.log(rows.K_over_F)))
        near = rows.sort_values(["distance","strike","option_type"]).iloc[0]
        proxy = implied_vol(float(near.price_mid_U),float(near.forward),float(near.strike),
                            float(near.T_years),near.option_type)
        assert proxy is not None and proxy > 0
        z = rows.distance / ((proxy/100)*np.sqrt(rows.T_years))
        rows = rows[z <= cutoff]
        info["display_quotes"] = len(rows)
        info["display_proxy_volatility"] = proxy/100
        info["display_proxy_strike"] = float(near.strike)
        info["proxy_source"] = "nearest ATM quote within this market-only slice; visualization only"
        points = []
        for _, r in rows.sort_values("K_over_F").iterrows():
            p = dict(x=float(r.K_over_F), k=float(r.strike), side=r.option_type,
                     f=float(r.forward), t=float(r.T_years), scale=float(r.price_scale_U),
                     mid=float(r.price_mid_U), bid=float(r.price_bid_U), ask=float(r.price_ask_U),
                     anchor=None, selected=None, ase=None, sse=None)
            for key in ("mid", "bid", "ask"):
                iv = implied_vol(p[key], p["f"], p["k"], p["t"], p["side"])
                assert iv is not None
                p["iv_" + key] = iv
                roundtrip = max(roundtrip, abs(black_price(iv/100, p["f"], p["k"], p["t"], p["side"])-p[key]))
            for key in ("anchor", "selected", "anchorlo", "anchorhi", "selectedlo", "selectedhi"):
                p["iv_" + key] = None
            points.append(p)
        data["facets"].append(dict(instrument="SPX", days=int(days), role_label=ROLES[role]+" · cutoff 3", points=points))
    data["facets"].sort(key=lambda f: (f["instrument"], f["days"]))
    data["coverage"] = coverage
    assert len(coverage) == 7 and sum(r["model_overlay"] for r in coverage) == 4
    assert sum(r["objective_quotes"] for r in coverage) == 707
    assert [f["days"] for f in data["facets"] if f["instrument"] == "VIX"] == [27,55,118,153]
    assert roundtrip < 1e-7

    # Derive a new literal template from the preserved original. Guard every
    # replacement so later template drift cannot silently change the view.
    template = (PREVIOUS / "joint-fit-template.html").read_text()
    replacements = {
        "m2-joint-fit-plot": "m2-joint-focus-plot",
        "m2-plot-view": "m2-focus-view",
        "m2-iv-status": "m2-focus-status",
        "Joint SPX / VIX fit · fresh checks": "SPX / VIX · cutoff 3 · selected maturities",
        "${facet.days} days · ${facet.points.length} options": "${facet.days} days · ${facet.role_label || 'fit'} · ${facet.points.length} options",
        "}return r;": "}for(const key of ['anchor','selected'])if(p[key]===null)r[key]=[null,null,null];return r;",
        "if(!u[key]||!v[key])return null;": "if(!u[key]||!v[key]||!rows.some(r=>Number.isFinite(r[key]?.[0])))return null;",
        "`Uninvertible IV: ${means} model means; ${bounds} band endpoints`": "`7 SPX expiries · 4 model overlays · 3 market only · ${bounds} unavailable IV-band endpoints`",
        "'SPX cutoff 3 · 707 SPX + 175 VIX options'": "'707 SPX + 175 VIX fit options; extra SPX panels are market only'",
    }
    for old, new in replacements.items():
        assert old in template, old
        template = template.replace(old, new)
    (OUT / "context-template.html").write_text(template)
    raw = json.dumps(data, separators=(",", ":"), allow_nan=False)
    (OUT / "plot-data.json").write_text(raw + "\n")
    chart = template.replace("__PLOT_DATA__", raw)
    assert len(chart.encode()) < 1_000_000
    (OUT / "joint-curve-focus.html").write_text(chart)
    pd.DataFrame(coverage).to_csv(OUT / "spx-maturity-coverage.csv", index=False)
    pngs = render_static(data)
    sources = [MARKET, SCOPE, PREVIOUS / "plot-data.json", PREVIOUS / "joint-fit-template.html", Path(__file__), ROOT / "plot_joint_curve_fit.py"]
    manifest = dict(source_sha256={str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                    output_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [OUT / "joint-curve-focus.html", OUT / "context-template.html", OUT / "plot-data.json", OUT / "spx-maturity-coverage.csv", *pngs]},
                    additional_market_iv_repricing_error_U=roundtrip,
                    model_overlay_expiries=[29,57,211,302],
                    market_only_expiries=[r["days"] for r in coverage if not r["model_overlay"]],
                    excluded_display_expiries=maturity_scope["deferred_expiry_days"],
                    cutoff=3,
                    cutoff_diagnostic_policy="Existing model curves/proxies unchanged. Additional market-only slices use nearest-ATM mid-IV from their own maturity for display only.",
                    reserved_119_day_status="Excluded from all fitting; market curve visually exposed at user request. No model repricing or held-out score computed.",
                    changes="No model, loss, numerical settings, original run artifacts or fitted membership changed.")
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps({"spx_maturities":len(coverage),"model_overlays":4,"market_only":3,"max_added_iv_roundtrip_error":roundtrip,"html_bytes":len(chart.encode())}))


def render_static(data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    plt.rcParams.update({"font.size":10,"svg.fonttype":"none","axes.spines.top":False,"axes.spines.right":False})
    facets = [f for f in data["facets"] if f["instrument"] == "SPX"]
    fig, axes = plt.subplots(4,2,figsize=(13,14),layout="constrained")
    for ax, facet in zip(axes.flat, facets):
        rows = facet["points"]
        for j, side in enumerate(("put", "call")):
            group = [p for p in rows if p["side"] == side]
            x = np.array([p["x"] for p in group])
            values = lambda k: np.array([p[k] if p[k] is not None else np.nan for p in group])
            ax.fill_between(x,values("iv_bid"),values("iv_ask"),color="0.3",alpha=.15,label="Market bid–ask" if j == 0 else None)
            ax.plot(x,values("iv_mid"),color="0.2",lw=1)
            for k,style,color,label in (("anchor","--","#c17621","Anchor"),("selected","-","#246cb1","Joint fit")):
                if any(p[k] is not None for p in group):
                    ax.fill_between(x,values("iv_"+k+"lo"),values("iv_"+k+"hi"),color=color,alpha=.12)
                    ax.plot(x,values("iv_"+k),style,color=color,lw=1.5,label=label if j == 0 else None)
        ax.set_xscale("log")
        x = np.array([p["x"] for p in rows])
        ticks = np.linspace(x.min(),x.max(),4)
        ax.set_xticks(ticks,[f"{v:.2f}" for v in ticks]);ax.set_xticks([],minor=True)
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.set_xlabel("Strike / market forward (log scale)");ax.set_ylabel("Black IV (%)")
        ax.set_title(f"{facet['days']} days · {facet['role_label']}\n{len(rows)} displayed OTM quotes",loc="left",fontsize=11)
        ax.axvline(1,color="0.6",lw=.7);ax.grid(alpha=.15)
    axes.flat[-1].axis("off")
    handles,labels=axes.flat[1].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncols=3,frameon=False)
    fig.suptitle("SPX: seven selected expiries · cutoff 3\nMarket-only panels at 119/666/1030 days; model shading is ±2 sampling SE only.",fontsize=14)
    outputs=[]
    for suffix in ("png","svg"):
        p=OUT/f"spx-focused-maturities.{suffix}";fig.savefig(p,dpi=140,bbox_inches="tight");outputs.append(p)
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(10,4.6),layout="constrained")
    ends=data["curves"]["ends_days"]
    x=[0,ends[0],ends[1],365]
    for key,ls,color,name in (("anchor","--","#c17621","Anchor"),("selected","-","#246cb1","Joint fit")):
        v=100*np.sqrt(data["curves"][key]);ax.step(x,[*v,v[-1]],where="post",ls=ls,color=color,lw=2,label=name)
        for a,b,y in zip(x[:-1],x[1:],v):
            ax.text((a+b)/2,y+(0.32 if key=="selected" else -0.62),f"{y:.2f}%",ha="center",fontsize=11)
    ax.set_xticks([0,ends[0],ends[1],365],["0","91.25","273.75","365"])
    ax.set(xlabel="Time from the option snapshot (calendar days)",ylabel="Square root of annual variance (%)",title="Three deterministic variance levels · first year")
    ax.legend();ax.grid(alpha=.2)
    fig.supxlabel("The final level continues to day 2129; data used in this fit extend only to day 302.\nVariance levels themselves are the squares of the displayed decimal volatilities.",fontsize=10)
    for suffix in ("png","svg"):
        p=OUT/f"variance-three-levels.{suffix}";fig.savefig(p,dpi=150,bbox_inches="tight");outputs.append(p)
    plt.close(fig)
    return outputs


if __name__ == "__main__":
    main()
