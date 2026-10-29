# Focused SPX maturity view — 23 October 2026

The user reaffirmed cutoff **3** and asked to avoid SPX maturities **1, 1401,
1765 and 2129 days**. Current displayed/active SPX maturities are therefore
**29, 57, 119, 211, 302, 666 and 1030 days**. All four VIX maturities and forwards
remain present. The previous all-11 view and all raw data are retained.

[SPX figure](spx-focused-maturities.png) and [variance levels](variance-three-levels.png)
also have SVG versions. The conversation view is `joint-curve-focus.html`.

The existing four fitted SPX slices reuse their exact 707 cutoff-3 training
rows, market quotes and pooled model estimates. There is no new optimization,
model pricing, numerical setting change, or change to the fitted objective.
The 119-day slice remains outside fitting; 666/1030 days remain diagnostics.
Their curves are **market only**, without model overlays.

Cutoff 3 is also applied to the three additional market-only displays. These
expiries have no training split, so their display-only volatility proxy is the
Black IV of their own nearest-to-ATM market midpoint quote, using the same
absolute log-moneyness rule and deterministic tie-breaking. This changes neither
the frozen training-only proxies at the fitted expiries nor the fitted membership.
`spx-maturity-coverage.csv` records original eligible counts and display counts,
proxy strikes, and proxy values. The earlier exposure of the reserved 119-day
market curve remains documented; no held-out model score is computed.

## Reading the two sides of the line

In each option panel, the vertical line is K/F = 1. Below it the chosen quotes
are OTM puts; above it they are OTM calls. Puts price protection against falls,
and calls price exposure to rises. They probe different sides of **the same**
underlying terminal distribution, so a smooth asymmetric smile is meaningful.

Changing from puts to calls does not by itself justify an abrupt jump in Black
IV. At a fixed strike and consistent forward, put-call parity links the prices
and gives the same Black IV. The visible VIX model call/put separation can be
caused by using the fixed market forward to invert prices from a model whose
own mean differs; separately controlled finite-sample price estimates can also
have parity errors. Do not automatically diagnose the gap as a discontinuous
model law or explain it away as an inevitable OTM switch.

The third variance level still continues to the original 2129-day horizon in
the saved numerical model. Omitting late maturity panels does not refit that
level or shorten the historical curve. The prior qualification limits remain:
model shading is sampling uncertainty only, with timestep, normalizer and
quadrature error outside the bands; no qualified joint fit has been established.

## Rebuild

```sh
MPLCONFIGDIR=/private/tmp/pdv-matplotlib XDG_CACHE_HOME=/private/tmp/pdv-font-cache \
 /private/tmp/pdv-options-venv/bin/python \
 matteo/model_comparison/option_pricing/plot_joint_curve_focus.py
```

The manifest hashes inputs, scope, source and outputs. PNG layout was visually
inspected, data membership and unchanged original slices were checked, and
JavaScript syntax was checked. Browser interaction testing remains unavailable
as recorded in the preceding visualization guides.
