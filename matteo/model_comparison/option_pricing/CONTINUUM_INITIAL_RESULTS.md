# Continuous M2/stress response results — 27 October 2026

The user authorized a continuous probability measure of normalized geometric M2/stress components, with four VIX forward equality constraints. SPX remains paused; the original baseline is preserved.

## What was fitted

For u in [0,1], R_t(u)=G_t^(1-u)/E_Q[G_t^(1-u)] * exp(b*u*Z_t - (b*u)^2*q(t)/2). The variance is f+(xi(t)-f)*integral R_t(u) mu(du), with mu=Beta(m*kappa,(1-m)*kappa). One old reservoir and one independent OU are shared across the spectrum. No independent bank per quadrature node. All mixing occurs before the conditional variance integral, VIX square root and option payoff.

217 native quotes, four maturities (27,55,118,153 days), independently fitted seven active levels of the 16-piece forward-variance curve, fitted floor5–12% volatility. The original readout/cap and positive spike_shift are fixed. No clipping/smoothing or discarded signed estimates.

## Calibration and fresh numerical check

Training8192points/expiry (eight nets, previous checking streams explicitly reused as training); fresh16384points/expiry (eight new independent scrambled joint importance nets after selection). All market targets are calibration inputs; fresh refers to numerical randomness, not market out-of-sample.

| Candidate | Training option loss | Fresh option loss | Fresh total comparison score | Quotes in bid/ask | Max fresh F error |
|---|---:|---:|---:|---:|---:|
| continuous | 2.13781815 | 2.04817044 | 2.07945475 | 84/217 | 0.16560750 |
| boundary | 2.13296679 | 2.04123832 | 2.07217480 | 84/217 | 0.16486879 |
| previous | 2.17492493 | 1.82846550 | 1.83139891 | 80/217 | 0.04097156 |

The previous candidate is the last fitted-floor model, repriced on the same fresh streams. Boundary is its two-atom family refitted on this training sample. Continuous is selected solely by forward-feasible training option loss. The total comparison score retains the old forward-error penalty for comparability; the new optimizer minimizes only option loss subject to equality constraints.

## Forward constraints

The four equalities constrain the frozen training estimator, rather than adding a soft penalty. They do not certify literal population equality. Fresh residuals are retained without shifting model prices or forwards.

| Expiry | Target forward | Continuous model | Error | Eight-net SE |
|---|---:|---:|---:|---:|
| 27 | 20.10000000 | 20.08557726 | -0.01442274 | 0.02209122 |
| 55 | 20.73000000 | 20.56439250 | -0.16560750 | 0.02480476 |
| 118 | 21.19000000 | 21.13026535 | -0.05973465 | 0.02771286 |
| 153 | 21.72000000 | 21.71093223 | -0.00906777 | 0.02931952 |

## Selected measure and remaining bounds

```json
{
  "b": 3.2160615999801028,
  "measure_mean": 0.6594439218853119,
  "kappa": 0.010000000000000005,
  "alpha": 0.006594439218853122,
  "beta": 0.003405560781146883,
  "half_life_days": 110.78072386526819,
  "floor_variance": 0.014399999999999996,
  "floor_volatility": 0.11999999999999998,
  "active_levels": [
    0.014399999999999996,
    0.04759646218097119,
    0.0633915725766426,
    0.054813343940472645,
    0.19999999999999996,
    0.06451340657861795,
    0.07472918587730684
  ]
}
```

All bounded multistart results are in optimizer_ledger.json; no global optimum is claimed. In particular, a small kappa concentrates the continuous density near the endpoints and approximates the old two-atom model. A boundary solution is not evidence of a useful spread of intermediate components. Floor and curve-bound saturation must remain visible.

## Paired comparison

```json
[
  {
    "candidate": "continuous",
    "reference": "boundary",
    "loss_change": 0.007279954676977063,
    "paired_jackknife_se": 0.0008714843593290897,
    "leave_one_net_changes": [
      0.0069836567029191166,
      0.006743298447632107,
      0.007581509716299006,
      0.007442684212100037,
      0.007267738378292421,
      0.007825529008817167,
      0.007385979747768889,
      0.00702352205330703
    ]
  },
  {
    "candidate": "continuous",
    "reference": "previous",
    "loss_change": 0.2480558455104891,
    "paired_jackknife_se": 0.03374661173863103,
    "leave_one_net_changes": [
      0.2559975848899043,
      0.22754489661581467,
      0.2620290741330411,
      0.25141474297363486,
      0.26028267377829883,
      0.25116632967164465,
      0.24645969824324587,
      0.22683258223085256
    ]
  }
]
```

The paired jackknife SE measures numerical checking variability conditional on the selected models, not fitting variability or generalization to another market date.

## Native low-strike puts and tails

Signed controlled estimates are preserved. Separate positive-payoff raw native put diagnostics appear in audit/raw_native_puts.csv; negative estimates and their SEs appear in audit/negative_native_prices.csv. Neither raw estimates nor Black-IV omissions are substituted into calibration. VIX second-moment and highest-strike payoff ESS/concentration appear in fresh_check/score_continuous.json.

## Verification

- Original endpoints replayed (u=0 original M2, u=1 pure stress). Independent analytic linear-Gaussian power oracle and Beta moment checks passed.
- Twelve price/forward derivative directions checked against central differences.
- Independent NumPy powered smooth-cap oracle and off-grid log-power interpolation checks passed.
- Independent physical-parameter VIX integration replay passed. All32 fresh Sobol/mixture scrambles regenerated exactly.
- Source/input hashes checked, floor support and all nine inactive curve levels checked. No new Lean theorem or population certificate.

Measure quadrature refinement24→48 on every fresh state:

```json
{
  "coarse": 24,
  "fine": 48,
  "max_price_difference": 2.395561082835229e-09,
  "max_forward_difference": 4.1418246610192e-09,
  "by_day": [
    {
      "day": 27,
      "max_VIX_difference": 1.867460014182143e-06
    },
    {
      "day": 55,
      "max_VIX_difference": 1.969140430446714e-05
    },
    {
      "day": 118,
      "max_VIX_difference": 0.00016425305511802435
    },
    {
      "day": 153,
      "max_VIX_difference": 0.0001810321118682623
    }
  ]
}
```

Combined u33/measure24/GH32/time8→u65/measure48/GH48/time16 refinements use three deliberately selected stressed states per expiry (largest weighted far-right payoff, minimum VIX, largest old-state norm); they are sampled diagnostics, not uniform error bounds.

```json
[
  {
    "day": 27,
    "probe_count": 3,
    "max_VIX_difference": 1.7728220029766817e-06,
    "ids": [
      6527,
      13407,
      13747
    ],
    "coarse": "u33/measure24/GH32/time8",
    "fine": "u65/measure48/GH48/time16"
  },
  {
    "day": 55,
    "probe_count": 3,
    "max_VIX_difference": 3.456615445429634e-06,
    "ids": [
      7535,
      10926,
      15528
    ],
    "coarse": "u33/measure24/GH32/time8",
    "fine": "u65/measure48/GH48/time16"
  },
  {
    "day": 118,
    "probe_count": 3,
    "max_VIX_difference": 3.3984466085712484e-06,
    "ids": [
      6378,
      12464,
      13844
    ],
    "coarse": "u33/measure24/GH32/time8",
    "fine": "u65/measure48/GH48/time16"
  },
  {
    "day": 153,
    "probe_count": 3,
    "max_VIX_difference": 2.398836841166485e-06,
    "ids": [
      3028,
      5406,
      11100
    ],
    "coarse": "u33/measure24/GH32/time8",
    "fine": "u65/measure48/GH48/time16"
  }
]
```

## Files and reproducibility

- Sources: continuum_response_core.py, continuum_response_study.py, continuum_power_kernel.cpp.
- Independent audit: audit_continuum_response.py.
- Plots: visualizations/continuum-response022 (PNG and SVG), including smiles, native prices, forwards/curve and fitted measure.
- Experimental model: runs/m2_continuum_response_002/experimental_model_continuous.json.
- Current selected original production baseline is unchanged.

The initial run001 stopped at an overly strict Beta moment tolerance: observed floating point error1.18e-11 versus1e-11. Run002 raises the oracle tolerance to1e-9 without changing the integrator/model; exact training power tables were reused with hashes and source/input matches. The interrupted numerical test is preserved.

This report is numerical evidence for the specified bounded family; it neither proves structural impossibility nor establishes a globally optimal or jointly calibrated SPX/VIX model.
