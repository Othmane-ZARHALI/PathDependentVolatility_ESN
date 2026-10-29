# Continuous M2/stress measure — completed VIX-only test, 27 October 2026

The continuous construction is implemented and tested with four VIX forward equality constraints. In the bounded searches it has not improved the two-component reference: the best Beta measure concentrates near the endpoints, while a broader continuous measure fits less well and reaches the stress-strength bound. Keep the previous working VIX fit; this experiment is not adopted. SPX remains paused.

## Model

For u in [0,1],

    R_t(u)=G_t^(1-u)/E_Q[G_t^(1-u)]
           * exp(b*u*Z_t - .5*(b*u)^2*q(t)),
    v_t=f+(xi(t)-f)*integral R_t(u) mu(du),
    mu=Beta(m*kappa,(1-m)*kappa), q(t)=1-exp(-2*a*t).

G is the unchanged capped M2 readout and Z is one independent OU. Components share the original reservoir and Z; no bank per node. Every power has its own normalizer and conditional expectation. The integral is formed before the VIX square root and option payoff. mu is a deterministic probability measure under Q, whose mass describes expected non-floor variance contributions, not regime probabilities.

All217 native VIX quotes and all four forwards (27,55,118,153days). Seven independent active xi levels of the16-piece curve; nine inactive levels unchanged. Floor volatility5–12%; b0–5; m.01–.99; kappa.01–30; half-life1–365days. Original readout/cap, positive beta_J, architecture and premiums frozen. No clipping, quote removal, curve smoothing, jumps, SPX repricing or baseline replacement.

## Forward condition and estimator repair

The optimizer minimizes native bid/ask scaled option-price distance subject to four forward equalities, not a forward penalty. Equality is enforced on the calibration estimator to numerical tolerance; a finite simulation does not certify exact population means.

The initial experiment002 used8192 powered paths/expiry and had a55-day checking gap-.1656 versus.0248SE. The corresponding two-atom refit also failed. Its training forward estimates were noisy, partly because legacy frozen controls were poorly matched to the new response. That calibration is rejected as reliable fulfillment of the forward condition; its outputs and independent model audit are preserved. See CONTINUUM_INITIAL_RESULTS.md and run002/forward_estimator_diagnostics.json.

Run003 keeps the exact model and uses the multilevel identity

    E[P_cont]=E[P_boundary]+E[P_cont-P_boundary].

The reference boundary uses the same b,m,half-life,floor and curve. Calibration estimates its cheap first term on32768 existing paths, plus the paired difference on8192 powered paths. Matching eight-net blocks retain the covariance of the nested samples. All model prices still come from the continuous variance law; subtracting the numerical reference is variance reduction, not a second model component or a mixture of finished option laws.

Checking uses65536 NEW paths/expiry for the boundary term (eight8192-point scrambles, seed20261027117), plus the existing16384 paired-difference checking paths. The difference stream is deliberately reused after002; this is not a wholly independent fresh study. Models and selection are frozen before generating the new large reference stream. Checking uncertainty and the reuse limitation are explicit.

## Fit comparison

Smaller price loss is better. Total comparison scores also include the historical forward-error penalty, retained for comparability; this penalty is not the new optimization objective.

| Model | Training option loss | Checking option loss | Checking total score | Within bid/ask | Max checking forward gap |
|---|---:|---:|---:|---:|---:|
| continuous | 1.73765986 | 1.73067537 | 1.73159404 | 83/217 | 0.02634676 |
| boundary | 1.73326449 | 1.72555497 | 1.72649073 | 83/217 | 0.02677712 |
| previous | 1.78148906 | 1.76584390 | 1.76686274 | 82/217 | 0.02444121 |
| interior | 1.80515466 | 1.86362964 | 1.87188206 | 80/217 | 0.06234895 |

Previous is the last fitted-floor model repriced on the same checking streams. Boundary is refitted using the same large calibration estimator. Continuous is selected by forward-feasible training loss; interior is the best completed feasible start with kappa>=.1, using training only. All starts and local minima are retained; no global optimum is claimed.

## Four forward checks

| Expiry | Observed | Continuous checking estimate | Gap | Eight-net SE |
|---|---:|---:|---:|---:|
| 27 | 20.10000000 | 20.12634676 | +0.02634676 | 0.00357915 |
| 55 | 20.73000000 | 20.71736705 | -0.01263295 | 0.00539449 |
| 118 | 21.19000000 | 21.18680012 | -0.00319988 | 0.00490941 |
| 153 | 21.72000000 | 21.72739517 | +0.00739517 | 0.00702323 |

Calibration gaps are in optimizer_ledger.json; calibration SEs are reported in training/forwards_*.csv. Numerical equalities should be read with both calibration and checking uncertainty. No forward shifting is applied to plotted prices. Canonical call IVs use the market forward for every curve; native IVs use each model forward and preserve unavailable inversions.

## Fitted continuous measures

continuous: m=0.65944476, kappa=0.01000000, b=3.18767493, half-life=105.287565days; mass in u=[.1,.9] is 0.978%. Floor=12.000000%.

interior: m=0.57415911, kappa=0.99014537, b=5.00000000, half-life=220.819719days; mass in u=[.1,.9] is 57.354%. Floor=12.000000%.

The main kappa lower-bound hit approximates the previous two-atom construction. The broader measure remains genuinely spread, but its b=5 ceiling and12% floor hit limit conclusions about its best possible fit. Forward-variance bounds are also visible in curve_levels.csv. These results do not exclude richer continuous measures, a different interpolation path or wider stress-strength bounds.

## Paired checking differences

```json
[
  {
    "candidate": "continuous",
    "reference": "boundary",
    "loss_change": 0.005103307291186798,
    "paired_jackknife_se": 0.0008166849073044103,
    "leave_one_net_changes": [
      0.004938941103217465,
      0.00465765091796988,
      0.00534605144668876,
      0.004983925049481419,
      0.004838957928025245,
      0.005722966715689726,
      0.0051926659066872105,
      0.005143440344829386
    ]
  },
  {
    "candidate": "continuous",
    "reference": "previous",
    "loss_change": -0.03526870486380895,
    "paired_jackknife_se": 0.0014288813348167977,
    "leave_one_net_changes": [
      -0.035677936139764554,
      -0.0362059099929497,
      -0.03433834663964119,
      -0.03553852885462949,
      -0.034770622620115654,
      -0.03496736338274431,
      -0.03527332043337217,
      -0.0353700541936155
    ]
  },
  {
    "candidate": "interior",
    "reference": "boundary",
    "loss_change": 0.14539132822639433,
    "paired_jackknife_se": 0.0405822763988726,
    "leave_one_net_changes": [
      0.13893515251197686,
      0.12312909931823457,
      0.1560997050629216,
      0.14260698685264228,
      0.13355806373734502,
      0.176819036473296,
      0.14401768468621756,
      0.1553783550992347
    ]
  },
  {
    "candidate": "interior",
    "reference": "previous",
    "loss_change": 0.10501931607139858,
    "paired_jackknife_se": 0.04059826158964984,
    "leave_one_net_changes": [
      0.09831827526899484,
      0.08226553840731499,
      0.11641530697659164,
      0.10208453294853137,
      0.09394848318920412,
      0.13612870637486196,
      0.10355169834615818,
      0.11486486056078982
    ]
  },
  {
    "candidate": "interior",
    "reference": "continuous",
    "loss_change": 0.14028802093520754,
    "paired_jackknife_se": 0.039789370256759146,
    "leave_one_net_changes": [
      0.1339962114087594,
      0.11847144840026469,
      0.15075365361623283,
      0.13762306180316086,
      0.12871910580931978,
      0.17109606975760627,
      0.13882501877953035,
      0.1502349147544053
    ]
  }
]
```

Eight-net paired jackknife errors measure conditional numerical checking variability, not fitting variability, global optimality or market out-of-sample performance.

## Verification and unresolved limits

- Original M2 and pure-stress endpoints, Beta moments, independent analytic Gaussian power oracle,12-direction continuum and12-direction multilevel price/forward Jacobians passed.
- Independent NumPy smooth-cap/power oracle, off-grid power checks and combined u33/measure24/GH32/time8 versus u65/measure48/GH48/time16 probes passed.
- Independent direct model-VIX replay is scale-aware (absolute1e-9 OR relative1e-11): rare tiny-weight proposal paths have very large VIX values; absolute, relative and weighted-mean errors are retained. No unweighted maximum is used as tail evidence.
- Independent reconstruction of base-plus-paired-difference prices/forwards, all eight blocks and both candidate models passed.
- All32 new large reference scrambles and32 reused difference scrambles regenerated exactly. Source/input hashes, floor support and all nine inactive xi levels checked.
- Signed native estimates and unavailable IVs retained. Separate raw positive-payoff put diagnostics and negative-price/SE diagnostics are saved, without substituting them into calibration.
- Second-moment targets are exact model identities; numerical moment estimates and SEs appear in multilevel_audit/second_moments.json. Full bid/ask feasibility is not established.
- Sampled deterministic quadrature checks are not uniform error bounds. No new Lean theorem or exact population-forward certificate.

Combined sampled refinement:

```json
[
  {
    "day": 27,
    "probe_count": 3,
    "max_VIX_difference": 1.8128710053133545e-06,
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
    "max_VIX_difference": 3.565720589904231e-06,
    "ids": [
      6605,
      7535,
      15528
    ],
    "coarse": "u33/measure24/GH32/time8",
    "fine": "u65/measure48/GH48/time16"
  },
  {
    "day": 118,
    "probe_count": 3,
    "max_VIX_difference": 3.464787027951388e-06,
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
    "max_VIX_difference": 2.6254806471115444e-06,
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

## Reproducible artifacts

- Main run: runs/m2_continuum_response_003. Initial rejected calibration: run002; initial overly strict Beta-moment assertion: run001, preserved.
- Mathematical and numerical core: continuum_response_core.py, continuum_power_kernel.cpp, continuum_multilevel_core.py, continuum_multilevel_study.py.
- Independent audits: audit_continuum_scaled_fast.py, audit_continuum_response_v2.py, audit_multilevel_continuum.py.
- Experimental models: experimental_model_continuous.json and experimental_model_interior.json.
- Plots (PNG/SVG): visualizations/continuum-response023, including all strikes, right wing, native puts, forward errors/curve, and measure comparisons.
- Frozen calibration sources and postprocessing sources are saved; artifact_manifest.json records hashes.

The original selected production baseline and previous working VIX candidate are preserved.

## Subsequent explanatory notes

See [the empirical reduction and its mathematical meaning](CONTINUUM_EMPIRICAL_REDUCTION.md) and [the27day wing diagnosis](VIX_27_DAY_DIAGNOSIS.md). These add documentation and read-only diagnostics; the archived calibration, model parameters and original audit manifests are unchanged.
