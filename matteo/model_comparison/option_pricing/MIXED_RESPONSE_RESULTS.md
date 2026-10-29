# Two variance responses: VIX-only calibration with a free forward-variance curve

Completed 26 October 2026. Run `runs/m2_mixed_response_001`.
Plan: [MIXED_RESPONSE_PLAN.md](MIXED_RESPONSE_PLAN.md).
Implementation: `mixed_response_core.py`, `mixed_response_study.py`.
Independent replay: `audit_mixed_response.py` (PASS).

The additive M2/stress mixture produces a substantial steepening, with smaller
fresh price and IV-level errors. The longer-dated wing slopes are now close to
market slopes; the short-dated slopes overshoot. This is an improved exploratory
VIX fit, not a satisfactory full-smile calibration or an adopted joint model.

## Equations and preserved scope

Let R(t)=G_old(X_t)/M_old(t), and introduce the same independent Q OU as before:

    dZ=-a Z dt+sqrt(2a)dB, Z0=0, q(t)=1-exp(-2at),
    L_b(t)=exp(b Z_t-b²q(t)/2).

Two variants implement the two-component proposal from the research note:

    Product:  v=f+(xi(t)-f) R(t) [(1-w)L_b1(t)+w L_b2(t)].
    Additive: v=f+(xi(t)-f) [(1-w)R(t)+w L_b(t)].

The mixture is inside the pathwise variance integral. We do not average completed
option prices or replace variance by a mixture of volatilities. Mean-one Gaussian
components and independence preserve E[v(t)]=xi(t), so the old Gaussian normalizer
can be retained. Future conditional factors are evaluated analytically before
integrating over the 30-day VIX period. The stock martingale reasoning, under
independence and the old capped response, is documented in the plan; it is not a
new Lean proof or a sampling-convergence claim.

All **16 original forward-variance intervals** remain in each experimental model.
The seven VIX-sensitive levels are independently fitted alongside mixture
parameters; the other nine have zero sensitivity to these VIX contracts and
retain their previous values. No smoothing or monotonicity ties them together.
Bounds are xi in[0.00264,0.20], b1 in[0,3], b2 in[b1,5], w in[0,1], half-life
in[1,365] days. The old readout/cap, positive beta_J=0.369729973, banks, premiums,
initial state and floor0.00264 stay fixed. All217 native quotes and four forwards
remain; SPX is paused. The same native bid–ask price objective, expiry weighting
and0.5-point forward scale are used, without new wing weighting.

## Results on fresh simulations

Each quoted strike uses all16384 outer states per expiry in eight independent
2048-point scrambles. These are new numerical draws on the same calibration
market quotes, not market out-of-sample results. Baseline and previous fitted
single-factor parameters are replayed on exactly the same fresh proposal.
Loss values therefore differ slightly from earlier runs' estimates.

| Candidate | Training loss | Fresh loss | Fresh wing IV MAE (percentage points) | Quotes inside bid–ask |
|---|---:|---:|---:|---:|
| baseline | 25.5184 | 27.3338 | 33.89 | 17/217 |
| curve_only | 19.1813 | 20.5171 | 27.42 | 22/217 |
| previous_single_factor | 19.6616 | 20.0052 | 30.81 | 22/217 |
| product_price | 17.5555 | 18.2300 | 26.60 | 19/217 |
| product_shape | 17.5555 | 18.2302 | 26.60 | 19/217 |
| additive_price | 9.7697 | 10.2960 | 8.51 | 49/217 |
| additive_shape | 13.9628 | 14.8127 | 24.14 | 20/217 |

The additive price selection reduces fresh native VIX+forward loss from20.0052
to10.2960 relative to the prior fitted single factor: paired difference-9.7093,
eight-net jackknife SE0.2195. The product price selection changes it by-1.7752,
SE0.1549. These are conditional numerical sensitivity estimates, not uncertainty
about the complete model-selection process. Every optimizer outcome is retained.

Measured IV slopes versus log(K/F), using every retained quote above2F:

| Expiry (days) | Market | Prior single factor | Product mixture | Additive mixture |
|---|---:|---:|---:|---:|
| 27 | 52.38 | 27.73 | 25.26 | 75.52 |
| 55 | 38.31 | 27.03 | 26.86 | 56.19 |
| 118 | 31.23 | 15.95 | 21.41 | 35.81 |
| 153 | 27.32 | 10.31 | 20.80 | 29.35 |

The additive fit has enough slope at118/153days and too much at27/55days.
Its wing IV-level error falls from30.81 to8.51percentage points, and ATM IV MAE
from14.34 to9.66. It still misses the low-strike/near-ATM shape and remains below
the longer-dated market wing levels. Only49/217 native prices are within bid–ask;
fresh forward errors are +0.289,+0.176,+0.468,+0.886 VIX points. This does not
solve the whole calibration problem.

The descriptive shape selection minimizes the earlier slope/curvature score
among completed local fits, independently of native-price selection. It is not
necessarily closest in IV levels: the additive shape selection has fresh wing
MAE24.14 despite its better shape-score ranking. Both selections were frozen
before fresh checking and both are retained, without checking-based retuning.
Curvature here is measured against **log strike/forward**; convexity depends on
the horizontal coordinate. Paired slope/curvature jackknife results are saved in
`fresh_check/paired_wing_shape.json`.

## What the optimizer chose

The additive price winner has b=4.9999998 (at the declared strength bound),
w=0.9957389, half-life6.588days. About99.6% of its expected non-floor variance
contribution comes from the new stress component and0.4% from the old M2 response.
It therefore gains mainly by replacing the old response within this exploratory
readout, rather than making a small adjustment to M2. The mixture weight is not
the probability of entering a stress regime: both components act on every path.

The product winner has b1 approximately0, b2=4.3780, w=0.11877 and half-life365days
(at its upper bound). Several starts find nearly zero added response instead;
one reaches its120-evaluation budget. No global optimum is claimed.

The additive curve also increases the mean-variance budgets: exact E[VIX²] is
511.33,707.12,710.98,733.01. These exceed all four archived necessary bounds for
matching the native quotes plus exact forward. Thus the old second-moment
obstruction no longer excludes this candidate, although the full prices/forwards
remain imperfect. Two narrow curve pieces,55–57 and118–119days, reach the0.20
variance upper bound. The curve is not adopted as a production term structure.
The complete16 levels are in `curve_levels.csv` and each experimental model JSON.
The earlier necessary bounds and their exact-rational checks are recorded in
[runs/m2_vix_design_research_001/moment_certificate_audit.json](runs/m2_vix_design_research_001/moment_certificate_audit.json).

## Sampling and verification

The stronger exponential needs explicit stress sampling. Training reuses old
4096-state conditional ratios, with two separately scrambled importance samples
for Z per old state:8192 points/expiry, clustered into the original eight nets.
The fixed Z proposal has shifts0,2,4,6 with probabilities.5,.25,.2,.05. Exact
weights are w_old*w_Z, without self-normalization. Checking instead uses fresh
37-coordinate states from the68-component product proposal. Original independent
pilot control coefficients remain frozen; signed estimates are saved unchanged.
No selected candidate has a negative fresh native price or unavailable fresh IV.

For the additive price fit, likelihood-weighted last-strike payoff ESS is
1663–1819 and the top ten raw payoff contributions account for roughly0.9–1.0%
per expiry. This improves tail visibility, without proving population convergence.
Fresh sampled second moments agree with the exact budgets to about1.3reported
scramble SE or less. These diagnostics concern weighted sampling, not raw maxima.

72 implementation checks cover old/single-factor nesting,21Jacobian directions,
and Gaussian conditional/mean-one integral oracles. Maximum Jacobian discrepancy
is4.47e-4relative; Gaussian conditional discrepancies are below2.0e-15. Independent
nodal replay checks14candidate/stream price sets, exact interval budgets, all16
levels, unchanged off-window entries and positive beta_J; maximum native price
difference is2.3e-15. All32 fresh scrambles and product-proposal weights replay.
285source hashes and five input hashes pass; a new C++ library is compiled.

Sampled GH32/time8 to GH48/time16 refinement uses44–46states/expiry, including
large weighted payoffs and extreme Gaussian coordinates. Maximum VIX difference
is0.000118point for the additive price winner. Across all candidates it is0.0332
point, in the product-mixture extreme probes. These are sampled value differences,
not uniform option-price error bounds. No population convergence, SPX fit,
physical-measure compatibility or new Lean theorem is asserted.

## Plots and continuation

- [All four smiles](visualizations/mixed-response020/vix-all-strikes.png)
- [Right-wing comparison](visualizations/mixed-response020/vix-right-wing.png)
- [Alternative shape selections](visualizations/mixed-response020/vix-shape-selections.png)
- [All16 curve intervals](visualizations/mixed-response020/variance-curve.png)
- [Forwards and exact moments](visualizations/mixed-response020/forwards-and-moments.png)

PNG and SVG exports are saved. Pointwise sampling bands exclude systematic
quadrature error and uncertainty from frozen pilot controls. The old selected
model is unchanged. No SPX paths were read or priced.

This experiment identifies a response that can supply substantial right-wing
steepness. The remaining research question is how to adjust the short/long
maturity balance and the near-ATM shape while retaining that gain, with forward
variance still fitted. A second response timescale is a possible next experiment;
it has not been implemented by this run. No further model change is adopted.
