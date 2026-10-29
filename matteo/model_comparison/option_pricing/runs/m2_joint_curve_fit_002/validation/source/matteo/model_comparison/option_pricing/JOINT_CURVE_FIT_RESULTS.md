# Actual joint SPX/VIX curve fitting with Lean checks — 23 October 2026

**Actual joint optimization and all three fresh checks are complete.** The
selected curve lowers pooled fresh loss from **597.90 to 347.10** (41.95%),
but the fit remains poor: only **50/882 objective options** are inside bid/ask,
and VIX forwards worsen substantially. This is a persistent sampled
improvement and a restricted-family trade-off, **not a qualified joint fit**.
Ten scoped Lean theorems and 170 Python tests pass; numerical accuracy and
physical/model-comparison qualification remain unresolved.

## What is being fitted

This checkpoint implements the user's decision to attempt a joint fit now,
with numerical qualification following candidate selection. It supersedes
the previously proposed nine-point/two-level/24-SPX profile, which was never
launched. It does not start another variance-reduction method search.

Optimize **all three levels of the declared nuisance variance curve**, with
knots 0.25/0.75/5.832876712328767 years and bounds [0.00010001,0.25] per level.
Keep Stage-3 candidate 18, readout/architecture, zero projected initial state
and Q scenario 0 fixed. There is no latent-state, structural or Q-risk-price
calibration in this run. A restricted-family failure cannot establish M2
infeasibility or physical incompatibility.

The searching objective contains **all 707 cutoff-3 SPX training quotes** at
29/57/211/302 days, **all 175 VIX training options** at 27/55/118/153 days, and
**all four VIX forwards**. The other **273 SPX training quotes** are priced as
wing diagnostics after selection and do not enter this provisional priority
objective. All 980 SPX training rows remain in scope. Exact membership/weights
are saved in `spx_membership.csv` and `objective_membership.csv`.

Use the existing squared distance to bid/ask intervals divided by frozen quote
scales, expiry-balanced within SPX and VIX. Both instrument weights are one;
the forward family has weight one and fixed scale 0.5 VIX points. The existing
curve penalty is `0.1*sum(diff(xi)^2)` and is reported separately. No model-SE
denominator or midpoint tie-breaker is added. The score is not the full-980-SPX
comparison score, and does not establish an ESN/M2 ranking.

The anchor is the old selected VIX-only curve
`[0.0542496112, 0.0600318990, 0.0642080964]`, retained as a frozen numerical
comparison point. Its old train/check reversal remains recorded; using it as
an anchor does not adopt it as a validated fit.

## Efficient repricing without changing the stock target

At fixed state/readout/Q, the exogenous OU path and
`a_k=G(eta_k)/M0(t_k)` do not depend on deterministic xi. For segment j(k),

`V_k(xi)=floor+(xi_j(k)-floor)*a_k`,

`log Y_T(xi)=sum_k[-V_k(xi)*dt_k/2+sqrt(V_k(xi))*dW_k]`.

The stock is repriced for every candidate. Only a_k and the primitive equity
increments are cached. Analytical pathwise derivatives use

`d log(Y_T)/d xi_j = sum_{k:j(k)=j} a_k*(-dt_k/2+dW_k/(2*sqrt(V_k)))`.

This targets the same adapted exponential-Euler discrete stock with GH48
normalization and 512 steps per trading day. It does not make that scheme
continuous-time exact. Chunked vector summation is checked against the existing
pricer at every saved anchor batch, not claimed bitwise identical.

SPX sampling uses the original horizon-specific vol7 proposal. To keep fitting
all strikes economical, the fit freezes **anchor** stock/reference/two adapted
hedge controls plus mass. Each target strike gets the nearest same-side
representative's two hedge gains. The existing ridge5 penalty remains one,
with 256 independent anchor-pilot paths per expiry. The anchor controls and
coefficients stay fixed across curves; they are not recentered on search or
checking data. This preserves the target expectation under a fixed primitive
law, but does not assert equal variance to per-strike production controls.
The original production reference remains unchanged.

VIX uses the unchanged targeted Gaussian mixture and pooled-28-pilot Gaussian
control coefficients from `m2_vix_gaussian_stability_001`. Fresh full states
generate fresh segment exposures. Conditional GH32, normalizer GH48 and time
order 8 stay fixed. Anchor control coefficients are frozen across curve
candidates; known Gaussian means do not force VIX-squared to its model target.

`ExploratoryJointLoss` extends the same interval-distance formula to finite
**signed Monte Carlo estimates**, recording negative values. The production
`JointPriceObjective` still rejects negative model prices. This avoids turning
negative numerical estimates into either acceptable prices or silent search
failures. It is an exploratory signed score, not price qualification. On
nonnegative inputs all family/penalty components match production in tests.

## Bounded search and independent checking

Run: `runs/m2_joint_curve_fit_002`. Searching uses 512 SPX paths per expiry,
common across all candidate curves, and eight independently scrambled
64-point VIX nets per expiry. There are eight declared optimizer starts,
at most 25 function evaluations per start, and a 60-minute searching stop.
Analytical Jacobians avoid repeated finite-difference stock simulations.
Every start/evaluation and signed failure is retained. These budgets are
below the plan's 4096-path screening minimum and carry unresolved status.

Selection is made only on searching data and written to `selected.json`
before generating any fresh checking paths. Three checking trials each use
512 new SPX paths and eight new 64-point VIX nets per expiry. The selected
curve is never reselected using those checks. They reuse **training market
contracts**, so they are numerical checks, not held-out market performance.

SPX SE uses iid proposal paths conditional on the frozen pilot. VIX SE uses
independent complete scrambles, not individual Sobol points. Paired price and
loss differences compare selected and anchor on the same checking draws.
Eight paired blocks per trial, or 24 pooled, give a descriptive leave-one-block
jackknife for the nonlinear interval loss; it is not a calibrated tail-error
certificate. All raw price/path/block data needed for score reconstruction
are retained.

## Searching result

All eight starts terminate successfully near the same searching minimum.
The selected annual variance levels are

`[0.0600941554, 0.0748505719, 0.0383634764]`.

Searching loss falls from **637.3158** to **333.6764**. This was an apparent
search improvement before the fresh checks below. The first two levels
rise and the last falls. The terminal level affects only the final portion
of the 302-day SPX horizon; it does not affect the VIX windows. No boundary
is hit. Agreement of optimizer starts is fixed-sample numerical evidence,
not proof of population identification or a global M2 optimum.

Saved-price replay reproduces all **882 selected searching option estimates**.
Two poor starting points contain **4 and 7 negative objective estimates**;
their exact signed quote vectors are retained in
`search_signed_evaluation_013.csv` and `search_signed_evaluation_090.csv`.
None is clipped. The selected searching estimates are nonnegative.

The searching residual Jacobian has singular values approximately
1406.43/782.52/47.81, or 88.49/55.89/1.84 after scaling each column by its
selected curve level. These describe the conditional numerical curve-only
problem, not statistical/full-model identification. Selected-curve stock
derivatives on 128 saved paths at every expiry agree with finite differences;
the maximum absolute discrepancy is **1.36e-9**.

## Fresh-check results and disposition

The run takes **2015.24 seconds (33.6 minutes)**, including path/exposure setup,
optimization and all three fresh checks. The optimizer itself takes **131.67
seconds**, with **104 evaluated curves** across eight successful starts.
Thus the costly stage was pricing/data setup, not the three-parameter search.
No searching stop or numerical accuracy threshold was relaxed to finish.

| Sample | Anchor loss | Selected loss | Selected minus anchor | Descriptive paired block-jackknife SE |
| --- | ---: | ---: | ---: | ---: |
| Searching | 637.3158 | 333.6764 | -303.6394 | 51.2652 |
| Fresh trial 0 | 565.2016 | 369.0914 | -196.1102 | 90.4198 |
| Fresh trial 1 | 644.1134 | 363.7140 | -280.3994 | 82.5129 |
| Fresh trial 2 | 662.7354 | 416.3827 | -246.3527 | 51.7908 |
| Pooled fresh prices | **597.9019** | **347.0984** | **-250.8035** | **40.5001** |

The pooled score applies the loss to the pooled **price estimates**; it is not
the average of the three losses. The leave-one-block-out pooled difference
ranges from -263.24 to -232.27. Every fresh trial and every observed deletion
favors the selected curve. Unlike the old VIX-only fit, its sampled total-loss
improvement does not reverse. However, tail coverage and deterministic bias
are unresolved, so these diagnostic SEs are not confidence certificates or
proof of a validated market improvement.

### The gain is a compromise across markets

| Pooled fresh component | Anchor | Selected |
| --- | ---: | ---: |
| SPX option loss | 529.1980 | 249.0437 |
| VIX option loss | 67.8792 | 84.9146 |
| VIX forward loss | 0.8247 | 13.1399 |
| Curve penalty | 0.0000051 | 0.0001549 |

SPX 57/211/302-day losses improve; the 29-day loss slightly worsens. VIX
55/153-day option blocks improve, while 27/118-day blocks worsen. Both
instrument-family weights remain one. The optimizer sacrifices VIX to reduce
the much larger SPX discrepancy under the declared fixed-scale objective;
no weights were retuned after inspecting checking results.

| VIX days | Market forward | Anchor pooled forward +/- SE | Selected pooled forward +/- SE | Selected error |
| --- | ---: | ---: | ---: | ---: |
| 27 | 20.10 | 20.7177 +/- 0.0886 | 21.8035 +/- 0.0962 | +1.7035 |
| 55 | 20.73 | 20.3868 +/- 0.0605 | 21.4570 +/- 0.0658 | +0.7270 |
| 118 | 21.19 | 21.1940 +/- 0.0818 | 23.6686 +/- 0.0976 | +2.4786 |
| 153 | 21.72 | 21.1496 +/- 0.1076 | 23.6084 +/- 0.1259 | +1.8884 |

These are large forward misses relative to the *reported sampling SEs*,
especially at 118 days. They are useful fitting evidence even though the
SEs do not include quadrature bias or guarantee unseen-tail coverage.

Pooled interval coverage improves only from **41/882 to 50/882**:
SPX **23/707 to 29/707**, VIX **18/175 to 21/175**. There are no negative
primary anchor/selected estimates in searching or the three fresh checks.
This positivity observation does not rescue the many outside-interval prices.
Secondary SPX wings still contain **10/7/13 negative selected estimates** in
fresh trials 0/1/2 (anchor: 10/9/13). These are saved in the full-training
diagnostic tables without clipping and do not enter the priority objective.

### Sampling and identity warnings remain

Only **27/882 selected pooled nominal SE screens** pass, all VIX; **0/707
SPX** passes. Median selected SE/scale is **2.3853**, against the target 0.125.
All four selected pooled forward SEs remain above 0.0625. These are nominal
SE comparisons only, with no extra claim about exercise hits or tail clearance.
Concentration and the complete signed data remain in the saved reports.

The numerical identity audit retains **four stock/mass monitor warning rows**
beyond two empirical SEs; these are dependent diagnostics, not four independent
tests or proof of a law bug. All four concern the shared 211-day trial-2 sample:
mass is 0.9545 +/- 0.0224 for both candidates; weighted stock is about
0.9514–0.9515 +/- 0.0234–0.0236. The pooled controlled VIX-squared estimates have
no two-SE warnings at this budget. That does not certify the conditional
Gaussian tower, normalization or nonlinear option prices.

Necessary market-option second-moment margins are rechecked after changing
the curve. At the pooled sampled model mean, the selected 153-day margin is
still **-4.8424 squared VIX points**; mean uncertainty and pricing bias prevent
a certified incompatibility verdict. The anchor retains **-26.7946/-46.2175**
at 118/153 days. The selected 118-day margin is positive (+16.5194). No
identity is imposed and no price is rescaled to eliminate these warnings.

### Concrete continuation decision

Retain this curve as a **provisional joint candidate**, not an adopted fit.
Stop the standalone variance-reduction feature search. The next fitting stage
should test the declared **restricted C/Y/F initial-state coordinates**, with
the existing prior penalty and Q scenario 0, profiling curve/state effects
without freeing readout and risk prices simultaneously. This tests whether
the fixed zero-state distribution is responsible for part of the SPX/VIX
compromise. It is a hypothesis, not a demonstrated remedy.

Changing state requires updated dated normalization and candidate state paths;
the current G/M0/exposure caches cannot be carried over unchanged. Preserve
common primitive randomness where mathematically valid, use a bounded fit,
and freshly reprice finalists. Keep the already identified stock-grid,
normalizer and VIX quadrature errors separate. A curve-plus-state failure
would still be a restricted-family/numerically unresolved result, not proof
of M2 infeasibility. Physical checks, reserved market evaluation and the
common-protocol ESN comparison remain future qualification work.

## What Lean checked

The [correspondence record](formalization/CORRESPONDENCE.md) maps ten theorems
to numerical roles, assumptions and remaining representation bridges. The
[proof source](formalization/CalibrationChecks.lean) compiles under pinned
Lean 4.35.0-rc4 / Mathlib `1f414401f69059aa7eead47b53ee40bd38455eeb` with warnings
treated as errors. The final axiom audit contains only `propext`,
`Classical.choice` and `Quot.sound`, with **no `sorryAx` or custom axioms**.

The checks establish, within their stated exact-real scope:

- A Gaussian exponential stock step has the current stock as its mean, and a
  fixed start-of-step hedge delta has mean-zero gain.
- Integrable known-zero anchor controls remain valid for another integrable
  candidate payoff under the **same measure** with fixed coefficients.
- A finite change-of-measure identity places likelihood weights on both
  payoffs and controls. The continuum Gaussian-mixture representation remains
  a separate source/oracle obligation.
- Cached integrated VIX variance is affine in curve levels and positive under
  the stated positive-floor/nonnegative-exposure conditions. VIX itself is
  nonlinear, and exact conditional Gaussian integration is not certified.
- Weighted residual and curve-penalty vectors have the claimed squared loss.
- Mean-zero numerical noise adds a variance contribution to squared error
  outside a fixed interval edge in a finite model. An optimizer's searching
  loss therefore need not rank the underlying expected prices correctly.

These are genuine scoped proofs, not a formal verification of the entire
Python simulator, the Monte Carlo sample, quadrature accuracy or the fitted
market parameters. The correspondence record explains the ordinary
integrability/adaptation argument and the assumptions still outside Lean.

All **170 Python tests pass**: 124 package and 46 study-driver tests. Four new
tests cover stock-cache/original-pricer correspondence, stock/VIX/loss
derivatives, signed-loss/production correspondence and market-only hedge
assignment. Finite-state/data replay checks complement these tests.

## Separate unresolved errors

Sampling and bias are not combined into one standard error. The new fit uses
the existing 512/day stock scheme and GH48/conditional32/time8 rules; it does
not perform new stock-grid or quadrature refinements. Historical matched-grid,
normalizer and VIX refinements remain configuration-limited evidence. The
changed curve requires its own eventual accuracy assessment. Rare-tail
concentration, signed estimates, moment/forward checks and necessary market
moment margins remain visible. No normalizer reset or VIX moment rescaling
is applied.

## Reproduction and retained failures

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
PYTHONPATH=matteo/exogenous_reservoir_eq_m2/src \
/private/tmp/pdv-options-venv/bin/python \
  -m matteo.model_comparison.option_pricing.joint_curve_fit \
  --output /private/tmp/NEW_JOINT_CURVE_FIT
```

Use `audit_joint_curve_fit`, `replay_joint_curve_cache` and
`audit_joint_fit_identities` with `--run PATH` for independent score
reconstruction, same-path derivative/negative-start replay, and numerical
stock/VIX identity diagnostics respectively. The identity audit follows the
score audit because it reads pooled forward reports.

Large fixed-shock repricing caches live in the unique `/private/tmp/` directory
recorded by `design.json`. They are acceleration/replay artifacts, not a
new data source. Source snapshots, seeds, coefficients, exposures, signed
prices and comparison paths/blocks remain under the run directory. If the
temporary cache disappears, regenerate in a new output directory using the
frozen source and seeds; do not claim the cache replay passed without it.

`m2_joint_curve_fit_001` is retained as an interrupted startup. A new test found
that the exploratory score report omitted its zero midpoint component; the
numerical loss was unchanged. The attempt was stopped before optimization,
the reporting component was corrected, all four new tests then passed, and
the complete run uses new directory 002. Partial pilot draws in 001 are not
fresh checking evidence. No prior unsuccessful experiment is overwritten.

## Undergraduate interpretation

A joint optimizer chooses one curve for both markets. It may reduce the total
error by improving a large SPX discrepancy while worsening VIX. That is a
useful diagnosis of this restricted family, not a satisfactory joint fit.
Fresh paths ask whether the apparent gain survives different simulated
scenarios. Lean checks the mathematical accounting behind the construction;
it cannot supply rare scenarios that the simulations missed or remove a
finite-grid stock error.
