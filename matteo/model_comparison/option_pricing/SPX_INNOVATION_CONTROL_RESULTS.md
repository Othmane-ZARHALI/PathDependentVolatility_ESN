# SPX OU-innovation martingale controls — 22 October 2026

This tests the next hypothesis recorded after the
[Gaussian reflection study](SPX_BRIDGE_COUPLING_RESULTS.md): use cheaper
volatility-factor martingale controls rather than simulate an additional fine
path. Model parameters, dated initial state, selected nuisance curve, vol7
proposal and GH48 normalizers remain fixed. All VIX maturities
27/55/118/153 days, 175 training options and four forwards remain active;
SPX maturity/split scope is unchanged. Held-out quotes are unused.

## Target-zero controls and full-state correspondence

The exact target-Q state transition is

X_next = D(dt) X + b_Q(dt) + eta,
E_Q[eta | past] = 0,    Cov_Q(eta | past) = Q(dt).

For projection matrix P selecting (C,Y,F), define

eps = P eta,    R = P Q P',    a = Cov_Q(eps,dW | past).

Then eps_i, eps_i²-R_ii and eps_i*dW-a_i all have conditional mean zero.
The last centering matters because equity/OU shocks are correlated. R is the
**full OU innovation covariance**, not its covariance conditional on dW.
With predictable coefficients, sums of these centered quantities are target
martingales. Under the path proposal q the complete gain is multiplied by
w=f/q; no proposal-drift recentering or sample-mass normalization occurs.

`InnovationRecorder` observes read-only full states during the existing exact
OU simulation. It computes

P X_next - P [D X + b_Q]

at each requested nested-grid boundary, retaining only projected innovations
and covariance arrays. Full states are retained only for currently open
intervals. **The three projections alone are not treated as Markov:** states
with identical (C,Y,F) can have different next projection means when mode rates
differ. Centering uses the target Q risk-price drift, including under proposal
sampling. The ordinary state simulation and its random stream are unchanged.

The predictable coefficient is Black-inspired log-variance vega

v_t = (1/2) S_t sqrt(W(t,T)) phi(d1),

using the current normalized stock and deterministic remaining curve variance.
The three linear and mixed gains multiply this by the start-of-step gradient
of log G; the three square gains use half the diagonal Hessian G/G. There are
**nine new gains**. The readout derivative predictors concern the smooth exponential formula
before numerical exponent guards; their target-zero property does not require
them to equal the derivatives of a guard-clipped readout. These are exposure proxies, not exact stochastic-volatility
option Greeks. In particular, the variance floor and future mean reversion
are not claimed to be fully reflected in those predictors. The zero-mean
argument holds for predictable integrable coefficients regardless of proxy
accuracy.

Each stock step retains the existing adapted exponential-Euler update. Its
payoff and existing stock-hedging controls match the previous implementation;
controls change the estimator, not the stock model or time discretization.

## Independent fits and checking experiment

Twelve predeclared training representatives at 29/57 days are simulated on one
512-step/day state grid, with coupled 128/256 coarsenings. Separate conditional
innovations are recorded from full-state means on **each** grid. This supports
three direct-price targets and 128->256 / 256->512 correction targets.
Correction controls subtract the two grid gains; the target-zero mass control
remains w-1, rather than subtracting two identical mass columns into zero.

A 512-path independent pilot fits every target's coefficients. The primary
comparison is **baseline ridge5 versus augmented ridge14**, both using penalty
one on centered, pilot-standardized columns. Five baseline controls are the
weighted stock, Black-reference and two adapted stock gains plus mass control.
The additional nine weighted OU gains make 14. The original OLS5 estimator is
retained separately to distinguish regularization effects. Disjoint 256-path
pilot halves provide sensitivity fits for both ridge dictionaries. No checking
result selects a dictionary, regularizer, fit or proposal.

Two fresh 512-path checking trials share paths across dictionaries, grids and
targets, conditional on frozen pilot fits. Standard errors use iid paths;
method differences are paired. The ordinary sampler/pricer is rerun with the
same seeds and alternating execution order for actual timing and numerical
correspondence. The augmented timing includes full-state innovation recording,
all three grid prices/controls and innovation diagnostics. It is a **joint
three-grid diagnostic cost**, not standalone fine-price cost. Setup/quadrature
and pilot overheads are separate. No new MLMC allocation is constructed.

## Validation and reproducibility

**146 tests pass: 112 package and 34 comparison-driver tests.** New independent
checks cover full-state conditional means at fine/coarse events; Gaussian
innovation covariances and equity cross-covariances; null-projection states
with distinct conditional means; sampler backward correspondence; smooth
readout derivatives against finite differences; target/weighted-proposal
martingale gains; independent Black payoffs; identical stock/legacy hedges;
and end-to-end ridge/correction mass bookkeeping. Constant readout produces
exactly zero new gains. These are algebraic/implementation/numerical checks,
not Lean formal verification or price convergence certification.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.spx_innovation_study --output matteo/model_comparison/option_pricing/runs/m2_spx_innovation_controls_001
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.verify_spx_innovation_study matteo/model_comparison/option_pricing/runs/m2_spx_innovation_controls_001
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.audit_spx_innovation_study matteo/model_comparison/option_pricing/runs/m2_spx_innovation_controls_001
```

Sources, scope, curve input and model/target hashes are frozen before execution.
The pilot and checking archives retain weighted payoffs, all gain dictionaries,
Gaussian scores, exact weights, per-path compensated innovation sums and all
estimator draws. Deterministic normalizer arrays/times are saved as numerical
artifacts, without declaring their population error resolved. The verifier
reconstructs dictionaries, coefficients, likelihoods, means/SEs and paired
method differences. Negative adjusted prices and signed corrections remain
unclipped. Variance reduction cannot remove stock-integration or normalization
bias; sparse payoff hits remain an independent limitation.

## Fresh results: useful direct-price variance reduction, nearly break-even cost

The [run](runs/m2_spx_innovation_controls_001/run.json) takes **396.1 seconds**,
including **174.1 seconds** of setup. The
[audit](runs/m2_spx_innovation_controls_001/audit.json) preserves all comparisons.
The replay verifies **960 estimates, 120 paired method differences, 64 bare
innovation/forward monitors and 1,680 contract-specific control means**.

The primary ridge14/ridge5 comparison uses the fixed central contracts 2/3/8/9:

| Target | Trial 0 variance gain | Trial 1 variance gain | Trial 0 joint cost-efficiency gain | Trial 1 joint cost-efficiency gain |
|---|---:|---:|---:|---:|
| Direct 128-step price | 1.8974 | 1.7851 | 1.0611 | 1.0082 |
| Direct 256-step price | 1.8743 | 1.8125 | 1.0482 | 1.0236 |
| Direct 512-step price | **1.8652** | **1.8421** | **1.0431** | **1.0404** |
| 128->256 correction | 1.2028 | 1.1304 | 0.6726 | 0.6384 |
| 256->512 correction | 1.2175 | 1.0956 | 0.6808 | 0.6188 |

A gain above one favors the new dictionary. The direct-price variance reduction
repeats on both fresh trials, but **measured runtime nearly absorbs it**:
ordinary/augmented joint collection costs are 26.69/47.72 and 26.58/47.07 seconds,
so the runtime multipliers are **1.7882/1.7707**. A roughly 4% measured efficiency
gain at 512 steps is small and not a robust speedup certificate. These are
costs for the whole shared three-grid diagnostic, including innovation recording
and diagnostics; they do not establish standalone fine-price efficiency.

At 512 steps, observed variance improves in **9/12 and 7/12** contracts, with
all-contract median gains **1.5473/1.4957**. The central aggregate gain and the
all-contract median describe different objectives. At 128->256 only 9/12 and
4/12 improve, and at 256->512 only 7/12 and 4/12 improve. **None of the four
correction cost comparisons supports replacing baseline correction controls.**
No dictionary is selected retrospectively from checking results.

All three primary price estimators (legacy OLS5, ridge5, ridge14) pass **0/24
positive-hit sampling screens at each grid**. At 512 steps, median observed
SE/scale falls from **2.8458 to 2.0331**, still far above the 0.125 screen.
Nonzero payoff counts can be as low as **one**. A finite sample-variance gain
therefore does not establish accurate far-out prices or convergence.

## Warnings that limit adoption

### Repeated far-call paired mean shifts

The complete primary paired audit has **17 related flags above two SEs**.
Several far-call shifts have the same sign across independent checking trials.
At 512 steps, augmented minus baseline ridge estimates differ by:

| Contract | Trial 0 difference ± paired SE | Trial 1 difference ± paired SE |
|---|---:|---:|
| 29-day 7900 call | **+0.14339 ± 0.03930** (3.65 SE) | **+0.16172 ± 0.03926** (4.12 SE) |
| 29-day 7950 call | **+0.11384 ± 0.03299** (3.45 SE) | **+0.11893 ± 0.03386** (3.51 SE) |

The two estimators have the same algebraic target if the stated martingale
conditions hold. These repeated finite warnings are nevertheless **not cleared
by citing that algebra or by dismissing them as ordinary independent multiple
comparisons**. Adjacent strikes and time grids share paths and are strongly
related. Dedicated control-tail checking remains necessary; the current run
does not establish the cause or population accuracy of these wing estimates.

All **64 bare compensated-innovation and forward** monitors are within two SEs.
That is not a test that every predictable contract-weighted gain is well
sampled. The post-hoc [contract control means](runs/m2_spx_innovation_controls_001/control_means.csv)
retain **126 correlated flags out of 1,680 means**, including existing controls:
reference-payoff means for those two calls are -3.39/-4.05 SE in trial 0;
existing deterministic-delta gains are -2.80/-2.92 SE in trial 1; new linear
cluster gains are -2.53/-2.50 SE in trial 1. Thus the warning is not isolated
to a new Gaussian covariance compensator. No control mean or price is forced
to agree with its known population identity.

A [roundoff audit](runs/m2_spx_innovation_controls_001/covariance_roundoff_audit.json)
composes the implemented fine-step OU root covariances onto each coarser grid.
Maximum projected covariance discrepancy from the target formula is
**2.78e-17**, maximum relative diagonal discrepancy **3.84e-16**. This check
finds only roundoff-scale compensation differences; it does not identify a
mean-centering implementation error or explain away the paired wing warnings.
The weighted Gaussian/Black and predictable-gain oracles pass. Rare-event
sampling, the distribution of control gains and coefficient sensitivity remain
open explanations to investigate, not an established diagnosis.

### Negative estimators and pilot instability

Preserved negative direct-price estimates include trial-0 baseline ridge
**29-day 3400 put at 128 steps: -3.69264 ± 6.62777**, and augmented
**57-day 3450 put at 256 steps: -2.12469 ± 2.59609**. These are noisy controlled
estimators of nonnegative payoffs, not an arbitrage result. No primary adjusted
512-step price is negative in this run; none is clipped at any grid.

Half-pilot sensitivity remains serious in some tails. Maximum SE inflation
relative to its pooled fit is **9.42** for augmented half-A at 128->256 and
**6.41** for augmented half-B at the direct 512-step price. Medians are much
closer to one; the maxima are retained rather than hidden by those medians.
No half-fit is chosen after evaluation.

### Timestep warnings survive the new controls

Augmented 256->512 changes exceed scale/8 in **16/24** cases, with median
absolute change **0.3276 quote scales**. In trial 0 the 57-day 6775 put changes
**-1.17242 ± 0.45488**, against scale **0.65303**. Trial 1 has a 29-day 7900-call
change **-0.00370 ± 0.00185**, with just **one nonzero raw correction**; that
thin-tail flag must not be treated as a reliable asymptotic significance test.
Both examples and all other corrections remain in the record. Controls can
reduce sampling noise; they have not removed 512-step discretization bias,
settled rare tails or certified GH48 normalization.

## Post-hoc pilot ablation: smaller dictionary is worth a fresh test

[innovation_design_audit.py](innovation_design_audit.py) also examines **only the
original training paths**, after the completed checking experiment. Two-fold
pilot analysis trains on one 256-path half and evaluates on the other, without
using checking data to refit. It compares baseline, linear-only, linear+square,
linear+mixed and all-nine dictionaries with the same ridge penalty.

For the direct 512-step central objective, the two pilot folds give variance
gains **1.6346/1.5209** for **three linear gains**, versus **1.6682/1.5558** for
all nine. Thus the three-column design captures most of the old-pilot benefit.
This is **post-hoc design evidence, not independent confirmation** of a new
method or a measured subset-cost advantage. No old checking result is rewritten
as evidence for a newly selected dictionary, and the original full-nine results
remain the primary completed experiment.

Reproduce the supplementary design/roundoff audit with:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.innovation_design_audit matteo/model_comparison/option_pricing/runs/m2_spx_innovation_controls_001
```

## Decision and next step

Retain the full-nine dictionary as an **experimental variance-reduction signal**,
not a default replacement or accuracy clearance. It improves central direct
variance, is nearly break-even after measured joint cost, loses on correction
costs and has unresolved far-wing/control-mean warnings. Direct vol7 with the
existing controls remains the reference; no model calibration is launched.

Next implement an actually cheaper **linear-only** control kernel, avoiding
unused quadratic/mixed calculations and profiling predictable-coefficient
computation. Freeze its design before a new independent checking sample.
Keep all representative strikes, audit contract-specific gain tails and paired
mean shifts, and measure standalone direct-price cost rather than assume a
saving from fewer columns. Check correction performance separately and keep
baseline corrections unless fresh cost evidence supports a change. This step
will not replace the still-open timestep, pilot-stability or normalization gates.
All VIX maturities and original scope remain included.
