# Frozen SPX multilevel Monte Carlo prototype — 22 October 2026

Following the user-authorized SPX investigation, this experiment targets time
corrections and cost allocation after
[volatility-path sampling](SPX_VOLATILITY_SAMPLING_RESULTS.md). The model,
selected nuisance variance curve and zero projected initial state remain frozen.
The prior design-selected vol7 proposal is retained without reselection.
All four VIX maturities, 175 training options and four forwards remain active;
this experiment changes no SPX maturity/split scope and uses no held-out quotes.

## Estimator and its target

For normalized payoff P_h from the adapted stock scheme with h steps/day,

E[P_512] = E[P_64] + E[P_128-P_64] + E[P_256-P_128] + E[P_512-P_256].

Each correction uses the same exact fine OU path and summed equity increments
for its coarse stock simulation. The four level sample sets are mutually
independent. Each level is generated **directly at its own fine resolution**;
the base does not pay for a 512-step state cache. Curve jumps, expiries and
7/21-day events are exact grid points.

The vol7 Gaussian path basis is rebuilt on each level grid. Therefore proposal
laws differ across levels, but their exact full-path target/proposal weights
recover the same target Gaussian coarse marginal. This makes the telescope
valid in expectation. Independently drawn level paths do not telescope
sample by sample. Weights are not self-normalized.

Each level fits target-zero controls on 512 independent pilot paths. The base
uses the weighted stock, Black-reference and two adapted gains plus mass control.
Corrections use the **difference** of fine/coarse controls, with the same path
weight, plus the target-zero mass control. The reference-payoff column cancels.
Separate level coefficients are valid because every such control has zero
target expectation. Three frozen coefficient fits use all 512 paths and two
disjoint 256-path halves. No half is selected after checking.

The direct comparator estimates the same 512-step payoff using its own
independent checking paths and coefficients trained on the fine payoff of the
top-level training pilot. It does not price an unnecessary coarse path during
checking. Pilot sharing is explicit; checking uncertainty is conditional on all
frozen fits.

If X_l denotes the controlled contribution at level l, the reported estimate
and sampling variance are

mu_hat = sum_l mean(X_l),    SE² = sum_l sample_variance(X_l)/N_l.

The variance formula requires independent level sets. It estimates sampling
uncertainty for the **discretized 512-step model**, not continuous-time bias,
normalizer error or unresolved rare tails.

## Cost allocation and independence

Separate 256-path allocation pilots evaluate the frozen coefficients and
estimate costs/variances. The design objective is the average sampling variance
divided by quote_scale² over fixed central contracts **2/3/8/9**. It is not a
market calibration objective and does not optimize every tail option.

For variance proxies V_l and measured costs C_l, the unconstrained continuous
solution has N_l proportional to sqrt(V_l/C_l). An active-set step enforces a
64-path minimum per level; counts round to 64-path blocks, then affordable
blocks are added by marginal variance reduction per cost. This is a bounded
integer approximation, not a guarantee of globally optimal integer allocation.
Zero pilot variance still receives its minimum sample count and is never
accepted as proof of accurate tail pricing.

The checking budget is the pilot-predicted cost of **1,024 direct 512-step
paths per trial**. Counts, coefficients and predicted costs are frozen before
any checking draw. Two fresh checking trials use that same frozen design.
Actual cost matching must be assessed from the recorded checking times;
pilot timing noise and different batch sizes can change realized costs.

Deterministic GH48 normalizers are warmed separately on all four grids.
Setup and pilot overheads are excluded from the checking budget and reported
separately; this does not assume they are free or already amortized. The direct
cost proxy uses fine-state/reference/fine-price timing from top-level pilots;
its ancillary weighting/bookkeeping costs are measured in actual checking.

## Validation and reproducibility

**138 tests pass: 106 package and 32 comparison-driver tests.** New tests check
allocation against an enumerated tractable optimum, minimum/zero-variance and
invalid-budget handling, unequal independent-level variance addition, and an
end-to-end constant-volatility Black oracle. In that oracle coupled timestep
payoff corrections vanish to floating-point precision. These are mathematical
correspondence and numerical tests, not Lean formal verification or a general
convergence theorem.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.spx_mlmc --output matteo/model_comparison/option_pricing/runs/m2_spx_mlmc_001
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.verify_spx_mlmc matteo/model_comparison/option_pricing/runs/m2_spx_mlmc_001
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.audit_spx_mlmc matteo/model_comparison/option_pricing/runs/m2_spx_mlmc_001
```

Captured sources include the model package, runner dependencies, curve-fit
input and full scope. Independent seeds identify training, allocation, each
checking level, direct checking and trials. Full raw/controlled level draws,
weights, Gaussian scores, fine/coarse hit counts and costs are retained.
The verifier reconstructs likelihoods, coefficients, allocation and estimator
means/SEs. Signed corrections are not mistaken for negative option prices;
negative final estimates and zero-hit levels are preserved without clipping.

## Results: this hierarchy is not an efficiency improvement

The [completed run](runs/m2_spx_mlmc_001/run.json) took **383.1 seconds**:
196.3 seconds of setup and 54.0 seconds of recorded design-path computation,
with remaining time including checking and artifact work. Grids contain
2,521/5,042/10,084/20,168 intervals, with max(lambda*dt) approximately
0.30982/0.15491/0.07746/0.03873. These are trading-day resolutions; contract
year fractions remain calendar ACT/365.

The frozen [allocation](runs/m2_spx_mlmc_001/allocation.json) is
**[2752, 640, 384, 192]** paths for the base and three corrections.
Predicted checking cost is 31.06 seconds versus a 31.22-second direct budget.
The [audit](runs/m2_spx_mlmc_001/audit.json) records the realized comparison:

| Pooled-fit checking diagnostic | Trial 0 | Trial 1 |
|---|---:|---:|
| Actual MLMC / direct seconds | 32.48 / 31.68 | 37.42 / 31.14 |
| Actual MLMC / direct cost ratio | 1.0253 | 1.2017 |
| Median direct/MLMC estimator-variance ratio over 12 contracts | 0.8550 | 0.7300 |
| Contracts with lower MLMC estimator variance | 5/12 | 4/12 |
| Direct/MLMC mean scaled variance on the four central contracts | **0.6223** | **0.6598** |
| Central variance-cost efficiency gain, direct/MLMC | **0.6069** | **0.5491** |

For the predeclared central objective, values below one favor direct simulation.
Thus this multilevel design is **worse in both fresh checking trials**, even
before charging its extra design cost. Trial 1 also overspends the predicted
checking cost because realized timings differ from the pilot. The budget was
matched by prediction, not exactly by realized runtime; cost-adjusted gains
above account for that difference. No checking-based reallocation is applied.

Both methods pass **0/24 nominal quote-relative sampling screens**, with no
negative pooled final price estimates. Separate cross-contract medians of
SE/scale (2.3985 MLMC versus 2.4918 direct) must not be mistaken for an efficiency
win: they do not represent the contract-matched ratios or the predeclared
central mean-variance objective. All raw and pilot-half results are retained.

### Why the naive multilevel hierarchy is expensive

The allocation-pilot central variance proxies for the corrections are
**638.70, 381.46 and 217.12**, while their measured per-path costs are
**0.00944, 0.01889 and 0.03815 seconds**. Correction variance drops about
40–43% per refinement while cost roughly doubles. The variance-cost terms
therefore do not become negligible on the fine levels. These are empirical
proxies from small pilots, not an asymptotic convergence-order theorem.

An explanatory replay of the existing design data gives predicted direct
central estimator variance **3.4457**, versus **4.6696** for MLMC: a direct/MLMC
ratio **0.7379**, already favoring direct. The experiment froze the allocation
without selecting between estimators; fresh checks confirm the negative design
signal rather than supplying grounds to scale this hierarchy up.

The correction controls themselves do help: on the central checking contracts,
aggregate uncontrolled/controlled path-variance ratios are approximately
**6.82, 10.35 and 8.03** for the three corrections (26.08 for the base).
This compares controlled and raw contributions on identical paths. It does not
mean a 6–10-fold improvement over the previous controlled pricer. After those
reductions, the remaining correction variance is still too expensive for this
hierarchy and budget. More independent level samples alone is not the next
justified change.

### Retained consistency, tail and pilot warnings

Two trial-0 MLMC/direct price differences exceed two independent combined SEs:

- **29-day 3400 put:** MLMC 1.51204 ± 0.39403 versus direct
  0.13051 ± 0.36820; difference +1.38152 ± 0.53928 (2.56 SE).
- **29-day 7030 call:** MLMC 47.70949 ± 1.08386 versus direct
  51.38283 ± 0.99150; difference -3.67334 ± 1.46895 (-2.50 SE).

They are finite, multiple-comparison warnings, not proofs that two estimators
have different population targets. Trial 1 has no differences above two SEs.
No estimator, fit or trial is removed to recover agreement.

The top correction has only 192 paths. In trial 0 its mean importance weight
is **1.09478 ± 0.03365** (2.82 SE above one), with correlated stock/reference
warnings. Trial 1 mass is **1.07058 ± 0.03476** (2.03 SE). A direct trial-1
reference-forward warning also remains. No likelihood weights are rescaled to
force these identities. Exact Gaussian weight/transition oracle tests pass;
these finite warnings remain unresolved checking evidence.

Six pooled level/contract rows have zero observed fine and coarse payoffs:
the 29-day 3400 put at levels 2 and 3, and the 29-day 7900/7950 calls,
57-day 3450 put and 8800 call at level 3 (all trial 0). Controlled corrections
can nevertheless be nonzero because controls have nonzero sample means.
Their raw zero variance and small controlled SEs do not establish negligible
population tail corrections. They weaken the apparent finest-level check.

Pilot halves remain unstable in some tails: relative to pooled fits, maximum
SE inflation is **4.00** for direct half-A and **2.53** for MLMC half-B.
All final full/half estimates in this run are nonnegative; signed level
corrections can legitimately be negative. No half is promoted after checking.

### 512 steps/day still does not establish convergence

For the pooled 256->512 correction, observed changes exceed scale/8 in
**8/12** contracts in trial 0 and **9/12** in trial 1. Median absolute changes
are **0.2333** and **0.8359** quote scales. In particular:

- Trial 0 **57-day 6775 put:** **-2.57035 ± 0.87963**, scale 0.65303.
- Trial 1 **57-day 5550 put:** **+1.69703 ± 0.69914**, scale 0.20024.

Each exceeds two observed correction SEs. The signs concern different contracts
and trials; they are not a uniform bias direction. Sparse fine-level paths and
zero-hit tails also prevent treating small other corrections as convergence.
These checks do not bound the remaining 512->continuous-model error. GH48
normalization is unchanged and its population error remains a separate gate.

## Verification and next direction

The artifact verifier reconstructs the frozen allocation, likelihood weights,
full/half pilot coefficients, **192 final estimates, 384 level statistics and
50 expectation monitors**, including the independent-level variance sum. All
138 tests pass. No model refit, maturity exclusion, forced normalization or
price clipping has been introduced.

Keep this MLMC implementation as a verified exploratory comparator; **do not
replace the direct vol7 pricer with it**. The next useful design question is
whether a better coarse/fine Gaussian coupling can reduce correction noise
before another allocation study. A candidate is conditional Gaussian bridge
reflection that preserves the coarse equity/OU endpoints; its target-law
invariance, adapted stock update and importance-weight accounting must be
validated before a market experiment. This remains a proposal, not implemented
or established to improve variance. Any alternative stock integration rule must preserve adaptedness and the
martingale law and needs its own correspondence checks. Pilot stability and
rare tails remain open too.
