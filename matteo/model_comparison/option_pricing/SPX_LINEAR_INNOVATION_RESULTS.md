# Three linear OU innovation controls: direct512 experiment

## Frozen experiment design

Declared before running `m2_spx_linear_innovation_001`: retain the Stage3
candidate18/scenario0 model, premium, curve, initial state, defensive volatility
path proposal, all scopes and 12 previously declared SPX training representatives.
All four VIX maturities remain active; this experiment runs SPX only.

Compare ridge8 (the five existing controls plus three linear OU gains) to ridge5,
both with ridge penalty 1. Retain legacy OLS5 and fits on each half of the pilot.
Use a fresh 512-path pilot and two independent 1024-path checking trials, batches
of 128, seed 20261030. Freeze coefficients before checking. No checking-based
control selection or parameter refitting. Both methods use identical weighted
payoffs, existing controls and original full-path likelihood weights.

Compute direct512 only, without coarser pricing/correction calculations. Alternate
baseline/augmented execution order between batches. Measured warm-normalizer
pipeline costs include path preparation, recording, reference calculation,
pricing, weights and bare-innovation monitors. Setup quadrature is recorded
separately. The linear branch skips Hessians, square and mixed noises and the
nine-column gain tensor; recording still uses the exact full-state OU transition
and retains its covariance metadata. This is an optimized three-control branch,
not an optimization of every common simulation operation.

Retain all contract-control mean checks, paired adjusted-mean warnings,
zero/low-hit cases, signed negative estimates, pilot sensitivity and SE/scale
screens. Previous repeated far-call mean warnings remain unresolved regardless
of whether this smaller sample happens to pass. GH48 and stock discretization
accuracy are separate open gates. No new MLMC or full-surface clearance is implied.

## Mathematical and implementation correspondence

At each step define `eps = P (X_next - D X_start - b_Q)` using the **full**
OU state. Each new gain is a sum of `vega_start * grad_i(log G)_start * eps_i`.
The multiplier is known at the start of the step and the innovation has zero
conditional target-Q mean. Under the proposal the full gain is multiplied by
`w = dQ/dProposal`. The final estimator subtracts pilot-frozen coefficients
times these weighted gains; no sample mean is subtracted from the controls.
This zero-mean argument does not certify finite-sample tail accuracy.

The source oracle verifies the linear branch is elementwise identical to the
first three full-nine gains on nonconstant-volatility paths and that stock
hedges and terminal stock paths are identical. Mocking the quadratic derivative
function to fail confirms the linear branch never calls it. Existing target
Gaussian, weighted-proposal and constant-volatility Black oracles remain active.
All 112 package tests and 35 study-driver tests pass. The previous full-nine
artifact is also replayed with the generalized verifier.

## Results

Run: `runs/m2_spx_linear_innovation_001`. It completed in 363.20 seconds,
including 116.25 seconds of setup. Coefficients were frozen before checking;
source snapshots and hashes are retained. The verifier replays 192 estimates,
24 paired differences, 16 bare monitors and 192 contract-control mean checks.

| Direct512 statistic | Trial 0 | Trial 1 |
| --- | ---: | ---: |
| Central normalized variance gain (ridge8 / ridge5 benefit) | 1.7953 | 1.7105 |
| Median contract variance gain | 1.4786 | 1.3821 |
| Contracts with lower variance | 10/12 | 11/12 |
| Baseline pipeline seconds | 39.44 | 38.51 |
| Linear-control pipeline seconds | 61.99 | 57.87 |
| Runtime multiplier | 1.5718 | 1.5030 |
| Central variance-cost efficiency gain | 1.1422 | 1.1381 |

The central objective is the mean variance divided by squared quote scale on
contracts 2,3,8,9; efficiency is variance gain divided by measured cost multiplier.
These are observed sample ratios from two trials, not confidence bounds or a
measured 14% wall-clock speedup at a guaranteed error tolerance. Pricing-kernel
seconds alone were 20.16/19.38 baseline versus 37.25/35.36 augmented; common path
simulation amortizes that overhead in the pipeline. Neither timing includes
pilot learning or quadrature setup. A full end-to-end accuracy comparison must
account for those fixed costs too. Do not compare these standalone timings
directly to the previous joint three-grid timings as if the workloads matched.
No same-path fresh nine-control comparison was run in this experiment.

## Warnings retained

- All three primary methods (legacy OLS5, ridge5, ridge8) pass **0/24** positive-hit
  `SE <= scale/8` screens. Median SE/scale falls from 2.1507 to 1.6617, still far
  above 0.125. The smallest raw payoff count is five; zero observed payoff cases
  do not arise in these trials, but sparse tails remain.
- Trial 0, 29-day 3400 put: ridge5 is **2.10347 +/- 0.36940**, ridge8 is
  **-4.59003 +/- 4.76450**. The negative estimate is retained. Path 162 alone
  contributes **-4.77923** to the augmented-minus-base mean and **93.99%** of its
  uncentered squared differences. This is a descriptive leverage audit with no
  path removal, clipping or refitting. It demonstrates tail sensitivity despite
  the central variance improvement.
- Half-pilot sensitivity remains severe: augmented-A SE inflation reaches
  **30.52x** on trial 0's 57-day 3450 put; augmented-B reaches 3.31x. The baseline
  itself reaches 10.32x on a far call. The pooled fit does not establish stable
  learning from a small pilot.
- One of 24 paired mean comparisons exceeds two empirical SEs: 29-day 5525 put,
  trial 0, **+0.36946 +/- 0.18002** (2.05 SE). The 29-day 7900/7950 call shifts are
  -0.10727 +/- 0.08400 and -0.09725 +/- 0.07071 in trial 0; -0.18000 +/- 0.11009
  and -0.12226 +/- 0.08191 in trial 1. These do not reproduce the previous
  full-nine positive 3–4 SE shifts. Changed dictionary, pilot and seeds prevent
  attributing that difference to a single cause; the earlier warning is **not
  cleared**.
- All 16 bare innovation/forward monitors fall within two empirical SEs, but
  **19/192 correlated contract-control means** exceed two SEs. Four are extremely
  small reference expectations with no observed reference payoff (means near
  -1.08e-29 and -1.86e-16): their huge standardized flags diagnose unobserved
  reference tails, not economically large discrepancies. Other flags include
  repeated stock means, adapted stock hedges and one linear-C gain. Inspect
  `control_means.csv`; the flags are not independent tests.

## Decision and next bounded step

Keep the linear branch available as an experimental option; keep the ordinary
vol7/stock-control estimator as reference. The central efficiency gain is modest
and does not offset the unresolved wing and pilot risks for general adoption.
Do not use these results to refit the market model. All VIX maturities remain.

Before increasing path counts or adding more control features, diagnose pilot
coefficient instability and contract-control tail leverage. A bounded candidate
is a pilot-only coefficient-stability rule (stronger shrinkage or omitting
unstable gains by a predeclared pilot criterion), assessed against ridge5 and
this fixed ridge8 on new checking seeds. First report coefficient/scaling and
largest-contribution diagnostics on saved pilot/checking data as exploratory
evidence; choose no new rule using checking performance. Any resulting rule
must be frozen before a fresh independent experiment. Keep paired timestep and
normalizer refinement as separate pending accuracy work.
