# ESN option pricing: code and report correspondence

Audit started 21 October 2026. Scope: `CrossSectionalData/`, relative to the
repository root. The original files are preserved. This is a source review
and an initial numerical diagnostic, not a reproduced optimization or an M2 fit.
The historical calibration audit lives separately in
[`../othmane_replication/OTHMANE_APPROACH_ISSUES.md`](../othmane_replication/OTHMANE_APPROACH_ISSUES.md).

## Recorded discrepancies

1. **VIX fitting range differs from the report's wing discussion.**
   `calibrate_smile.py`, `VIX_MONEYNESS_RANGE` and `load_vix_points()`, retain
   only `0.40 <= K/F <= 3.00`. The report's final VIX-smile discussion describes
   calibration and weights at moneyness 5–9. Increasing the weighting bandwidth
   cannot include quotes already removed by the loader. The runner records
   actual retained targets and raw rows outside the code's range. The figures
   may come from another code revision; this has not been established.

2. **Failed implied-volatility inversions disappear from the loss.**
   `calibrate_smile.py`, `residuals()`, sets NaN differences to zero;
   `diagnostics()` removes NaNs before computing errors. This is acknowledged
   in the report's final VIX figure caption, but limits any claim to fit all
   quoted strikes. The runner retains failures and their original weight,
   labels surviving-point RMSE explicitly, and does not optimize this loss.

3. **The plotted/labeled objective is not the executed optimizer objective.**
   The report writes a weighted IV sum of squares and describes a convergence
   plot as that sum divided by 858 points. Executed `residuals()` additionally
   divides by fixed reference-fit scales and `sqrt(number of expiry specs)`,
   applies expiry multipliers (`vix_near=3`, `vix_next=1.5`), and the entry point
   uses `least_squares(loss='soft_l1', f_scale=1)`. Multipliers in the residual
   are squared in the ordinary sum of squares. The history logger writes the
   raw residual sum of squares; the checkpoint after fitting stores twice the
   robust optimizer cost under the same `sumsq` key. These quantities must not
   be treated as interchangeable. Exact original scales/history are not in
   either supplied ZIP archive.

4. **Several headers describe an older calibration.**
   `calibrate_smile.py`'s opening documentation says seven parameters and four
   expiries; `unpack()` and `CombinedCalibrator.__init__()` implement eight
   parameters and six expiries. Old bandwidth descriptions also differ from
   current constants. The report's eight-parameter/six-expiry description
   agrees with the executable class on those two points.

5. **Printed final parameters cannot reproduce the exact optimizer output.**
   The report prints approximate final parameters. The entry point contains
   an exact reference `theta0`, not the printed final solution. Neither ZIP
   archive contains final parameter, scale, or history JSONs. The runner
   labels its `report-rounded` and `reference` presets distinctly. Agreement
   or disagreement at rounded parameters is not an exact reproduction claim.

6. **Standalone execution depends on missing original filesystem paths.**
   Several scripts hard-code `/mnt/user-data/uploads/...`; output paths expect
   sibling `results/` directories. The runner substitutes explicit local CSV
   paths in memory and writes only to a fresh output directory. It does not
   change numerical equations or target-processing rules.

7. **Archive and checkout are not wholly byte-identical.**
   All ten top-level Python files were compared with `esn_market_test_code.zip`.
   Nine match; `joint_calibration.py` differs in historical target paths and a
   leverage-target index correction. The option engine only imports its fixed
   architecture and parameter constants, which match. These differences do not
   explain the option code/report discrepancies above.

## Pricing conventions and limitations to keep separate from bugs

- The daily Gaussian log-return correction makes the discrete stock ratio a
  martingale under the simulator's chosen law. It does not identify that law
  with an empirically estimated change from the historical measure.
- Burn-in averages over an ensemble of possible current states. It does not
  condition prices on a state filtered for 19 February 2026.
- VIX strikes are evaluated at market moneyness times the **model's** VIX
  forward. A shape fit therefore does not enforce market forward levels or
  prices at the actual absolute strikes. The report's model means around
  15–16 versus market forwards around 20–22 remain material level differences.
- The conditional VIX uses a finite SPX strike range (0.60–1.60) and finite
  inner simulation. Wing truncation, integration error, inner bias and outer
  sampling error must be measured separately before comparing model families.
  A deterministic 20% constant-volatility oracle at 30/365 years gives variance
  0.0383673 on the calibration's 29-point grid versus the exact 0.04; refining
  to 284 points gives 0.0399595. This isolates a real grid effect without any
  Monte Carlo noise, but is not a bias estimate for the calibrated ESN.
- Maturities are rounded to trading days for diffusion, but exact calendar-day
  fractions enter Black inversion. The runner records both clocks.
- The original validation script compares independent samples. The new check
  also compares scalar and batched single-path outputs with the same seed,
  across smoothing factors and bias offsets, and checks conditional restart.
- Put/call mid-IV pooling is reproduced for the legacy benchmark only. Our
  common comparison will use the existing audited bid/ask price pipeline.

## Reproduction status

`benchmark.py` records source and input hashes, archive byte correspondence,
versions, parameters, seed, budgets, target counts, every pricing point,
failed-IV weight, and VIX forward estimates. Its mean standard errors do not
include finite-inner or strike-truncation bias. It does not reconstruct the
legacy robust/scaled optimization objective or establish global optimality.

Next: assess independent-seed and sampling/grid refinement; recover the exact
final parameters and scales if available; then freeze the common ESN/M2 target
contract before optimization. Missing provenance does not block a new,
explicitly versioned benchmark fit.

## Reporting correction in our own VIX diagnostics

The first M2 VIX runner counted `2*sample_SE <= price_scale/4` as an MC gate
without checking whether any positive payoff had been observed. A deep OTM
option with zero payoffs in the outer sample has zero empirical SE; that does
not bound its rare-event expectation. This is a diagnostic reporting problem
in **our runner**, not a newly identified defect in Othmane's source.

The revised runner reports payoff-hit counts, separates the nominal SE screen,
and excludes fewer-than-two-hit cases from even the minimal screen. Two hits
still do not certify tail accuracy. Earlier `run.json` files are preserved;
`audit_tail_samples.py` adds `tail_diagnostics.json` with input/source hashes.
In the 256-outer-state run, all 21 nominal screen passes are zero-hit cases
(1/6/7/7 across the four expiries). No quote passes the corrected minimal
screen. Model prices, residuals and original records are not rewritten.

The first development runs recorded source hashes at completion. That is not
an immutable record of loaded runtime sources if files change during a long
calculation. The new runners capture sources before calculation, reject
mid-run changes, and save the complete captured files under `source/` in each
run. The `m2_vix_frozen_*` and `m2_spx_frozen_controls_001` runs replay the main
diagnostics with this stronger provenance. Initial runs remain preliminary
records; no old hash or numerical result is silently updated.

## Existing time-step decision not carried into initial pricing grids

On 21 October the user pointed back to the model's high mean-reversion rates.
The pre-existing `experiments/timestep_convergence_001/RESULTS.md` in the M2
package explicitly rejects 2/4/8 steps per day for final evidence, sets 64 as
the working baseline and requests a targeted 64-versus-128 confirmation. The
model source's later fine-grid appendix also records 64-step experiments.
The new SPX diagnostics nevertheless began at 2/8 steps per day. They were
labeled numerically unresolved, but the initial implementation work did not
carry forward that existing resolution decision. Those runs are coarse engine
and variance-reduction diagnostics only; they cannot establish accurate
continuous-time option prices.

At the frozen maximum rate 5000/year, lambda_max*dt is approximately
9.92/2.48/0.310/0.155/0.0775 at 2/8/64/128/256 steps per trading day. Exact OU
endpoints do not make the nonlinear stock integral exact between endpoints.
The next option-pricing convergence study must keep the model/state fixed and
compare coupled 64/128 grids, extending to 256 if price-relative tolerances
remain unresolved. The old physical-statistic convergence study motivates this
starting resolution; it does not itself establish option-price convergence.
