# Return to M2 market fitting — 23 October 2026

## Subsequent user decision: actual joint optimization

The user subsequently directed an actual exploratory joint fit now, with Lean
checks, instead of requiring another standalone profile first. This supersedes
the unlaunched nine-point/two-level/24-SPX proposal below. The current run
`m2_joint_curve_fit_002` frees all three declared curve levels and scores all
707 cutoff-3 SPX priority quotes, all 175 VIX options and four forwards; the
other 273 SPX training quotes remain diagnostics. Structural/state/Q objects
remain fixed. See [JOINT_CURVE_FIT_RESULTS.md](JOINT_CURVE_FIT_RESULTS.md) and
[formal correspondence](formalization/CORRESPONDENCE.md) for current results,
scoped Lean checks and remaining numerical limitations. The earlier decision
below is retained as historical rationale; it is not the current next task.

## Decision

**No-go for a qualified market fit or fit ranking.** The current evidence can
support an explicitly provisional, curve-only joint exploration, but cannot
validate its selected minimum. Do not free readout, state and Q risk prices
together, and do not launch another variance-reduction feature search by
default. The bounded conditional-control result and final estimator disposition
are recorded in [SPX_CONDITIONAL_CONTROL_RESULTS.md](SPX_CONDITIONAL_CONTROL_RESULTS.md).

This is a numerical readiness decision, not a claim that M2 cannot fit the
market. The user's authorization covers further exploratory investigation;
the lack of numerical accuracy limits its conclusions, not that authorization.
No new market parameters were fitted in the conditional-control experiment.

## Separate readiness gates

The table is an evidence assessment, not a sum of unrelated observed changes.
An SE at fixed grid/GH order excludes deterministic approximation error.

| Quantity | Evidence actually available | Consequence |
| --- | --- | --- |
| SPX sampling | Fresh controlled checks on 24 representatives of cutoff 3; see the new conditional-control report and its signed estimates, pilot halves and tail concentration. Earlier original central-panel run passed 0/48 primary screens. | No full-priority-population accuracy clearance; positivity and an efficiency gain do not supply one. |
| SPX stock timestep | At the same frozen model/selected curve, `m2_spx_mlmc_001` records 256-to-512 corrections above scale/8 for 8/12 and 9/12 short-expiry representatives. Median absolute observed changes are 0.2333 and 0.8359 scales. The 57-day 6775 put changes -2.57035 +/- 0.87963 points in trial 0. | Bias is unresolved even at 512/day. This historical paired-grid evidence is not a new run on all 24 representatives, and does not bound the 512-to-continuous-time remainder. |
| Date-zero SPX normalizer | `m2_spx_path_importance_001` holds paths fixed while changing GH48 to GH64 at 256/day; observed contract changes are below 0.0013 scales. The curve, state and frozen model agree with the current records. | Smaller than sampling/grid effects in that bounded check. Not a population error bound, full-maturity certificate or license to reuse normalization after changing state/Q/readout. |
| VIX sampling | `m2_vix_gaussian_stability_001`: pooled 28-pilot controls pass 45/175 option sampling screens. All four pooled forward SEs pass 0.0625, but the 153-day SE is 0.06157 and one balanced-pilot alternative gives 0.06666. | Most option prices remain unresolved. A forward-SE pass is conditional on the pilot and fixed quadrature and is fragile at the last expiry. |
| VIX conditional Gaussian quadrature | Larger-pilot study selected-state GH32-to-48 changes have maxima 0.00055/0.00041/0.00326/0.00120 VIX points. Earlier matched-state refinements did not explain the original large moment warning. | Selected-state sensitivity only; no global conditional-price bias bound or universal clearance under parameter changes. |
| VIX time quadrature | The recorded 55-day, three-state Gaussian time-order 8-to-16 change is at most 1.75e-9 VIX points; all other rules held fixed. | Validates that limited comparison. Does not measure stock timestep error or establish a full-distribution/all-candidate VIX time bound. |
| VIX moment / rare tails | Independent direct-Gaussian marginal checks support the exact normalized-model second moment. The larger-pilot study still records a pooled 118-day residual -4.55861 +/- 2.24254 squared VIX points. | Keep the warning; do not enforce the identity by rescaling quotes, or interpret it as a proved normalization bug or model infeasibility. Moment agreement alone cannot price nonlinear VIX payoffs. |

Sources: [SPX MLMC](SPX_MLMC_RESULTS.md),
[SPX path importance](SPX_PATH_IMPORTANCE_RESULTS.md),
[VIX stability](VIX_GAUSSIAN_STABILITY_RESULTS.md),
[VIX moment diagnosis](VIX_MOMENT_DIAGNOSIS.md),
[time integration](CALIBRATION_PROGRESS.md).
The actual saved run configurations were checked for common curve and dated
state; the old SPX grid/normalizer experiments cover only the short expiries.
No old short-expiry result is promoted to the long-expiry surface.

The common requirement remains `SE <= scale/8` for every scored quote, plus
separately assessed total timestep/normalizer/conditional/time bias at most
`scale/4`. A small observed refinement with uncertainty larger than its budget
does not clear bias. Multiple contract warnings are dependent diagnostics, not
independent significance tests. Repeated sampling checks do not replace a
refinement assessment, and vice versa.

## Smallest useful restricted fitting experiment

The next fitting question is **whether joint SPX constraints change the curve
direction favored by VIX alone**. Repeating the completed VIX-only two-level
optimizer would not answer it. The previous train/check reversal remains a
failed validation result, not a preferred fitted curve.

Use the declared three-segment nuisance family, initially freeing **only its
first two variance levels**. Keep the third at 0.0642080964154938, with its fixed
status visible; it is not identified by the VIX windows, but it does affect the
longest SPX expiry and is not globally unidentifiable. Freeze the Stage-3
readout, architecture, zero projected state and Q scenario 0. Report the
penalty `0.1*sum(diff(levels)^2)` separately. Varying a deterministic curve does
not allow terminal-stock caches to be reused without repricing the stock;
the OU projections/equity shocks and compatible VIX segment exposures may be
reused. A changed state, readout or Q scenario invalidates those exposures and
date-zero normalization.

For a first bounded joint diagnostic, explicitly score the **24 fixed SPX
representatives plus all 175 VIX training options and all four VIX forwards**.
SPX membership is the immutable `contracts.csv` in
`runs/m2_spx_conditional_controls_001`; the other **683 cutoff-3 priority quotes
remain unpriced by this provisional subset objective**. Retain all 980 SPX
training quotes and wings in the overall research scope. Do not report this
subset score as the full 707-priority or 980-SPX common score.

Use `JointPriceObjective` with expiry-balanced SPX/VIX family weights one,
forward-family weight one and fixed 0.5-point forward scale. Keep the existing
market quote scales and interval loss; no noisy-model-SE denominator, midpoint
tie-breaker or tail exclusion. Save the exact subset weights. Use common
randomness across curve candidates, independent control pilots, and frozen
coefficients during each comparison. With independent frozen coefficients,
reuse is legitimate for target-zero controls under candidate-specific prices;
recompute their actual values/reference Black means for each candidate.

A concrete low-budget first step is a **nine-point 3-by-3 local profile**:
multiply each of the first two frozen selected variance levels by
`{0.95,1,1.05}`, retaining every point and applying the declared family bounds.
The frozen selected curve is a numerical comparison anchor, not an adopted
market fit. Use 256 pilot paths and 512 common searching paths per SPX expiry;
for VIX use eight independent 64-point nets per expiry with the existing
conditional/time rules. This is below the plan's 4096-path screening budget:
call it a **numerically unresolved local direction diagnostic**, not a qualified
screen or an optimizer-converged fit. A 60-minute wall budget is a stop boundary,
not an instruction to drop contracts or declare unfinished points losers.
Retain interrupted/failed points explicitly. Do not run wider parameter searches
automatically when the profile is noisy or its boundary is best.

Freeze the selected point using searching data only. A distinct numerical
checking stage compares it with the anchor using **three fresh seeds** and at
least the searching budget, without reselecting on those checks. Save paired
price/loss differences and replication-level diagnostics, including the
uncertainty of the nonlinear interval loss. More checking paths may be needed;
minimum counts do not establish accuracy. A reversed/unstable result remains
unsuccessful and does not trigger another profile on the checking draws.

**This profile has not been launched in this checkpoint.** The work completed
here is the bounded conditional-control investigation and readiness decision.
The next checkpoint can execute this declared profile without reopening scope
or repeating the old VIX-only fit. A low-budget local direction is useful only
as a proposal for later accurate numerical checking. Its apparent gain cannot
establish market improvement, feasibility, or uniqueness.

After a persistent curve direction is established with sufficient numerical
accuracy, add the third curve level, then restricted C/Y/F state coordinates
with the existing stationary-prior penalty (weight 0.1), then the predeclared
Q scenarios one at a time. Recompute dated normalization, prices and appropriate
P diagnostics after shared parameter/state/Q changes. Do not identify state or
risk prices uniquely from a noisy single-date fit. P-feasibility gates, the
reserved market splits, complete quote coverage and the common ESN protocol
remain future qualification requirements.

## Undergraduate interpretation

There are two distinct questions: “Which parameters look best on these random
paths?” and “Which parameters give the best expected market prices?” An
optimizer can answer the first very accurately while missing the second. The
previous VIX reversal is an example. Averaging out a little randomness may help,
but it does not recover rare volatility paths that the retained sample misses,
and it does not correct a stock update made on a finite time grid.

The restricted profile therefore asks a small, useful economic question with
its numerical uncertainty visible. It does not turn a lower searching score
into a defensible joint market fit.
