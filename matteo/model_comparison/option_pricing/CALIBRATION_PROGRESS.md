# Common loss and faster conditional VIX — 21 October 2026

This checkpoint advances the implementation plan. No M2 market calibration,
P-feasibility decision or ESN/M2 ranking has been performed. The current branch
is `codex/m2-option-pricing-benchmark`; these results are local.

## Common objective

`esn_eq/calibration_loss.py` implements the planned common price objective,
independent of either pricing engine. It averages squared, frozen-scale
interval distances within each expiry, then averages expiries within each
instrument. SPX and VIX each have total weight one. The separate forward
family averages squared VIX level errors divided by 0.5 points, also with
weight one. State and curve penalties are separate. The 0.5-point scale is a
declared research convention, not measured quote uncertainty. Midpoint loss
defaults to zero and can be enabled as a separately reported sensitivity.

The objective copies targets, rejects mixed/implicit held-out splits, checks
contract uniqueness and VIX forward consistency, and returns a residual vector
whose squared norm equals the scalar objective. Invalid prices fail explicitly.
Scales never depend on model sampling errors. It reports family components,
contract residuals, put/call expiry blocks and forward errors separately.
It supplies no numerical acceptance decision.

Tests cover weighting with unequal quote counts, frozen targets, held-out
rejection, missing/failed predictions, configurable weights and penalties.
A synthetic constant-variance surface with shared SPX/VIX volatility 22% is
recovered from starts at 10% and 40%. This tests optimizer/loss wiring with an
analytic oracle; it does not establish M2 parameter recovery or a market fit.

## Conditional VIX time quadrature

`ConditionalVIXPricer.gaussian_window()` uses composite Gauss-Legendre time
integration. Lag panels start at 1/lambda_max and double in size; curve knots
always split panels. This resolves rapid initial conditional mean reversion
without imposing an equally fine grid across the whole window. The full-state
conditional law, date-zero normalizer, ACT/365 clock and 30/365 window are
unchanged. Gaussian time integration and the three-dimensional Gauss-Hermite
expectation are distinct approximations. Stock-path integration retains its
64/128/256 refinement requirements.

`runs/m2_vix_gaussian_time_001` compares time orders 4/8/16 with the archived
128/day midpoint reference at the same three diagnostic states and 55-day
expiry. Conditional GH order is 32 and normalizer order 48 throughout.

| Time order | Nodes | Largest difference from 128/day midpoint | Runtime for three states |
| --- | ---: | ---: | ---: |
| 4 | 44 | 0.0000241644 | 0.43 s |
| 8 | 88 | 0.0000237389 | 0.74 s |
| 16 | 176 | 0.0000237390 | 1.47 s |

The reference has 2652 nodes. Order 8 uses about 30 times fewer time nodes;
this is a node-count comparison, not a measured end-to-end speedup. The largest
8-to-16 change is 1.75e-9 VIX points. The remaining difference from midpoint
is compatible with the previously observed midpoint grid error, but no rigorous
error bound or guarantee over all states is claimed. Constant-variance and
nonconstant-state tests cover curve jumps and fine-midpoint correspondence.

## First combined anchor diagnostic

`runs/m2_joint_anchor_001` uses frozen Stage-3 candidate 18/scenario 0,
zero projected initial state and `training_mid_piecewise_v1`. It reuses the
4096-path, independent-2048-pilot, 256/day SPX prices from
`m2_spx_128_256_001`, after checking model/state/target hashes, core pricing
dependencies and contract data. It freshly prices VIX with 64 iid exact OU
outer states per expiry. Conditional GH orders 20/32 and time orders 8/16 are
compared on the same states. Normalizer order is fixed at 48.

Only SPX 29/57 days are included; planned 211/302-day training maturities are
absent from this diagnostic. All four VIX training maturities and level targets
are included. No validation or final-test quotes are evaluated.

| VIX expiry | Time nodes, order 8/16 | Max option time-order change | Max option conditional GH change | Max option MC SE | Zero-hit quotes |
| --- | ---: | ---: | ---: | ---: | ---: |
| 27d | 88/176 | 6.08e-10 | 0.0020692 | 0.53704 | 10 |
| 55d | 88/176 | 3.12e-10 | 0.0046703 | 1.58169 | 6 |
| 118d | 80/160 | 1.01e-9 | 0.0045863 | 2.05408 | 7 |
| 153d | 80/160 | 2.43e-10 | 0.0046922 | 1.13969 | 12 |

Price changes/SEs are undiscounted VIX-point premiums. Across these states,
the largest change in the conditional VIX level is 5.80e-8 points. Conditional
GH error and outer sampling remain substantially larger. Zero-hit quotes
remain unresolved despite zero empirical SE. Outer SE excludes shared
normalizer and quadrature error; fixed order 48 is not an accuracy certificate.

The descriptive loss is SPX 963.9281 + VIX 218.2625 + forwards 29.7125 =
1211.9030. State/curve penalties are zero because those quantities are fixed.
Do not compare this restricted, noisy score to a full-maturity fit or interpret
its components as sampling-standardized statistical evidence. The original
variance-curve compatibility flag remains unresolved. Outputs include block,
contract and forward residuals, sampled VIX levels, configurations, input and
dependency hashes and frozen sources. Runtime was about 86 seconds.

## Validation and next work

68 package tests and eight legacy/driver checks pass. Earlier reports saying
60 package tests describe their older source state. The plan's stale 2/8/32
starting ladder was corrected and logged in `AUDIT.md`; original numerical
records remain unchanged.

Next priorities:

1. Benchmark exact outer Gaussian sampling with randomized Sobol/PCA against
   iid sampling. Use independent scrambles and estimate uncertainty across
   scrambles; Sobol points are not iid observations. Retain tail diagnostics
   and refine conditional GH/normalizer orders.
2. Improve SPX sampling with stronger conditional/control methods, retaining
   high-resolution coupled stock-grid checks.
3. Implement the declared restricted variance-curve/state profiles and run
   a labelled exploratory search before any accuracy-qualified fit claim.
4. Complete all planned SPX training maturities and the common ESN evaluator.
   Preserve differences in dated conditioning, clock and VIX definition.
5. Require fresh-seed numerical gates, physical checks and untouched held-out
   evaluation before ranking fits.

Reproduce from the repository root in the project Python environment:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -m matteo.model_comparison.option_pricing.vix_time_rule \
  --output /tmp/m2_vix_time_new
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -m matteo.model_comparison.option_pricing.joint_diagnostic \
  --outer 64 --output /tmp/m2_joint_anchor_new
```

The environment used here is `/private/tmp/pdv-options-venv/bin/python`.
Commands refuse existing output directories. The archived input runs are
required for these reproducibility checks; their hashes are recorded in the
new outputs. A new conversation can resume from this file plus
`documentation/PRICING_STATUS.md` without reconstructing chat history.
