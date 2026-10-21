# Conditional M2 VIX and SPX controls — 21 October 2026

## Implemented scope

Conditional VIX valuation is now implemented for baseline Gaussian M2, using
two engines: three-dimensional conditional Gauss-Hermite integration and
nested full-state OU continuations. Both use the same date-zero normalizer,
average future variance before taking its square root, and integrate the
continuous 30/365-year window with curve-jump-aware midpoint rules. The full
conditional state is retained; the three projections are not assumed Markov.

The VIX market diagnostic prices **actual absolute strikes** and compares the
model forward separately with the quoted forward. It does not rescale the
model's VIX draws to market moneyness. Only training quotes are evaluated.

SPX now has a stock-martingale control and a deterministic-volatility reference
option control. Coefficients are fitted on independent pilot paths and frozen
for evaluation; controlled estimates are not clipped to nonnegative values.

There are **60 passing M2 package tests and six passing legacy benchmark tests**.
The original physical-stability edit and all legacy ESN sources are preserved.
No calibration, physical-fit update or held-out market scoring was performed.

## Run ledger

All run directories are under `runs/`, with source/configuration hashes and
quote-level results. The frozen Stage-3 anchor, projected reference state
(C,Y,F)=0 and diagnostic training-mid curve are inherited from the SPX study.
The preprocessing flag requiring nuisance-curve profiling remains unresolved.

| Run | Outer states | Conditional GH orders | Nested inner counts | Window nodes |
| --- | ---: | --- | --- | --- |
| `m2_vix_reference_001` | 64 | 12,20 | 128,512 | 42 |
| `m2_vix_outer_refinement_001` | 256 | 20,32 | 128,512 | 42 |
| `m2_vix_window_inner_001` (55d only) | 64 | 20,32 | 512,2048 | 84 |

Date-zero normalization uses GH order 48 throughout. Every maturity samples
its outer state directly with one exact OU transition. Consequently the first
64 outer states in the 256-state run match those in the 64-state run. Nested
streams are keyed by seed, outer index and integration node; increasing inner
count at a fixed window preserves the earlier continuations as prefixes.

The main diagnostics are replayed in `m2_vix_frozen_outer_001`,
`m2_vix_frozen_window_001` and `m2_spx_frozen_controls_001`. These save the
captured executable sources under `source/`, check for mid-run changes and
apply the corrected payoff-hit screen directly. The initial development
records above remain intact; they recorded hashes at completion without a
frozen source snapshot. Prefer the frozen replays for resuming this work.
All shared quote-level pricing columns in the frozen replays exactly match
the corresponding development records. Every stored source file was checked
against its recorded hash; `replay_comparison.json` retains these checks.

The two engines share the specified OU law/readout, while evaluating conditional
expectations by different numerical methods. Their agreement checks numerical
correspondence; it is not an independent verification of the model specification.

## VIX numerical findings

For the 256-state, conditional-order-32 run:

| Expiry | Training quotes | Model forward ± outer SE | Market forward | Model second moment ± outer SE | Curve second moment |
| --- | ---: | ---: | ---: | ---: | ---: |
| 27d | 40 | 19.38 ± 0.69 | 20.10 | 495.92 ± 61.72 | 476.52 |
| 55d | 47 | 20.56 ± 0.69 | 20.73 | 542.43 ± 55.30 | 601.08 |
| 118d | 44 | 21.43 ± 0.74 | 21.19 | 600.47 ± 72.08 | 609.50 |
| 153d | 44 | 20.67 ± 0.76 | 21.72 | 573.14 ± 68.95 | 609.50 |

The second moment is checked against
`10000 * integral_T^(T+30/365) xi(u) du / (30/365)`. Differences are consistent
with the broad sampling errors here. These SEs exclude normalizer, conditional
quadrature and integration bias; agreement at this budget is not a certificate.
The initial 64-state 27d moment was substantially lower and moved upward when
more outer states were added, illustrating sensitivity to the sampled tail.

Conditional GH order 20 to 32 changes the estimated forward by at most 0.00187
VIX points and option prices by at most 0.00408 points in this outer sample.
The nested 512-minus-128 inner-count forward changes are approximately
0.11/0.46/0.26/0.50 points, with paired outer SEs 0.27/0.26/0.29/0.30. This
supports further inner refinement rather than treating fixed-inner SE as total
error. Finite-inner futures bias is generally downward for an unbiased positive
inner variance estimate; no common option-bias sign is assumed.

At 55d with 2048 inner paths and 84 integration nodes, nested minus conditional
order-32 quadrature forward is -0.0545 points (paired outer SE 0.1805).
This is closer than the initial smaller-inner diagnostic, but neither result
establishes a converged nested estimator.

An isolated **time-integration** comparison uses the same first 64 outer states,
conditional order 32 and normalizer order 48 at 42 versus 84 midpoint nodes.
The forward changes by -0.00341 (paired SE 0.00292); the largest individual-state
VIX change is 0.1806. Details and source hashes are retained in
`m2_vix_window_inner_001/window_integration_comparison.json`. Other settings in
the window/inner run change together, so raw cross-run nested differences must
not be attributed to the time rule alone.

## Tail sampling and reporting correction

Original VIX `run.json` files include nominal sample-SE screen fractions.
An additive audit found those passes were all zero-hit payoffs: none of the
sampled outer VIX values generated a positive payoff at those strikes. The
empirical option SE is then zero without bounding the true tail expectation.

Each run now has `tail_diagnostics.json`; original outputs remain intact. In
the 256-state run, 1/6/7/7 quotes across the four expiries have zero hits, exactly
matching the 21 nominal screen passes. **No quote passes even the corrected
minimal nonzero-hit screen at this budget.** The revised runner records hits
and nominal screens separately. Requiring two hits is only a minimal diagnostic,
not a confidence guarantee. The correction is also recorded in [AUDIT.md](AUDIT.md).

Neither the apparently close forward levels nor the small conditional-order
changes establish a joint option fit. Outer tail sampling and the flagged
input curve remain material unresolved issues.

## SPX variance-reduction diagnostic

`m2_spx_controls_001` uses 4096 evaluation paths and 2048 independent pilot
paths, a coupled 2-to-8-step/day grid, and normalizer orders 32/48. The two
controls have known expectations under the implemented discretized stock law.
No evaluation-path coefficient fitting or stock-path rescaling is used.

| Expiry | Plain maximum price SE | Controlled maximum price SE | Median plain/controlled variance ratio |
| --- | ---: | ---: | ---: |
| 29d | 4.52 | 2.97 | 1.46 |
| 57d | 6.46 | 3.11 | 1.80 |

All controlled price estimates in this run are nonnegative, but zero quotes
pass the planned sampling gate. Controlled paired time-step differences also
remain large relative to quote scales. The controls improve this configuration;
they do not yet make the pricer calibration-ready at the chosen budget.

## Next work

Improve outer/tail integration and SPX stock-integration/sampling accuracy,
including independent final seeds. Then add the common-protocol ESN price
evaluator, nuisance-curve/state profiles and joint objective. Keep the original
diagnostic curve conflict visible throughout. Do not compare these price errors
with the legacy pooled-IV RMSE as if the target definitions were identical.

Reproduction commands (choose fresh output directories):

```sh
python -m matteo.model_comparison.option_pricing.m2_vix --outer 256 \
  --order 20 --check-order 32 --output /tmp/m2_vix_outer_new
python -m matteo.model_comparison.option_pricing.m2_vix --expiry-days 55 \
  --inner 512 --refined-inner 2048 --order 20 --check-order 32 \
  --window-steps-per-day 4 --output /tmp/m2_vix_window_new
python -m matteo.model_comparison.option_pricing.audit_tail_samples /tmp/m2_vix_outer_new
```
