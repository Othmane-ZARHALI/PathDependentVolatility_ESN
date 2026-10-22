# SPX instability by moneyness and positivity — 22 October 2026

User direction: check whether current instability is confined to OTM options;
far-OTM difficulties can be tolerated temporarily, but negative prices are not
usable option prices. This is a retrospective audit of saved results, not a
new fit, scope reduction, estimator replacement or fresh convergence experiment.
All SPX representatives and all four VIX maturities remain recorded and active.

## Scope and definition

Source: `runs/m2_spx_linear_innovation_001`, independent pilot and two 1024-path
checking trials, direct512. All 12 representatives are technically OTM relative
to their forward: puts have K/F < 1 and calls K/F > 1. Thus an OTM/ITM split alone
cannot answer the question. No ITM or exactly ATM contract was evaluated here.

Use distance `abs(K/F - 1)` for descriptive reporting: near ATM <=5%, OTM 5–25%,
far OTM >25%. These bands are reporting conventions introduced for this audit;
they are not learned estimator selection rules or exclusions from calibration.
The four near-ATM representatives are actually about 2% OTM.

## Findings

Counts below are contract/trial observations for the pooled-pilot linear estimator.

| Group | Observations | Negative estimates | Raw payoff hits / 1024 | Estimated SE / absolute price |
| --- | ---: | ---: | ---: | ---: |
| Near ATM (four contracts) | 8 | 0 | 229–476 | 0.95–1.43% |
| OTM, 5–25% from forward | 10 | 0 | 6–46 | 5.59–29.34% |
| Far OTM, >25% from forward | 6 | 1 | 5–16 | 16.55–166.31% |

The severe new-control instability is concentrated in the wings in this sample:

- Every one of the **eight negative estimates across all methods**, including
  half-pilot sensitivity fits, is a put with K/F about **0.493 or 0.499**.
  The primary pooled-pilot negative remains the 29-day 3400 put, -4.5900 +/-
  4.7645. Raw, legacy OLS5 and baseline ridge5 have no negative price estimates
  in these trials, but the latter two are not mathematically positivity-preserving.
- For near-ATM options, augmented half-pilot SE/full-pilot SE ranges from
  **0.9931 to 1.0264**. The severe **30.52x** case is a far-OTM put. The middle
  wing band still reaches **9.75x**, so instability is not exclusive to strikes
  more than 25% from the forward.
- The largest single path's share of centered squared augmented residuals is
  **1.34–9.12%** near ATM, **7.04–76.71%** in the middle band and **21.96–94.49%**
  in the far band. This is a descriptive concentration statistic, not proof of
  population variance or reliable tail coverage. It differs from the earlier
  93.99% statistic, which concerned squared *paired estimator differences*.
- Near-ATM estimates are positive and substantially more stable, but still have
  **SE/quote-scale = 1.22–2.05**, versus the earlier target **0.125**. None of the
  eight near-ATM observations passes that sampling screen. Quote-relative
  accuracy, stock timestep and normalizer convergence remain unresolved even
  if wing difficulties are temporarily accepted. No market refit is cleared.

## Why a controlled estimate can be negative

Raw importance-sampled pricing averages `w_i * payoff_i` with positive weights,
positive forward scaling and nonnegative positive-part payoffs. Each path
contribution is nonnegative, so the raw average is nonnegative. This pathwise
property is checked on all saved samples.

The controlled estimator instead averages `w_i * payoff_i - beta dot C_i`.
Even if the controls have exactly zero population mean and beta is fitted on an
independent pilot, an individual contribution or the whole finite-sample average
can be negative. Unbiasedness and positivity are different properties. All actual
payoffs are nonnegative; it is the correction that produces the negative estimate.

A negative controlled result must be flagged as **unusable as a price**, with
its signed estimate/uncertainty retained for research. Clipping it to zero does
not fix sampling and changes the estimator's expectation. Selecting another
estimator only when this checking sample looks unfavorable also needs its own
selection-bias analysis. No clipping or sample-dependent fallback was applied.

## Practical interpretation and next step

The data support focusing the next stability work on the four near-ATM contracts
while retaining all wing diagnostics. That is a research priority, not a claim
that the near-ATM prices meet quote tolerances or an authorization to drop quotes.

For wings, the existing raw weighted-payoff estimator is a **nonnegative diagnostic**
with larger uncertainty; do not present its positivity as accurate pricing. To
make a controlled production price nonnegative by construction, a separate
positivity-preserving estimator (for example, a suitable nonnegative conditional
payoff integration) needs implementation/cost/accuracy validation. Ordinary
unconstrained linear controls, including the baseline, cannot guarantee positivity.
Any fixed choice of estimator by moneyness should be declared before new checking
samples, rather than selected on these results. No new estimator was adopted here.

Reproduce the saved-data audit with:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/pdv-options-venv/bin/python -m matteo.model_comparison.option_pricing.spx_moneyness_audit \
  matteo/model_comparison/option_pricing/runs/m2_spx_linear_innovation_001 \
  /private/tmp/fresh-spx-moneyness-audit
```

Outputs in `runs/m2_spx_moneyness_audit_001`: contract-level estimates and bands,
pilot sensitivity, path concentration and JSON summary with input/script hashes.
