# Recorded next steps — 22 October 2026

User direction: record the VIX tail-control next step, move forward for now,
and keep all maturities. The working interpretation is to return to SPX
numerical diagnostics while leaving the VIX extension pending.

## Scope retained

[active_scope.json](active_scope.json) remains authoritative and unchanged:
**27/55/118/153-day VIX maturities, all 175 training options and four forwards**.
No VIX maturity is deferred or removed. The existing SPX maturity/split scope
also remains unchanged; representative pilot contracts are diagnostics, not
an exclusion of the other targets. Held-out market observations stay reserved.

## Pending VIX extension — recorded, not launched

Screen payoff-shaped Gaussian tail controls with analytically known means,
using the existing nine-control/ridge design and pooled 28-pilot fit as the
reference. Learn and compare candidates on design-only pilot replications,
freeze the dictionary and coefficients, then validate on new shared checking
paths at all four maturities. Retain raw/legacy estimates, payoff-hit counts,
negative adjusted estimates and paired sample-count checks. Keep conditional
quadrature, time integration and normalization checks separate.

This extension remains pending for a later bounded experiment. It is not a
completed result, an automatic background job or a convergence clearance.
The current checkpoint is
[VIX_GAUSSIAN_STABILITY_RESULTS.md](VIX_GAUSSIAN_STABILITY_RESULTS.md): all four
pooled forward SE screens pass, but the 153-day pass is borderline and not
supported by every pilot fit; only 45/175 pooled option SE screens pass, and
118-day second-moment warnings remain recorded. No market refit is authorized
by those numerical screens alone.

## SPX work to advance next

Resume the unresolved sampling, stock-integration and normalization work in
[NUMERICAL_STABILITY_RESULTS.md](NUMERICAL_STABILITY_RESULTS.md):

1. Design a small path-level Gaussian importance-sampling pilot against the
   current adapted hedging-control baseline. SPX depends on the whole stock/
   volatility path; the VIX terminal-state proposal cannot simply be reused.
   Preserve the Q model, use exact path likelihood weights and separate pilot
   learning from checking. Check martingale/reference expectations and rare
   option payoffs with tractable Gaussian/Black oracles before market runs.
2. Keep coupled 64/128/256 steps-per-trading-day comparisons. Exact OU state
   transitions do not remove nonlinear stock-integration bias; the old 2/8-step
   diagnostics remain coarse historical evidence.
3. Measure sampling error, paired-grid correction variance/cost and normalizer
   refinement separately. Preserve zero-hit and unsuccessful cases; do not
   infer accurate tail prices from a zero empirical SE or fit a new market
   model before the numerical gates are resolved.

The latest source review confirms that the SPX runner already supplies a
12-contract diagnostic, independent pilot coefficients, adapted hedging controls
and coupled 64/128/256 grids. Its recorded results do not qualify the full
surface. No new SPX simulation has been launched by this queue update.
