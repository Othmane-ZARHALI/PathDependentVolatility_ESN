# Recorded next steps — 22 October 2026

User direction: record the VIX tail-control next step, move forward for now,
and keep all maturities. The working interpretation is to return to SPX
numerical diagnostics while leaving the VIX extension pending.

The user subsequently confirmed: **"Yes, let's focus now on that for a while."**
SPX sampling, timestep and normalization diagnostics are therefore the current
explicit focus. The first bounded path-sampling experiment is recorded in
[SPX_PATH_IMPORTANCE_RESULTS.md](SPX_PATH_IMPORTANCE_RESULTS.md); the VIX queue
and full maturity scope below remain in force.

SPX follow-up: [SPX_VOLATILITY_SAMPLING_RESULTS.md](SPX_VOLATILITY_SAMPLING_RESULTS.md)
records independently selected volatility-path sampling that improves fresh
variance in 20/24 cases, while no quote-relative sampling screen passes.
Additional Gaussian-score controls help only modestly. Next prioritize paired
time corrections: standalone level costs, independent correction-pilot fits,
and a bounded matched-cost MLMC comparison. The 128->256 warnings, pilot
instability and tail precision remain open; no model refit is cleared.

The [first MLMC prototype](SPX_MLMC_RESULTS.md) now checks independent
64/128/256/512 levels with frozen allocation and correction-specific controls.
It loses on the central variance-cost objective in both fresh trials and
**does not replace direct vol7 simulation**. Sparse top-level tails and
256->512 warnings remain. Before another allocation study, investigate better
coarse/fine Gaussian coupling on tractable oracles; conditional bridge
reflection is a candidate only, with target-law/likelihood/adaptedness checks
still to be established. Do not interpret the completed prototype as a
completed improved coupling or continuous-model accuracy clearance.

The [endpoint reflection study](SPX_BRIDGE_COUPLING_RESULTS.md) has now completed
the algebraic/Gaussian/Black and independent checking step. Endpoint and target
law checks pass, but all four central variance-cost comparisons lose (gains
0.67–0.77). Direct vol7 and ordinary corrections remain the references.
The same-coefficient averaged residual retains 86–96% of controlled variance.

**Next SPX hypothesis, pending:** cheap projected-OU innovation martingale
controls, including centered quadratic innovations with known conditional
covariance, as complements to the existing stock hedges. Preserve full-state
conditional means; the three readout projections alone are not Markov. Prove
the zero-mean/likelihood correspondence, validate Gaussian/Black and forward
oracles, train on separate pilots, then check direct-price and correction
variance/cost separately. This is not an implemented control or a bias fix;
512-step, tail/pilot and normalization gates stay open. VIX scope is unchanged.

The [full OU-innovation study](SPX_INNOVATION_CONTROL_RESULTS.md) now implements
and independently checks that hypothesis. Central direct-price variance gains
around 1.8–1.9 are nearly consumed by a 1.8 runtime multiplier; corrections
lose on cost. Repeated far-call paired mean shifts, contract-control mean
warnings and pilot sensitivity remain unresolved, despite Gaussian/covariance
and bare-innovation checks passing. The full-nine dictionary remains experimental.

**Next candidate:** an actually cheaper three-linear-control kernel. Old-pilot
two-fold ablations retain most direct-price benefit, but are post-hoc design
evidence only. Avoid unused square/mixed calculations; profile standalone
pricing; freeze before fresh checking and retain every representative strike.
Audit contract-specific control tails/paired mean shifts as well as bare
innovation sums. Keep correction performance and timestep bias separate;
baseline controls remain the reference and no market refit is cleared.

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
## Linear-only follow-up completed; next diagnostic

The [direct512 linear-only experiment](SPX_LINEAR_INNOVATION_RESULTS.md) retains
a modest central variance-cost gain (1.14x in both fresh trials), but far-put
outlier leverage and 30.52x half-pilot SE inflation prevent general adoption.
All primary sampling screens still fail. Keep the baseline as reference.
Next examine saved-data coefficient/scaling and tail-contribution diagnostics;
any shrinkage/stability rule is exploratory until frozen using pilot data and
validated on new checking samples. Do not select a rule on existing checking
performance. Stock-timestep and normalization accuracy remain separate open gates.
