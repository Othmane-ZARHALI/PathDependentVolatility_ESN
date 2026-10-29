# Research decisions and open questions

Snapshot: 29 October 2026. This is a synthesis of recorded decisions and outcomes, with source links. It adds no modelling authorization or candidate adoption.

## Constraints and evidence standards to retain

| Decision | Consequence | Source |
|---|---|---|
| Retain the exogenous M2 structure as far as possible | A new readout is different from feedback into OU dynamics; establish the version before comparing models | [Ledger](../../matteo/model_comparison/option_pricing/RESEARCH_LEDGER.md) |
| Strictly positive spike shift and q_J | Historical zero-shift runs remain historical references | [Calibration plan](../../matteo/exogenous_reservoir_eq_m2/documentation/CALIBRATION_PLAN.md) |
| All retained maturities, cutoff 3, native quotes and signed estimates | No hidden quote exclusions, clipped negatives or smoothed smiles | [Current scope](../../matteo/model_comparison/option_pricing/active_scope.json) |
| Current feasibility work is in sample | Fresh random draws do not create market holdout evidence | [All-maturity result](../../matteo/model_comparison/option_pricing/ALL_MATURITY_FIT_RESULTS.md) |
| Preserve the experimental record | Keep unsuccessful tests, learnings, pictures and raw-source links | User direction in the consolidation conversation |
| Quintic and jumps remain deferred | Their mention in literature or brainstorming is not a completed M2 remedy | [Ledger](../../matteo/model_comparison/option_pricing/RESEARCH_LEDGER.md) |
| Keep mathematical, numerical and implementation evidence distinct | Scoped Lean checks and finite probes retain their original assumptions | [Correspondence](../../matteo/model_comparison/option_pricing/formalization/CORRESPONDENCE.md) |

## Tested approaches and reasons to avoid repeating them unchanged

| Approach | Recorded learning | A reason to revisit would be |
|---|---|---|
| Higher quadratic curvature everywhere | Stage 2 failed the joint physical proxy; fast-scale variance energy was a problem | A different scale allocation and untouched validation design |
| Coarse-grid physical validation | Exact OU transitions did not ensure accurate nonlinear stock/readout integration | Resolution evidence for the proposed candidate |
| More curve levels alone | Small/inconsistent extra gains; protected curve lowered price loss but worsened IV error | A distinct strike-shape sensitivity, rather than another blind level fit |
| Cap-only explanation of VIX wings | Cap changes redistributed variance and could hurt SPX; no single-cause finding | A controlled normalized comparison with an adequate variance budget |
| Continuous or four-speed stress extension | No demonstrated advantage over tested two-speed reference | A distinct mechanism or better forward precision, not branch count alone |
| Constant stock/stress rotation | Large SPX benefit, but missing short-expiry downside IV remains | A distinct dependence mechanism, beyond repeating the theta scan |
| Local or forward-variance adaptive directions | Local covariance maximization did not optimize finite-expiry smile shape | A finite-expiry rationale with testable new sensitivity |
| Global or conditional-mean downside multipliers | SPX improvements could cost VIX shape; conditional means are weaker than VIX protection | Explicit distributional protection or a quantified acceptable tradeoff |
| Slower existing-bank downside signals | Modest changes; stronger settings cost more VIX fit | A distinct nonlinear mechanism with ablations |
| Positive linear memory filters | Implemented; no material improvement; averaging reduces instantaneous loading | A different hypothesis, such as nonlinear history response, with its own derivation |
| First and richer bounded generator corrections | Valid quoted-date protection; modest SPX improvement and a large residual gap | A more effective shape under an explicit global positivity bound |
| Existing-stress correction budget | Larger budget was exercised, but no statistically established gain over baseline | A mechanism addressing the missing downside level, not an assumption that amplitude was never enlarged |

The corresponding reports are linked in [the experiment history](EXPERIMENT_HISTORY.md). These are reasons to retain knowledge, not universal exclusions of model families.

## Open questions

- Which exogenous mechanism can raise the missing short-expiry SPX downside-IV level while retaining the successful VIX shape?
- Can a protected variance family provide a materially more useful sensitivity direction under a practical positivity budget?
- Is the endpoint/path distinction empirically useful enough for the eventual M2 design? Existing tests establish feasibility more strongly than fitting success.
- How much of each candidate's apparent gain survives tail, stock-grid, normalizer and conditional-integration refinement?
- Can any promising Q candidate also retain acceptable physical-measure behaviour? Synthetic P proxy success and the later Q experiments remain separate evidence.

The VIX-first inverse construction discussed in the ledger is literature/design work. It is not a fitted M2 replacement. Endogenous external October architectures also remain distinct from the legacy A2 implementation and our tested exogenous modifications.
