# VIX maturity scope — 21 October 2026

**Current decision, subsequent user instruction:** restore **27, 55, 118 and
153 days**, all 175 training options and all four associated forwards, for new
analyses. No maturities are deferred. SPX scope is unchanged. This supersedes
the earlier same-day temporary exclusion below. Restoration provides the full
picture before the user chooses the next fitting scope; it does not clear any
numerical accuracy, moment compatibility or feasibility gate.

`active_scope.json` is authoritative. `near_scope_historical.json` preserves
the former policy. The near-only driver now requires that historical policy
to be selected explicitly; it cannot silently discard the restored maturities.

## Superseded two-maturity decision


**Historical decision, explicitly authorized by the user on 21 October 2026:**
temporarily remove the 118- and 153-day VIX expiries from exploratory fitting.
The remaining VIX maturities are **27 and 55 days**. Fit their two associated
forwards and defer the two longer forwards as well. SPX scope is unchanged.
The machine-readable historical decision is `near_scope_historical.json`.

This is a deliberate reduced-scope research assumption, not a conclusion that
the deferred quotes are wrong or that the model cannot fit them. Original data,
conflict flags and four-maturity results remain preserved. No reduced-scope
result establishes full-maturity SPX/VIX feasibility. Deferred training
maturities are future diagnostics, not pristine new held-out tests.

## Historical diagnostic before the explicit two-maturity decision

The last-maturity-only explanation is not supported. The earlier audit also
identified the 118-day block. `runs/m2_vix_maturity_scope_001` now evaluates
the same two frozen curves under three scopes, with no refitting or candidate
reselection. Both options and the associated forward are removed for each
excluded expiry. Expiry weights are re-averaged over included maturities;
individual market scales and the curve penalty remain unchanged.

Let D = loss(selected) - loss(starting). Negative favors selection. These
are point diagnostics at unresolved numerical accuracy.

| Included VIX expiries | Fitting-sample D | Fresh-sample D |
| --- | ---: | ---: |
| 27/55/118/153 days | -10.3056 | +5.3688 |
| 27/55/118 days | -3.0358 | +3.4569 |
| 27/55 days | -0.1470 | -0.4645 |

Dropping 153 days alone does not remove the reversal for this candidate pair.
The near maturities have consistent ranking, but their small gains have not
been shown reliable and their losses/pricing errors remain substantial.
Total losses across different scopes are not directly comparable fit-quality
measures. This fixed-candidate check does not predict a new reduced-scope fit.

At that earlier checkpoint, the four-maturity protocol was retained because
removing only 153 days did not resolve the reversal. The user's subsequent
explicit two-maturity instruction supersedes that scope decision. The
four-maturity protocol was then a preserved baseline; the latest instruction
above has now reinstated it for analysis.

## Historical records and reintroduction gate for the earlier deferral

A reduced-scope run must explicitly record:

- Included/excluded expiries and treatment of each associated forward. The
  default sensitivity removes both options and forward at an excluded expiry.
- That its fit is conditional on the reduced scope, with no claim about the
  deferred maturity, full VIX surface or full SPX/VIX feasibility.
- Preservation of original quotes, failures, conflict flags and the four-expiry
  baseline. Deferred training expiries are diagnostics, not pristine new tests.
- That scope choices used these training/numerical findings. Existing checking
  draws cannot become a fresh validation set through relabeling.
- A return gate: once shorter-expiry sampling, tails and quadrature are
  controlled, reinstate deferred expiries, recheck moment bounds and forwards,
  and evaluate the same candidate with fresh draws before any full-scope
  ranking. Final-test market quotes remain reserved.

Of the four next numerical steps, tail-aware sampling and separate error gates
apply to both VIX and SPX. Adapted hedging controls and proposed multilevel
stock-grid allocation primarily address SPX. Current quadrature-based VIX
pricing does not require a stock-path multilevel method.

The first historical near-only refit is recorded in [VIX_NEAR_RESULTS.md](VIX_NEAR_RESULTS.md).
It uses a single identifiable near variance coefficient and fresh numerical
draws. The large prior ranking reversal is absent in the point comparisons,
but the new gain is smaller than its uncertainty and pricing tolerances remain
unmet. This simplifies the research scope; it does not certify convergence.
