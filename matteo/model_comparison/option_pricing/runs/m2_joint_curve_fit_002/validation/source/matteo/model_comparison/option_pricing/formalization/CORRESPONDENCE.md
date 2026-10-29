# Lean checks for the joint curve fit

## Checked scope and environment

`CalibrationChecks.lean` contains ten checked theorems, compiled with
warnings treated as errors under **Lean 4.35.0-rc4** and Mathlib commit
`1f414401f69059aa7eead47b53ee40bd38455eeb`. The successful audit is
`audits/attempt_003/audit.json`; its compiler log includes axiom reports and
explicit elaboration of the principal assumptions. All ten declarations use
only `propext`, `Classical.choice`, and `Quot.sound`. There are no admitted
proofs, custom axioms, or `sorryAx` dependencies in the successful audit.

The first two compiler attempts remain as failed development evidence. Their
automatically generated error placeholders do not count as checked proofs.
The final source and final axiom audit supersede them. The proof statements
were not weakened to remove the Gaussian law, integrability, or zero-mean
control obligations.

This project has no separate downloaded dependency tree. `check_lean.py`
reuses the installed pinned Mathlib environment **read-only**, checks its
actual Git revision and tracked-source cleanliness, then writes proof outputs
only here. It neither invokes Lake nor upgrades dependencies. Reproduce with:

```sh
/private/tmp/pdv-options-venv/bin/python \
  matteo/model_comparison/option_pricing/formalization/check_lean.py
```

For a different installation, pass `--mathlib-project PATH` with the same
pinned toolchain, manifest revision and compiled dependencies. Each invocation
creates a new audit directory. An old passing audit does not certify edited
source: compare the recorded SHA-256 before reusing it.

## Theorems, numerical roles, and remaining obligations

| Lean declaration | Ordinary mathematical statement | Corresponding numerical object | Limits of correspondence |
| --- | --- | --- | --- |
| `gaussian_exponential_mean_one` | For Z~N(0,v), E[exp(Z-v/2)]=1, including v=0. | Exponential-Euler stock innovation, conditional on start-of-step variance. | Uses the exact real Gaussian distribution. Does not prove that floating-point draws, roots or stock updates realize that distribution without error. |
| `gaussian_stock_step_mean` | For a fixed current stock S, E[S exp(Z-v/2)]=S. | Adapted normalized stock update in `pricing.py` and `hedging_controls.py`. | Relating this kernel to the full correlated OU filtration and iterating conditional expectations remains an ordinary mathematical/source argument, not a Lean formalization of the entire simulator. |
| `gaussian_hedge_step_mean_zero` | For fixed start-of-step S and delta, E[delta*(S_next-S)]=0. | Each existing adapted hedge increment. | Future-dependent delta would invalidate the interpretation. The source uses start-of-step information; full path integrability/tower composition is not encoded here. |
| `frozen_anchor_control_reuse` | Under one fixed measure, integrable known-zero controls may be subtracted from any integrable candidate payoff with arbitrary fixed coefficients. | Anchor stock/reference/hedge controls reused across deterministic curve candidates, and fixed Gaussian VIX controls. | Zero means are explicit premises, not inferred from observed sample averages. This theorem does not itself prove the Black reference mean or all Gaussian-control means, nor justify checking-dependent coefficients. Apply it conditional on an independently frozen pilot. It cannot justify reusing a changed-law state/Q/readout cache. |
| `finite_importance_identity` | Sum q_i*(p_i/q_i)*H_i equals sum p_i*H_i for nonzero q_i. | Finite-law oracle for likelihood placement. | Algebraic finite identity only. It is not a formal proof of the continuum Gaussian-mixture density, score basis or f/q implementation. Positivity is needed for the probability interpretation, although the algebra holds more generally. |
| `integrated_variance_affine` | I(xi)=c+sum A_j*(xi_j-floor) is affine in xi. | `vix_curve_profile.py`, `VIXData.price`, and reuse of segment exposures. | Concerns integrated variance, not VIX or option prices. The derivation of A from exact conditional Gaussian expectations is not formalized; the implemented finite-GH/time cache has independent source/tests. |
| `integrated_variance_positive` | If c>0, A_j>=0 and xi_j>=floor, then I(xi)>0. | Positivity of cached VIX integrated variance and square-root input. | Does not establish quadrature accuracy or an independently correct normalizer. |
| `weighted_residual_norm_identity` | Sum (sqrt(w_i)*r_i)^2 = sum w_i*r_i^2 when w_i>=0. | Residual vector supplied to least-squares with frozen expiry/family weights. | Checks the real arithmetic identity. Correct row membership, weights and finite-precision evaluation are checked in Python, not proved by this theorem. |
| `penalty_residual_identity` | Two sqrt(0.1)-scaled adjacent curve differences have squared norm 0.1*sum(diff(xi)^2). | Existing declared nuisance-curve penalty represented as two residuals. | Exact real 1/10 differs from the finite binary representation of 0.1; numerical correspondence is tested, not bitwise formally verified. |
| `finite_noise_square_decomposition` | With normalized weights and mean-zero error, E[(a+error)^2]=a^2+E[error^2] in a finite model. | Explanation of noise in the fixed-scale squared interval loss outside a fixed interval edge. | Does not prove finite-sample fit ranking. With changing active interval edges or adaptive parameter selection, this fixed-edge formula alone is insufficient. |

## Why the anchor controls may be reused

The current fit changes only the deterministic variance curve. The Gaussian OU
state law, its Q drift, primitive equity shocks, importance proposal, state and
readout stay fixed. An anchor hedge is therefore still the same known-zero
function of the same underlying random experiment when another curve is
priced. Its control value and its reference Black centering stay anchored;
the **candidate stock payoff is repriced**. This differs from recomputing
candidate-specific controls and must be named explicitly.

To cover all SPX strikes economically, the numerical fit assigns each strike
the nearest same-side representative's two anchor hedge gains. These gains
retain their zero-mean interpretation; no claim of equal variance or efficiency
to the per-strike production dictionary follows. The pilot fits all five
coefficients independently at the anchor and freezes them before search/checking.
This is a control-reuse choice inside fitting, not a new variance-reduction
study or a change to the M2 price target.

The integrability premises also have a concrete ordinary-mathematical argument
for this fixed finite-grid model: the tempered readout is bounded, each dated
normalizer is deterministic and strictly positive, curve levels are bounded,
and the grid is finite. Thus each step variance has a finite deterministic upper
bound. Gaussian exponential moments then give finite stock moments recursively.
Black-inspired deltas are bounded by one, so finite sums of their stock-increment
gains are integrable; reference lognormal payoffs are integrable too. The
defensive importance weight is at most two. This supports the required
integrability but is not itself a Lean proof of the Python readout bounds,
full recursive process or Gaussian-mixture change of measure. Poor tail
precision can persist despite the existence of all these moments.

Changing readout, initial state or Q risk prices changes the law or normalization
and requires revisiting these conditions. Lean's generic fixed-measure theorem
does not remove that requirement.

## Signed estimates and numerical selection

The production `JointPriceObjective` rejects negative model prices. The new
`ExploratoryJointLoss` explicitly extends the **same interval-distance formula**
to finite signed Monte Carlo estimates. Negative estimates remain in the data
and the search residual; none is clipped or turned into an admissible quote.
Python tests compare all family components against production on nonnegative
inputs and check that production still rejects negative inputs. A lower signed
search score supplies a candidate for checking, not price or fit clearance.

Common searching paths permit efficient candidate comparisons. Because search
selects a candidate using those paths, the fixed-candidate expectation identities
do not make its searching score unbiased evidence of improvement. The selected
point freezes before three fresh numerical checks. Market validation/final-test
observations remain reserved.

## Separate assessments

- **Formal proof:** ten identities/conditional statements checked as above.
- **Semantic correspondence:** fixed-measure, adapted-information and curve-only
  scope reviewed; remaining representation bridges explicitly listed above.
- **Implementation correspondence:** stock-cache/original-pricer and direct/cache
  VIX comparisons, analytical derivative checks, source hashes and independent
  score reconstruction. These tests are finite/configuration-specific.
- **Numerical evidence:** the joint optimizer and its fresh checking runs determine
  observed fit and stability. Lean does not certify rare-tail sampling, finite
  stock-grid bias, GH/time errors, physical compatibility or M2 feasibility.

At undergraduate level: Lean checks why the algebra and probability identities
are valid under their stated assumptions. Numerical checks ask whether the
computer calculation is accurate enough to apply them to market prices. Both
are needed, and neither substitutes for the other.
