# Lean checks for the joint curve fit

## 29 October 2026: normalized linear memory

`memory_filter.cpp/.py` filter the realized nonlinear normalized M2 and stress signals while retaining the exogenous Gaussian factor dynamics. Constant-input and mean-one filter identities are algebraic; they do not imply VIX-law invariance or certified population-forward matching. Conditional pricing uses realized filter histories plus the full-state conditional future signal law. Numerical scalar stock/history recurrences, zero-share replay, whitened Gaussian bridge covariance/terminal probes and likelihood-weighted payoff refinement checks are in `../runs/m2_memory_filter_001`. The corrected bridge uses a rank-truncated Gaussian conditioning approximation; rejected high-rank conditioning remains separate. None of these numerical checks is a new Lean certificate, continuous-stock martingale theorem or model-family impossibility result.

## 28 October 2026: slower existing-bank downside multiplier

`slow_readout.cpp` and `slow_readout.py` retain the frozen exogenous Gaussian bank and fast M2 function while adding a selected existing slow OU mode to a one-sided readout. The conditional identity E[Hnew|stress]=a+wA follows by integrating correlated R,D jointly over the independent old bank. E[RD] must not be factored. Unit first moment is not invariance of the conditional VIX distribution. Full-state four-projection Gaussian means/covariances, conditional-X-on-Y/F quadrature, scalar tensor comparisons, alpha-zero replay and scalar stock recurrence are numerically audited in `../runs/m2_slow_readout_001` and `002`. These are finite numerical checks, not a new Lean certificate, continuous-stock martingale theorem, population-forward certificate or model-family impossibility result.

## 28 October 2026: future-integrated-variance sensitivity

`forward_leverage.cpp` computes the Brownian loading of conditional next30day integrated variance using full24mode decays, conditional Gaussian readout derivatives, analytic stresses and exact curve-break hat weights. Integration/recurrence numerical audits are in `../runs/m2_forward_leverage_001`. Primitive finite differences independently reprice the finite-rule conditional variance target; this checks the derivative correspondence at declared states, not an exact infinite-dimensional quadrature theorem. Unit predictable stock loading preserves the unchanged exogenous factor/variance law, so VIX is inherited unchanged. No new Lean proof, option-skew optimum or continuous-stock martingale result is claimed. Fixed projection memory is separately diagnosed from stationary Gaussian covariance; it is not proof of calibration infeasibility or of the proposed slower-readout remedy.

## 28 October 2026: two targeted skew numerical implementations

`targeted_skew.cpp` adds13-coordinate state-dependent unit stock exposure and a stress-regime-normalized variance readout. Scalar recursion, readout gradient, primitive/OU covariance, unit norm, stress quadrature, interpolated conditional means and new forward derivatives have numerical audits in `../runs/m2_targeted_skew_001`. No new Lean theorem is claimed. The variance/VIX-law invariance argument applies only to the stock-loading experiment; the new readout is repriced. The conditional-mean cancellation is not preservation of the full future conditional variance law. Maximizing instantaneous covariance is not a theorem of optimal finite-maturity implied skew. Normalized H0 sensitivity uses its continuous extension at xi=floor; the original rotation is the fallback only at zero H0 gradient. Existing formal sources/audits remain unchanged.

## 28 October 2026: asymmetric readout numerical correspondence

`asymmetric_readout.cpp` implements H0*D/E[H0*D] using the existing exogenous feedback projection. Because old-bank/stress factors are independent, numerator and normalizer are decomposed with E[R*D], not E[R]E[D]. Direct tensor quadrature, scalar stock recursion, zero-strength replay and forward-derivative numerical audits pass; see `../runs/m2_asymmetric_readout_001/implementation_audit.json`. No new Lean theorem was added. Existing normalization identities have their original scoped premises; they do not certify this compiled sampler, new conditional quadrature or continuous-time stock martingality. The new VIX law is repriced rather than claimed invariant.

## 28 October 2026: numerical rotation implementation

`leverage_scan.cpp` now implements the proposed rotation through four statistics (I,U0,Us,Uf), so each direction selects the third sum in the abstract three-statistic theorem. `audit_leverage_scan.py` independently reconstructs target/reference statistics and checks exact finite-step stress OU/Brownian covariance. Theta=0 full replay matches the frozen protected-candidate stock output. These floating-point checks address implementation correspondence, without constituting a formal proof of compiled code or convergence. Numerical response and signed-price/invalid-IV evidence are in `../LEVERAGE_SCAN_RESULTS.md`; the existing Lean statements and audits are unchanged.

## 28 October 2026: proposed exogenous leverage rotation

`LeverageRotation.lean` is a separate, newly checked source. Its successful audit
is `leverage_audits/attempt_001/audit.json`, using the pinned environment below,
warnings as errors and an axiom audit. It proves three restricted claims:

- `rotation_coefficient_norm`: the three Brownian-combination coefficients have
  squared norm one for theta,r in[0,1]. The assertion that the inputs are
  independent Brownian drivers and the resulting process has the intended joint
  law is a separate stochastic correspondence obligation.
- `stress_covariance_contribution_nonpositive`: the proposed stress contribution
  has nonpositive sign given nonnegative loadings and theta. It does not identify
  abstract loadings with a compiled SDE or prove finite-wing implied skew.
- `three_statistics_stock_reconstruction`: the explicitly declared discrete
  log-stock recurrence is exactly reconstructed by three finite path sums. The
  proof derives the reconstruction by induction; correct numerical increments,
  variance paths, floating-point arithmetic and the recurrence's implementation
  remain separate obligations.

No new VIX-law invariance theorem, continuous-time stock martingale theorem,
pricing-control identity or calibration claim is included. The proposal keeps
all OU dynamics exogenous by user instruction. See
`../SPX_SKEW_RESEARCH_AND_FAST_ABLATIONS.md` for the mechanism and proposed scan.

## Checked scope and environment

**Latest numerical-scope continuation:** audit `attempt_006/audit.json`
rechecks the unchanged thirteen-theorem proof source for the all-retained-
maturity, all-observation in-sample feasibility run. Warnings/axiom checks pass.
The user has superseded reserved market roles for that run; this changes data
selection, not these probability/algebra statements. The new explicit split
adapter and same-side fixed-control transfer are implementation tests, not Lean
proofs. The finite recursion and fixed-control statements do not impose a
four-maturity restriction. Longer grids still require separate numerical bias
checks. Earlier audits remain evidence for their historical checkpoints.

`CalibrationChecks.lean` now contains thirteen checked theorems, compiled with
warnings treated as errors under **Lean 4.35.0-rc4** and Mathlib commit
`1f414401f69059aa7eead47b53ee40bd38455eeb`. The state/refinement checkpoint's audit is
`audits/attempt_005/audit.json`; current audit 006 rechecks the same proof source.
The compiler logs include axiom reports and
explicit elaboration of the principal assumptions. All thirteen declarations use
only `propext`, `Classical.choice`, and `Quot.sound`. There are no admitted
proofs, custom axioms, or `sorryAx` dependencies in the successful audit.

The first two compiler attempts and attempt 004 remain as failed development evidence. Their
automatically generated error placeholders do not count as checked proofs.
The final source and final axiom audit supersede them. Attempt 003 remains the
successful ten-theorem checkpoint for the original curve-only experiment.
Attempt 004 failed while developing call/put integrability (the available
library operation is `Integrable.sup`, not `.max`); its compiler-generated
`sorryAx` dependency is a failure, not accepted evidence. The proof statements
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
| `affine_state_path_difference` | Two affine recursions with the same coefficients and forcing differ by the product of coefficients times their initial difference. | Per-mode OU initial-state translation on the same primitive Gaussian path. | The OU exponential coefficients, grid accounting, Gaussian sampler and projection mapping are source/implementation obligations; not a formal proof of the full simulator. |
| `normalized_readout_mean` | For integrable G under a probability law with nonzero exact mean, the normalized variance readout has mean xi. | Candidate-specific date-zero normalization after an initial-state change. | The denominator is the exact expectation in the theorem. Finite GH, interpolation and time-grid approximations are not certified. |
| `option_put_call_parity` | For an integrable terminal value X, expected call payoff minus expected put payoff is E[X]-K. | Same-maturity VIX model forward is E[VIX_T]; a supplied market forward is a separate target. | Does not certify finite-sample controlled prices, Black inversion, settlement conventions or an option-implied market forward reconstructed from IV-derived premiums. |

## Why the anchor controls may be reused

The original curve-only fitting stage changes only the deterministic variance curve. The Gaussian OU
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

- **Formal proof:** thirteen identities/conditional statements checked as above,
  including the restricted-state and forward extension below.
- **Semantic correspondence:** fixed-measure and adapted-information scope
  reviewed for the original curve stage and the primitive-space state extension;
  remaining representation bridges explicitly listed above and below.
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

## Restricted-state extension and forward diagnostics

`joint_state_fit.py` keeps the primitive Gaussian law, Q risk prices and readout
fixed. Candidate full states are minimum-prior-norm lifts of C/Y/F projections;
the OU path translates by exp(-rates*t)*X0. Candidate G, dated M0, stock prices
and conditional VIX exposures are rebuilt. Old **zero-state anchor controls**
remain fixed functions of the original primitive shocks, whose law and importance
likelihood do not change. This is a new primitive-space correspondence argument;
it does not justify evaluating old G/M0 or old candidate VIX exposures as if
they belonged to the changed state. Four implementation checks cover translated
projections/equity shocks, direct candidate stock repricing, prior coordinates
and normalizer ownership. A new log-GH normalizer interpolation adds a separately
screened deterministic approximation and no Lean numerical error certificate.

The generic `integrated_variance_affine`, positivity and weighted-residual results
already quantify over arbitrary finite dimension, so their algebraic scope
also applies to a five-piece curve. New segment membership, derivatives,
time-rule splitting and approximation errors still require implementation tests.

For Black IV display, using the model's own SAME-expiry mean is meaningful for
diagnosing smile shape. It does not change any option price and cannot eliminate
the market-forward calibration residual. Independent market-forward inference
from call/put parity requires independently observed option premiums; these
source quotes are IVs converted using the supplied forward, so that inference
would otherwise be circular. Any nearby-strike IV comparison is not itself a
same-strike put/call parity test.

## 29 October 2026: quoted-VIX-protected variance correction

`m2_protected_variance.cpp/.py` use a bounded additive generator correction on existing centered OU modes. Conditional Dynkin cancellation over four quoted windows is an analytic argument under the unchanged filtration and integrability assumptions. It protects model values, not exact market forwards. `runs/m2_protected_variance_001` contains compiled/Python formula checks, conditional quadrature probes, zero-amplitude replay, independent signed-correction stock recurrence and two-grid conditional left-sum audits. These finite diagnostics are not new Lean proofs, a uniform numerical error bound, a stock martingale theorem or an impossibility claim. Positivity uses a lower bound f/2 rather than the original exact f.

## 29 October 2026: interval-scaled protected corrections

`m2_rich_protected_variance.cpp/.py` preserve four conditional VIX-window integrals using piecewise C1 time potentials and bounded tanh/erf state potentials. Local analytic normalization bounds each generator correction by1 and the same L1 budget proves variance>=f/2. Endpoint zeros of both potential and time derivative justify joining interval scales without distributional jump terms. Numerical formula, endpoint, conditional quadrature, replay, independent stock recurrence, finite-difference generator and two-grid left-sum audits are in `../runs/m2_rich_protected_variance_001`. These checks are not new Lean proofs, a uniform pricing bound, exact market-forward calibration or a continuous-stock martingale theorem.

## 29 October 2026: existing-stress variance budget

`m2_available_variance.cpp/.py` apply the full generator to stress-scaled, Gaussian-damped potentials, including all time/stress terms. Analytic bounds give |h_j|<=c_j A_j and retain original floor f for |eta|<=1. Zero potentials/time derivatives at all fitted curve breaks avoid jump terms and retain four conditional VIX-window cancellations. Numerical formula, full-generator differences, analytic/GH conditional oracle, old replay, scalar stock recurrence and conditional discrete sums are in `../runs/m2_available_variance_001`. These are scoped numerical correspondence checks, not new Lean proofs, a uniform error bound, a continuous-stock martingale theorem or exact market-forward calibration.
