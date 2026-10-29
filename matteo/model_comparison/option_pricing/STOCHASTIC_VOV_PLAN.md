# Dedicated continuous sensitivity state: bounded experiment

Authorized 25 October 2026 following the user's approval of the two-cap,
two-strength, three-persistence comparison. No new production model is adopted.

## Equations and frozen objects

Under Q, add one independent scalar state, initially zero:

    dZ = -aZ Z dt + sqrt(2 aZ) dBZ,
    A(Z) = exp(s tanh Z),
    eta = cC C + qY (Y-betaY)^2 + qJ A(Z) (F-betaJ)^2,
    G = exp(2 h(eta)), M(t) = E_Q[G(t)],
    v = f + (xi(t)-f) G/M(t).

The first cap is the selected logistic cap, L=5.502780056785304, k=2.
The second adds .5 log(1+softplus(k(eta-L))/k) to that same h.
No cap-location or sharpness changes. The bounded A preserves the original
bounded-cap argument; for the logarithmic continuation, at-most-quadratic
variance growth in the enlarged Gaussian state gives the same nonexplosion
martingale argument as CAP_TAIL_MATHEMATICS.md. These are mathematical arguments,
not new Lean proofs or certificates of finite stock-payoff/estimator variance.

Strengths s=log(1.5), log(2); half-lives 7,30,90 calendar days,
aZ=log(2)*365/half_life. Z0=0 and zero Q mean are explicit new-layer assumptions;
the existing reservoir/risk-price scenario stays unchanged. No P calibration is
claimed. No new correlations or feedback in old state dynamics. s=0 nests baseline.

Freeze ALL 16 xi levels and floor .00264, old score coefficients, original OU
rates/geometric grids/state/correlations, premiums and quote conventions.
All 4096 observations at every strike; seven SPX/four VIX dates/four forwards;
original OTM quotes, cutoff3, common calls and native model-forward puts.
No jumps, quintic proposals, clipping, smoothing or hidden curve refits.

## Computation and verification

Reuse original full terminal OU states and likelihoods. Independently scrambled
one-dimensional Gaussian draws add Z_T in training, sharing draws across all
six strengths/persistences. Original likelihoods remain exact, with new Z density
ratio one. Full-scramble uncertainty remains based on eight 512-point nets.
This independent pairing is unbiased under scramble randomization; it is not a
claim of optimal combined-dimensional stratification or adequate rare-tail sampling.

Conditional VIX integrates over the original full Markov state AND Z_T.
Original conditional GH32 and normalizer GH48 remain; the extra Z integral uses
GH128 conditionally and GH192 unconditionally, after lower orders showed sampled
errors up to 1e-4. Check refinements at representative ordinary/extreme states.
Time panels/order8 and original stock512-step mesh stay fixed.

A Chebyshev interpolation of the original-state expectation as a function of
A on [.5,2] shares work among the six cases. Start32 nodes, adapt64/128 based
on six off-node comparisons with DIRECT finite Gaussian quadrature, including
extreme projected states. This is numerical interpolation of pricing integrals,
not a polynomial volatility model. Record every degree/error and retain failures.
Normalizer log-time cubic interpolation uses midpoint1e-7/independent2e-7 checks.

Stock volatility and terminal stocks MUST be recomputed. Stream exact independent
OU updates at the existing mesh, with deterministic counter-based Gaussian draws;
seeds depend on expiry/batch, not parameters, preserving common random numbers.
Use Z at the left endpoint; never use the next Z innovation in current volatility.
Freeze existing reference/hedge controls; their zero means survive independent
augmentation. Check streamed calculations against an independent Python replay.

## Selection, checking and interpretation

1. Price both zero-strength cap references and all12 fixed-curve layer cases for
   all four VIX dates. Save signed prices, IVs, forwards, weighted tails, second
   moments, slopes/curvatures, interpolation/quadrature diagnostics and seeds.
2. Before SPX or checking results, select minimum equal-maturity right-wing IV MAE
   (K/F>2) per cap among complete-wing cases; ties by native VIX+forward loss.
   Selection is descriptive, with no assumed aggregate/ATM/forward acceptance.
   Retain even poor winners. Take the union of chosen strength/half-life pairs
   and cross with both caps, plus two no-layer references: at most6 full reprices.
   This matched comparison separates layer and cap interactions. It does not
   search every grid point for the best joint fit.
3. Reprice all seven SPX dates for those cases, all4096 paths each. Native/common
   joint losses, bid/ask counts, signed errors, shape flags and paired block
   sensitivities are diagnostics. Do not adopt a candidate automatically.
4. After selection, generate32 new full-dimensional VIX scrambles using old
   proposal shifts augmented by one zero component for Z. Joint Gaussian draws
   and exact likelihoods are saved; no self-normalization. Recheck selected
   cases and both references without reselection. These are new numerical
   streams on the SAME in-sample quotes, not market out-of-sample evidence.
5. Independent saved-draw replay and source/seed checks, all-maturity figures,
   report and research-record updates. No global-fit/convergence/P claims.

Code: stochastic_vov_core.py, stochastic_vov_kernel.cpp,
stochastic_vov_study.py, test_stochastic_vov.py. Historical sources/runs preserved.
