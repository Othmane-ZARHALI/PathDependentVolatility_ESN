# Full-maturity VIX importance-sampling diagnostic — 21 October 2026

The latest user instruction restores all four VIX maturities, 27/55/118/153
days, all 175 training options and their four forwards. This supersedes the
earlier same-day near-only scope. No market observations are removed or
reclassified as held out. SPX scope is unchanged.

## Design and interpretation

The experiment fixes the existing Stage-3 model, zero initial projected state,
risk premium and the original four-maturity fitted curve
(0.0542496112, 0.0600318990, 0.0642080964). It does not optimize parameters.
The machine-readable run supplies full-precision inputs and frozen sources.

Write the exact terminal OU state as X = mean + L z, with z standard normal
and every covariance-factor column retained. Ordinary scrambled Sobol draws
z from its original Gaussian density f. The defensive proposal is

    g(z) = 0.5 f(z) + 0.125 [f(z-2e1)+f(z+2e1)+f(z-2e2)+f(z+2e2)].

The first two coordinates are ordered by deterministic decayed-projection
variance contribution, without tuning to these evaluation payoffs. A further
Sobol coordinate selects the mixture component. The exact likelihood weight is

    w(z) = f(z)/g(z)
         = 1 / sum_k pi_k exp(m_k dot z - ||m_k||^2/2).

Each forward, option payoff and second moment is averaged after multiplication
by w. No sample-weight normalization or forward rescaling is applied. The
unshifted half implies 0 < w <= 2. This bounds likelihood weights, not payoffs
or pricing error. Leading Gaussian modes may fail to cover nonlinear VIX tails.

Both methods evaluate 64 conditional states in each of eight independently
scrambled replications, at each maturity, with nested 32-point prefixes.
There are 4096 conditional-state evaluations overall. Standard errors use
complete independent scrambles, rather than treating Sobol points as iid.
The comparison matches conditional-pricing counts; mixture generation adds a
small proposal overhead. Eight-scramble variance comparisons are exploratory.

Conditional Gauss-Hermite order 32, global normalizer order 48 and panelwise
Gaussian time order 8 remain fixed across methods. Thus this experiment isolates
outer sampling behavior at those approximations. It neither clears their biases
nor replaces the separate refinement studies. VIX terminal-state sampling uses
the exact Gaussian OU law; there is no SPX stock time step in this experiment.

Implementation checks separately recover Gaussian mass, first/second moments
and a three-standard-deviation call expectation with a known analytic answer.
Weight ESS and payoff concentration are diagnostics, not acceptance tests.

Reproduce from the repository root (with the package dependencies installed):

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m matteo.model_comparison.option_pricing.vix_importance_sampling --output /path/to/new/run
```

The driver refuses an existing output directory, freezes its sources before
computation and rejects source changes during execution. Its default scope
check requires the restored four maturities; it cannot inherit the former
near-only restriction silently.

## Results

Run: `runs/m2_vix_importance_sampling_001` (454.8 seconds). Frozen code, seeds,
raw VIX values and exact weights are saved alongside the summaries.

| Days | Sobol forward ± SE | Mixture forward ± SE | Median option variance ratio, Sobol / mixture | Zero-hit options, Sobol / mixture |
| --- | ---: | ---: | ---: | ---: |
| 27 | 20.67454 ± 0.05950 | 20.96321 ± 0.12930 | 0.448 | 5 / 0 |
| 55 | 20.51722 ± 0.24754 | 20.52601 ± 0.10266 | 3.663 | 2 / 0 |
| 118 | 21.21253 ± 0.17144 | 21.14958 ± 0.08428 | 2.798 | 2 / 0 |
| 153 | 21.23292 ± 0.11219 | 21.31299 ± 0.16393 | 1.620 | 8 / 0 |

Ratios exclude contracts with zero observed SE under either method; those
are unresolved, rather than perfectly priced. Ratios above one favor the
mixture. The eight-scramble estimates are noisy and do not certify superiority.
The proposal improves many 55/118-day estimates but is not uniformly better.
For example, the 27-day forward SE increases from 0.0595 to 0.1293 and the
153-day forward SE increases from 0.1122 to 0.1639.

Median option SE / market scale, Sobol → mixture, is respectively
1.050 → 1.550 (27d), 4.026 → 1.799 (55d), 2.370 → 0.713 (118d),
and 0.770 → 0.825 (153d). The sampling target is at most 0.125.
Only six mixture option estimates satisfy that observed SE threshold, and
even those are not accepted without tail and bias checks. No full-maturity
sampling gate passes.

Zero-hit contracts fall from 17/175 to 0/175 at 8x64 points per expiry.
This improves coverage, but the maximum single-observation payoff share
within a scramble still reaches 100% at every maturity for both methods.
Observing a payoff once does not establish its expectation accurately.

Per-scramble mixture weight means range from 0.97324 to 1.03536; the
largest observed weight is 1.76154, below the analytic upper bound 2.
These checks describe weight behavior and cannot establish payoff accuracy.

## Next decision supported by this pilot

Keep all four maturities in subsequent analyses. Do not replace ordinary
Sobol with this untuned mixture as a universal default, or refit against its
small noisy sample. The next sampling experiment should freeze a proposal
using independent pilot information about nonlinear payoff-tail directions,
then test it on fresh complete scrambles and larger nested budgets at all
four expiries. Preserve unshifted defensive mass and exact weights.

SPX hedging/importance sampling and coupled grid corrections remain separate
work items. Conditional GH, time integration and global normalizer refinement
must accompany any later larger VIX run; this outer-only experiment does not
clear those existing error gates. No final held-out market quotes were used.
