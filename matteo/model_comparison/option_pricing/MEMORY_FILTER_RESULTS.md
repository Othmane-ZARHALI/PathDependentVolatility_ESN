# Linear memory filters of normalized variance signals — 29 October 2026

This experiment implements the linear-memory idea inspired by the two user-supplied Guyon slides. It tests M2-only filtering first, alongside filtering the whole multiplier. It does not implement the nonlinear GLE paper, feedback from actual returns, or a new Gaussian bank. The original fast M2 function and both stress factors remain. This report records a bounded numerical experiment, not a model-family impossibility or a new Lean certificate.

## Construction and interpretation

The frozen baseline is H0=(1-w)R+wA, where R is the original capped and Q-normalized M2 readout and A=(1-r)As+rAf is the existing two-speed lognormal stress mixture. Strictly positive spike shift and qJ are unchanged. Write b=1-w. For each normalized signal J in {R,As,Af}, add deterministic filtering of its realized path:

$$dM^J_t=a(J_t-M^J_t)dt,\qquad M^J_0=1,\qquad a=\log(2)/h.$$

$$M^J_t=e^{-at}+a\int_0^t e^{-a(t-s)}J_sds.$$

M2-only filtering uses

$$H^{M2}_t=b[(1-\chi)R_t+\chi M^R_t]+wA_t.$$

Whole-multiplier filtering uses

$$H^{all}_t=(1-\chi)H^0_t+\chi[bM^R_t+w((1-r)M^{As}_t+rM^{Af}_t)].$$

In both cases V=f+(xi(t)-f)H. Since the inputs and initial filters have mean one, E[H]=1 and E[V]=xi(t), assuming exact Q normalization. Positivity is preserved. The dated Monte Carlo implementation inherits the existing normalizer's approximation error. No new cap or sample normalization is introduced. Share0 recovers the baseline. Share chi controls how much current signal is replaced by its history; h is the effective filter half-life, not the original feedback-kernel decay constant in the slides.

This is a normalized-signal adaptation of linear hysteresis. Literal feedback on level variance with a time-varying xi would require a different deterministic input curve. Filtering a nonlinear signal also differs from passing a slower Gaussian state through a nonlinear readout: average(g(Z)) generally differs from g(average(Z)). Gaussian factors remain exogenous and reusable, but histories, stock paths and VIX exposures must be recomputed. Filters are ordinary path integrals; they add no Brownian motion or volatility feedback to the OU equations.

## Conditional VIX calculation

At expiry T the current OU state alone is insufficient: M_T is a realized random state. Retain that history. For u>=0,

$$E_T[M^J_{T+u}]=e^{-au}M^J_T+a\int_0^u e^{-a(u-s)}E_T[J_{T+s}]ds.$$

Future R means come from the existing full-state conditional Gaussian law, not from treating its three projections as Markov. Stress means are analytic. Integrate conditional filter means and the deterministic xi curve over the same 30-day VIX window. Four VIX forwards are constrained numerically in training; preserving E[V] alone is not a population-forward certificate or invariance of VIX smiles.

To obtain M_T jointly with each importance-sampled terminal state, simulate a Gaussian bridge history, then filter the actual nonlinear R and stress signals. Bridge residuals are sampled independently of the terminal state, with its exact target/proposal likelihood weight retained. There is no replacement of M_T by its conditional mean and no sample-weight normalization.

The terminal covariance has nearly redundant OU modes. Whitening and a relative eigenvalue cutoff1e-11 retain21,21,22,22 directions at the four expiries. This is a documented finite-precision Gaussian conditioning approximation. It is checked against the unconditional covariance of the implemented history, as well as cross-covariances and terminal projections. An earlier cutoff1e-14 run was rejected because its projected history covariance error reached .02146 despite tiny endpoint errors; it is preserved separately in `m2_memory_filter_001.conditioning_diagnostic001` and excluded from all tables and plots here. Initial preparation failures are retained separately. Calibration and checking were repeated after correction, using new checking seeds. The algebra did not change.

## Experiment

Generic half-lives 3, 10 and 30 days, and shares .25 and .5, applied uniformly across maturities; M2-only and whole filtering give12 nonzero settings. Each has a fixed-curve control and a fit of the same seven active variance levels under four training-forward equalities. Nine off-window levels, M2/stress parameters, stock theta=.8 and all Gaussian dynamics remain frozen. There are26 pilot stock configurations including baseline and baseline-refit controls. All 1,751 eligible SPX quotes (1,294 in the primary set), seven SPX maturities, 217 VIX quotes and four VIX forwards are retained.

The pilot uses 512 stock paths in eight independent antithetic batches, and 512 training VIX terminal states per expiry, augmented with independently simulated conditional histories. The terminal training stream is reused research data. Coupled historical grids with 64/128 steps per trading day and stock grids with 128/256 steps are saved. Every nonzero fit is attempted; unsuccessful or uncompetitive settings are retained.

Selection is frozen before new checking draws: training-forward gap<1e-5, training VIX loss<=alpha-zero-refit+.1, then smallest pilot fine-grid native-IV error among refits. The global selection is M2_3d_share0.25_refit. Additional checking compares baseline, M2 fixed, M2 refit, the same refitted curve with no filter, and a whole-filter refit. M2 candidate M2_3d_share0.25_refit; whole candidate whole_3d_share0.25_refit. If a scope has no feasible candidate, its lowest-VIX-loss diagnostic is retained with its feasibility status. The tolerance is an exploratory gate, not a satisfactory-fit certificate.

New checking uses 1,024 stock paths and 1,024 joint VIX terminal states plus independent conditional histories per expiry. All four training equalities are numerical sample statements; checking gaps below remain visible.

| Configuration | SPX interval loss | SPX IV error, pp | VIX interval loss | Max VIX forward gap, points | Invalid primary SPX IVs |
|---|---:|---:|---:|---:|---:|
| Baseline | 1903.59 | 3.238 | 0.050 | 0.0552 | 5 |
| M2 filter, fixed curve | 1899.76 | 3.236 | 0.052 | 0.0500 | 5 |
| M2 filter + curve refit | 1901.94 | 3.237 | 0.049 | 0.0076 | 5 |
| Curve refit alone | 1905.77 | 3.239 | 0.046 | 0.0108 | 5 |
| Whole filter + curve refit | 1912.77 | 3.241 | 0.051 | 0.0238 | 3 |

M2 refit SPX price-loss reduction versus baseline is 0.09 percent. IV error equally averages available native IVs across seven maturities; unavailable IVs and signed negative prices remain in the data. At 29 days, quoted right-minus-left IV gaps near strike/forward .9 and 1.1 are -5.57, -5.56, -5.56, -5.57, -5.29 percentage points for the five rows, versus -14.05 in the market. Both wing levels and skew must be assessed; aggregate improvement alone does not establish a satisfactory fit.

Approximate paired 1,024-resample eight-batch bootstrap intervals for candidate-minus-baseline price loss are [{"candidate": 1, "loss_delta": -3.8323652867031797, "approximate_batch_bootstrap_CI95": [-20.471795931506954, 12.830074381970723]}, {"candidate": 2, "loss_delta": -1.6527795849751783, "approximate_batch_bootstrap_CI95": [-18.280641196284982, 15.918271836982116]}, {"candidate": 3, "loss_delta": 2.176509747429691, "approximate_batch_bootstrap_CI95": [0.23315609320969208, 3.7702190281389787]}, {"candidate": 4, "loss_delta": 9.184866288836929, "approximate_batch_bootstrap_CI95": [-93.96125695762511, 91.4845001067481]}]. Candidate and baseline share expensive draws; cheap reference batches are resampled independently and shared across candidates. These intervals exclude deterministic errors, fitting variability and market uncertainty. The existing 65,536-path stress-only auxiliary reference is reused, with valid conditional-lognormal call prices and fixed deterministic delta controls. Full checking is therefore not wholly independent of previous research even though the target draws are new; no market out-of-sample claim is made.

## Numerical audits and limits

Independent 97-step scalar stock recurrence error 0.0; full-horizon zero-share replay error 0.0; independent scalar bridge-history/filter recurrence error 2.220446049250313e-15; constant-input identity error 0.0. These are implementation checks, not formal simulator certification. No new Lean theorem is claimed.

Corrected bridge probe maximum unconditional projected covariance error 5.52e-07, maximum actual-terminal projection error 9.49e-09, and stable cutoff1e-10vs1e-11 matched-history change 0. The bridge's dropped covariance directions and residual cross-covariance errors remain in `bridge_audit.json`; this is not exact full-rank conditioning or a uniform tail bound.

Historical 64/128 grid comparison: max VIX change 116.4 points, max absolute weighted-mean change 9.021e-05. Future conditional GH32/time8/32lag knots versus GH48/time12/64lag knots on 64 matched states per expiry: max VIX change 32.64, max absolute weighted-mean change 1.283e-05. Coupled stock 128/256 median price changes in quote-scale units are 0.127, 0.097, 0.099, 0.126, 0.098. These are probes, not rigorous error bounds.


Large unweighted maxima above occur in extreme importance-sampled states and must be interpreted with their likelihood weights. Across these declared refinement probes, the maximum likelihood-weighted absolute VIX change is 0.000955629 points. The maximum absolute native option-price change is 0.0012763 quote scales. Future quadrature changes include the paired baseline control; history changes use all 1,024 checking histories, while future refinement uses 64 matched states per expiry. No extreme states are clipped. Detailed paired payoff changes and SEs are in `fresh_validation/weighted_numerical_payoff_audit.json`; these are finite probes, not uniform tail bounds.

Raw stock payoff concentration, importance-weighted VIX tail concentration, likelihood-weight sample mass, memory means/SEs, stock means/SEs, signed prices, primary membership, invalid IVs, and separate anchor/difference forward SEs are saved. No clipping, smoothing, sample-forward renormalization or hidden quote deletion is used. Rare-history sampling can remain uncertain even when the terminal covariance and deterministic quadrature agree.

## Assessment

The M2 filter produces no material improvement in this checking sample. The missing SPX skew remains. Whole filtering and curve adjustments must be read together with the VIX costs and paired uncertainty above; no satisfactory double calibration has been demonstrated.

The mechanism is persistent averaging. In particular dM has no direct Brownian shock, so replacing a fraction of a signal by its filtered history reduces that signal's instantaneous Brownian loading. Longer persistence does not automatically strengthen short-maturity downside leverage. The numerical result is specific to these positive linear filters, shares, memories and frozen base parameters; it does not exclude nonlinear hysteresis, other readouts, or exogenous Gaussian models in general.

No candidate is adopted. The production model and previous experiments remain unchanged. Do not repeat this same filter grid as an untested hypothesis. Figures: [seven SPX smiles](visualizations/memory-filter029/SPX-memory-filter.png), [four VIX smiles](visualizations/memory-filter029/VIX-memory-filter.png). Raw run: `runs/m2_memory_filter_001`.
