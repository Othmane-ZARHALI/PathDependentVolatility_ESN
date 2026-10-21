# Pilot-designed VIX sampling — 21 October 2026

## Scope and question

Keep all four VIX maturities (27/55/118/153 days), 175 training option quotes
and four forwards. SPX scope is unchanged. At the same frozen M2 model, state,
risk premium and selected variance curve as the preceding experiment, does a
proposal designed from independent pilot data improve fresh sampling precision?
This is a numerical diagnostic, not another model calibration.

Run: `runs/m2_vix_pilot_sampling_001`. Sources and proposals are frozen before
any evaluation draws. Each method uses eight independent scrambles with 128
points each per maturity, with nested 64-point prefixes. There are 8192 base
conditional-state evaluations, plus twelve selected-state GH refinement checks.
Compare ordinary Sobol against the pilot mixture on fresh, independent seed
streams. Prior untuned-mixture results are historical comparisons, not a fresh
validation baseline. No held-out market quotes enter this experiment.

## What the proposal learns

The 512 mixture draws per maturity from `m2_vix_importance_sampling_001` are
explicitly promoted to pilot data. Regenerate their whitened normals from saved
seeds, and verify their exact likelihood weights against the saved observations.
Use their VIX values without any new conditional-pricing approximation. This
reuse is not new validation evidence; every checking draw is generated afresh.

Three nonnegative pilot targets are the means of payoff / market price scale
for (a) puts, (b) calls and (c) calls in the upper strike quartile. For each,
estimate its likelihood-corrected payoff-weighted mean and covariance in the
first eight deterministic PCA coordinates. The first moment alone could miss
a symmetric two-sided tail, so retain its largest covariance direction too.

For pilot weights w_i and nonnegative target h_i, let

    p_i = w_i h_i / sum_j w_j h_j,
    a = sum_i p_i z_i,
    C = sum_i p_i (z_i-a)(z_i-a)^T.

Take two proposal centers a +/- sqrt(max(lambda_max(C)-1, 0.25)) v_max(C),
cap their Euclidean norms at 3, and set every unused coordinate to zero. This
is a declared heuristic; it is neither a price fit nor an optimal proposal.
Normalize pilot masses only for learning these moments. Evaluation weights
are never self-normalized.

The final mixture retains probability 0.5 at the original Gaussian and gives
each of the six learned centers probability 1/12. Every checking payoff,
forward and second moment uses the exact original-density / mixture-density
weight, bounded above by 2. No forward rescaling is applied. All four proposals
are serialized and hashed before the checking loop; the hash is checked again
at completion. Raw pilot arrays and fresh weighted observations are saved.

## Pilot limitations, recorded before evaluation

| Days | Positive far-call pilot contributions | Contribution ESS | Largest share |
| --- | ---: | ---: | ---: |
| 27 | 29 | 11.54 | 19.6% |
| 55 | 11 | 5.42 | 29.0% |
| 118 | 12 | 5.87 | 30.5% |
| 153 | 13 | 4.65 | 38.0% |

Those far-call summaries are based on few influential observations. In contrast,
all-put target contribution ESS is about 150–181 and all-call ESS about 62–79.
The pilot proposal may therefore overfit a tail observation. The first-eight-
coordinate restriction and norm cap reduce proposal complexity but can miss
relevant directions. Fresh evaluation must retain negative results. Only one
pilot is fitted; variability across independent fitted pilots is not measured.

## Error accounting

Standard errors use eight complete independent scrambles, conditional on the
frozen proposal. Prefix changes use paired scramble differences. Neither point
independence within a Sobol net nor independence of nested prefixes is assumed.
Report payoff contribution ESS and maximum within-scramble and pooled shares;
these are diagnostics, not accuracy certificates. Zero observed payoffs or SE
are not evidence of an accurately priced rare-event expectation.

Conditional GH order 32, date-zero normalizer order 48 and Gaussian time order
8 remain fixed. The first pilot-mixture scramble at each maturity receives a
conditional 32/48 comparison at its minimum, median-rank and maximum VIX state.
This selected-state diagnostic cannot bound full-price quadrature error.
Time integration and normalizer biases remain separate unresolved gates.
Exact terminal Gaussian OU sampling requires no SPX stock time step. Existing
market moment conflicts remain flagged regardless of sampling results.

Implementation checks distinguish a symmetric quadratic tail (whose weighted
mean alone is zero), payoff-scale invariance, norm caps, density correction and
zero-pilot-mass rejection. Current tests: 90 package plus 11 driver/legacy tests.

## Predeclared fixed-candidate follow-up

Once the full fresh run is complete, compare the original starting and selected
four-maturity curves on these same weighted draws, without fitting or choosing
new candidates. Both curves were fixed before this experiment. Because each of
these four VIX windows lies wholly within one curve segment and the date-zero
normalizer is independent of the variance curve, the alternative VIX value is
exactly reconstructed at the same quadrature approximation from

    VIX_new^2 = 10000 floor
                + (xi_new-floor)/(xi_old-floor) (VIX_old^2-10000 floor).

The audit refuses a window crossing a curve knot and independently checks this
identity against the original full conditional-exposure cache at every expiry.
The proposal remains unchanged. Measure the nonlinear loss change with paired
leave-one-scramble-out sensitivity; this is descriptive, not a confidence bound
or final ranking. This follow-up avoids additional conditional-state simulations.

## Exploratory control added after the first maturity

The completed 27-day comparison showed reduced option SE and concentration,
but worse forward SE. After inspecting that block, add a separate diagnostic:
for a payoff H and a center c learned solely from old independent pilot data,
use c + mean[w(H-c)] = mean[wH] - c(mean[w]-1). Since E_g[w]=1, this preserves
the expectation without a stochastic ratio. It is not self-normalization or
substitution of the market forward. Center each option, forward and second
moment separately using its old pilot estimate. Keep every raw result.

This addition was prompted by a fresh checking result. Its comparisons are
therefore explicitly exploratory, rather than validation of an entirely
predeclared estimation method. No center uses checking payoffs. Do not clip
negative controlled prices or interpret a small observed SE as tail convergence.
`vix_mass_control` saves this diagnostic separately from the main study.

## Completed results

The run completed in 1068.5 seconds (17.8 minutes). The independent artifact
check regenerated all likelihood weights, verified 35 frozen source hashes,
confirmed 64 unique checking seeds disjoint from the old pilot streams, and
recomputed all 732 summary rows and paired prefix changes from 8192 saved
observations. These are correspondence checks, not accuracy certificates.

| Days | Median option SE / scale: Sobol | Pilot mixture | Median option variance ratio: Sobol / mixture | Zero-hit options: Sobol / mixture |
| --- | ---: | ---: | ---: | ---: |
| 27 | 1.281 | 0.871 | 2.59 | 0 / 0 |
| 55 | 1.386 | 0.699 | 1.93 | 5 / 0 |
| 118 | 0.878 | 0.624 | 1.75 | 4 / 0 |
| 153 | 0.622 | 0.302 | 2.42 | 6 / 0 |

Variance ratios exclude options with zero observed SE under either method.
A zero-hit option is unresolved rather than perfectly priced. The proposal
improves median option SE at all four maturities, but only 27/175 estimates
meet the observed SE / scale <= 0.125 threshold. Those 27 are not certified
without tail and deterministic error checks; the full sampling gate fails.

The largest pooled single-observation share falls to 5.61%, 12.15%, 8.47%
and 16.59% at 27/55/118/153 days, versus 71.54%, 100%, 100% and 100%
under ordinary Sobol. Minimum payoff-contribution ESS becomes 47.65, 23.78,
23.67 and 13.55, versus 1.69, 1, 1 and 1. These diagnostics indicate less
concentration, not iid equivalent sample counts or error guarantees.

Raw forward estimates are noisier under the learned mixture at every expiry:

| Days | Sobol forward ± scramble SE | Pilot mixture forward ± scramble SE |
| --- | ---: | ---: |
| 27 | 20.73109 ± 0.06267 | 21.00368 ± 0.12669 |
| 55 | 20.43127 ± 0.09099 | 20.53520 ± 0.15442 |
| 118 | 21.22352 ± 0.05303 | 21.18901 ± 0.09911 |
| 153 | 21.12030 ± 0.06922 | 21.00893 ± 0.18754 |

The learned mixture is not a universal replacement for ordinary Sobol. The
forward 64→128 prefix changes have paired SEs large relative to their changes;
this budget does not establish stable forward convergence. Per-scramble weight
means range from 0.97043 to 1.02865; the largest weight is 1.84942 (bound 2).

### Separate exploratory mass control

`runs/m2_vix_mass_control_001` preserves the raw estimator and centers each
payoff using the old independent pilot only. Forward SEs become 0.09673,
0.15864, 0.12339 and 0.14997. This improves 27/153 relative to the raw mixture,
worsens 55/118, and remains worse than ordinary Sobol at all four expiries.
Ten of 175 controlled option means at 128 points are negative; they are saved
without clipping and are unsuitable as accepted prices or inputs to the
nonnegative-price calibration objective. Thus simple mass centering does not
resolve the trade-off, and it is not adopted as the default estimator.

### Independent fixed-candidate curve ranking

`runs/m2_vix_pilot_rank_audit_001` applies the original common objective,
including all four forwards and the original curve penalty. Negative selected
minus starting loss would favor the old selection. Both new methods give a
positive change instead:

| Method | Starting loss | Selected loss | Change | Descriptive jackknife SE | Leave-one-scramble-out change range |
| --- | ---: | ---: | ---: | ---: | ---: |
| Sobol | 63.98618 | 66.64123 | +2.65505 | 1.31297 | +1.76048 to +3.41010 |
| Pilot mixture | 63.25632 | 66.13988 | +2.88356 | 1.41654 | +2.44454 to +4.18613 |

The same point ranking under two fresh sampling schemes, with positive changes
under every one-scramble omission, adds evidence against the reliability of
the original training gain. It does not establish a statistically certified
population ranking. The selected curve improves the forward loss but worsens
option loss enough to worsen the total. Neither curve is refitted here.

### Conditional integration and normalization remain separate

At three states from the first learned-mixture scramble, the maximum conditional
GH 32→48 changes are 0.000599, 0.003380, 0.001884 and 0.003077 VIX points.
The selected extreme states reach VIX values of about 346–503, so these checks
include some proposal-driven tail states. They do not bound full-price errors,
and no new time-order or global-normalizer convergence claim follows.

The exact normalized-model second moment for these single-segment windows is
542.49611 (27/55d) and 600.31899 (118/153d). Mixture estimates are respectively
548.19 ± 4.08, 546.18 ± 9.21, 595.09 ± 5.27 and 593.00 ± 7.21. Their agreement
at this noisy resolution does not exclude normalizer or conditional quadrature
bias. No moment identity is imposed on the prices by rescaling.

## Decision and next bounded work

Keep all four maturities and forwards. Preserve both ordinary Sobol and the
learned proposal: the former currently has better forward precision, the latter
better option-tail coverage. Do not combine them by selecting whichever gives
the most attractive fitted loss, and do not start a new optimizer against
unresolved noisy prices.

A useful next numerical comparison is component-stratified defensive sampling
with fixed allocations and exact mixture weights, designed before new draws.
This tests whether controlling mixture allocation can retain option-tail gains
without degrading central expectations. Pilot choices and controls need a new
independent check, with the previous findings treated as design evidence.
Continue quote-relative conditional/time/normalizer checks and the separate SPX
hedging, tail sampling and coupled-grid work. No joint-fit or physical-feasibility
claim has been established.
