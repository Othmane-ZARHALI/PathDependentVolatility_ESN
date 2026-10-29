# Joint SPX/VIX refit of the two-speed stress model

28 October 2026. The user authorized returning to the two-speed stress model and resuming SPX fitting.

**This first bounded joint fit is not satisfactory overall.** It substantially improves the sampled SPX fit, but loses most of the VIX smile fit. SPX numerical accuracy remains unresolved, so the result is not evidence that the entire two-speed model family cannot fit both markets.

## What was fitted

The model is unchanged:

$$v_t=f+(\xi(t)-f)[(1-w)R_t+w\{(1-r)A_t^s+rA_t^f\}],$$

$$A_t^i=\exp(b_iZ_t^i-\tfrac12b_i^2q_i(t)),\qquad dZ_t^i=-a_iZ_t^i\,dt+\sqrt{2a_i}\,dB_t^i.$$

The two stress innovations remain independent of each other and the old M2 bank/equity driver. The original M2 dynamics, initial state, risk premia, capped readout and positive spike shift remain fixed. This is a restricted joint calibration, not a refit of the seven original M2 readout coefficients.

All16 variance-curve levels are now free, including the nine intervals invisible to the VIX-only fit. Also free are the floor, both stress sensitivities, both ordered half-lives, the total stress weight and the fast fraction:23parameters. Bounds retain floor volatility5–12%, variance levels between the floor and0.2, sensitivities0–5 and half-lives1–365days. No new factor, correlation, jump, cap modification or expiry-specific rule was introduced.

The fit includes1294primary SPX quotes at29,57,119,211,302,666,1030days under the agreed market-only cutoff3 membership. All1751eligible SPX quotes, including457diagnostic-wing quotes, are repriced in the timestep check. All217native VIX quotes and all four forwards are retained. Original source splits are metadata; no market observations are held out.

The objective gives SPX and VIX families equal total weight, equal weight to their respective maturities, and bid–ask-scaled squared interval distances within each maturity. There is no extra VIX-price preservation constraint or new smoothing penalty. The four VIX forwards are equality constraints on the calibration estimator. Quote signs and unavailable IVs are preserved.

The first stage fits16curve levels and the floor for40SLSQP iterations. The second allows all23parameters for60iterations. Both stop at their iteration limits; their four forward residuals meet the eligibility threshold1e-5. The best eligible training candidate is selected before checking. There is one local staged search, not a global optimization or systematic multistart study.

## Price-fit trade-off

Smaller losses are better. Losses are scaled pricing errors, not monetary P&L or implied-volatility errors.

| Candidate | Training SPX loss | Checking SPX loss | Training VIX loss | Checking VIX loss |
|---|---:|---:|---:|---:|
| Original VIX-only parameters | 5371.218 | 5467.104 | 0.043630 | 0.050081 |
| Curve/floor joint stage | 2970.576 | 5134.115 | 95.952267 | 94.293497 |
| All23 joint stage | 67.689 | 618.742 | 38.147245 | 37.971099 |

| Candidate | SPX inside bid–ask | VIX inside bid–ask | Negative SPX estimates |
|---|---:|---:|---:|
| VIX-only | 9/1294 | 153/217 | 0 |
| Curve/floor joint | 9/1294 | 0/217 | 0 |
| All23 joint | 9/1294 | 20/217 | 20 |

The negative values are signed statistical estimates retained as unusable price estimates, not negative theoretical option prices. Nothing was clipped or removed.

## What the optimizer changed

| Parameter | VIX-only fit | Joint fit |
|---|---:|---:|
| Total stress weight $w$ | 0.922838 | 0.019435 |
| Old M2 weight $1-w$ | 0.077162 | 0.980565 |
| Fast fraction $r$ | 0.658093 | 1.000000 |
| Fast sensitivity | 3.535495 | 5.000000 |
| Slow sensitivity | 2.285440 | 0.009585 |
| Fast half-life | 71.390 days | 8.430 days |
| Slow half-life | 365.000 days | 292.624 days |
| Floor volatility | 12.000% | 5.000% |

In plain language, the VIX-only fit gets most of its variance response from the new stress factors. This joint search moves almost all that weight back to the original M2 response. Its remaining stress contribution is almost entirely fast; the slow branch is effectively inactive, so its fitted speed and sensitivity are not meaningfully identified. The floor hits its5%lower bound, the fast sensitivity hits5, and one curve level hits0.2. All curve levels are saved in the model artifact.

This is a finding about this bounded, sample-dependent fit with the original M2 readout frozen. It does not establish that SPX universally requires these weights, that stress independence is the unique cause, or that a joint fit is mathematically impossible. The original M2 and stress response remain defined in the same variance process for every maturity.

## All seven SPX maturities

| Days | VIX-only IV MAE | Joint IV MAE | Joint checking loss | Median joint SE / price scale |
|---:|---:|---:|---:|---:|
| 29 | 5.035 pp | 2.242 pp | 670.333 | 14.558 |
| 57 | 4.317 pp | 1.666 pp | 708.328 | 18.348 |
| 119 | 3.663 pp | 1.593 pp | 838.278 | 17.205 |
| 211 | 4.150 pp | 1.413 pp | 731.739 | 15.347 |
| 302 | 3.591 pp | 1.577 pp | 816.634 | 21.732 |
| 666 | 4.266 pp | 3.227 pp | 533.095 | 7.060 |
| 1030 | 3.423 pp | 1.253 pp | 32.784 | 2.914 |

These native IV averages omit only values for which Black inversion is undefined; all signed prices remain in the price loss. The CSVs retain every quote and flag negative estimates. The tiny bid–ask price scales make IV differences of a few points correspond to large interval losses.

## All four VIX maturities and forwards

| Days | VIX-only IV MAE | Joint IV MAE | Joint checking loss |
|---:|---:|---:|---:|
| 27 | 3.543 pp | 22.539 pp | 32.883 |
| 55 | 2.933 pp | 26.254 pp | 35.151 |
| 118 | 2.639 pp | 24.688 pp | 44.111 |
| 153 | 2.164 pp | 20.996 pp | 39.739 |

| Days | Target forward | Joint checking error ± block SE |
|---:|---:|---:|
| 27 | 20.1000 | +0.012528 ± 0.002980 |
| 55 | 20.7300 | -0.009374 ± 0.002938 |
| 118 | 21.1900 | +0.013672 ± 0.002625 |
| 153 | 21.7200 | +0.017822 ± 0.004549 |

The selected joint training maximum forward residual is 1.24e-08 VIX points. Checking maximum error is 0.017822. These constraints control sampled calibration means; exact population-forward matching is not certified. The much worse VIX option fit persists even though its checking forwards are close.

## Numerical accuracy and independent checks

The old exogenous M2 projections, equity increments and their full-path importance weights are reused from separate512-path-per-expiry search/check caches. The model initial state, architecture and premium match those caches exactly. Independent stress Gaussian innovations use separate declared seeds and are cached for common-random-number repricing. Every candidate rebuilds stress OU paths, variance and stock paths. The deterministic reference-stock control has a fixed known expectation; it is not normalized by a sample forward.

The initial calibration uses64steps per trading day. Both fine-grid points bracketing each curve knot are retained; any remaining knot delay is at most one original512-step interval (the sampled missing-knot gaps are around9–106seconds). VIX conditional/time integration and its validated multilevel estimator are retained. The nine non-window curve pieces correctly have zero VIX sensitivity.

The implementation audit verifies197source/input hashes, all56stress seed replays, independent NumPy stock reconstruction (maximum error1.33e-15), and all23analytic stock derivative directions. A separate scalar implementation reproduces the new coupled64/128/512stock checker (maximum error0). The positive spike shift remains0.369729972517.

**SPX accuracy gates fail:** none of1294primary SPX checking quotes passes the pre-existing sampling screen $2SE\le price\_scale/4$. This holds for the original and fitted models. For the joint candidate, median per-maturity SEs are roughly3–22price-scale units. Its training SPX loss67.7 versus checking618.7 is therefore not a stable calibrated-price result. Eight-block SEs exclude fitting variability and normalizer error.

As a read-only diagnostic, the existing frozen hedge-control estimator was compared on the same paths. It reduces SEs for several joint-fit maturities but does not close the quote-scale accuracy gap. It was not substituted into the fit or selected using checking prices. See control_variance_diagnostic.csv. Raw likelihood-weighted payoff hits, ESS, concentration and controlled-estimator variance concentration for every primary SPX quote are retained in SPX_weighted_tail_diagnostics.csv; zero sampled payoffs do not certify a zero theoretical price.

## Coupled stock timestep and additional stress-stream check

The selected model and original VIX-only model are repriced with new fine-grid stress innovations on the same stored old checking paths. Exact stress transitions are simulated on the512-step grid; coarse64/128states are the corresponding fine states, and all stock grids share the same equity Brownian path. This separates paired timestep differences from changes of stress sampling.

| Model | Steps/trading day | SPX primary loss | Inside bid–ask | Median SE / price scale |
|---|---:|---:|---:|---:|
| vix_only | 64 | 6363.928 | 12/1294 | 11.929 |
| vix_only | 128 | 6353.152 | 9/1294 | 11.830 |
| vix_only | 512 | 6381.396 | 11/1294 | 12.046 |
| selected | 64 | 643.546 | 5/1294 | 16.370 |
| selected | 128 | 637.987 | 27/1294 | 15.900 |
| selected | 512 | 606.332 | 27/1294 | 15.727 |

The finer-grid check preserves the broad SPX improvement and its sampling limitation. It does not turn the overall fit into a satisfactory result. The paired128→512changes can themselves exceed quote scales, so timestep precision is also unresolved. Full paired changes and block SEs are saved; these are sampled diagnostics, not uniform bounds. All1751eligible SPX quotes are present in this check, including the457outside the primary band.

## Conclusion and next decision

Retain the two-speed VIX-only candidate and this joint candidate separately. The tested joint compromise is unsatisfactory for VIX and is not numerically qualified for SPX. First improve SPX pricing precision before interpreting a new joint optimizer as a reliable market calibration. Further searches could also revisit the frozen original M2 readout or stress/equity dependence, but those are separate experiments; neither is implemented here. This run supplies no global infeasibility theorem and no new Lean certificate. Production and previously selected model artifacts are unchanged.

## Artifacts

- [Seven SPX smile panels](visualizations/two-speed-joint026/spx-joint-smiles.png)
- [Four VIX smile panels](visualizations/two-speed-joint026/vix-joint-smiles.png)
- [Variance-curve comparison](visualizations/two-speed-joint026/joint-variance-curves.png)
- [Checking forward errors](visualizations/two-speed-joint026/joint-forward-errors.png)
- [Experiment plan](TWO_SPEED_JOINT_PLAN.md)
- Run: `runs/m2_two_speed_joint_001`.
- Reusable SPX primitive/ratio cache: `runs/m2_two_speed_joint_cache_001`.
