# ESN and M2 experiment history and learnings

Research through 29 October 2026, ordered by the progression of questions. Early physical-measure stages are ordered by experiment stage rather than inferred calendar dates. This record retains both progress and failed or partial results. The [complete catalogue](REPORT_CATALOG.md) links the individual numerical studies and plans, including studies grouped together here.

Price losses below are comparable within their declared experiment, not automatically across sections. Different scopes, price estimators, numerical grids, auxiliary anchors and fitted parameters can change their scale. Fresh simulation checks generally reuse the same market snapshot; some also reuse auxiliary samples. Source reports specify those dependencies.

## 1 Legacy ESN A2 and historical comparison

**Question and setup.** Establish the executable ESN reference and compare an M2 synthetic anchor with the saved historical nine-asset A2 fit.

**Result and learning.** The 30 September diagnostics compare a fitted historical A2 model with an unfitted M2 anchor, so they are not a fair final family ranking. Raw versus used variance and market-volatility proxies affect physical-statistic comparisons. Audits found starting-vector, residual-length, fitted-alpha and serializer inconsistencies in the checked-in entry point; these do not establish how every archived fit was produced.

**Consequence.** Keep architecture, implementation provenance and calibration status separate. The later externally described October model is another version, with endogenous feedback; do not infer its implementation or full path reuse from legacy A2 or M2.

Sources: [historical comparison](../../matteo/model_comparison/othmane_replication/RESULTS.md), [code issues](../../matteo/model_comparison/othmane_replication/OTHMANE_APPROACH_ISSUES.md), [legacy ESN pricing audit](../../matteo/model_comparison/option_pricing/AUDIT.md). The historical comparison includes offline interactive figures.

## 2 Physical capacity from pilot through Stage 3

**Question and setup.** Can exponential-quadratic exogenous M2 jointly reproduce the declared stylised-fact proxies on synthetic paths, and do optional orthogonal-energy or echo channels help?

**Result.** The aggressive Stage-2 curvature box had zero strict joint acceptance. Slower spike memory, slow-mode weighting and a low-curvature/high-shift feedback ridge led to 62.0% M2 and 64.7% M2+O joint acceptance in untouched Stage 3. Echo benefit was negligible; Stage 3 disabled it in the recommended M2 baseline.

**Learning and consequence.** Allocate variance across timescales rather than simply raising every curvature. Orthogonal energy is an optional roughness tradeoff. These are synthetic physical-measure proxy results, not an empirical SPX fit or a Q option calibration.

Sources: [stage comparison](../../matteo/exogenous_reservoir_eq_m2/experiments/EXPERIMENT_COMPARISON.md), [Stage 2](../../matteo/exogenous_reservoir_eq_m2/experiments/stage2_001/RESULTS.md), [Stage 3](../../matteo/exogenous_reservoir_eq_m2/experiments/stage3_001/RESULTS.md).

## 3 Positive asymmetry and finer physical grids

**Question and setup.** Strengthen downside asymmetry with positive spike shift, then test frequency robustness using coupled grids rather than changing the model at each resolution.

**Result.** Positive shifts strengthened downside response but could weaken rank-Zumbach acceptance. The nested 2-to-64-step study showed that exact OU transitions did not remove errors in nonlinear variance and stock integration. The Stage-4 primary passed four search seeds but only two of five untouched validation seeds under the unchanged seven-gate protocol.

**Learning and consequence.** Preserve the later strict positive-shift constraint, while retaining its measured tradeoffs. Stage-3 coarse-grid success is historical evidence, not fine-grid validation. Stage-4 failed the declared validation; Pearson asymmetry cannot substitute for the rank-based gate.

Sources: [asymmetry](../../matteo/exogenous_reservoir_eq_m2/experiments/asymmetry_001/RESULTS.md), [time-grid audit](../../matteo/exogenous_reservoir_eq_m2/experiments/timestep_convergence_001/RESULTS.md), [Stage 4](../../matteo/exogenous_reservoir_eq_m2/experiments/stage4_fine_grid_001/RESULTS.md), [Pearson/rank discrepancy](../../matteo/exogenous_reservoir_eq_m2/experiments/zumbach_discrepancy_001/RESULTS.md).

## 4 Empirical physical statistics and uncertainty

**Question and setup.** Determine how stable SPX physical targets are across disjoint five-year periods, using moving-block bootstrap intervals and block-length sensitivity.

**Result and learning.** Roughness, kurtosis and asymmetry vary across periods; Zumbach intervals are often broad. A magnitude-sensitive Pearson statistic and a rank statistic answer different questions. Sampling uncertainty and volatility-proxy definitions belong with a proposed calibration target.

**Consequence.** Do not turn synthetic acceptance bands or one historical point estimate into a universal market requirement.

Source: [physical stability study](../../matteo/exogenous_reservoir_eq_m2/experiments/physical_stability_001/RESULTS.md).

## 5 First option pricers and a fit that reversed on checking

**Question and setup.** Build dated M2 SPX pricing, conditional full-state VIX pricing, common price losses and exploratory variance-curve fitting.

**Result and learning.** A fixed-sample VIX curve improvement reversed on fresh checking scrambles. Cheap optimizer improvement and narrow internal errors were insufficient to rank expected prices. The full factor state, rather than a few non-Markov projections, matters for conditional VIX.

**Consequence.** Preserve the failed fit, matched-state refinements and independent pilot design. Keep a candidate's model change separate from estimator changes.

Sources: [initial SPX](../../matteo/model_comparison/option_pricing/M2_SPX_RESULTS.md), [initial VIX](../../matteo/model_comparison/option_pricing/M2_VIX_RESULTS.md), [curve profile](../../matteo/model_comparison/option_pricing/CURVE_PROFILE_RESULTS.md), [reversal and numerical diagnosis](../../matteo/model_comparison/option_pricing/NUMERICAL_STABILITY_RESULTS.md).

## 6 VIX tails and second-moment diagnostics

**Question and setup.** Investigate the reversal and moment inconsistency using direct Gaussian expectations, conditional quadrature refinement, defensive importance sampling, component allocation and independent-pilot controls.

**Result and learning.** Direct Gaussian checks supported the normalized-model second moment at the frozen model; rare outer-state contributions remained sensitive. Score proposals, payoff-targeted pilots, stratification and Gaussian controls helped selected quantities but did not establish full-surface accuracy. Zero payoff hits and zero empirical SE do not certify a tiny price.

**Consequence.** Retain weighted payoff concentration, ESS, likelihood mass, pilot stability and per-maturity SEs. A passed forward screen is not a convergence certificate for options or the entire model.

Sources: [moment diagnosis](../../matteo/model_comparison/option_pricing/VIX_MOMENT_DIAGNOSIS.md), [pilot sampling](../../matteo/model_comparison/option_pricing/VIX_PILOT_SAMPLING_RESULTS.md), [targeted pilots](../../matteo/model_comparison/option_pricing/VIX_TARGETED_PILOT_RESULTS.md), [allocation refinement](../../matteo/model_comparison/option_pricing/VIX_ALLOCATION_REFINEMENT_RESULTS.md), [Gaussian-control stability](../../matteo/model_comparison/option_pricing/VIX_GAUSSIAN_STABILITY_RESULTS.md). Additional stages remain in the catalogue.

## 7 SPX sampling controls and stock integration

**Question and setup.** Improve pricing efficiency with path importance sampling, volatility tilts, adapted controls, multilevel allocation, Gaussian bridge reflection and full/linear OU-innovation dictionaries.

**Result and learning.** Some direct-price variance reductions were consumed by runtime. The linear-only follow-up retained about 1.14x central variance-cost gain, while tail concentration, pilot instability and failed quote-scale screens prevented general adoption. Reflection did not pay for its runtime. Conditional integration of independent equity noise offered limited benefit because that residual variance was very small. Negative adjusted estimates remained numerical failures to record, not prices to clip.

**Consequence.** Compare variance per unit cost, freeze pilots before checks, and retain unsuccessful tails and grid corrections. Representative near-ATM diagnostics do not qualify the whole surface.

Sources: [path importance](../../matteo/model_comparison/option_pricing/SPX_PATH_IMPORTANCE_RESULTS.md), [volatility sampling](../../matteo/model_comparison/option_pricing/SPX_VOLATILITY_SAMPLING_RESULTS.md), [MLMC](../../matteo/model_comparison/option_pricing/SPX_MLMC_RESULTS.md), [bridge coupling](../../matteo/model_comparison/option_pricing/SPX_BRIDGE_COUPLING_RESULTS.md), [innovation controls](../../matteo/model_comparison/option_pricing/SPX_INNOVATION_CONTROL_RESULTS.md), [linear-only controls](../../matteo/model_comparison/option_pricing/SPX_LINEAR_INNOVATION_RESULTS.md), [conditional controls](../../matteo/model_comparison/option_pricing/SPX_CONDITIONAL_CONTROL_RESULTS.md).

## 8 Scope expansion and renewed joint fitting

**Question and setup.** Move from restricted training diagnostics to actual exploratory fitting, state/curve refinement and all retained SPX maturities under cutoff 3.

**Result and learning.** Restricted state-plus-three-piece fitting helped; five pieces added small inconsistent gains. The all-maturity study expanded the fitting population rather than validating the old near-only fit. Later feasibility work uses all eligible original splits as calibration inputs.

**Consequence.** The current seven-SPX/four-VIX scope supersedes older training-only counts. Preserve the old scope with each old result, and do not call new simulations market holdout validation.

Sources: [reentry decision and subsequent override](../../matteo/model_comparison/option_pricing/FITTING_REENTRY_DECISION.md), [state refinement](../../matteo/model_comparison/option_pricing/JOINT_STATE_REFINEMENT_RESULTS.md), [all maturities](../../matteo/model_comparison/option_pricing/ALL_MATURITY_FIT_RESULTS.md), [active scope](../../matteo/model_comparison/option_pricing/active_scope.json).

## 9 Original readout fitting and maturity-aligned variance curve

**Question and setup.** Fit eleven original readout/curve parameters with fixed reservoirs, then test a maturity-aligned 16-piece curve/floor search while enforcing positive spike shift.

**Result.** The eleven-parameter search improved aggregate SPX fitting but incurred VIX and maturity-specific costs. The later accepted maturity-curve proposal changed curve/floor; original-OTM VIX-plus-forward training loss moved 35.43 to 27.17, while the VIX right wing remained too shallow.

**Learning and consequence.** Level improvements do not supply independent smile curvature. Preserve exactly what the accepted proposal changed, rather than treating every searched coordinate as a successful readout refit.

Sources: [eleven-parameter fit](../../matteo/model_comparison/option_pricing/FULL_READOUT_FIT_RESULTS.md), [maturity curve](../../matteo/model_comparison/option_pricing/MATURITY_CURVE_FIT_RESULTS.md). [Visual comparison](VISUAL_EVIDENCE.md#curve-and-cap).

## 10 Cap and stochastic response experiments

**Question and setup.** Diagnose cap compression, then compare gentler capped responses with limited curve adjustments, inside-score stochastic sensitivity and independent after-cap modulation.

**Result and learning.** Gentler responses could raise VIX wings but move the centre and harm SPX. Fresh aggregate cap gains were smaller than one reported paired SE. Inside-score modulation changed the full normalizer and could flatten measured right-wing slopes despite bringing levels closer. Placement, functional form and persistence changed together in some comparisons.

**Consequence.** Keep normalized comparisons, weighted tail diagnostics and SPX costs together. Cap influence is established more strongly than a cap-only causal explanation or remedy.

Sources: [cap-tail](../../matteo/model_comparison/option_pricing/CAP_TAIL_RESULTS.md), [controlled cap](../../matteo/model_comparison/option_pricing/CAP_CONTROLLED_RESULTS.md), [inside-score response](../../matteo/model_comparison/option_pricing/STOCHASTIC_VOV_RESULTS.md), [post-cap response](../../matteo/model_comparison/option_pricing/POSTCAP_VOV_RESULTS.md). Figures 1–2 in [the gallery](VISUAL_EVIDENCE.md).

## 11 Frozen variance budget and mixed-response mechanism

**Question and setup.** Separate insufficient VIX second-moment budget from inadequate smile shape, then compare curve-only, single-response, product and additive mixtures.

**Result.** Static quote-payoff moment certificates excluded the old frozen budgets at 55/118/153 days. A normalized marginal mixed-response benchmark demonstrated additional wing curvature at matched moments, without being an M2 fit. In the later M2 comparison, additive response fresh loss 10.2960 was below product response 18.2300 and curve-only 20.5171, with residual errors.

**Learning and consequence.** Free enough shared variance budget before interpreting a shape failure. The marginal mechanism, scoped Lean statements and implemented pricing experiment have separate evidential roles.

Sources: [wing findings](../../matteo/model_comparison/option_pricing/WING_MECHANISM_FINDINGS.md), [mixed response](../../matteo/model_comparison/option_pricing/MIXED_RESPONSE_RESULTS.md), [formal correspondence](../../matteo/model_comparison/option_pricing/formalization/CORRESPONDENCE.md).

## 12 Forward constraints and fitted floor

**Question and setup.** Fit the four VIX forwards tightly and allow the floor to move in the additive model.

**Result.** Tightening forwards alone worsened the option fit; freeing the floor reduced fresh comparable total loss from 12.7780 to 1.7830. It remained a partial VIX fit, with some boundary and unavailable-IV issues.

**Learning and consequence.** The floor materially changes response dispersion and left-tail options. Training forward equalities remain estimator constraints, not exact population identities.

Source: [forward and floor](../../matteo/model_comparison/option_pricing/FORWARD_FLOOR_RESULTS.md).

## 13 Continuous mixture and shared two-speed stress

**Question and setup.** Fit a continuous Beta mixing measure, then use a simpler shared slow/fast stress representation across all VIX expiries.

**Result.** The continuous measure concentrated near endpoints and checking option loss 1.7307 did not beat the boundary's 1.7256. The shared two-speed experiment reduced checking option loss 1.725555 to 0.050081 and increased inside-bid–ask coverage 83 to 153 of 217 quotes. Its stress weight was about 92.3%, floor volatility 12%, and half-lives 365/71.4 days, with boundary hits. Checking forward gaps reached about 0.055 points.

**Learning and consequence.** This was the strongest partial VIX fit, using one common model across maturities. Checking reused old/slow coordinates with new fast coordinates; it was not wholly fresh. The experiment did not establish SPX fit or production adoption.

Sources: [continuum](../../matteo/model_comparison/option_pricing/CONTINUUM_RESPONSE_RESULTS.md), [two speeds](../../matteo/model_comparison/option_pricing/TWO_SPEED_RESULTS.md). [Preserved figure](VISUAL_EVIDENCE.md#two-stress-speeds).

## 14 Four fixed stress speeds

**Question and setup.** Add a shared positive stress bank at 7/28/112/448-day half-lives rather than a special short-maturity model.

**Result.** Checking option loss was 0.062160 versus 0.050081 for two speeds; paired difference +0.012079 had jackknife SE 0.024117. Checking forward errors reached 0.216 points despite nearly exact sampled training equalities.

**Learning and consequence.** No demonstrated superiority. Preserve the completed bank experiment separately; improve population-forward precision before interpreting small option-loss rankings or adding branches.

Source: [stress-bank study](../../matteo/model_comparison/option_pricing/STRESS_BANK_RESULTS.md), including full native-IV and forward figures.

## 15 Joint refit and protected curve

**Question and setup.** Refit SPX/VIX together from the two-speed candidate, then fit only the nine curve levels outside quoted VIX windows to isolate VIX-preserving SPX flexibility.

**Result.** The bounded joint fit reduced SPX price loss 5467.10 to 618.74 while VIX loss rose 0.0501 to 37.97 and stress weight collapsed to 1.94%. The separate protected-curve study reduced fresh SPX price loss 5979.32 to 3895.07 but worsened IV MAE 4.345 to 4.710 pp, leaving downward skew missing.

**Learning and consequence.** Shared variance-budget optimization can destroy a successful VIX distribution. Off-window curve freedom protects VIX in the exogenous construction, but lower price loss is not automatically better smile shape. Neither local search proves impossibility.

Sources: [joint fit](../../matteo/model_comparison/option_pricing/TWO_SPEED_JOINT_RESULTS.md), [protected curve](../../matteo/model_comparison/option_pricing/PROTECTED_CURVE_RESULTS.md). Figures 4–5 in [the gallery](VISUAL_EVIDENCE.md).

## 16 Stock and stress correlation

**Question and setup.** Keep the full variance-factor law fixed while rotating the stock Brownian direction towards negative exposure to stress shocks.

**Result.** Theta=0.8 improved fresh SPX price loss 3788.48 to 1861.38 and IV MAE 4.493 to 3.199 pp. VIX law was unchanged. The 29-day right-minus-left IV gap moved about -1.52 to -5.58 pp versus market -14.05 pp; only 18/1294 prices were inside bid–ask.

**Learning and consequence.** This supplies a useful SPX control with mathematical VIX protection under the fixed-factor/information assumptions. State-dependent exposure opposite instantaneous variance sensitivity or expected next-month variance sensitivity did not produce a convincing further benefit. Maximizing local covariance does not maximize finite-expiry option skew.

Sources: [rotation](../../matteo/model_comparison/option_pricing/LEVERAGE_SCAN_RESULTS.md), [targeted skew](../../matteo/model_comparison/option_pricing/TARGETED_SKEW_RESULTS.md), [forward sensitivity](../../matteo/model_comparison/option_pricing/FORWARD_LEVERAGE_RESULTS.md). [Preserved figure](VISUAL_EVIDENCE.md#stock-and-stress-correlation).

## 17 Asymmetric readouts and slower existing signals

**Question and setup.** Multiply by a downside response of an existing factor; compare global versus stress-regime normalization, curve refits and slower existing-bank signals.

**Result.** Stronger asymmetry could improve SPX while badly damaging VIX. Global normalization altered effective stress allocation; preserving its conditional mean reduced damage without protecting the full VIX law. A weak slower signal plus refit reduced SPX loss 1866.88 to 1811.15 but moved VIX loss about 0.050 to 0.085; the 29-day skew gap changed only -5.98 to -6.17 pp versus market -14.05. Stronger settings cost more VIX fit.

**Learning and consequence.** Use ablations: some striking changes required both the readout and changed curve. Neither mean-one normalization nor a slower signal guarantees distributional protection or enough SPX downside level. No candidate adopted.

Sources: [asymmetric readout](../../matteo/model_comparison/option_pricing/ASYMMETRIC_READOUT_RESULTS.md), [targeted normalization](../../matteo/model_comparison/option_pricing/TARGETED_SKEW_RESULTS.md), [slow signal](../../matteo/model_comparison/option_pricing/SLOW_READOUT_RESULTS.md). Figures 7–8 in [the gallery](VISUAL_EVIDENCE.md).

## 18 Linear memory filters

**Question and setup.** Test positive normalized-signal history filters on M2 alone and on the whole multiplier, retaining exogenous factors and realized filter histories in conditional VIX pricing.

**Result.** On 1024 new stock paths at 256 steps/day, baseline SPX loss 1903.59 became 1899.76 for M2 filtering and 1901.94 with refit. Whole filtering plus refit gave 1912.77. The report concludes no material improvement; missing skew remains. Numerical bridge and likelihood-weighted refinement limitations are recorded.

**Learning and consequence.** A persistent average has no direct Brownian shock; replacing part of an instantaneous signal by history can reduce instantaneous loading. This tested positive linear filtering, not every nonlinear hysteresis model. It is completed, not pending, and no candidate was adopted.

Source: [memory-filter study](../../matteo/model_comparison/option_pricing/MEMORY_FILTER_RESULTS.md). [SPX and VIX figures](VISUAL_EVIDENCE.md#linear-memory-filters).

## 19 First and richer protected-generator corrections

**Question and setup.** Change variance paths while preserving conditional integrated variance at all four quoted VIX dates. Use h=(partial_t+L)P with endpoint-zero potentials built from existing OU modes; derive global positivity bounds before fitting.

**Result.** The first bounded correction reduced SPX loss 1841.07 to 1801.18 and IV MAE 3.169 to 3.121 pp; its 29-day skew remained -6.59 pp versus market -14.05. Richer local sine-squared/erf shapes on a separate fresh sample reduced loss 1945.52 to 1884.31, with IV MAE 3.209 to 3.137 pp and gap -5.85 pp. Bounds guarantee v>=f/2=0.0072, not the original floor f.

**Learning and consequence.** Conditional Dynkin protection can be translated into finite implementation audits, but this useful feasibility result gave only modest empirical progress. It protects the four model VIX values, inheriting baseline errors, not every unquoted date. No new Lean proof or adopted candidate.

Sources: [first correction](../../matteo/model_comparison/option_pricing/PROTECTED_VARIANCE_RESULTS.md), [richer shapes](../../matteo/model_comparison/option_pricing/RICH_PROTECTED_VARIANCE_RESULTS.md). [Figures](VISUAL_EVIDENCE.md#protected-variance-corrections).

## 20 Existing stress variance as correction budget

**Question and setup.** Enlarge the usable correction budget by scaling the potential with existing stress variance, applying the full generator and retaining the original floor f=0.0144. This also changes calendar zeros and damping, so it is not a pure amplitude ablation.

**Result.** A 38-configuration pilot selected the new candidate before 1024 fresh stock paths. At 256 steps/day, losses were baseline 1948.66, previous rich correction 1888.60 and new correction 1925.31. New-minus-baseline loss -23.35 had paired approximate 95% interval [-120.12,76.59]; new-minus-previous also included zero. Downside IV was 19.26% versus market 26.90% at 29 days near K/F=0.9.

**Learning and consequence.** The enlarged budget was exercised in probes, but no improvement was established. The missing downside level persists. Preserve both the endpoint/path insight and these discouraging results; neither mathematical feasibility nor larger amplitude supplies an empirical success. No adoption.

Source: [available-variance study](../../matteo/model_comparison/option_pricing/AVAILABLE_VARIANCE_RESULTS.md). [Full seven-maturity figure](VISUAL_EVIDENCE.md#existing-stress-variance-budget).

## What remains exploratory

The literature-inspired VIX-first inverse construction in the ledger is a design study, not a completed calibrated M2 replacement. Quintic and jump remedies are outside the completed M2 joint-fit tests. Existing scoped Lean results preserve their original premises; the latest filters and protected corrections rely on analytic derivations plus finite numerical correspondence checks. Reassessing the physical fit of promising Q variants and resolving numerical accuracy remain open.

Use [the decision record](DECISIONS.md) before selecting another test. The [catalogue](REPORT_CATALOG.md) retains additional parity, acceleration, conditional-pricing, readout and diagnostic stages that support this sequence.
