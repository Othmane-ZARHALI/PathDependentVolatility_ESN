# Fixed-reservoir joint readout calibration — 24 October 2026

User authorizes recalibration of all coefficients compatible with unchanged
exogenous reservoir paths and new plots at every agreed maturity.

Free (11): cluster_loading, feedback_curvature, feedback_shift,
spike_curvature, spike_shift, cap_level, cap_sharpness, three forward-variance
levels and variance_floor. Enforce spike_shift >= 0.01 (strictly positive)
and strictly positive curvatures. The numerical lower bound is a declared
experimental choice, not an inferred economic threshold.

Fixed: complete Stage-3 architecture (OU counts/rates, correlations, projection
powers, orthogonal/echo settings), Q risk prices and projected initial state
C0=1,Y0=F0=0. Optional orthogonal/echo channels stay disabled: activating them
changes the model cell. Stock paths and Q normalization are recomputed from
cached features for every candidate. No historic artifacts are overwritten.

Scope: same 1294 OTM SPX cutoff-3 quotes at 29/57/119/211/302/666/1030 days,
217 OTM VIX quotes at 27/55/118/153 days, four VIX forwards. All market targets
are in sample. Keep original source labels as metadata. Same interval loss,
expiry balance and forward scale as existing all-maturity objective; retain
curve penalty and constant state penalty. No MC-SE normalization.

Estimator: positive raw importance-weighted payoffs with frozen proposal, no
strike-dependent controls, no clipping or smoothing. Reprice old zero-shift
parameters and a positive-shift start with the same estimator. This changes
the numerical estimator, so old saved loss values are not direct comparators.
One fixed master stock path grid at 512 steps/trading day, exact event dates,
128 searching paths, 256 separate checking paths. Eight VIX outer scrambles
of 16 points searching, 32 checking. Frozen existing targeted VIX proposal.
Search conditional GH12, normalizer GH48, VIX time order 8; frozen selection
checked at conditional GH24 and compared on matching states with GH12.
These finite budgets/quadrature choices are exploratory, not convergence.
Candidate normalizers use log interpolation on a fixed transformed-time grid;
independent same-order probes report its error separately.

Bounded optimization: two starts, at most 35 residual/Jacobian evaluations per
start, 30 minutes total optimization wall budget (checked between evaluations).
Analytic derivatives, finite-difference checks before search. Preserve every
objective evaluation, selection, raw values, seeds, source/input hashes and
failures. Selection freezes before separate numerical ensemble checking.
Plot baseline and selected model on all seven SPX/four VIX slices plus forwards,
variance curve and parameter changes. Report missing IVs, sampling bands,
boundary hits, numerical reversals and unresolved P compatibility explicitly.
No new theorem or formal certification is asserted.
