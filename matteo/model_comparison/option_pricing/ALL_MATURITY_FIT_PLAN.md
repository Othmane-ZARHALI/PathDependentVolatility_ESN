# All retained maturities: in-sample feasibility run, 23 October 2026

The user explicitly requests all SPX maturities except days 1, 1401, 1765 and
2129, cutoff 3, and no out-of-sample results. This supersedes the earlier
four-maturity training-only objective. Use **29/57/119/211/302/666/1030 days**,
including all eligible original quote splits. Retain all four VIX maturities
and use all their eligible original quote splits as in-sample calibration data.
Historical source labels and numerical runs are retained without rewriting.

Market-only selection before pricing: nearest-to-ATM undiscounted midpoint
Black IV per retained SPX expiry, tie by strike/option type, using all eligible
quotes. Include `abs(log(K/F))/(sigma_proxy*sqrt(T)) <= 3`. Counts are
277/298/228/141/167/112/71, **1294 SPX options** out of 1751 eligible. All
**217 VIX options** and four market forwards enter the joint objective. The
new adapters label copied rows `in_sample_fit`, preserve `original_split` and
`original_maturity_role`, and explicitly pass that split to the production
objective validator. No raw/preprocessed market split or default is edited.

First price the existing C0=1/Y0=F0=0 three-piece working reference on this
expanded population. Then refit only its three variance levels. Keep structure,
readout and Q scenario 0 fixed; do not mix a risk-price change into this scope
expansion. The earlier planned risk-price comparison is deferred until this
requested expanded baseline is established. The curve keeps its existing
quarter-year/three-quarter-year/final boundaries and original 2129-day horizon.

The fixed-scale price-interval loss remains expiry-balanced within each family,
with SPX/VIX/forward weights one, forward scale 0.5 points, weak adjacent-curve
penalty 0.1 and the constant selected-state penalty. Absolute objective values
are not directly comparable to the previous four-SPX-maturity objective.

Reuse compatible candidate-state ratio/equity caches and fixed known-zero
controls at the original four SPX expiries, and state-specific conditional VIX
exposures at all four VIX expiries. New strikes use nearest same-side frozen
controls and coefficients; exact old keys must replay their previous prices.
For SPX 119/666/1030, generate new zero-state primitive paths, rebuild C0=1
candidate G/M0 and stock repricing caches, and construct the original five
anchor-control family. Coefficients transfer from independent frozen 302-day
pilots with forward scaling; known-zero means justify reuse, while efficiency
at new maturities remains unqualified. This is not a new control-feature search.

Budget: 512 searching paths per expiry and a separate 512-path numerical-check
ensemble at the SAME market quotes. Existing searching and prior fresh-t0
ensembles are reused at original dates; added SPX dates receive distinct new
search/check streams. The check is not a new independent complete calibration
validation and has no market out-of-sample interpretation. New maturity jobs
run with at most two workers. Three variance starts, at most 30 optimizer
evaluations each, 30-minute optimization budget after cache preparation. Freeze
the selected curve before checking scores are evaluated.

Keep 512 stock steps per trading day, normalizer GH48 with the declared 2e-7
sampled interpolation screen, VIX conditional GH32 and time order 8. Added
maturity grids preserve exact maturity/curve events; these grids differ from
old grids and do not supply a stock-grid convergence test. Sampling SE excludes
stock-grid, normalizer/GH/time bias and pilot uncertainty. Preserve signed
estimates and unsuccessful starts. No parameter fit proves global feasibility
or infeasibility. The current Lean identities apply under their original
assumptions; scope selection, control transfer and repricing remain tested
implementation bridges, not formally verified code.

All seven SPX and all four VIX panels must show actual model estimates, with
market forwards used for all IV displays, and explicit **in-sample fit** labels.
Show the old working parameters and expanded-scope refit separately. Save quote
membership, original labels, fixed scales, source hashes, seeds, optimization
records and independent score replay. Do not label a new maturity market-only,
reserved, diagnostic or out-of-sample in the new figures.

Changing active_scope.json is an explicit new user decision. Previous run source
snapshots preserve the old scope; old strict root-hash replay scripts now see
this intentional scope-version change. Their completed historical audits remain
evidence for their captured versions; do not refresh historical hashes.
