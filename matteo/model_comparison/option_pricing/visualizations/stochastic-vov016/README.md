# Dedicated stochastic sensitivity: reading the figures

The new state randomly changes spike sensitivity. Its half-life controls how long
an episode lasts. The entire expected variance curve, floor and original reservoir
are fixed. The two outer transformations are the current cap and a logarithmic
continuation with strength .5. This is a bounded diagnostic, not a new accepted fit.

- `vix-grid-cap0` / `vix-grid-cap05`: every proposed strength/persistence case.
- `strength-persistence-grid`: wing errors; lower is closer to market.
- `vix-all-maturities-iv`: matched training selections and cap-only controls.
- `spx-all-maturities-iv`: actual full stock repricing at all seven dates.
- `*-signed-errors`: signed call errors in quote scales, including failures.
- `vix-new-checking`: new joint scrambles generated after selection, with bands.
- `vix-forwards` / `vix-weighted-tails`: level and rare-tail diagnostics.

Solid lines add the new state; dashed lines are no-layer references. Market inputs
remain original OTM quotes converted to common calls. IVs use market forwards.
All 4096 observations enter every strike. Curves connect raw price-derived IVs;
there is no smoothing or price clipping. Missing IVs are marked. Sampling bands
exclude integration, stock timestep, pilot and parameter uncertainty. New numerical
streams use the same in-sample market targets, with no market out-of-sample claim.
See STOCHASTIC_VOV_RESULTS.md and the run's audit for evidence and limitations.
