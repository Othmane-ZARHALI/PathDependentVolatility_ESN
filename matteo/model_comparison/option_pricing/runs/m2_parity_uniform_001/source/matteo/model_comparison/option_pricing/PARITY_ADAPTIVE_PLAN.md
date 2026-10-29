# One payoff and larger wing budgets — 24 October 2026

Implement the user's two requested estimator changes before further calibration.
Keep all 1294 SPX and 217 VIX original OTM observations and four forwards, cutoff
3 and all fixed reservoir/model settings. No additional market observations.
Convert observed OTM puts to undiscounted calls by C=P+F_market-K, including
bid/ask bounds. Price a single call payoff at every strike. For SPX the model
forward is the given market forward (adapted exponential Euler has mean-one
normalized stock). For VIX it is an independently estimated model expectation:
model puts are C-F_model+K. Never substitute the market forward for that model
expectation. The canonical-call objective therefore differs from the former OTM
objective for VIX when forwards differ; preserve original-OTM residuals too.
Display canonical-call IVs with the observed market forward, and show VIX forward
errors explicitly. Invalid IVs, signed prices and shape failures remain visible.

Freeze estimator choices before examining new samples. SPX: one smooth analytic
lognormal-reference control at every strike; both actual stock and reference use
the same cached equity shocks. Formula (normalized forward 1):
C_BS(K)+w[(R-K)+-(R_ref-K)+-delta_BS(K)(R-R_ref)]. Both stocks have target mean 1;
reference call mean is analytic. No coefficients fitted on evaluation paths.
This can still give negative estimates; comparison with old controls is required.
VIX: fixed Gaussian polynomial controls fitted on an independent pilot at the
latest selected parameters; one common call-control family across the forward.

Path tiers 1024/2048/4096. Use 1024 up to standardized distance 1, continuously
blend 1024 and 2048 estimates between distances 1 and 2, and 2048/4096 between
2 and 3; 4096 beyond 3. Smoothstep weights depend only on frozen market inputs.
Nested samples share their covariance: do not treat prefix estimates as independent.
SPX SE uses iid corrected draws and the exact deterministic sample coefficients;
VIX uses eight independent scrambles, prefix sizes 128/256/512 per scramble, with
SE across the eight complete blended estimates. Report actual and effective counts.
No payoff-dependent stopping. Blending avoids sample-count jumps but does not
ensure convexity or uniform precision. Fixing precision in the wings remains
an empirical question, especially for small option prices.

Use fresh fixed exogenous OU/equity paths per expiry, saved in reusable chunks;
retain grid 512 stock steps/trading day, GH32 VIX, GH48 normalizer, time order 8.
Preserve previous caches/runs. First reprice old and latest selected vectors on
identical new draws; do not optimize while diagnosing the estimator change.
Save sources, seeds, counts, prices, controls, SEs, forward/parity/shape diagnostics
and all eleven maturity plots. Tests must check parity, known-price oracle,
shared-prefix SE covariance, analytic derivatives and allocation continuity.
