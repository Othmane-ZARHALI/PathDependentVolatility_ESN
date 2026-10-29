# Existing-bank slow downside readout — experiment record, 28 October 2026

Authorized question: can slower signals already present in the frozen Gaussian feedback bank supply persistent SPX downside leverage, while preserving the fast M2 component and an acceptable four-maturity VIX fit?

Select one existing mode at a time, nearest in log half-life to 2, 8, 32 days. Actual modes: indices 7, 6, 5, half-lives 3.329125, 9.481468, 27.003562 days. These are generic memory scales, not expiry-specific parameters. Test alpha .5 and 1 first. Their forward-constrained VIX training losses were too large; a second exploratory stage tests .125 and .25. Both stages use identical training draws and stock seeds so cross-stage differences are paired. The second grid is reactive, not a preregistered independent validation.

For centered unit-stationary-standard-deviation mode X, set D=exp(alpha*softplus(-X)). Original fast M2 R, both independent stresses, their parameters and fixed mixed stock loading theta=.8 stay intact. Use stress-regime conditional mean normalization:

H=(aR+wA)D*(a+wA)/(a*E[RD]+wA*E[D]), a=1-w.

E[H|A]=a+wA, but the VIX distribution need not be preserved. Compute correlated R/D product moments under the full frozen Gaussian law. No independent-product approximation. Rebuild every downstream normalization, stock path and conditional VIX exposure. Full Gaussian OU dynamics are exogenous and unchanged; no endogenous dependence on returns or volatility.

For each slow mode/strength compare a fixed 16-level variance curve and a refit of the same seven active levels under four numerical training-forward equalities. Nine other levels remain frozen. Keep all 1751 eligible/1294 primary SPX quotes and seven expiries, 217 VIX quotes/four forwards, signed Monte Carlo estimates and invalid-IV counts. Never clip or smooth prices. Production baseline remains unchanged.

Initial stock comparison: 512 fresh paths/eight independent antithetic batches, coupled 128/256 steps per trading day, fixed deterministic-delta controls and existing 65536-path conditional-lognormal auxiliary reference. VIX training 512 states per expiry, separate previously examined checking stream 4096. Training equality is not exact population forward matching.

Numerical audits: alpha-zero replay against unchanged stock kernel, independent scalar conditional tensor quadrature of E[GD], GH48/64 normalizer check, interpolation checks, and a matched-state VIX quadrature refinement. Freeze selection before any additional fresh stock/outer-state validation. Report approximate paired intervals and finite-grid limits. This experiment does not introduce a Lean theorem certifying the simulator or establish a model-family impossibility.
