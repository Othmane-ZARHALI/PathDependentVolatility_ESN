# Uniform 4096 paths — user steering during preparation

The user agrees that stock paths are shared across strikes and explicitly chooses
ALL paths at EVERY strike. Supersedes adaptive allocation in V1/V2. Preparation
of the same 4096-path cache and independent frozen-control pilots continues
unchanged. Production uses a uniform 4096 iid SPX paths and eight independent
512-point VIX nets per expiry. SPX SE is sample SD/sqrt(4096); VIX SE is SD of
eight complete net prices/sqrt(8). There are no moneyness budget transitions.
Retain uniform nested 512/1024/2048/4096 results as precision diagnostics only.

The call-payoff/forward conventions, fixed model parameter comparison, controls,
OTM input population, quadrature/time grids, auditing and signed-price reporting
remain as in V1/V2. The initially queued runners 001/002 were interrupted while
waiting for paths/controls, before any objective evaluation. Completed comparison
will be runs/m2_parity_uniform_001. No recalibration is claimed.

The user's question about starting-point reuse does not authorize changing the
agreed initial state in this comparison: it remains C0=1,Y0=F0=0. Implement/test
an exact OU translation helper for future reuse, explain that G/normalizer/stock/
VIX must be recomputed, and keep initial-state calibration separate from the
previously agreed eleven parameters.
