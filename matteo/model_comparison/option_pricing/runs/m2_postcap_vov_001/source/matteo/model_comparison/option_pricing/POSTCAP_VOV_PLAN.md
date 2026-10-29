# Independent normalized amplitude after the original cap

Authorized 2026-10-26: first alternative from the preceding discussion. No mixture.

Let R_t=G(X_t)/M(t) be the frozen baseline's capped normalized component.
Add independent Q OU dZ=-aZ dt+sqrt(2a)dB, Z_0=0, q(t)=1-exp(-2at).
L_t=exp(gamma Z_t-gamma^2 q(t)/2), v_t=f+(xi(t)-f)R_t L_t.
All baseline parameters, 16 xi levels, floor, cap, OU banks and correlations fixed.
The old normalizer is unchanged: independence and E[L_t]=1 imply E[R_tL_t]=1.
For u>=T, E[L_u|Z_T]=exp(gamma exp(-a(u-T))Z_T
- gamma^2 exp(-2a(u-T))q(T)/2). Conditional expectations factor by independence.
The added variance diffusion coefficient is gamma sqrt(2a)(v-f); there is no
additional tanh or cap on this channel. The old cap bounds only the old G.

For any finite horizon, conditioning on the independent Z path makes L a bounded
continuous time function. The original capped G is bounded and its deterministic
M is positive and continuous on that horizon. Conditional stock volatility is
therefore bounded, giving the continuous stock true-martingale identity, provided
the original equity/old-factor Brownian vector is independent of Z. All fixed-time
moments of L and the resulting variance are finite. This does not establish finite
stock-payoff estimator variance or sampling convergence. No new Lean claim.

## Bounded comparison, frozen before prices

- gamma 0.25, 0.5, 1.0; half-lives 7, 30, 90 calendar days; gamma=0 reference.
  gamma measures log-VARIANCE amplitude, not a bounded multiplier range.
- Same original training Gaussian states/weights, with common independent Z draws;
  4096 points/strike. Preserve old independent pilot controls and signed prices.
- Original 4 VIX maturities, 7 SPX dates, all1511 OTM quotes, cutoff3; common calls
  use market forwards. Native OTM/forward losses reported separately.
- Two descriptive VIX selections: minimum native VIX+forward price loss, and
  minimum combined relative wing slope/curvature discrepancy. Require complete
  wing IVs; shape score averages squared errors divided by max(10,abs(market
  coefficient))^2 over the four maturities and two coefficients. Also report ATM,
  level errors and unavailable IVs; these selections are not acceptance criteria.
- Full SPX reprice baseline and union of these at most2 selections using streamed
  exact OU updates on the old stock mesh. All downstream stock paths recomputed.
- Freeze selection before fresh joint37D VIX scrambles. Use eight512-point nets
  at each expiry, exact defensive-mixture weights, no self-normalization.
- Check analytic conditional factor against independent integration, actual
  stock kernel against Python and zero-factor original kernel; compare sampled
  old-state GH32->48 and time8->16. Report unresolved deterministic/sampling errors.
- Preserve existing selected_model.json and previous study. No automatic adoption.
