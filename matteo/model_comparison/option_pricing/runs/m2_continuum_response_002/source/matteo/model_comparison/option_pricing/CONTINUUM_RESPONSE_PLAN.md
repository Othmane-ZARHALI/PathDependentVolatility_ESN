# Continuous response VIX-only experiment, 27 October 2026

User authorized testing the proposed continuous probability measure with exact forward constraints.

## Mathematical model

v_t=f+(xi(t)-f) integral R_t(u) mu(du), u in [0,1],
R_t(u)=G_t^(1-u)/E_Q[G_t^(1-u)] exp(b u Z_t - (b u)^2 q(t)/2).
G is the unchanged capped M2 readout. Z is one independent OU, initially zero,
q(t)=1-exp(-2at). mu=Beta(m kappa,(1-m) kappa), deterministic under Q.
All components are combined before the VIX conditional variance integral and square root.
The two-atom reference is the kappa->0 boundary, not a second separately weighted integral.

## Scope and design frozen before fitting

All 217 native VIX quotes and four forwards, expiry days 27/55/118/153.
16 xi intervals, same seven active independently fitted levels. Other nine unchanged.
Old readout/cap, architecture, premiums, initial old states, controls and quote definitions frozen.
No clipping, smoothing or dropping signed prices. No SPX pricing or baseline replacement.
Bounds: b 0..5; m .01...99; half-life 1..365 calendar days;
kappa .01..30 (plus separately evaluated two-atom boundary); floor vol 5..12%;
active xi between fitted floor variance and .2. Each component normalized separately.

Training uses the first 1024 points from each of eight nets of the previous
forward/floor experiment's checking draws: 8192/expiry, explicitly reused as training.
Positive Gauss-Jacobi measure rule order24; cubic interpolation of log conditional
power ratios on 33 uniformly spaced u points; old conditional GH32/time8,
unconditional normalizers GH48 at exact time nodes. Refinement checks compare
measure48, power65, conditional GH48/time16 on sampled states.

Native bid/ask scaled option distance, equal expiry weighting, is minimized using
SLSQP with four equality constraints on the training forward estimator, analytic
chain-rule price/forward Jacobians (quadrature node derivatives by small central
finite differences). Selection requires max forward residual <=1e-5 VIX point.
Fit a two-atom boundary on the same training draws, then continuous starts near
that boundary and with dispersed/uniform/concentrated measures. Preserve all starts.
Fresh check uses eight new scrambled joint importance nets, 2048 points/net;
selection is frozen beforehand. Exact importance weights and frozen pilot controls.
Report fresh forward residuals with eight-net SE; equality on a training estimator
is not exact population equality. Same market targets, no market out-of-sample claim.

Audit Gaussian powered-readout oracles, original endpoints, Beta moments,
price/forward Jacobians, sampled deterministic quadrature refinement, tail weights,
second moments, negatives, artifact/source hashes and paired loss differences.
No new formal theorem or global optimization claim.
