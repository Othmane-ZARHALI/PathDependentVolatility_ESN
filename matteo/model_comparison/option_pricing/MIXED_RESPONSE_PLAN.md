# Two variance responses, VIX-only calibration

Authorized 26 October 2026: point 2 of the ranked experiments; forward variance
levels remain calibration parameters. Original 16 intervals, seven independent
VIX-sensitive levels (indices1,2,3,4,7,8,10); nine off-window entries unchanged
because sensitivity is zero. Retain all217 native OTM quotes and four forwards.
SPX remains paused. Old readout/cap, beta_J>0, OU banks and floor0.00264 fixed.

R=G_old/M_old. Independent Q OU dZ=-aZdt+sqrt(2a)dB, Z0=0,
q(t)=1-exp(-2at); L_b(t)=exp(b Z_t-b²q(t)/2). Two pathwise variance sums:

1. Product mixture: v=f+(xi-f) R [(1-w)L_b1+w L_b2].
2. Additive mixture: v=f+(xi-f) [(1-w)R+w L_b].

These are variance contributions on each path, never mixtures of completed
option prices. The first uses the explicit two-exponential response proposed in
the research note; the second is its stated M2/stress-component adaptation.
Each has exact expected variance xi(t), since components have mean1 and Z is
independent of old factors. Original M(t) is unchanged. Conditional factors use
E[L_b(u)|Z_T]=exp(b exp(-a(u-T))Z_T-b²exp(-2a(u-T))q(T)/2).
Conditional expectations factor only where independence requires it.

Both sums are positive, with no additional cap on L. On a finite horizon,
conditioning on the independent continuous Z path bounds L_b1 and L_b2 there;
the old capped G/M is bounded. The resulting conditional stock volatility is
bounded, giving the same continuous-model martingale argument as the prior
post-cap model. This is mathematical reasoning under the stated independence,
not a new Lean theorem or estimator convergence guarantee.

Product b1 in[0,3], b2=b1+(5-b1)d with d in[0,1]; w in[0,1]. Additive b in[0,5],
w in[0,1]. One half-life in[1,365] calendar days per family. xi in[f,0.20],
independent with no smoothing/monotonicity penalty. Endpoints nest original M2
and the single-factor mechanism. No new factor correlation or bank dynamics.

Objective: unchanged native OTM bid–ask distances /existing price scales, equal
total expiry weights, plus four forward errors /0.5 point. No SPX term or extra
wing weighting. Compare old baseline, refitted curve-only, prior fitted single
factor, native-loss winner and descriptive wing-shape winner for each family.
Finite multistart least squares, at most120 evaluations/start. Freeze selections
before fresh simulations. Save every optimizer outcome, including failures.

Training reuses prior old-factor4096 conditional ratios with two independently
scrambled Z importance samples per old state,8192 total/expiry in eight clusters.
Z proposal is a fixed mixture with shifts0,2,4,6 and probabilities.5,.25,.2,.05.
Exact combined weights w_old*w_Z, no self-normalization. Repeated old states are
clustered by their original net, not treated as independent observations.
This explicitly samples stress tails instead of relying on Gaussian maxima.

Fresh check: after selection, generate eight2048-point scrambles/expiry from
the product of original17-component old-state proposal and that Z proposal,
in37 coordinates (68 mixture components). Fixed training seed base20261026071;
checking20261026072. Same market targets, not market out-of-sample validation.
Frozen independently fitted old pilot controls remain unchanged.

Checks: analytic Jacobian in every free direction; nesting old/single-factor
values/prices; Gaussian integral oracle for conditional mixtures and mean1;
fresh paired losses and wing diagnostics; weighted tail concentration, moments,
ESS; GH32/time8 versus GH48/time16 sampled refinement including states with
large selected-candidate payoff contributions. Preserve signed estimates,
unavailable IVs, old baseline and all prior runs. No automatic model adoption.
