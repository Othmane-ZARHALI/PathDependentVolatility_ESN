# Continuous stochastic vol-of-vol after the cap profiles

25 October 2026. This distinguishes published results, proposed M2 adaptations and actual fixed-curve experiments. No jumps, quintic proposal or modification to the existing OU rates/geometric grids. See CAP_TAIL_RESULTS.md for numerical outcomes and CAP_TAIL_MATHEMATICS.md for baseline equations and stock-martingale arguments.

## What stochastic vol-of-vol means here

If v=F(t,X), Ito's rule makes the Brownian part of variance

    dv_t = ... dt + grad_X F(t,X_t) B dW_t.

Its instantaneous random-motion strength is ||B^T grad_X F||. Thus the current nonlinear M2 already has state-dependent instantaneous vol-of-vol. The new question is whether to give that response extra independently adjustable variation, rather than merely making variance bigger. A cluster state can alter the spike/feedback sensitivities; a dedicated latent state can alter their strength and persistence independently of the cluster's current role. Either changes the distribution of conditional future variance, hence potentially the VIX smile.

## Primary literature checked live

[Fouque and Saporito, Heston Stochastic Vol-of-Vol Model for Joint Calibration of VIX and S&P 500 Options, arXiv1706.00873v1](https://arxiv.org/html/1706.00873v1): their variance diffusion coefficient is eta(Y_fast,Z_slow)*sqrt(V). Fast/slow continuous factors govern eta. Under stated regularity and scale assumptions they derive first-order joint pricing approximations around Heston. This motivates modulating response strength with another state. Their perturbation formulas and accuracy theorem do not apply directly to M2's normalized Gaussian readout.

[Alos, Garcia-Lorite and Muguruza, On smile properties of volatility derivatives and exotic products: understanding the VIX skew, arXiv1808.03610v1](https://arxiv.org/html/1808.03610v1): they derive short-time ATM level/skew expressions using Malliavin derivatives. Their mixed generalized rough-volatility example adds two exponential components with different sensitivities; they discuss mixing as a way to overcome flat VIX smiles. This supports examining sensitivity mixtures. It does not prove any particular M2 cap or cluster modulation will fit finite-strike, finite-maturity wings.

[Jacquier, Muguruza and Pannier, Rough multifactor volatility for SPX and VIX options, arXiv2112.14310v2](https://arxiv.org/html/2112.14310v2): they study ATM level, skew and curvature, including a two-factor rough Bergomi model. Their calibration discussion separates the sensitivity/mixing parameters used for VIX from equity correlations used for SPX. The relevant lesson is to assess shape and cross-market coupling separately. We preserve M2 correlations; these asymptotic equations are not M2 pricing formulas.

Additional primary empirical evidence: [Romer, Empirical analysis of rough and classical stochastic volatility models to the SPX and VIX markets, Quantitative Finance2022](https://researchprofiles.ku.dk/en/publications/empirical-analysis-of-rough-and-classical-stochastic-volatility-m/). The author's university record describes a two-OU model with a nonstandard transformation and joint calibration results without jumps or roughness. This supports continuous transformation/multifactor exploration. It is an empirical result for another model and market exercise, not evidence that the present frozen reservoir is sufficient.

## Existing-cluster experiment actually implemented

Let eta0 be the original score and eta1=eta0+b*qJ*(F-betaJ)^2. Set

    w(C)=0.5+0.45*tanh(0.5*C),
    G_b=(1-w(C))*exp(2h(eta0)) + w(C)*exp(2h(eta1)),
    M_b(t)=E_Q[G_b(t)],
    v_b=f+(xi-f)*G_b/M_b.

As C rises, more weight moves to a component with stronger spike sensitivity. b=.5 and1 were profiled, holding curve/floor and all original coefficients fixed. b controls added spike response, the tanh scale controls how quickly cluster states alter that response, and the0.45 amplitude limits mixture weights to[.05,.95]. Those two tanh choices were frozen diagnostic design choices, not calibrated values. b=0 recovers the original model exactly. Two spike components can respond differently before the common cap; after both saturate, modulation cannot restore unlimited tail response. This experiment therefore tests finite-state response redistribution, not elimination of the ceiling.

The mixture remains positive and bounded by exp(2L); Gaussian moments and the bounded-volatility martingale argument remain available. No new random variables or OU parameters are introduced. Existing full outer states, conditional Gaussian laws, stock C/Y/F features, equity increments, likelihood weights and frozen zero-mean controls remain valid primitives. G, Q normalizers, VIX exposures, volatility and stock paths must all be recomputed. Actual cluster experiment reprices all four VIX dates only; it has no claimed SPX tradeoff score or joint calibration. It is separate from the two cap candidates receiving full SPX checks.

The scalar cap convolution is insufficient for state-dependent mixture weights. The experimental kernel adds a second interpolation coordinate for conditional C means (spacing.05), with Hermite score interpolation (spacing.02); direct same-GH checks and independent Python tests assess it. It changes no model equation. Measured full VIX evaluation times include this extra integration cost; compare saved timings with the roughly8-second baseline profile.

## Dedicated continuous latent layer: explicit proposal, not implemented

A minimal additional state under Q is

    dZ_t = -a_Z*Z_t dt + nu_Z dB_Z(t), Z_0=0,
    A(Z)=exp(s*tanh(Z)),
    eta_Z=c_C*C+qY*(Y-betaY)^2+qJ*A(Z)*(F-betaJ)^2,
    G_Z=exp(2h(eta_Z)), M_Z(t)=E_Q[G_Z(t)],
    v_Z=f+(xi-f)*G_Z/M_Z(t).

This adds one state to the unchanged reservoir. a_Z sets persistence, nu_Z controls how far Z moves, s sets the range of spike sensitivities, qJ*A in[qJ*exp(-s),qJ*exp(s)]. To avoid poor identification between nu_Z and s, start with Var_stationary(Z)=1, so nu_Z=sqrt(2a_Z), and fit only a_Z and s. s=0 exactly nests baseline. All states and stock paths remain continuous. A parallel feedback modulation is possible, but adding both simultaneously would obscure which mechanism helped.

For independent B_Z, existing OU/equity draws and their importance likelihoods can be reused unchanged; simulate independent Z values at the SAME stock times and outer VIX times. The new Z likelihood is1 when sampled directly under its Q law. With nonzero correlation, shared driving randomness and covariance/conditional laws must be reconstructed consistently; assigning arbitrary correlation to saved paths afterwards is invalid. Begin independent to isolate the new response mechanism. Correlation is not an authorized baseline change.

Conditional pricing must condition on Z_T as well as the original full OU state. Z_u conditional on Z_T is Gaussian with mean exp(-a_Z*(u-T))*Z_T and variance nu_Z^2*(1-exp(-2a_Z*(u-T)))/(2a_Z). Baseline three-projection marginal/conditional integrals become FOUR-dimensional, although C or Z may be integrated using reusable tables. A direct tensor rule adds approximately one quadrature-order factor of cost (e.g.32x GH nodes conditionally); this is an operation count, not a measured runtime. Outer importance sampling may need to target joint high-response states even with4096 points; existing proposal efficacy cannot be assumed.

Stock cache contains about3.496 billion step/path entries across all seven maturities. Caching one extra float64 Z value per entry adds26.05GiB; streaming Z avoids that storage but requires reproducible seeds/bridges and full stock recomputation for every changed response. New randomness does not require resimulating the unchanged exogenous reservoir, but neither volatility nor terminal stocks are reusable unchanged. Reference-control and hedge/Gaussian zero means survive independent augmentation, although variance reduction may deteriorate. Rebuild candidate-specific VIX exposures/forwards/option payoffs/normalizers and check new sampling, GH, time and stock errors.

G_Z remains bounded by exp(2L), so integrability and the bounded-volatility stock-martingale argument survive this augmentation. Simply using exp(Z)*(F-betaJ)^2 WITHOUT tempering could produce nonintegrable Gaussian/lognormal mixtures; the bounded sensitivity and outer cap avoid that issue. A future log-continuation adaptation can also preserve at-most-quadratic growth when sensitivity A is bounded.

## A second dedicated-layer option

A smooth mixture weight w(Z) can replace w(C) in the implemented two-component mixture. This allows changing mixture persistence independently of the cluster reservoir. It keeps bounded G, nests a constant mixture at zero modulation, and needs the same added state/pricing work. It differs from direct qJ*A(Z): it interpolates between two response functions instead of continuously changing the curvature. Compare these two mechanisms at matched mean variance and forwards, not by assigning the label 'stochastic vol-of-vol' and assuming success.

No dedicated-state simulation or calibration has been launched. The existing-cluster diagnostic is a bounded test, not an adopted architecture. A poor cluster profile does not reject the dedicated layer, whose new persistence and response dependence differ.
