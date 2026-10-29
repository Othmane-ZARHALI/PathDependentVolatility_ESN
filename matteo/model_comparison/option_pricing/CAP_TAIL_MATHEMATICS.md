# What the cap does, and continuous alternatives

25 October 2026. Equations below are analytical arguments, not new Lean certificates. Numerical results are in CAP_TAIL_RESULTS.md. The original selected model remains unchanged.

## Verified production equations

Let C,Y,F be the fixed OU reservoir projections (cluster, feedback, spike). The active channels are

    eta = c_C C + q_Y (Y-beta_Y)^2 + q_J (F-beta_J)^2,
    h(eta) = L - log(1+exp(k*(L-eta)))/k,
    G = exp(2h), M(t) = E_Q[G_t],
    v_t = f + (xi(t)-f) G_t/M(t).

Source correspondence: package `esn_eq/readout.py` and `pricing.py`; current fast `controlled_readout_convolution_kernel.cpp`; `maturity_curve.py`, `maturity_curve_stock.cpp`, `maturity_curve_predict.py`. The direct `full_readout_kernel.cpp` has the same mathematical cap. Production convolution exploits independent C to accelerate the SAME tensor GH32 rule; GH48 normalizes, spline interpolation evaluates date-zero M between time knots. Optional echo/orthogonal readout coefficients are zero. The package has exponent clipping at +/-700; it is inactive for these bounded positive-cap parameters. Experimental kernels have no such clipping.

Selected L=5.5027800568, k=2; beta_J=0.3697299725>0; f=.00264. The 16 positive curve levels satisfy xi_i>=f without monotonicity. Hence sqrt(v)>=sqrt(f)=5.1381%. This is an instantaneous floor, not an option-IV floor.

VIX is a conditional average of future VARIANCE, followed by a square root:

    VIX_T^2 = 10000/delta * integral_T^(T+delta) E_Q[v_u | full state at T] du.

The full OU state is Markov; C/Y/F alone need not be. Conditional means use full states. Forward-curve identity: E[VIX_T^2]=10000/delta * integral xi(u)du. No extra fit parameter changes this identity.

## Range, response and weights answer different questions

Because h<=L, G<=exp(2L). Therefore

    VIX_T <= 100 sqrt( f + 1/delta * integral (xi(u)-f) exp(2L)/M(u) du ).

This is a deterministic upper bound (numerically integrated in this study), not a proof the bound is attained: correlated OU paths and mean reversion can prevent simultaneous maximizing states. An upper bound above a strike alone does not establish attainability. Here saved conditional-pricing states actually produce VIX above the highest strikes; continuity and Gaussian support give the corresponding reachability argument within the implemented approximation. This says nothing by itself about option prices.

An option call price is E[(VIX-K)+]. For proposal draws with likelihood w=dQ/dPproposal, its raw estimate is mean[w*(VIX-K)+]. Count hits, sum likelihoods of hits, average excess above K, and distribution of WEIGHTED payoff contributions separately. A hundred extreme proposal draws can still represent a very small risk-neutral probability. `payoff_ESS=(sum contributions)^2/sum squares` is only a concentration descriptor; eight scramble means give the reported sampling SE. Neither certifies unobserved tails. Control corrections are separate signed zero-mean terms and are kept in actual prices.

The local response to the score is

    h'(eta)=1/(1+exp(k*(eta-L))),
    d log G/d eta = 2h'(eta).

At eta=L the response is half its uncapped value. Under a parameter theta,

    d log(G/M)/dtheta = d log G/dtheta - E_Q[G*d log G/dtheta]/E_Q[G].

The second term is essential. Raising L increases every raw G, but normalization redistributes relative variance toward states where G increases more. It need not increase every option price or VIX forward. In particular dlogG/dL=2(1-h'), so the marginal G-weighted response directly measures how strongly normalization reacts. Jensen also matters: fixed E[VIX^2] does not fix E[VIX]; redistributing variance toward rare extremes can lower the forward.

## Cap candidates actually profiled

All original score/curve coefficients are frozen in the initial profiles. Smoothness means smooth response, not smoothed plotted IVs.

1. Original logistic cap with L shifted +/-0.75, or k changed from2 to1/4. L controls transition location and plateau; k controls transition width.
2. Algebraic saturation, with d=2log(2)/k:

       h_A(eta)=L - ((L-eta)+sqrt((L-eta)^2+d^2))/2.

   Same h(L)=L-log(2)/k and h'(L)=1/2; convergence to the plateau is algebraic rather than exponential. Still bounded above by L.
3. Logarithmic continuation, u=softplus(k*(eta-L))/k:

       h_a(eta)=h_original(eta)+a*log(1+u), a=.25 or .5.

   Beyond the transition, G grows like exp(2L)*eta^(2a), rather than saturating. This is NOT an uncapped exponential-quadratic tail. All Gaussian moments of this G are finite. Since eta is at most quadratic in the finite Gaussian state, a<=.5 gives v<=K_T*(1+||X||^2) on each finite horizon. This growth control is useful for stock martingality as explained below. Negative eta has decaying exponential response.
4. Separate channels:

       G=exp(2c_C*C +2h_LY(qY*(Y-betaY)^2)+2h_LJ(qJ*(F-betaJ)^2)),
       LY=3, LJ=2.5, k=2.

   Quadratic channels are individually capped; the linear Gaussian cluster channel is outside those caps. All Gaussian multiplier moments are finite. This particular allocation is one diagnostic, not evidence that every channel-wise design fails.

These are proposed M2 adaptations. They are not formulas asserted by the literature below. The new value-only kernels do not supply readout gradients and must not enter gradient-based readout calibration. Curve/floor sensitivities remain algebraically valid. Cubic interpolation in numerical quadrature tables is an integration approximation, not a polynomial model proposal.

## Why removing the cap needs a calculation

For X=m+RZ, Z standard normal, let A=diag(0,qY,qJ). Integrability of exp(2p*eta) requires

    I - 4p R^T A R positive definite.

The quadratic exponent from the Gaussian density must still point down in every direction; otherwise the tail integral diverges. Shifts and linear terms change the finite value but not this positive-definiteness threshold. At the selected coefficients, the stationary largest eigenvalue of 4R^TAR is0.8970204882. Covariance at a finite date is below stationary covariance, so the FIRST uncapped moment is finite at every date. At long maturities only multiplier powers p<1.114801739 remain integrable. Its second moment is infinite there, including tested VIX-window midpoints. Thus 'uncapping always makes M infinite' would be false for THIS selected readout. Nearby larger curvatures can cross the first-moment threshold, and first-moment existence alone does not establish numerical stability or a stock martingale.

## Stock martingality and numerical stability

For the bounded original/algebraic caps (and bounded mixtures below), M(t)>0 is continuous and bounded away from zero on compact time intervals. With bounded xi and xi>=f, volatility has a deterministic finite bound on such intervals. Novikov's sufficient condition holds; the discounted continuous stock exponential is a true martingale.

For logarithmic continuation a<=.5, sqrt(v) has at most linear growth in the unchanged finite OU state. Under the stopped stock-density measure the state's drift is OU drift plus its covariance with equity noise times sqrt(v), hence still at most linear. A stopped second-moment/Gronwall bound gives nonexplosion and uniform exit control on finite horizons, establishing the usual density-martingale argument. This uses local regularity, nonnegative xi-f, strictly positive floor, fixed finite state and positive continuous exact M. It is an analytical argument, not a Lean-checked theorem or a proof that finite-GH/spline pricing is exact.

For separately capped channels with uncapped linear C, C is independent of the equity/feedback/spike driving noises in the verified baseline. Conditional on a continuous C path, the two other channels and hence variance have a finite deterministic time-dependent bound on a finite horizon. The conditional stock exponential has mean one; integrating over C preserves it. This argument would need reconsideration if cluster-equity correlation changed.

For uncapped exponential-quadratic G, finite E[G] alone provides none of those bounds; stock martingality remains unproved here. It is not adopted or priced as a solution.

Every adapted discrete stock step also has conditional mean-one exponential update under the exact Q shocks, since left-endpoint v is known before the step's Gaussian equity innovation. This finite-grid identity does not bound timestep bias. Weighted sample stock means can be noisy.
