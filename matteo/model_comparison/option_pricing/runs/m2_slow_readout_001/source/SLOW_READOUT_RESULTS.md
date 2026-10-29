# Slower downside signals from the existing OU bank — 28 October 2026

The slower signals produce a modest SPX improvement, but this particular multiplicative readout does not solve the joint fit. Strong signals improve SPX more while damaging VIX. A weak signal avoids most of that damage but leaves much of the SPX skew missing. Production model and previous results are unchanged.

## Model in plain terms

The existing fast signal continues to react to recent shocks. Add a separate slower signal, so that a negative shock can raise volatility for longer. No OU process is added or given feedback from the simulated price or volatility. Select existing feedback modes nearest in log half-life to generic targets 2, 8 and 32 days. Their actual half-lives are 3.329, 9.481 and 27.004 days; the same construction applies at every option maturity.

Let U_j be the selected existing mode and X_j=U_j-E_Q[U_j]. Each mode has stationary variance one; its dated conditional variance is not assumed to be one. For one mode at a time,

$$dX_j=-\kappa_jX_jdt+\sqrt{2\kappa_j}\,dB_j,\qquad h_j=\log(2)/\kappa_j,$$

$$D_j(X_j)=\exp\{\alpha\log(1+e^{-X_j})\}.$$

D increases as the state becomes more negative. Alpha controls its strength. The current fast M2 component, including its asymmetric positive spike shift, stays exactly in place:

$$\eta=c_CC+q_Y(Y-\beta_Y)^2+q_J(F-\beta_J)^2,\quad h=\ell-k^{-1}\log(1+e^{k(\ell-\eta)}),\quad R=e^{2h}/E_Q[e^{2h}].$$

Keep both lognormal stress factors and their frozen parameters, with A=(1-r)A_s+rA_f and a=1-w. The old multiplier is H0=aR+wA. The new one is

$$m_D(t)=E_Q[D_j],\qquad m_{RD}(t)=E_Q[RD_j],$$

$$H_j=(aR+wA)D_j\frac{a+wA}{a\,m_{RD}(t)+wA\,m_D(t)},\qquad v_t=f+(\xi(t)-f)H_j.$$

R and D are correlated: E[RD] cannot be replaced by E[R]E[D]. The whole old Gaussian bank remains independent of the stress bank. Consequently,

$$E_Q[H_j\mid A_s,A_f]=a+wA,\qquad E_Q[H_j]=1,\qquad E_Q[v_t]=\xi(t).$$

This preserves the conditional mean variance allocated to each stress regime. It does not preserve separate internal component budgets or the entire VIX law. VIX depends on the conditional future variance and a square root, so matching E[v] does not by itself match VIX forwards or options.

Stock loading stays fixed at the previous mixed theta=.8 value:

$$dW^S=.6\,dW^0-.8\{\sqrt{1-r}\,dB_s+\sqrt r\,dB_f\}.$$

The original M2 cap is unchanged. No new cap or maturity-specific exception is introduced. Reusable Gaussian factors are distinguished from rebuilt readouts, normalizers, stock paths and VIX exposures.

## Experiment and controls

First test alpha .5 and 1 for all three memories. Their VIX training losses were too large, so add a reactive exploratory grid .125 and .25. These grids were not jointly preregistered. Both stages share identical 512-path stock draws, eight independent antithetic batches, and identical 512-state-per-expiry training samples. Evaluate the fixed variance curve and refit the same seven active variance levels under four training-forward equalities; nine off-window levels remain frozen. All seven SPX maturities/1751 eligible/1294 primary quotes and all 217 VIX quotes/four forwards remain.

Separate previously examined 4096-state-per-expiry VIX checking is saved for every candidate. A weak candidate is frozen before additional new stock and outer-state draws: require training-forward equality and training VIX loss within .1 of the alpha-zero refit, then minimize the pilot's fine-grid SPX IV MAE among weak refits. This selects the 9.481-day mode, alpha=.125. The .1 threshold is a diagnostic tolerance, not evidence that the resulting VIX smile is satisfactory.

Additional checking uses 1024 newly generated stock paths and 2048 new full-factor importance-sampled VIX states per expiry. The five comparisons are baseline, selected slow signal with fixed curve, selected slow signal with refitted curve, that same refitted curve without the slow signal, and a stronger 3.329-day alpha=.5 refitted contrast. The strong contrast is illustrative, not a second newly selected winner. All four training forwards satisfy numerical equalities; their separate checking gaps below are retained.

## Independent checking results

| Configuration | SPX interval loss | SPX IV error, pp | VIX interval loss | Max VIX forward gap, points | Invalid primary SPX IVs |
|---|---:|---:|---:|---:|---:|
| Baseline | 1866.88 | 3.176 | 0.050 | 0.0552 | 1 |
| Slow signal, fixed curve | 1816.32 | 3.130 | 0.089 | 0.0433 | 2 |
| Slow signal + curve refit | 1811.15 | 3.128 | 0.085 | 0.0174 | 1 |
| Curve refit alone | 1862.10 | 3.174 | 0.041 | 0.0501 | 0 |
| Stronger signal + refit | 1749.73 | 3.030 | 0.553 | 0.0616 | 3 |

IV error equally averages available native-IV errors over seven maturities. Invalid IVs and signed negative price estimates are retained; native-IV metrics alone are not a fit certificate. Forward gaps include the inherited baseline estimator gap. Training equalities are not certified exact population forwards.

For the weak refit, SPX price loss falls 2.99 percent in this checking sample. The curve-only control and fixed-curve signal show how much comes from each change. At 29 days, the quoted right-minus-left IV gap near strike/forward .9 and 1.1 is -5.98, -6.17, -6.17, -5.98, -7.05 percentage points, respectively, for baseline, weak fixed, weak refit, curve alone and stronger refit, versus -14.05 in the market. The required downward slope remains missing.

Paired sampling intervals: [{"candidate": 1, "loss_delta": -50.56199836400674, "approximate_batch_bootstrap_CI95": [-85.63739866719007, -5.20013727585913]}, {"candidate": 2, "loss_delta": -55.73565137293258, "approximate_batch_bootstrap_CI95": [-95.70683002007192, -7.822502998978674]}, {"candidate": 3, "loss_delta": -4.784891297016429, "approximate_batch_bootstrap_CI95": [-9.1857322818127, -0.3899292965022488]}, {"candidate": 4, "loss_delta": -117.15274754054349, "approximate_batch_bootstrap_CI95": [-341.2378229603223, 162.7598610977011]}]. These approximate 1024-resample bootstrap intervals use eight antithetic batch means with paired candidate/baseline resampling and independent resampling of the reused cheap reference batches. They exclude deterministic quadrature/grid errors and market uncertainty. New full-model draws do not make the already reused 65536-path auxiliary reference independent of prior research. This is the same market snapshot, not market out-of-sample validation.

## Numerical correspondence and accuracy

Conditional moments use the full Gaussian state law projected onto C,Y,F and the chosen slower mode. Integrate independent C by a 4096-interval cubic-Hermite convolution of the capped readout, Y/F by tensor Gaussian quadrature, and X conditional on Y/F by one-dimensional Gaussian quadrature and an 8192-interval interpolation table. Direct one-dimensional marginal quadrature computes E[D]. This preserves the R/D correlation. Product normalizers use GH48 and adaptive logarithmic-time cubic splines; conditional VIX checking uses GH32 and time order8 with positive stress quadrature32.

Alpha-zero full-horizon stock replay error is 2.22e-16. Independent scalar conditional tensor quadrature relative errors for the three memories are [5.969134038164015e-07, 5.651782006843575e-07, 3.143547642192601e-07]. GH48/64 product-normalizer checks are below 6.66e-07. An independent 97-step scalar stock recurrence for five candidates, both grids and both antithetic partners has error 0.0. These are numerical implementation checks, not Lean certification or uniform error bounds.

Coupled stock grids use 128 and 256 steps per trading day; their median absolute price changes in quote-scale units are 0.117, 0.135, 0.128, 0.118, 0.183. A matched 64-state-per-expiry conditional VIX refinement, GH32/time8 to GH48/time12, gives maximum VIX change 7.321e-06 points and maximum absolute weighted mean change 1.442e-07 points. These probes do not bound unobserved tails. Raw stock payoff concentration, importance-weighted VIX tail concentration, sample mass, stock means and SEs, quote membership, forward SEs and negative/invalid estimates are preserved in the run. No clipping, smoothing, sample-weight normalization or sample-stock-forward renormalization.

Use the previously validated stress-only auxiliary stock reference with its conditional lognormal call prices and fixed deterministic delta controls. Adding the new readout changes the expensive target payoff correction, not the auxiliary model. Forward and option checking SEs combine the independent fresh difference and old anchor SEs; separate components remain available.

## Durable conclusion

The slow-mode hypothesis has now been tested for three generic existing memories and four strengths, with fixed-curve, constrained-curve and curve-only controls. This simple one-sided multiplier supplies some downside response, but stronger changes also alter the VIX distribution; protecting the first moment does not protect option smiles. Do not interpret the weak aggregate improvement as a solved SPX skew or repeat the same strength scan as the next default step. The result is specific to this readout and tested range, not a no-go theorem for slower projections or exogenous Gaussian banks.

No candidate is adopted. The production baseline and original fast/Zumbach-motivated readout remain intact. No endogenous return/volatility feedback, new OU bank, new constant-theta scan or quintic/jump extension has been introduced.

Figures: [all seven SPX smiles](visualizations/slow-readout028/SPX-slow-readout.png), [all four VIX smiles](visualizations/slow-readout028/VIX-slow-readout.png), [full exploratory trade-off](visualizations/slow-readout028/tradeoff.png). Data: `runs/m2_slow_readout_001` and `runs/m2_slow_readout_002`, with independent checking in the latter's `fresh_validation` folder.
