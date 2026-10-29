# Cap and tail investigation, 25 October 2026

Authorized sequence: cap/tail first, stochastic vol-of-vol afterwards. No jumps or quintic proposal. Preserve selected maturity-curve baseline and all historical runs. Experimental engines only; no production change or accepted calibration.

Fixed baseline: `m2_maturity_curve_fit_001/selected_model.json`, all 16 xi and floor frozen during profiles, beta_J positive. OU architecture, initial state, correlations, Q scenario, time mesh and importance proposals unchanged. All original 1511 OTM quotes, cutoff3, seven SPX/four VIX dates, four forwards, common calls, every4096 samples. Existing pilot controls frozen.

Bounded first pass: baseline; L +/-0.75; kappa=1 and4; algebraic saturation matched at transition; logarithmic continuation strengths0.25 and0.5; separate feedback/spike caps3/2.5 with linear cluster outside. Full4096 VIX pricing for all, rebuild GH48 Q normalizers, GH32 conditional/time8 with all curve splits. Score smile level/slope/curvature and forwards; raw likelihood-weighted tails and eight-net uncertainty separately from controlled prices. No clipping. Reprice all seven SPX maturities for at most two interesting candidates selected using TRAINING VIX; separate reused VIX streams afterwards. This is a profile, not optimizer convergence. Baseline replay first.

New scalar-cap convolution uses cubic Hermite interpolation of the SAME GH32 C sum, step0.01, checked against direct quadrature; numerical approximation is recorded. Separate-channel candidate factorizes the independent C Gaussian sum. No derivative-based fitting with these value-only experimental kernels.

After cap profiles: compare existing cluster response modulation with a dedicated continuous latent state mathematically and by a bounded additional diagnostic if supported by findings. Report cache/cost/integrability/martingality constraints. Dedicated state is an addition, never a change of existing lambda/grid.

Save configurations, source hashes/snapshots, raw draws, weighted contributions, replay tests, negative outcomes and plots in new001 run and cap-tail014 visualization directory. Used checking streams are numerical checks, not pristine holdout or market OOS. Finite-GH/time/normalizer/tail and SPX timestep uncertainties remain unresolved.
