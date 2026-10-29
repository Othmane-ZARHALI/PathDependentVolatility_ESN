# Current ESN and M2 research state

Snapshot: 29 October 2026. This brief supersedes the long continuation instructions for navigation; original reports retain their historical findings and limitations.

## Objective and present conclusion

Seek an interpretable common SPX/VIX model, retaining as much of exogenous M2 as possible. Two stress speeds gave the strongest partial VIX fit. Correlating stock with stress improved SPX while protecting that VIX law. Subsequent curve, asymmetric readout, slow-signal, linear-filter and protected-generator tests have left substantial SPX downside-IV errors. No satisfactory joint calibration or adoption of the latest candidates is established. See [the history](EXPERIMENT_HISTORY.md) and [figures](VISUAL_EVIDENCE.md).

## Model and reference roles

The original M2 uses an exogenous Gaussian OU bank and a capped exponential-quadratic response. In schematic notation,

\[
\eta=c_CC+q_Y(Y-\beta_Y)^2+q_J(F-\beta_J)^2,
\quad R_t=\frac{e^{2h(\eta_t)}}{E_Q[e^{2h(\eta_t)}]}.
\]

Here F is the spike factor. Positive `beta_J` corresponds to positive `spike_shift`; retain `q_J>0`. Initial states, rates and other fixed choices are experiment-specific and must come from the chosen model artifact.

The two-speed experimental variance model is

\[
v_t=f+(\xi(t)-f)H_t,\qquad
H_t=(1-w)R_t+w[(1-r)A_t^s+rA_t^f],
\]

where each independent Gaussian-OU exponential stress response A has mean one. The VIX-only fitted reference has stress weight about 0.922838, half-lives 365 and 71.3905 days, and volatility floor 12%, so f=0.0144. Its artifact is [experimental_model_two_speed.json](../../matteo/model_comparison/option_pricing/runs/m2_two_speed_001/experimental_model_two_speed.json). This is an experimental reference, not a newly adopted production model. Later SPX tests generally freeze this variance law and use stock/stress rotation theta=0.8; read the individual report for its precise baseline and reused auxiliary reference.

The earlier physical-measure Stage-3 baseline, original pricing anchor, VIX-only two-speed reference, joint-refit candidate and later SPX experiment baselines have different roles. Do not select one by recency or compare their losses across different estimators and scopes.

## Active market scope

Use [active_scope.json](../../matteo/model_comparison/option_pricing/active_scope.json), checked during consolidation:

- SPX: 29/57/119/211/302/666/1030 days; 1751 eligible OTM quotes, 1294 primary quotes under the user-selected cutoff 3. Other eligible wings remain diagnostics.
- SPX excluded maturities: 1/1401/1765/2129 days.
- VIX: 27/55/118/153 days; all 217 eligible OTM quotes and all four forwards.
- Current feasibility work uses all eligible original market splits as inputs. Old 175-VIX-training-quote and held-out-market instructions describe earlier phases. Fresh simulations check the same market snapshot, not market out-of-sample performance.

The market snapshot underlying the visual retrospective is 19 February 2026. Preserve native quotes, signed estimates, unavailable IVs and full maturity membership. No clipping or smoothing.

## Findings that should guide the next iteration

1. A frozen variance budget can be incompatible with VIX quote moment bounds. Extra response shape cannot repair that budget by itself.
2. Two stress speeds strongly improved VIX checking loss; adding a continuous mixing measure or four fixed speeds gave no demonstrated further advantage in the tested searches.
3. The bounded joint fit largely removed the successful stress contribution and damaged VIX. That is a local empirical tradeoff, not a model-family impossibility theorem.
4. Off-window curve changes preserve VIX but did not repair SPX shape. A lower price objective can accompany worse IV error.
5. Stock/stress rotation improved SPX without changing the fixed exogenous variance-factor law. Adaptive directions did not provide a convincing further improvement.
6. Asymmetric response changes and slower signals can affect SPX but also damage VIX. Unit means and fitted training forwards do not protect the VIX distribution.
7. Linear filtering has been implemented and tested; it produced no material improvement. It is no longer a pending experiment.
8. Three protected-generator correction studies have been tested. Their quoted-date VIX protection and positivity arguments are useful; empirical SPX progress is modest or unestablished. The stress-budget candidate gives 29-day downside IV about 19.26% versus market 26.90%. Its paired improvement interval includes zero.

## Mathematical protection and numerical limits

With unchanged factors and pricing information, a correction h=(partial_t+L)P whose potential vanishes at both ends of each quoted VIX window has zero conditional window integral. It can preserve model VIX values at those four dates state by state. This inherits baseline pricing errors and does not protect every unquoted date. The first two constructions guarantee v>=f/2; the latest stress-budget construction guarantees the original v>=f. These are different bounds.

The useful endpoint/path distinction is: VIX options constrain future VIX values, while SPX also depends on variance along stock paths and its interaction with shocks. Finite quotes do not determine the full distribution; the markets still share a variance process and budget. The tested constructions have not shown that this distinction supplies the missing empirical fit or is fundamental to an eventual M2 model.

Exogenous OU states and random shocks can be reused where valid. Readout, cap, curve, filtering or variance changes require rebuilding affected normalizers, conditional VIX exposures and stock paths. Stock-driver changes require rebuilding stock outcomes. Later endogenous-feedback external models have different reuse rules.

Keep sampling, tail coverage, stock timestep, normalizer and conditional-VIX integration errors separate. Numerical audits and scoped Lean identities do not certify the full simulator, continuous-stock martingality or convergence. Later Q variants have not re-established the original P fit.

## Continuation

No new experiment is launched or selected by this consolidation. Before proposing one, name the unresolved feature, the economic quantity it changes, what it preserves, and how it differs from the completed tests. Check [decisions and open questions](DECISIONS.md), then the relevant source report. Preserve the unsuccessful results and their figures. Quintic and jump remedies remain outside the completed M2 joint-fit programme unless explicitly reopened.
