# VIX-focused constrained search — supplement after rejected proposal

The first local joint-merit proposal in run001 had true common-call loss79.216
vs98.119, but native OTM VIX-plus-forward loss worsened by9.060. The acceptance
safeguard rejected it. Stop001 with all records intact, before checking prices;
it is not a completed calibration or an adopted candidate.

Run002 retains all original eleven parameters/bounds, beta_J>0, fixed OU/Q/state,
market scope,4096 samples, frozen controls and quadrature/time settings. Same
independent checking cache remains unseen for parameter selection.

Use exact nonlinear VIX and first-order SPX inside a trust box. Minimize the
average common-call/native-OTM VIX-option loss plus the forward loss. Constrain
both VIX-plus-forward criteria and both joint criteria to improve at the local model. Reserve2% of each anchor joint loss and0.5% of
each VIX-plus-forward loss at radius0.1, scaled linearly with radius, with
minimum1e-4. This leaves room for SPX approximation error before true acceptance. Analytic gradients, parameter units scaled by SCALE;
SLSQP with at most12 inner iterations and20 function evaluations per outer
proposal. Retain the best feasible inner point if its evaluation budget expires.
This is guidance for candidates, not a new reported price objective.

At most3 outer proposals, initial radius0.1*SCALE, shrink0.5 on true-price
rejection and grow1.5 up to0.3 on accepted agreement. Maximum75 optimization
minutes, checked between evaluations. Every acceptance still requires FULL true
joint4096 repricing and improvement in both joint and VIX-plus-forward criteria.
The prior fitted vector stays incumbent until those conditions pass.

Reuse the production Jacobian correspondence established in001, as the true
pricing/Jacobian code is unchanged. Add tests/production checks for the native
VIX gradients and constrained merit/criteria. Save full price/residual Jacobians
at every true anchor for interruption recovery. Same independent4096 checking,
all-maturity plots and uncertainty/tail/shape audits after final selection.
No model or objective-convention change is concealed, and no global optimum
or exact-market-fit claim is asserted.
