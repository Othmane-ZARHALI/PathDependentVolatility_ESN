# Frozen hedge controls for the parity/path-budget experiment

Supplement to PARITY_ADAPTIVE_PLAN.md, before fresh production results.
A replay on OLD selected 57-day paths showed that the simple lognormal-reference
control alone increases SE relative to the old richer controls, particularly near
strike 8000. Preserve that finding; do not call a new estimator better solely
because its sample count is larger.

Add known-zero stock, deterministic-lognormal-reference-call and two predictable
Black-inspired hedge controls, the same mathematical control class as before.
One CALL family at nine market-only anchors spanning both sides of the forward;
interpolate random correction draws continuously, with no separate put branch.
Hold the hedge delta between approximately four trades per trading day, evaluated
on the unchanged 512-step-per-trading-day stock path. Delta uses only adapted
state/stock at the trade time. Its gain therefore has exactly zero target mean
under the adapted discrete stock martingale. This is variance reduction, not a
new stock dynamics approximation.

Freeze control model at the previous selected vector, and fit all control
coefficients on NEW independent 256-path SPX pilots. Production primitives are
not used to fit or select coefficients. Save pilots, coefficients, all frozen
control draws and provenance. Reuse these controls unchanged for both old and
selected parameter repricing and for subsequent optimization. Analytic derivatives
must treat these frozen correction draws as constants.

Fresh path preparation already running under V1 is unchanged: that cache stores
only readout-independent OU paths/equity increments/weights. This V2 supplement
changes its subsequent pricing/control layer, not its captured cache inputs.
