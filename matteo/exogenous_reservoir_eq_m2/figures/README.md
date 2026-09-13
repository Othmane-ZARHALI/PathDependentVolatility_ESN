# Visual diagnostics

Figures `01`-`04` and `06` document the preserved zero-shift Stage-3 M2 baseline.
Figure `05` is now a paired event study comparing that control with the recommended
non-zero fast-factor shift. Figures `07` and `08` explain the asymmetric mechanism and
show its trade-off with the original joint stylised-fact gates.
The HTML dashboard combines the three new comparison figures with the five baseline views.

The comparison uses Stage-3 candidate 26, structural scenario 0, because it passes all
seven original gates at `spike_shift=0.10` on every one of the five fresh seeds. Control
and shifted simulations reuse identical latent paths, and event dates are selected once
from the control returns.

Regenerate the asymmetric figures from the package root with:

```bash
python figures/generate_asymmetry_diagnostics.py
```

The plotted volatility quantity is `100 * sqrt(daily average instantaneous variance)`.
It is not model-implied VIX and the figures are synthetic mechanism evidence, not an SPX
or VIX market calibration.
