# Prespecified five-piece refinement — 23 October 2026

The user asked whether more than three deterministic variance levels could help
while the seven-state joint search was running. This plan is written before its
fresh checks are inspected and before any five-piece optimization.

Finish the existing bounded state experiment. Then hold its search-selected
restricted state, structural/readout parameters and Q scenario fixed. Compare
the three-level curve with a nested **five-level** family, ending at
`[57/365, 0.25, 151/365, 0.75, 5.832876712328767]` years. The old three levels
`[a,b,c]` embed as `[a,a,b,b,c]` exactly. No old candidate is removed.

Numerical grid convention, declared before fitting: use the exact stored
57-day SPX event and snap day 151 to the nearest existing 512-step stock-grid
point. Record the actual dates and distance. This keeps the primitive stock
grid unchanged, respects every new knot exactly, and avoids treating a hidden
new stock discretization as a fitting improvement. The rounding uses dates
only, with no inspection of prices or scores.

Why these two additional dates: the first VIX window is days 27–57 and the
second is 55–85, so day 57 separates most of them. The later windows are
118–148 and 153–183, so day 151 separates them. These are date-based design
choices, not optimized knots. Existing 91.25/273.75-day boundaries remain.

All 707 cutoff-3 SPX options, all 175 VIX options and all four forwards enter
the unchanged fixed-scale price objective. No reserved/diagnostic SPX maturity
is added. Bounds remain [0.00010001,0.25] per variance level. Retain the existing
`0.1*sum(diff(xi)^2)` penalty over adjacent levels; it penalizes the new level
differences and is unchanged on the embedded old curve. This is weak smoothing,
not proof against overfitting. Report its magnitude separately.

Declare two curve starts (embedded three-piece selection and a flat curve),
at most 30 least-squares evaluations per start and a 30-minute optimization
budget after exposure construction. Select only on searching simulations.
Use independent checking draws with shared primitives for the embedded baseline
and new candidate, and explicitly record any checking draws reused across model
comparisons. Do not describe shared or previously generated checks as a new
independent experiment. Preserve all negative estimates and failed points.

At fixed state the multiplier ratio and primitive stock shocks remain reusable,
but candidate stock paths must be repriced. Segment membership changes. VIX
time windows must split at new knots; recompute the five-column exposures and
compare the embedded old curve with the prior three-piece exposure calculation.
That difference is a time-quadrature diagnostic, not an optimization gain.
Keep normalizer interpolation, GH, time quadrature, stock grid and sampling
errors separate. State calibration and this refinement remain exploratory.

A richer deterministic curve can improve term structure and VIX windows, but
cannot freely adjust every strike's smile shape. A lower searching loss is
expected from a nested family; fresh price/forward behaviour determines whether
the extra flexibility appears useful, without establishing a qualified fit.
