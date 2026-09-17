# Stage-4 targeted fine-grid M2 recalibration

This experiment keeps the Stage-3 OU rates and projection weights fixed, disables the rejected
echo and orthogonal blocks, and recalibrates only the nonlinear M2 readout and the equity
correlations of the feedback and spike banks. Every simulation uses 64 internal steps per
trading day.

The design has three explicit boundaries:

1. Four search seeds rank 65 readouts (the preceding candidate plus 64 Sobol points) under eight
   correlation scenarios.
2. The primary candidate and seven sensitivity candidates are frozen using search results only.
3. Five untouched seeds validate the frozen candidates. Validation may reject the primary, but
   it may not silently replace it with a runner-up.

The seven Stage-2 proxy gates are unchanged. Both the log-return and the Itô-relative-return
versions of the Zumbach statistic are reported; the latter is a robustness diagnostic and is not
an additional calibration gate. The gates remain synthetic capacity proxies rather than
empirical SPX confidence bands.

Run from the package root with:

```bash
python experiments/stage4_fine_grid_001/run.py
```

The next experiment is a nested 64-versus-128 time-step confirmation, but only if the primary
candidate passes the untouched validation.

The preceding M2 candidate is also run on the validation seeds as an explicitly post-selection
benchmark. It cannot replace the frozen Stage-4 primary.
