# Nested M2 time-step convergence audit

This experiment keeps the selected continuous-time M2 architecture and readout fixed while
changing only the internal integration resolution. It is a numerical audit, not a calibration.

The finest exact OU/equity path is generated at 64 steps per trading day. The 2-, 4-, 8-, 16-,
and 32-step paths use exact subsamples of those OU endpoints, while their equity Brownian
increments are sums of the same fine increments. Every level therefore shares the same
underlying path and retains 252 daily output observations per year.

The experiment uses five declared seeds, 12 paths per seed, six simulated years, and a one-year
burn-in. It reports the full stylised-fact vector, joint proxy acceptance, cap activity, and
pathwise errors in daily returns and daily average variance relative to the 64-step reference.

Run from the package root:

```bash
python experiments/timestep_convergence_001/run.py
```

Outputs:

- `per_seed.csv`: one row per seed and resolution;
- `summary.csv`: across-seed means and standard deviations;
- `manifest.json`: frozen architecture, readout, rate audit, and nesting convention.
