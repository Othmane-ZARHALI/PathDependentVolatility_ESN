# Exogenous-reservoir exponential-quadratic volatility: M2 Stage 3

This is the standalone, executable M2 research package. It started as a curated copy of the
`esn_eq_model_v0_3_0_stage3.zip` snapshot and is now the active development package. The
source archive remains unchanged in the sibling
`../m2_exponential_quadratic_stage3/` folder.

## Scope

The recommended baseline is M2: a structured exogenous Ornstein-Uhlenbeck reservoir with a
smoothly tempered exponential-quadratic variance readout. The package also retains the paired
M3, M2+O, and M4 ablations. The checked-in numerical evidence is synthetic,
daily-frequency, physical-measure capacity evidence; it is not an SPX/VIX option calibration.

Pilot-1 result files and its old archive are intentionally excluded. One small Pilot-1-derived
file is retained in `provenance/latent_calibration.json`, because it freezes the latent
architecture used by the Stage-2 and Stage-3 studies. This is provenance input, not the
Pilot-1 experiment ledger.

## Layout

- `src/esn_eq/` — version 0.4.0 model, simulation, diagnostics, calibration, and experiment CLIs.
- `tests/` — unit and integration tests included in the archived snapshot.
- `experiments/stage2_001/` — Stage-2 training/probe record.
- `experiments/stage3_001/` — Stage-3 held-out confirmation record.
- `experiments/asymmetry_001/` — paired fresh-seed test of the fast-factor shift.
- `experiments/EXPERIMENT_COMPARISON.md` — compact experimental conclusion.
- `documentation/` — the model specification (`.tex` and `.pdf`) and implementation guide (`.pdf`).
- `figures/` — M2 diagnostics figures, dashboard, summary, and reproduction script.
- `research_notes/` — technical review transcript.
- `provenance/` — the frozen latent-calibration input required to rerun Stage 3.
- `ARCHIVED_README.md` — verbatim README from the archived v0.3.0 package.

Version 0.4.0 adds a nested `spike_shift` readout parameter. A zero shift exactly recovers
the Stage-3 fast quadratic channel; a positive shift makes negative fast-factor shocks produce
larger volatility responses than equal positive shocks. The command `esn-eq-asymmetry` runs
the paired fresh-seed screen of this mechanism.

## Installation and checks

Use Python 3.10 or newer. From this directory:

```bash
python -m pip install -e .
python -m unittest discover -s tests -p 'test_*.py' -v
```

## Reproducing Stage 3

The archive's Stage-3 command expected the full Pilot-1 directory. This curated package keeps
only the frozen latent-calibration file, so explicitly point `--pilot` to `provenance`:

```bash
esn-eq-stage3 --pilot provenance --training experiments/stage2_001 \
  --output experiments/stage3_rerun
```

Run the focused downside-asymmetry study with:

```bash
esn-eq-asymmetry --output experiments/asymmetry_001
```

Do not overwrite `experiments/stage3_001/`; it is the preserved confirmation record. The
Stage-2 source remains included because Stage 3 imports shared analysis helpers and uses the
Stage-2 record as its fixed training reference. Reproducing Stage 2 itself requires the excluded
Pilot-1 full ledger and is therefore intentionally outside this package's scope.

## Known boundary

The package implements exact joint OU/equity transitions, common-random-number path reuse,
the deterministic P-to-Q state shift, and the cross-sectional Q forward-variance normalizer.
It does not yet implement conditional VIX valuation, option-surface pricing, or joint SPX/VIX
calibration.
