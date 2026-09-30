# TODO

- [ ] **Tests for option processing.** Extend the current synthetic and miniature-snapshot tests to cover malformed schemas, quote reason codes, strike splits, parity and convexity conflicts, training-only curve provenance, variance-floor enforcement, and deterministic output. Keep the supplied snapshot as a reconciliation check rather than a hard-coded golden fixture. See [`test_options.py`](../preprocessing/tests/test_options.py) and the [option preprocessing plan](documentation/OPTION_PREPROCESSING_PLAN.md).

- [ ] **Understand the cross-market checks better.** Review the VIX second-moment identity and finite-strike option lower bound, including the assumed model mean, omitted wings, settlement basis, and forward-source uncertainty. Investigate the 118- and 153-day conflicts in `data/processed/options/2026-02-19/cross_market_checks.csv`; profile a training-only nuisance curve before drawing a conclusion about joint feasibility.

- [ ] **Calibration plan: from theory to code.** Turn the [calibration plan](documentation/CALIBRATION_PLAN.md) into explicit module interfaces and an implementation sequence for the M2 pricing engines, dated latent state, objective and constraints, numerical accuracy gates, and train/validation/final-test workflow.
