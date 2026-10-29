# Research evidence before consolidation

Snapshot date: 29 October 2026.

- `CONTINUE_FITTING_PROMPT.md` and `NEXT_STEPS.md` are exact copies of the former active continuation files. Relative links inside these unchanged originals use their original repository context.
- `evidence_snapshot.zip` is a local duplicate snapshot, excluded from Git. The original reports and figures remain in their documented repository locations; EVIDENCE_MANIFEST.json verifies the reports directly without this ZIP.
- `missing_tracked_evidence_HEAD.zip` is a local recovery copy of 1934 tracked run files that were already absent from the working tree. These bytes come from the committed ancestor identified in RECOVERY_MANIFEST.json, not from missing uncommitted versions. It has not been restored over the current checkout. This redundant 99 MB ZIP is excluded from Git: the committed ancestor preserves the same bytes, and the consolidation verifier can check them directly from Git when the ZIP is absent.
- `RECOVERY_MANIFEST.json` identifies the commit and each recovered file. Availability gaps in the current checkout remain recorded in the main evidence manifest. An archived committed version does not establish that it was the last version used by a report.
- `git_status_before.bin` is an internal local Git-status dump, excluded from Git. The committed recovery manifest supplies the original missing-file list used by the preservation checker.

Existing research documents, figures, run outputs and executable sources remain in their original locations where present. The optional local ZIP files do not constitute a complete backup of every untracked file or simulation array.
