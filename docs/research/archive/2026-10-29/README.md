# Research evidence before consolidation

Snapshot date: 29 October 2026.

- `CONTINUE_FITTING_PROMPT.md` and `NEXT_STEPS.md` are exact copies of the former active continuation files. Relative links inside these unchanged originals use their original repository context.
- `evidence_snapshot.zip` contains 148 original Markdown research documents, the selected original figure files, the existing visual retrospective PDF and source metadata, and active_scope.json. Repository-relative paths are retained. The embedded EVIDENCE_MANIFEST.json records the original document hashes. Extract to a separate folder when inspecting the historical documentation; its instructions describe historical phases.
- `missing_tracked_evidence_HEAD.zip` is a local recovery copy of 1934 tracked run files that were already absent from the working tree. These bytes come from the committed ancestor identified in RECOVERY_MANIFEST.json, not from missing uncommitted versions. It has not been restored over the current checkout. This redundant 99 MB ZIP is excluded from Git: the committed ancestor preserves the same bytes, and the consolidation verifier can check them directly from Git when the ZIP is absent.
- `RECOVERY_MANIFEST.json` identifies the commit and each recovered file. Availability gaps in the current checkout remain recorded in the main evidence manifest. An archived committed version does not establish that it was the last version used by a report.
- `git_status_before.bin` records the pre-existing working-tree status as NUL-separated porcelain output. It is a preservation check, not a research result.

The source snapshot covers research documentation, selected figures and report inputs. Existing run outputs and executable sources remain in their original locations where present. The two ZIP files do not constitute a complete backup of every untracked file or simulation array.
