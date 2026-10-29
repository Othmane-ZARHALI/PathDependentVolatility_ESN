# Consolidation preservation checks

Checked on 29 October 2026. These are documentation and file-integrity checks, not new pricing simulations or scientific validation.

- The original Markdown reports, plans and records are checked against their original hashes directly. The two former continuation files have separate exact archived copies. One report in the locally retained formalization folder is unavailable in a fresh checkout and is reported separately by the verifier. The local duplicate documentation ZIP is not required for verification or included in the branch.
- The historical bodies of the pricing README, research ledger and M2 README are unchanged after their new navigation banners.
- Thirteen selected figures are exact copies of their originals. Hash checks also cover 215 original visual and report assets, including existing plot inputs and metadata.
- Current active_scope.json is unchanged. No model code or simulation output was edited by consolidation.
- Local links and heading anchors in the new active documents and replacement continuation files were checked. Machine-readable results are in [validation.json](validation.json).
- The local committed-HEAD recovery archive contains 1934 verified files that were already absent from the working tree. Their [manifest](archive/2026-10-29/RECOVERY_MANIFEST.json) records origin, sizes and hashes. The original deletions remain as they were. The redundant ZIP is excluded from Git; the same historical files remain in the recorded ancestor commit, and the verifier falls back to those Git objects.

## Existing availability gaps

The original reports reference 43 run-directory names absent from the checkout. Eight other referenced folders are present without a top-level run.json; report-specific manifests may use other layouts, so this is an availability observation rather than a claim of corruption. Exact paths are in [EVIDENCE_MANIFEST.json](EVIDENCE_MANIFEST.json).

Committed versions of already-deleted tracked run files were preserved separately. They may differ from the versions used by a later report and cannot recover missing uncommitted outputs. No claim of complete rerun reproducibility follows from this archive. Existing code and per-experiment audits retain their original scope.

The curated research history was cross-checked against current repository reports, active scope, the original visual retrospective source and the latest filtering/protected-correction reports. Memory was used to locate prior topics; historical memory counts were not treated as the current scope. No Lean compilation, calibration or Monte Carlo experiment was needed or rerun for this documentation-only task.

Recheck preservation and navigation from the repository root with `python3 docs/research/verify_consolidation.py`.
