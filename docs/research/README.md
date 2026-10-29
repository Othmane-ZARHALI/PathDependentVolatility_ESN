# ESN and M2 research record

Consolidated on 29 October 2026. Preserve the experiments, what they taught us, and the pictures that show the remaining problems. The strongest partial VIX result and the improvements in SPX dependence have not produced a satisfactory joint calibration.

Read only what the next task needs:

- [Current state and continuation brief](CURRENT_STATE.md): model, scope, constraints, completed tests and open decisions.
- [Experiment history and learnings](EXPERIMENT_HISTORY.md): the research sequence, including unsuccessful experiments and evidence limits.
- [Visual evidence](VISUAL_EVIDENCE.md): thirteen preserved figures with explanations and source links.
- [Decisions and questions](DECISIONS.md): what remains in force, which ideas have been tested, and what would justify revisiting them.
- [Complete report catalogue](REPORT_CATALOG.md): individual reports, plans, physical-measure experiments and formal correspondence.
- [Evidence manifest](EVIDENCE_MANIFEST.json): source hashes, available run folders and missing referenced run folders.
- [Consolidation checks](VALIDATION.md): preservation and navigation checks, plus existing evidence gaps.
- [Checkpoint and artifact policy](PUBLICATION.md): publication checks, local array retention and the approved cleanup boundary.
- [Cleanup receipt](CLEANUP_RECEIPT.json): temporary and compiled files removed after the verified push; simulation arrays retained.

The existing [ten-page visual retrospective](../../output/pdf/concise/spx_vix_visual_summary.pdf) remains useful for the earlier mechanism comparisons. Its source, plot inputs and provenance remain in [the original report folder](../../output/pdf/concise). The new history and gallery add the 29 October filtering and protected-variance experiments.

## Keeping future context small

Start a research session with CURRENT_STATE.md, then load the relevant experiment and original report. Read the full catalogue or historical ledger when investigating provenance or an older approach. Keep fitted parameter tables and per-contract diagnostics in their original evidence folders.

The former continuation prompt and next-steps file are preserved verbatim in [the dated archive](archive/2026-10-29). Their active paths now point to this record. The [historical ledger](../../matteo/model_comparison/option_pricing/RESEARCH_LEDGER.md) retains its accumulated content; dated proposals there may already have been tested by later entries.

## Preservation boundary

This consolidation changes research navigation and continuation instructions. It adopts no model, changes no pricing scope, and launches no calibration. Original result reports, source code, run artifacts and figures remain in place. Selected figures are copied without alteration and have [source hashes](FIGURE_PROVENANCE.json). The evidence manifest verifies the original reports directly, using the archived copies of the two replaced continuation files. Duplicate ZIP snapshots, raw test transcripts and the internal Git-status dump remain local and are excluded from the branch.

The `option_pricing/formalization/` and `option_pricing/literature/` folders are also kept locally and excluded from the branch at the user's request. Historical reports may refer to files in those folders; those references require the local copies. The research history preserves the conclusions and their evidence limits.
