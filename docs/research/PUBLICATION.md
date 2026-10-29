# Research checkpoint and local artifact policy

29 October 2026. The checkpoint preserves accumulated code, experiment reports, compact results, input/provenance records, source snapshots, figures, formal sources and the consolidated research record on the existing research branch. The 35 earlier local commits remain in its ancestry. No model adoption or market-fit certification follows from publication.

## What belongs in Git

Research source and tests, small run metadata, fitted candidate JSON files, CSV diagnostics, selected audit logs, figures and final report artifacts are versioned. The 148-document pre-consolidation evidence snapshot remains versioned. Existing tracked deletions are recorded; their committed historical bytes remain accessible through the ancestor and recovery manifest.

New generated NPY/NPZ simulation arrays remain local and ignored, following the user's explicit retention choice. Small fixtures already tracked by earlier commits retain their existing Git status. [LOCAL_ARRAY_INVENTORY.json](LOCAL_ARRAY_INVENTORY.json) records local array paths, sizes, modification times and symlink targets. The approximately 54 GB total counts symlink references; regular array files account for about 47 GB. This inventory is not a remote backup of arrays.

Raw email material, the prepared derivative package for the supplied October PDFs, temporary PDF review renders, compiled Python/C++/Lean caches and the redundant local recovery ZIP are excluded from Git. The original supplied PDFs and scientific result figures are retained. The recovery ZIP's contents can be verified from the recorded ancestor commit without uploading a duplicate 99 MB archive.

## Checks before cleanup

- Consolidation preservation: 148 original documents, 13 copied figures, 215 visual/source assets, 299 local links, and 1934 historical Git files verified. The check passes both with the local recovery ZIP and directly from ancestor Git objects.
- Physical model suite: 124 tests passed.
- PDF preprocessing tool: 11 tests passed.
- Pricing suite: 117 tests ran; 113 passed and four errored on previously absent fixtures. Two require the absent `m2_parity_adaptive_cache_001` stock/contract cache, and two require historical `spx_scaled_moneyness_003/training_membership.csv`. This is not a fully passing pricing suite. Logs are in [checks/2026-10-29](checks/2026-10-29).
- Tests used a temporary Python 3.12 environment with NumPy 2.5.4, SciPy 1.18.1, pandas 3.1.0 and matplotlib 3.11.2; the separate installed PDF tool used PyMuPDF 1.26.7. An initial incompatible Python 3.9 invocation was replaced by the supported environment before these results.
- New consolidation documents, tool files and ignore rules pass the scoped whitespace check. Historical source documents, extracted literature and immutable run snapshots retain their original bytes, including inherited whitespace.
- Staging was checked to exclude local simulation arrays, compiled libraries, temporary renders, raw email material and the duplicate recovery ZIP. This exclusion check applies to new content; already-deleted tracked arrays remain historical deletions.

## Cleanup boundary

After the remote checkpoint is verified, remove temporary PDF review renders, Python bytecode and test caches, and compiled `.dylib`, `.olean` and `.ilean` files. Keep every simulation array and its symlink. The compiled kernels and Lean outputs must be rebuilt before dependent tests or exact replays; the above test outcomes describe the pre-cleanup environment.

Keep final PDFs, report source, figures, price tables, fitted parameters, run source snapshots, audit records, original literature and the local recovery ZIP. Git object storage and unrelated local materials are outside this cleanup. A separate cleanup receipt records the files actually removed after publication.

## Cleanup completion

The research checkpoint `2aa44db` was pushed and verified on the approved public research branch before cleanup. The [cleanup receipt](CLEANUP_RECEIPT.json) records 595 removed temporary or compiled files, totalling 20,945,437 bytes. No tracked file or simulation array was deleted. All 2379 inventoried array files and symlink references retained their paths, sizes, modification times and symlink targets. Research evidence preservation was checked again after cleanup.
