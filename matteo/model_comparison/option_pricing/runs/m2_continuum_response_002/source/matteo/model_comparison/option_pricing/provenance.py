"""Freeze executable sources before a diagnostic and reject mid-run edits."""
import hashlib
from pathlib import Path


def capture_sources(paths):
    return {Path(p): Path(p).read_bytes() for p in [*paths, Path(__file__)]}


def source_hashes(captured, root):
    return {str(path.relative_to(root)): hashlib.sha256(raw).hexdigest()
            for path, raw in captured.items()}


def check_sources(captured):
    changed = [str(p) for p, raw in captured.items() if p.read_bytes() != raw]
    if changed:
        raise RuntimeError(f"Sources changed during calculation; rerun from frozen code: {changed}")


def save_sources(output, captured, root):
    for path, raw in captured.items():
        target = output/"source"/path.relative_to(root)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
