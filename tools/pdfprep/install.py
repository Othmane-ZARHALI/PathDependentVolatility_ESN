#!/usr/bin/env python3
"""Install the standalone CLI and its self-contained user skill."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, default=Path.home() / ".local/share/pdfprep")
    parser.add_argument("--bin-dir", type=Path, default=Path.home() / ".local/bin")
    parser.add_argument("--skill-root", type=Path, default=Path(os.environ.get("CODEX_HOME", str(Path.home() / ".codex"))) / "skills")
    parser.add_argument("--offline", action="store_true", help="Install using only uv's existing cache")
    args = parser.parse_args()
    uv = shutil.which("uv")
    if not uv:
        parser.error("uv is required for installation (https://docs.astral.sh/uv/).")
    if sys.version_info < (3, 10):
        parser.error("Run this installer with Python 3.10 or later.")
    source = Path(__file__).resolve().parent
    prefix = args.prefix.expanduser().resolve()
    skill = args.skill_root.expanduser().resolve() / "pdf-preprocess"
    command = args.bin_dir.expanduser().resolve() / "pdfprep"
    marker = prefix / "installation.json"
    if prefix.exists() and not marker.is_file():
        parser.error(f"Refusing an unrelated existing installation directory: {prefix}")
    if skill.exists() and not (skill / ".pdfprep-installed").is_file():
        parser.error(f"Refusing an unrelated existing skill: {skill}")
    venv = prefix / "venv"
    expected_command = venv / "bin/pdfprep"
    if command.exists() or command.is_symlink():
        if not command.is_symlink() or command.resolve() != expected_command.resolve():
            parser.error(f"Refusing to replace an existing command: {command}")
    prefix.mkdir(parents=True, exist_ok=True)
    # Mark partial installation too, so an interrupted install can be retried.
    marker.write_text(json.dumps({"tool": "local-pdfprep", "source": str(source)}, indent=2) + "\n")
    skill.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".pdf-preprocess-install-", dir=skill.parent) as scratch:
        bundle = Path(scratch) / "pdf-preprocess"
        (bundle / "scripts").mkdir(parents=True)
        for name in ("SKILL.md", "pyproject.toml"):
            shutil.copyfile(source / name, bundle / name)
        shutil.copyfile(source / "scripts/pdfprep.py", bundle / "scripts/pdfprep.py")
        (bundle / ".pdfprep-installed").write_text("local-pdfprep\n")
        if not (venv / "bin/python").exists():
            subprocess.run([uv, "venv", str(venv), "--python", sys.executable], check=True)
        install = [uv, "pip", "install", "--python", str(venv / "bin/python")]
        if args.offline:
            install.append("--offline")
        subprocess.run([*install, str(bundle)], check=True)
        # Keep only the files we own; an unrelated existing skill was rejected.
        skill.mkdir(exist_ok=True)
        (skill / "scripts").mkdir(exist_ok=True)
        for name in ("SKILL.md", "pyproject.toml", ".pdfprep-installed", "scripts/pdfprep.py"):
            shutil.copyfile(bundle / name, skill / name)
    command.parent.mkdir(parents=True, exist_ok=True)
    if not command.is_symlink():
        command.symlink_to(expected_command)
    subprocess.run([str(command), "--doctor"], check=True)
    print(f"CLI installed: {command}\nSkill installed: {skill / 'SKILL.md'}")
    if str(command.parent) not in os.environ.get("PATH", "").split(os.pathsep):
        print(f"Use the absolute command path or add {command.parent} to PATH.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
