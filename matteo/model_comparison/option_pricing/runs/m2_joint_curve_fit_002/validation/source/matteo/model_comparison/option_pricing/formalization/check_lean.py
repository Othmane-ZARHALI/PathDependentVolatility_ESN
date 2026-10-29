"""Read-only reuse of a pinned installed Mathlib environment; audit local proof.

No lake invocation, dependency update, or writes in the dependency checkout.
Each attempt has a new audit directory, including compiler failures.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--mathlib-project', type=Path, default=Path('/Users/matteogambara/Desktop/code/LLM handoff/HandoffResearch/Lean env'))
args = parser.parse_args()
toolchain = (HERE/'lean-toolchain').read_text().strip()
if toolchain != (args.mathlib_project/'lean-toolchain').read_text().strip():
    raise ValueError('Pinned toolchain mismatch')
manifest = json.loads((args.mathlib_project/'lake-manifest.json').read_text())
package = next(p for p in manifest['packages'] if p['name'] == 'mathlib')
expected = '1f414401f69059aa7eead47b53ee40bd38455eeb'
if package['rev'] != expected:
    raise ValueError('Pinned Mathlib mismatch')
mathlib = args.mathlib_project/'.lake/packages/mathlib'
revision = subprocess.run(['git', '-C', str(mathlib), 'rev-parse', 'HEAD'], text=True, capture_output=True, check=True).stdout.strip()
if revision != expected:
    raise ValueError('Actual Mathlib checkout revision mismatch')
dirty = subprocess.run(['git', '-C', str(mathlib), 'status', '--porcelain', '--untracked-files=no'],
    text=True, capture_output=True, check=True).stdout
if dirty.strip():
    raise ValueError('Mathlib tracked sources differ from pinned revision')
lean = Path.home()/'.elan/toolchains'/toolchain.replace('/', '--').replace(':', '---')/'bin/lean'
version = subprocess.run([str(lean), '--version'], text=True, capture_output=True, check=True).stdout.strip()
paths = sorted((args.mathlib_project/'.lake/packages').glob('*/.lake/build/lib/lean'))
env = dict(os.environ); env['LEAN_PATH'] = ':'.join(str(p) for p in paths)
source = HERE/'CalibrationChecks.lean'; raw = source.read_text()
if re.search(r'\b(sorry|admit|axiom|unsafe)\b', re.sub(r'/\-.*?\-/', '', raw, flags=re.S)):
    raise ValueError('Admitted, custom-axiom or unsafe declaration in proof source')
base = HERE/'audits'; base.mkdir(exist_ok=True)
index = 1
while (base/f'attempt_{index:03d}').exists(): index += 1
output = base/f'attempt_{index:03d}'; output.mkdir()
(output/'CalibrationChecks.lean').write_text(raw)
command = [str(lean), '-DwarningAsError=true', '-o', str(output/'CalibrationChecks.olean'), str(source)]
completed = subprocess.run(command, env=env, text=True, capture_output=True, timeout=180)
log = completed.stdout+completed.stderr
(output/'compiler_and_axioms.log').write_text(log)
theorems = ['gaussian_exponential_mean_one', 'gaussian_stock_step_mean', 'gaussian_hedge_step_mean_zero',
    'frozen_anchor_control_reuse', 'finite_importance_identity', 'integrated_variance_affine',
    'integrated_variance_positive', 'weighted_residual_norm_identity', 'penalty_residual_identity',
    'finite_noise_square_decomposition']
axioms = {name: [a.strip() for a in values.split(',') if a.strip()]
    for name, values in re.findall(r"'M2Calibration\.([^']+)' depends on axioms: \[([^]]*)\]", log)}
axiom_audit_passes = set(axioms) == set(theorems) and all(
    set(values) <= {'propext', 'Classical.choice', 'Quot.sound'} for values in axioms.values())
audit = dict(command=command, toolchain=toolchain, lean_version=version, mathlib_revision=revision,
    dependency_project=str(args.mathlib_project), lean_path=env['LEAN_PATH'], exit_code=completed.returncode,
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    sorryAx_present='sorryAx' in log, theorem_axioms=axioms, expected_theorems=theorems,
    axiom_audit_passes=axiom_audit_passes,
    status='checked' if completed.returncode == 0 and 'sorryAx' not in log and axiom_audit_passes else 'failed')
(output/'audit.json').write_text(json.dumps(audit, indent=2)+'\n')
print(log); print(json.dumps(audit, indent=2))
raise SystemExit(0 if audit['status'] == 'checked' else 1)
