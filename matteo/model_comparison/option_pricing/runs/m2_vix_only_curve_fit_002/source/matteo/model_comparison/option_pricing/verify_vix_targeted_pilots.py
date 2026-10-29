"""Deterministic pilot-artifact replay; no proposal-performance certificate."""
import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path
import tempfile
import numpy as np
from .vix_targeted_pilots import run as reproduce


def run(folder):
    record=json.loads((folder/'run.json').read_text())
    for relative,digest in record['implementation_sha256'].items():
        assert hashlib.sha256((folder/'source'/relative).read_bytes()).hexdigest()==digest,relative
    assert record['proposal_sha256']==hashlib.sha256((folder/'proposals.json').read_bytes()).hexdigest()
    with tempfile.TemporaryDirectory(prefix='vix-pilot-replay-') as scratch:
        check=Path(scratch)/'replay'
        with contextlib.redirect_stdout(io.StringIO()):reproduce(check)
        names=('proposals.json','pilot_diagnostics.json','screening_replications.csv','screening_metrics.csv','screening_families.csv','screening_ranks.csv','run.json')
        for name in names:assert (folder/name).read_bytes()==(check/name).read_bytes(),name
        a=np.load(folder/'design_arrays.npz');b=np.load(check/'design_arrays.npz')
        assert set(a.files)==set(b.files)
        for key in a.files:np.testing.assert_array_equal(a[key],b[key])
    result=dict(status='deterministic saved-artifact replay passed',source_snapshot_hashes=len(record['implementation_sha256']),replayed_files=list(names),array_blocks=len(a.files),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limitation='Regenerates design through the same tested implementation; not an independent numerical integration, fresh validation or performance certificate.')
    (folder/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path);run(parser.parse_args().folder)
