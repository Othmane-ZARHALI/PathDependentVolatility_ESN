"""Verify training/checking byte identities against previously frozen manifests."""
import argparse
import hashlib
import json
from pathlib import Path
from .full_readout_fit import HERE,dump


def verify(run):
    r=json.loads((run/'run.json').read_text())
    training=json.loads((HERE/'runs/m2_parity_uniform_001/immutable_inputs_before.json').read_text())
    fresh=json.loads((Path(r['check_cache'])/'immutable_inputs_before_prices.json').read_text());checking=fresh['sha256']
    for phase,items in [('training',training),('checking',checking)]:
        for path,digest in items.items():
            h=hashlib.sha256()
            with Path(path).open('rb') as f:
                for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
            assert h.hexdigest()==digest,('Immutable input changed',path)
        print('INPUTS VERIFIED',phase,len(items),flush=True)
    record=dict(status='all previously frozen inputs unchanged',training_arrays=len(training),checking_arrays=len(checking),
                checking_manifest_before_selection=bool((Path(r['check_cache'])/'immutable_inputs_before_prices.json').stat().st_mtime_ns<(run/'selected.json').stat().st_mtime_ns))
    dump(run/'immutable_inputs_verified.json',record)
    model=json.loads((run/'selected_model.json').read_text());p=r['selected'];names=r['parameter_names']
    for n,x in zip(names[:7],p[:7]):assert model['readout'][n]==x
    assert model['forward_variance']['levels']==p[7:10] and model['variance_floor']==p[10]
    assert model['initial_state']==r['initial_state']
    print('MODEL CONFIG VERIFIED',flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--run',type=Path,required=True);args=a.parse_args();verify(args.run.resolve())
