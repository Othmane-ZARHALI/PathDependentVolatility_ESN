"""Hash completed independent checking inputs without evaluating model prices."""
import argparse
import hashlib
import json
from pathlib import Path
import time


def run(cache,controls,output):
    begin=time.perf_counter()
    while not (cache/'run.json').exists() or not (controls/'run.json').exists():
        if time.perf_counter()-begin>14400:raise TimeoutError('Waiting for completed checking inputs')
        time.sleep(3)
    hashes={}
    for root in [cache,controls]:
        for p in sorted(root.rglob('*.np*')):
            h=hashlib.sha256()
            with p.open('rb') as f:
                for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
            hashes[str(p)]=h.hexdigest()
    output.write_text(json.dumps(dict(status='inputs frozen before checking-price evaluation',sha256=hashes,arrays=len(hashes)),indent=2)+'\n')
    print('FROZEN CHECK INPUTS',len(hashes),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--controls',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.cache.resolve(),a.controls.resolve(),a.output.resolve())
