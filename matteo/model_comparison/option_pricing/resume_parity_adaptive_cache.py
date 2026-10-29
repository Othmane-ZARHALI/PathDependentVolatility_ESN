"""Resume unchanged cached simulation after interruption, preserving provenance.

The original resume assertion compares Python tuples with JSON lists. Validate
its JSON-normalized design instead, then invoke the EXACT original workers.
The original source snapshot and design are not replaced.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time
from .m2_spx import ROOT,frozen_stage3
from .full_readout_fit import dump
from .parity_adaptive_cache import stock_expiry
from .provenance import capture_sources,source_hashes,save_sources,check_sources
from esn_eq.pricing import DatedM2Pricer


def run(out):
    begin=time.perf_counter();d=json.loads((out/'design.json').read_text())
    for name,digest in d['source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,('Original source changed',name)
    a,_,premium,_=frozen_stage3();p=DatedM2Pricer(a,premium);initial,_=p.initial_state([1,0,0])
    for key,value in dict(architecture=asdict(a),premium=asdict(premium),initial_state=initial.tolist()).items():
        assert json.loads(json.dumps(value))==d[key],('Actual model changed',key)
    assert d['paths']==4096 and d['seed']==2026102409 and d['stock_steps_per_trading_day']==512
    folder=out/'resume_001';folder.mkdir(exist_ok=True)
    captured=capture_sources([Path(__file__),out/'design.json']);save_sources(folder,captured,ROOT)
    record=dict(reason='Interrupted task; original resume assertion compares Python tuples to JSON lists. JSON-normalized settings and all original source hashes verify.',
                source_sha256=source_hashes(captured,ROOT),original_design_preserved=True)
    dump(folder/'design.json',record)
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(stock_expiry,str(out),days,d['paths'],d['seed']) for days in [1030,666,302,211,119,57,29]]
        results=[f.result() for f in futures]
    check_sources(captured);dump(out/'run.json',d|dict(status='complete reusable fresh reservoir',resumed=True,resume_seconds=time.perf_counter()-begin,expiries=results,resume_record=record))
    print('RESUMED CACHE COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',type=Path,required=True);args=a.parse_args();run(args.cache.resolve())
