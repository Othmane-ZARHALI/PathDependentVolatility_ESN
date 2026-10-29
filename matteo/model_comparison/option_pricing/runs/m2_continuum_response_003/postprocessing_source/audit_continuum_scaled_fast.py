"""Independent power recurrence and scale-aware VIX replay audit."""
import argparse,json,hashlib
from pathlib import Path
from . import audit_continuum_response_v2 as audit
from .audit_continuum_response_fast import tables
from .full_readout_fit import dump

def run(folder):
    audit.tables=tables;audit.run(folder);p=folder/'audit/audit.json';d=json.loads(p.read_text());d['replay_tolerance']='absolute1e-9 OR relative1e-11; report both and weighted mean error, including extremely large tiny-weight proposal VIX values';d['backend']='NumPy exact finite-GH power recurrence, with scalar crosschecks';d['postprocessing_source_hashes']={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() for name in ['audit_continuum_scaled_fast.py','audit_continuum_response_v2.py','audit_continuum_response_fast.py']};dump(p,d)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);run(ap.parse_args().run.resolve())
