"""Training-selected broad-spectrum diagnostic, with preset selection rule."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from .continuum_response_study import record,metadata
from .continuum_response_core import ContinuumFit
from .vix_only_curve_fit import PRIOR,ARCHIVE,DAYS,load_draws,paired
from .full_readout_fit import dump
from .forward_floor_core import unpack


def run(run):
    assert (run/'secondary_design.json').exists()
    ledger=json.loads((run/'optimizer_ledger.json').read_text());eligible=[z for z in ledger if z['label'].startswith('continuous') and z['eligible'] and z['model']['kappa']>=.1]
    if not eligible:
        dump(run/'secondary_selection.json',dict(status='no eligible broad-spectrum candidate'));return
    best=min(eligible,key=lambda z:z['option_loss']);x=np.array(best['parameters']);dump(run/'secondary_selection.json',dict(rule=json.loads((run/'secondary_design.json').read_text()),selected=best,selection_used_training_only=True))
    prior=json.loads((PRIOR/'run.json').read_text());p=np.array(prior['selected']);ends=np.array(prior['ends_years']);draws=load_draws(run/'fresh_cache');nodes={'u_grid':np.linspace(0,1,33)}
    for day in DAYS:
        with np.load(run/'fresh_power_tables'/f'd{day}.npz') as a:
            for key,v in a.items():nodes[f'd{day}_{key}']=v
    levels,floor,_=unpack(x,p);pp=p.copy();pp[7:-1]=levels;pp[-1]=floor;dump(run/'experimental_model_interior.json',dict(parameters=pp.tolist(),ends_years=ends.tolist(),response='normalized geometric M2/stress continuum',measure='Beta',readout=metadata(x,p),architecture=prior['architecture'],premium=prior['premium'],initial_state=prior['initial_state'],adopted=False))
    market=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/vix_contracts.csv');forwards=pd.read_csv(ARCHIVE/'m2_parity_adaptive_cache_001/forward_targets.csv');controls=json.loads((PRIOR/'vix_controls.json').read_text());e=ContinuumFit(nodes,draws,market,forwards,controls,p,ends);s=record(e,x,run/'fresh_check','interior');q=ContinuumFit(nodes,draws,market,forwards,controls,p,ends,measure_order=48);r=q.predict(x);s['measure24_to_48_max_price_change']=float(max(abs(r['prices']-pd.read_csv(run/'fresh_check/prices_otm_interior.csv').model_price_U)));s['paired_losses']=[paired(e,run/'fresh_check','interior',ref) for ref in ['boundary','continuous','previous']];dump(run/'fresh_check/interior_comparison.json',s);print(json.dumps(s),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);run(ap.parse_args().run.resolve())
