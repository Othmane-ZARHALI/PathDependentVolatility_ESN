"""Paired selected-state GH32/48 VIX probes, never used in calibration."""
import argparse
import hashlib
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
from .m2_spx import frozen_stage3
from .joint_state_core import StateDatedPricer,AdaptiveDatedNormalizer
from .full_readout_fit import NAMES,dump
from esn_eq.pricing import ForwardVarianceCurve
from esn_eq.vix_pricing import ConditionalVIXPricer


def check(run):
 r=json.loads((run/'run.json').read_text());cache=Path(r['cache']);x=np.array(r['selected_parameters']);a,p,premium,_=frozen_stage3();p=replace(p,**dict(zip(NAMES[:7],x[:7])));dated=StateDatedPricer(a,premium);initial=np.array(r['initial_state']);curve=ForwardVarianceCurve(tuple(r['ends_years']),tuple(x[7:10]));horizon=(153+30)/365
 norm=AdaptiveDatedNormalizer(dated,p,initial,horizon,order=48,tolerance=2e-7);dated.state_normalizer=norm;engine=ConditionalVIXPricer(dated,p,curve,initial,normalizer_order=48,variance_floor=x[10]);arrays=np.load(run/'paths_check_direct_selected.npz');rows=[]
 for days in [27,55,118,153]:
  source=cache/'check'/f'vix_d{days}.npz';states=np.load(source)['states'];v=arrays[f'vix_{days}'];ix=np.argsort(v)[[0,len(v)//2,-1]];window=engine.gaussian_window(days/365,time_order=8);low=engine.evaluate(states[ix],window,order=32).vix;high=engine.evaluate(states[ix],window,order=48).vix
  rows.append(dict(days=days,indices=ix.tolist(),saved_direct_GH32=v[ix].tolist(),paired_GH32=low.tolist(),paired_GH48=high.tolist(),GH48_minus_GH32=(high-low).tolist(),maximum_saved_vs_paired_GH32=float(max(abs(low-v[ix]))),state_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
 result=dict(status='completed selected-state paired quadrature probes, not surface convergence',rule='Low/median/high saved raw VIX states at each expiry; no parameter selection uses these probes.',normalizer=norm.audit,rows=rows,maximum_absolute_GH48_minus_GH32=max(max(abs(np.array(t['GH48_minus_GH32']))) for t in rows),limitations='Twelve observed states, not a uniform bound or coverage of unseen rare states. Both orders share the same GH48 normalizer and time rule. No original prices or parameters changed.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(run/'validation').mkdir(exist_ok=True);dump(run/'validation/selected_quadrature.json',result);return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();print(json.dumps(check(a.run),indent=2))
