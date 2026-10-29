"""Reconstruct matched-state refinements from saved node expectations and prices."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT
from esn_eq.market_targets import MarketTargets


def run(folder):
    record=json.loads((folder/'run.json').read_text())
    for relative,digest in record['implementation_sha256'].items():
        assert hashlib.sha256((folder/'source'/relative).read_bytes()).hexdigest()==digest,relative
    assert record['scope']['vix_fit_expiry_days']==[27,55,118,153] and not record['scope']['vix_deferred_expiry_days']
    a=np.load(folder/'nodes_and_rules.npz');old=Path(record['outer_source'])
    raw=pd.read_csv(old/'outer_values.csv');raw=raw[(raw.days==118)&(raw.method=='targeted_mean')].sort_values(['replicate','point'])
    np.testing.assert_array_equal(a['states'],np.load(old/'states.npz')['targeted_mean_118'])
    np.testing.assert_allclose(a['weight'],raw.weight,rtol=1e-12,atol=1e-14)
    np.testing.assert_allclose(a['GH32_M48_time8'],raw.vix,rtol=1e-12,atol=1e-14)
    floor=record['variance_floor'];delta=record['delta']
    for name,norm in [('GH48_M48_time8','normalizer48'),('GH48_M64_time8','normalizer64')]:
        integrated=(floor+(a['xi']-floor)*a['GH48']/a[norm])@a['time_weights']
        v=100*np.sqrt(integrated/delta)
        np.testing.assert_allclose(a[name],v,rtol=1e-13,atol=1e-13)
    target=MarketTargets.load(ROOT/'data/processed/options/2026-02-19');assert target.artifact_sha256==record['target_sha256']
    q=target.select('VIX',expiry_days=[118]);sign=np.where(q.option_type.to_numpy()=='call',1.,-1.)
    R,N=record['replications'],record['points_per_replication'];w=a['weight'];summary=pd.read_csv(folder/'summary.csv');replications=pd.read_csv(folder/'replications.csv')
    def payoffs(v):return np.column_stack((v,v*v,np.maximum(sign*(v[:,None]-q.strike.to_numpy()),0.)))
    baseline=payoffs(a['GH32_M48_time8'])
    for name in ('GH32_M48_time8','GH48_M48_time8','GH48_M64_time8'):
        h=payoffs(a[name]);values=(w[:,None]*h).reshape(R,N,-1).mean(axis=1);changes=(w[:,None]*(h-baseline)).reshape(R,N,-1).mean(axis=1)
        rows=summary[summary.rule==name]
        for expected,column in ((values.mean(axis=0),'estimate'),(values.std(axis=0,ddof=1)/np.sqrt(R),'replication_se'),(changes.mean(axis=0),'change_from_base'),(changes.std(axis=0,ddof=1)/np.sqrt(R),'paired_change_se')):
            np.testing.assert_allclose(rows[column],expected,rtol=2e-12,atol=1e-12)
        for rep in range(R):
            rows2=replications[(replications.rule==name)&(replications.replicate==rep)]
            np.testing.assert_allclose(rows2.estimate,values[rep],rtol=2e-12,atol=1e-12)
            np.testing.assert_allclose(rows2.change,changes[rep],rtol=2e-12,atol=1e-12)
    for check in record['checks']:
        label=check['check'];ix=a[label+'_indices'];v=a[label+'_vix'];ref=a['GH48_M64_time8'][ix]
        d=(w[ix]*(v*v-ref*ref)).reshape(R,check['points_per_replication']).mean(axis=1)
        np.testing.assert_allclose([check['weighted_moment_change'],check['paired_change_se'],check['max_pointwise_vix_change']],[d.mean(),d.std(ddof=1)/np.sqrt(R),np.max(abs(v-ref))],rtol=2e-12,atol=1e-12)
    result=dict(status='saved-state node integration and statistics replay passed',source_snapshot_hashes=len(record['implementation_sha256']),full_states_matched=True,original_weights_matched=True,node_integrals_reconstructed=True,summary_rows_recomputed=len(summary),subset_checks_recomputed=len(record['checks']),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limitation='Does not recompute Gaussian node expectations or certify quadrature or unseen-tail error.')
    (folder/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path);run(parser.parse_args().folder)
