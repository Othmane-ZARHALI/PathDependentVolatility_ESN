"""Independent scores and pooled fresh sampling checks for state/curve fits."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .m2_spx import ROOT
from .audit_joint_curve_fit import score,paired_loss


def audit(run):
    record=json.loads((run/'run.json').read_text())
    for name,digest in record['implementation_sha256'].items():
        assert hashlib.sha256((run/'source'/name).read_bytes()).hexdigest()==digest,name
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    state_run='selected' in record and 'states' in record
    if state_run:
        levels={'baseline':np.array(record['baseline_levels']),'selected':np.array(record['selected']['levels'])}
        penalties={'baseline':0.,'selected':.1*record['selected']['state']['prior_squared_norm']/3}
        assert record['selected']['fresh_checking_not_used']
        fresh=[s['seed'] for s in record['seed_records'] if not s['reused_old_search']]
        assert len(fresh)==len(set(fresh))
    else:
        levels={'baseline':np.array(record['embedded_baseline_levels']),'selected':np.array(record['selected_levels'])}
        penalties={k:.1*record['state']['prior_squared_norm']/3 for k in levels}
        assert json.loads((run/'selected.json').read_text())['checking_not_used']
    spx=pd.read_csv(run/'spx_membership.csv');vix=pd.read_csv(run/'vix_contracts.csv')
    allq=pd.concat([spx,vix],ignore_index=True);mask=np.r_[spx.priority_selected.to_numpy(bool),np.ones(len(vix),bool)]
    q=allq[mask].reset_index(drop=True);f=pd.read_csv(run/'forward_targets.csv')
    assert len(q)==882 and mask[:len(spx)].sum()==707 and len(vix)==175
    archive=np.load(run/'comparison_blocks.npz');paths=np.load(run/'comparison_paths.npz')
    scores=[]
    for row in record['scores']:
        role,trial,name=row['role'],row['trial'],row['candidate']
        prices=pd.read_csv(run/f'prices_{role}_t{trial}_{name}.csv');forwards=pd.read_csv(run/f'forwards_{role}_t{trial}_{name}.csv')
        pd.testing.assert_frame_equal(prices[['instrument','expiry','strike','option_type']],q[['instrument','expiry','strike','option_type']])
        total,components,res=score(q,f,prices.model_price_U.to_numpy(),forwards.model_forward.to_numpy(),levels[name])
        total+=penalties[name]
        np.testing.assert_allclose(total,row['total'],rtol=1e-12,atol=1e-9)
        components['state_penalty']=penalties[name]
        for k,value in components.items():np.testing.assert_allclose(value,row['components'][k],rtol=1e-12,atol=1e-9)
        b=archive[f'{role}_t{trial}_{name}_prices'];bf=archive[f'{role}_t{trial}_{name}_forwards']
        np.testing.assert_allclose(b[:,mask].mean(axis=0),prices.model_price_U,rtol=1e-12,atol=1e-10)
        np.testing.assert_allclose(bf.mean(axis=0),forwards.model_forward,rtol=1e-12,atol=1e-10)
        scores.append(dict(role=role,trial=trial,candidate=name,total=float(total),**components))
    pd.DataFrame(scores).to_csv(run/'score_replay.csv',index=False)
    pooled={};summary=[];expiry_rows=[]
    for name in ['baseline','selected']:
        b=np.concatenate([archive[f'fresh_t{j}_{name}_prices'] for j in range(3)])
        bf=np.concatenate([archive[f'fresh_t{j}_{name}_forwards'] for j in range(3)])
        means=b.mean(axis=0);se=b.std(axis=0,ddof=1)/np.sqrt(len(b))
        for days,group in spx.groupby('expiry_days'):
            individual=np.concatenate([paths[f'fresh_t{j}_{name}_spx_{days}'] for j in range(3)])
            np.testing.assert_allclose(individual.mean(axis=0),means[group.index],rtol=1e-12,atol=1e-10)
            se[group.index]=individual.std(axis=0,ddof=1)/np.sqrt(len(individual))
        fm=bf.mean(axis=0);fs=bf.std(axis=0,ddof=1)/np.sqrt(len(bf))
        result=q.copy();result['model_price_U']=means[mask];result['mc_se_U']=se[mask]
        total,comp,res=score(q,f,means[mask],fm,levels[name]);total+=penalties[name];comp['state_penalty']=penalties[name]
        result['interval_residual']=res;result['negative_estimate']=means[mask]<0
        result['nominal_se_screen']=(means[mask]>0)&(se[mask]<=q.price_scale_U.to_numpy()/8)
        result.to_csv(run/f'prices_fresh_pooled_{name}.csv',index=False)
        ft=f.copy();ft['model_forward']=fm;ft['mc_se']=fs;ft['error_points']=fm-ft.forward_vix_points
        ft.to_csv(run/f'forwards_fresh_pooled_{name}.csv',index=False)
        full=allq.copy();full['model_price_U']=means;full['mc_se_U']=se;full['priority_objective_member']=mask
        full.to_csv(run/f'full_training_fresh_pooled_{name}.csv',index=False)
        summary.append(dict(candidate=name,total=float(total),components=comp,inside_bidask=int((res==0).sum()),
                            negative_primary=int((means[mask]<0).sum()),nominal_se_passes=int(result.nominal_se_screen.sum()),
                            median_se_over_scale=float(np.median(se[mask]/q.price_scale_U))))
        for (instrument,days),group in result.groupby(['instrument','expiry_days']):
            expiry_rows.append(dict(candidate=name,instrument=instrument,days=int(days),quotes=len(group),
                                   mse=float(np.mean(group.interval_residual**2)),inside=int((group.interval_residual==0).sum()),
                                   negatives=int(group.negative_estimate.sum()),nominal_se_passes=int(group.nominal_se_screen.sum())))
        pooled[name]=(b,bf)
    a,fa=pooled['baseline'];b,fb=pooled['selected']
    comparison=paired_loss(q,f,a[:,mask],b[:,mask],fa,fb,levels['baseline'],levels['selected'])
    delta=penalties['selected']-penalties['baseline'];comparison['selected_minus_anchor']+=delta
    comparison['leave_one_block_range']=[v+delta for v in comparison['leave_one_block_range']]
    comparison['baseline_label']='baseline, not original VIX-only anchor'
    comparison['selected_minus_baseline']=comparison.pop('selected_minus_anchor')
    pd.DataFrame(expiry_rows).to_csv(run/'pooled_expiry_scores.csv',index=False)
    report=dict(status='scores, membership, archived mean estimates and source hashes replayed',pooled=summary,paired=comparison,
                uncertainty='SPX pooled iid path contributions; VIX pooled independent complete scrambles. Frozen-pilot conditional SE only. Separate normalizer interpolation, GH/time and stock-grid bias unresolved.',
                scores_replayed=len(scores),market_quotes_per_score=882,physical_fit_or_ranking=False)
    (run/'score_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True)
    args=parser.parse_args();audit(args.run)
