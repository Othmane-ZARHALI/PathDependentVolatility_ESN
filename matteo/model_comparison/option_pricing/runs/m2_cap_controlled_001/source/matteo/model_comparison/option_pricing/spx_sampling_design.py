"""Pilot-only component diagnosis and fixed call-tilt hypothesis; no check data."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .m2_spx import ROOT
from .spx_volatility_sampling import proposal
from .spx_path_score_controls import ridge_fit
from esn_eq.importance_sampling import DefensiveGaussianMixture


COMPONENTS=['unshifted','equity_plus','equity_minus','cluster_plus','cluster_minus','feedback_plus','feedback_minus','spike_plus','spike_minus']


def call_proposal(positive_shift=3.):
    if not np.isfinite(positive_shift) or positive_shift<=0:raise ValueError('Positive finite shift required')
    original=proposal(1);shifts=original.shifts.copy();shifts[1,0]=positive_shift
    return DefensiveGaussianMixture(shifts,original.probabilities)


def tilted_second_moment_ratio(y,b,beta,old_weights,z,new_proposal):
    # The mass column is w-1, so residual = w*h + beta_mass, not w*h.
    residual=y-np.einsum('njk,jk->nj',b,beta)-beta[:,-1]
    ratio=new_proposal.weights(z)/old_weights  # q_old/q_new, not its inverse
    return np.mean(residual**2*ratio[:,None],axis=0)/np.mean(residual**2,axis=0)


def diagnose(parent,output):
    output.mkdir(parents=True,exist_ok=False);components=[];screens=[];hashes={}
    for folder in sorted(parent.glob('expiry_*'),key=lambda p:int(p.name.split('_')[-1])):
        for name in ['pilot_draws.npz','frozen.json','contracts.csv']:
            path=folder/name;hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        quotes=pd.read_csv(folder/'contracts.csv');days=int(quotes.expiry_days.iloc[0]);fit=json.loads((folder/'frozen.json').read_text())
        beta=np.asarray(fit['coefficients']['base'])
        with np.load(folder/'pilot_draws.npz') as data:
            y,b,w,z=data['g3_y'],data['g3_base'],data['w'],data['z'];labels=data['components']
            h=y-np.einsum('njk,jk->nj',b,beta)-beta[:,-1];half=len(y)//2
            for j,row in enumerate(quotes.itertuples()):
                total=float(np.sum(h[:,j]**2))
                for k,label in enumerate(COMPONENTS):
                    mask=labels==k
                    components.append(dict(expiry_days=days,contract=j,strike=row.strike,option_type=row.option_type,band=row.band,
                        component=label,paths=int(mask.sum()),payoff_hits=int(((y[:,j]>0)&mask).sum()),
                        weighted_payoff_sum=float(y[mask,j].sum()),
                        weighted_h_second_moment_share=float(np.sum(h[mask,j]**2)/total) if total else None))
            for train,test,fold in [(np.arange(half),np.arange(half,len(y)),0),(np.arange(half,len(y)),np.arange(half),1)]:
                beta=np.stack([ridge_fit(y[train,j],b[train,j]) for j in range(y.shape[1])])
                h=y[test]-np.einsum('njk,jk->nj',b[test],beta)-beta[:,-1]
                base_squares=(h*quotes.forward.to_numpy()[None,:]/quotes.price_scale_U.to_numpy()[None,:])**2
                for shift in [2.,2.5,3.,3.5,4.]:
                    ratios=tilted_second_moment_ratio(y[test],b[test],beta,w[test],z[test],call_proposal(shift))
                    density_ratio=call_proposal(shift).weights(z[test])/w[test]
                    screens.append(dict(expiry_days=days,fold=fold,positive_shift=shift,
                        quote_scale_objective_ratio=float(np.mean(base_squares*density_ratio[:,None])/np.mean(base_squares)),
                        mean_contract_ratio=float(ratios.mean()),worst_contract_ratio=float(ratios.max()),outer_call_ratio=float(ratios[-1])))
    pd.DataFrame(components).to_csv(output/'pilot_components.csv',index=False);pd.DataFrame(screens).to_csv(output/'pilot_reweighting.csv',index=False)
    record=dict(inputs_sha256=hashes,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        selected_positive_equity_shift=3.,original_shift=2.,
        reason='Moderate targeted outer-call coverage hypothesis; all outer-call pilot hits used upward equity component. Not an aggregate-objective optimum or validated improvement.',
        math='Residual Z_old = w_old*h + beta_mass. For fixed beta, E_qnew[(w_new*h)^2] = E_qold[(Z_old-beta_mass)^2 * w_new/w_old]. Two-fold pilot-only proxies; not fresh variance measurements.',
        limitations=['No historical checking archive read; only previous pilots are design inputs.',
            'One-hit folds cannot identify unseen tails; several candidate aggregate proxies worsen. No robust pilot-wide win is claimed.',
            'Fresh comparison must retain original proposal, equal learning budgets, fixed old-pilot coefficients and paired fresh checking.',
            'Target Q model, cutoff3, all quote/maturity scopes and full-path likelihood remain unchanged.'])
    (output/'design.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print('Pilot-only diagnosis saved; fixed +3 equity-tilt hypothesis, not a validated sampler.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('parent',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();diagnose(args.parent,args.output)
