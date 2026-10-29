"""Separate the new multiplier from its VIX-driven variance-curve refit."""
import json,time
import numpy as np
import pandas as pd
from .asymmetric_readout import OUT,PARENT,PRIOR,LEVERAGE,Kernel,evaluate_stocks,VIX,arr,dump
from .vix_only_curve_fit import DAYS,ACTIVE

def main():
    begin=time.perf_counter();config=json.loads((OUT/'candidates.json').read_text());x=arr(config['models']);models=arr([x[0],x[0],x[5],x[5]]);alpha=arr([0,1.5,0,1.5]);d=dict(np.load(PARENT/'fresh_check/grid_and_transitions.npz'));p=arr(json.loads((PRIOR/'run.json').read_text())['selected']);z=np.load(OUT/'normalizers.npz')['Z'][:,4];Z=arr(np.column_stack([np.ones(len(z)),z,np.ones(len(z)),z]));k=Kernel(OUT/'asymmetric.dylib');scale=json.loads((OUT/'design.json').read_text())['scale'];old=np.load(OUT/'fresh_stocks.npz');pieces=[]
    from .leverage_scan import combine
    for seed in old['seeds']:pieces.append(k.stocks(d,64,p,models,alpha,Z,scale,int(seed)))
    stocks=combine(pieces);error=float(np.max(abs(stocks[:,:,3]-old['stocks'][:,:,5])));assert error<1e-12;folder=OUT/'ablation';folder.mkdir(exist_ok=False);np.savez_compressed(folder/'stocks.npz',stocks=stocks,seeds=old['seeds']);ref=np.load(LEVERAGE/'reference.npz')['stats'];q=pd.read_csv(OUT/'SPX_membership.csv');scores=evaluate_stocks(stocks,ref,models,alpha,q,folder/'SPX');draws={day:dict(np.load(OUT/f'fresh_VIX_draws/d{day}.npz')) for day in DAYS};H={day:np.load(OUT/f'fresh_VIX_exposures/d{day}.npz')['H'][[0,4]] for day in DAYS};qv=pd.read_csv(OUT/'VIX_membership.csv');f=pd.read_csv(OUT/'forward_targets.csv');v=VIX(H,draws,x[0],qv,f);rows=[];names=['original','multiplier_only','curve_only','multiplier_and_curve']
    for c,name in enumerate(names):
        row=next(s for s in scores if s['candidate']==c and s['steps']==128).copy();r=v.predict(1 if c in [1,3] else 0,models[c,ACTIVE]);row.update(name=name,VIX_loss=r['loss'],maximum_VIX_forward_gap=float(max(abs(r['forwards']-f.forward_vix_points.to_numpy()))));frame=pd.read_csv(folder/f'SPX/SPX_prices_c{c}_128.csv');a=frame[(frame.expiry_days==29)&frame.priority_selected];lo=a.loc[(a.K_over_F-.9).abs().idxmin()];hi=a.loc[(a.K_over_F-1.1).abs().idxmin()];row['day29_right_minus_left_IV_pp']=float(hi.model_IV-lo.model_IV);rows.append(row)
    dump(folder/'results.json',dict(rows=rows,full_candidate_replay_max_error=error,wall_seconds=time.perf_counter()-begin,scope='paired post-selection mechanism ablation, not new independent validation; all Gaussian draws match the fresh stock comparison'));pd.DataFrame(rows).to_csv(folder/'results.csv',index=False);print('ABLATION',rows,flush=True)

if __name__=='__main__':main()
