"""Replay smoke check, then frozen1024path128/256 future-leverage experiment."""
import json,time
import numpy as np
import pandas as pd
from .forward_leverage import OUT,TARGET,LEVERAGE,setup,relabel,arr,dump
from .leverage_scan import combine
from .asymmetric_readout import evaluate_stocks

def main():
    begin=time.perf_counter();audit=json.loads((OUT/'integration_audit.json').read_text());assert audit['production_rule_approved'];prior,dated,norm,_,p,x,d,k=setup(OUT);b=dict(np.load(OUT/'basis.npz'));ss=dict(np.load(OUT/'signals.npz'));old=np.load(TARGET/'grid_refinement/stocks.npz');tick=time.perf_counter();probe=k.stocks(d,b,ss,p,x,2,128,int(old['seeds'][0]));idx=np.array([0,1,1024,1025]);expected=old['stocks'][idx];error=float(np.max(abs(probe[:,:,:2]-expected[:,:,:2])));aux=float(np.max(abs(probe[:,:,-1]-expected[:,:,-1])));assert error<1e-11 and aux<1e-11;dump(OUT/'baseline_replay.json',dict(maximum_baseline_and_instant_stock_error=error,maximum_auxiliary_stock_error=aux,seconds=time.perf_counter()-tick,original_counter_stride=128,compared_paths=4));print('FORWARD REPLAY PASS',error,'seconds',time.perf_counter()-tick,flush=True)
    parts=[];seeds=[]
    for batch in range(8):
        seed=int(np.random.SeedSequence([20261028577,700,batch]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed);tick=time.perf_counter();parts.append(k.stocks(d,b,ss,p,x,64,64,seed));print('FORWARD FRESH',batch,'seconds',round(time.perf_counter()-tick,2),flush=True)
    stocks=combine(parts);np.savez_compressed(OUT/'stocks.npz',stocks=stocks,seeds=np.array(seeds,dtype=np.uint64));q=pd.read_csv(TARGET/'SPX_membership.csv');q.to_csv(OUT/'SPX_membership.csv',index=False);ref=np.load(LEVERAGE/'reference.npz')['stats'];evaluate_stocks(stocks,ref,arr([x]*4),arr([0]*4),q,OUT/'prices');scores=relabel(OUT/'prices');dump(OUT/'run.json',dict(status='completed four fixed stock couplings, no market-parameter selection',names=['fixed_theta0.8','instantaneous_variance','future30day_4updates','future30day_8updates'],scores=scores,stock_paths=len(stocks),batches=8,antithetic_pairs=512,factor_paths_new=True,variance_parameters_changed=False,VIX_law_unchanged=True,VIX_forward_numerical_accuracy_inherited=True,auxiliary_reference_reused=True,wall_seconds=time.perf_counter()-begin,adopted=False));print('FORWARD STOCK COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':main()
