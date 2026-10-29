"""Coupled128/256 grid check for baseline, adaptive and representative readouts."""
import json,time
import numpy as np
import pandas as pd
from .targeted_skew import OUT,Kernel,ComponentMeans,primitive_bridge,ALPHA
from .asymmetric_readout import PARENT,PRIOR,LEVERAGE,evaluate_stocks,arr,dump
from .protected_curve_check import context,transition_arrays
from .leverage_scan import combine

def main():
    begin=time.perf_counter();folder=OUT/'grid_refinement';folder.mkdir(exist_ok=False);prior,dated,norm,_=context();p=arr(prior['selected']);initial=arr(prior['initial_state']);old=dict(np.load(OUT/'augmented_grid.npz'));times=np.r_[np.column_stack([old['time'],old['time']+.5*old['dt']]).ravel(),old['time'][-1]+old['dt'][-1]];d,ta=transition_arrays(dated,initial,arr(times));n=len(times)-1;reference=np.zeros(n,np.int32);reference[::2]=old['reference'];expiry=np.full(n,-1,np.int32);expiry[1::2]=old['expiry'];d.update(time=arr(times[:-1]),dt=arr(np.diff(times)),seg=np.repeat(old['seg'],2).astype('int32'),coarse=(np.arange(n)%2==0).astype('int32'),reference=reference,expiry=expiry,lognorm=arr(norm(times[:-1])[:,0]),expiries=int(old['expiries']));ba=primitive_bridge(d,dated);dump(folder/'transition_audit.json',dict(gaussian=ta,primitive=ba));np.savez_compressed(folder/'grid.npz',**d);k=Kernel(OUT/'targeted.dylib');config=json.loads((OUT/'candidates.json').read_text());idx=[0,1,4,9];models=arr(config['models'])[idx];alpha=arr(config['alpha'])[idx];mode=np.array(config['mode'],np.int32)[idx];scale=json.loads((OUT/'design.json').read_text())['scale'];M=ComponentMeans(k,dated,initial,p,norm,ALPHA,scale).build();mom=arr(M(d['time'])[:,np.array(config['mean_indices'])[idx]]);parts=[];seeds=[]
    dump(folder/'design.json',dict(original_candidate_indices=idx,steps=[128,256],fresh_paths=2048,selection_not_changed=True,stage='extra grid/precision check motivated by64/128 changes; previously frozen illustrative alpha1.5'))
    for batch in range(8):
        seed=int(np.random.SeedSequence([20261028491,900,batch]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed);tick=time.perf_counter();parts.append(k.stocks(d,128,p,models,alpha,mode,mom,scale,seed));print('REFINED STOCK',batch,'seconds',round(time.perf_counter()-tick,2),flush=True)
    stocks=combine(parts);np.savez_compressed(folder/'stocks.npz',stocks=stocks,seeds=np.array(seeds,dtype=np.uint64));ref=np.load(LEVERAGE/'reference.npz')['stats'];q=pd.read_csv(OUT/'SPX_membership.csv');stage=folder/'prices';scores=evaluate_stocks(stocks,ref,models,alpha,q,stage)
    # The shared evaluator labels its two input grids64/128. Relabel every
    # stored artifact explicitly to the actual128/256 steps used here.
    for path in list(stage.glob('SPX_prices_*.csv')):
        frame=pd.read_csv(path);oldstep=int(frame.steps.iloc[0]);frame['steps']=oldstep*2;frame.to_csv(path.with_name(path.name.replace(f'_{oldstep}.csv',f'_{oldstep*2}.csv')+'.new'),index=False);path.unlink()
    for path in stage.glob('*.new'):path.rename(path.with_name(path.name[:-4]))
    scores=[{**s,'steps':2*s['steps'],'original_candidate':idx[s['candidate']]} for s in scores];dump(stage/'SPX_scores.json',scores);pd.DataFrame(scores).to_csv(stage/'SPX_scores.csv',index=False);blocks=dict(np.load(stage/'SPX_blocks.npz'));np.savez_compressed(stage/'SPX_blocks.npz',**{key if key=='cheap' else key.rsplit('_',1)[0]+'_'+str(2*int(key.rsplit('_',1)[1])):v for key,v in blocks.items()});dump(folder/'run.json',dict(status='completed coupled finer-grid check',wall_seconds=time.perf_counter()-begin,scores=scores));print('REFINEMENT COMPLETE',time.perf_counter()-begin,flush=True)

if __name__=='__main__':main()
