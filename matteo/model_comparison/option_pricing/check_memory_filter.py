"""Selection-frozen stock and complete terminal-plus-history VIX checking."""
from .memory_filter import *
from .vix_only_curve_fit import ARCHIVE
from .mixed_response_study import ZSHIFTS,ZPROB
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture
from .two_speed_joint_core import physical

def run():
    folder=OUT/'fresh_validation';folder.mkdir(exist_ok=False);cc=json.loads((OUT/'candidates.json').read_text());ledger=json.loads((OUT/'VIX_fit_ledger.json').read_text());pilot=json.loads((OUT/'SPX/SPX_scores.json').read_text());selection=json.loads((OUT/'selection.json').read_text());fine=[r for r in pilot if r['steps']==256 and r['candidate']>=14]
    def choose(kind):
        rows=[s for s in fine if cc['config'][s['candidate']][0]==kind];feasible=[s for s in rows if ledger[cc['exposure_indices'][s['candidate']]]['eligible'] and ledger[cc['exposure_indices'][s['candidate']]]['loss']<=ledger[0]['loss']+.1];return min(feasible,key=lambda s:s['IV_MAE_pp'])['candidate'] if feasible else min(rows,key=lambda s:ledger[cc['exposure_indices'][s['candidate']]]['loss'])['candidate']
    chosen=choose(1);whole=choose(2);fixed=chosen-13;models=arr([cc['models'][0],cc['models'][fixed],cc['models'][chosen],cc['models'][chosen],cc['models'][whole]]);cfg=[CONFIG[0],tuple(cc['config'][chosen]),tuple(cc['config'][chosen]),CONFIG[0],tuple(cc['config'][whole])];idx=[0,cc['exposure_indices'][chosen],cc['exposure_indices'][chosen],0,cc['exposure_indices'][whole]];names=['baseline','M2_filter_fixed_curve','M2_filter_refit','curve_only','whole_filter_refit'];dump(folder/'selection.json',dict(global_pilot_selection=selection,M2_candidate=chosen,M2_name=cc['names'][chosen],whole_candidate=whole,whole_name=cc['names'][whole],names=names,models=models.tolist(),config=cfg,exposure_indices=idx,selection_before_fresh_draws=True,adopted=False));print('FILTER FROZEN CHECK',cc['names'][chosen],cc['names'][whole],flush=True)
    prior,dated,norm,_=context();p=arr(prior['selected']);initial=arr(prior['initial_state']);x=arr(cc['models'][0]);k=Kernel(OUT/'memory.dylib');d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));pieces=[];seeds=[]
    for b in range(8):
        seed=int(np.random.SeedSequence([20261029831,777,b]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed);tick=time.perf_counter();pieces.append(k.stocks(d,64,p,models,cfg,seed));print('FILTER NEW STOCK',b,round(time.perf_counter()-tick,2),flush=True)
    stocks=combine(pieces);np.savez_compressed(folder/'stocks.npz',stocks=stocks,seeds=np.array(seeds,dtype=np.uint64));q=pd.read_csv(PARENT/'spx_membership.csv');evaluate_stocks(stocks,np.load(LEVERAGE/'reference.npz')['stats'],models,arr([c[1] for c in cfg]),q,folder/'SPX');relabel(folder/'SPX')
    draws={};props=json.loads((ARCHIVE/'m2_vix_targeted_check_001/proposals.json').read_text());(folder/'VIX_draws').mkdir();seeds=[]
    for day in DAYS:
        dec=props[str(day)];shifts=arr([np.r_[old,s,f] for old in dec['shifts'] for s in ZSHIFTS for f in ZSHIFTS]);prob=arr([po*ps*pf for po in dec['probabilities'] for ps in ZPROB for pf in ZPROB]);mix=DefensiveGaussianMixture(shifts,prob);sampler=GaussianOuterSampler(dated,day/365,initial);zs=[];ws=[]
        for batch in range(8):
            seed=int(np.random.SeedSequence([20261029831,888,day,batch]).generate_state(1,dtype=np.uint64)[0]);z,w,_=mix.draw(7,seed);zs.append(z);ws.append(w);seeds.append(dict(day=day,batch=batch,seed=seed))
        z=np.concatenate(zs);draws[day]=dict(states=sampler.mean+z[:,:-2]@sampler.root.T,z_normals=z[:,-2],fast_normals=z[:,-1],weights=np.concatenate(ws),nets=np.repeat(np.arange(8),128),source_indices=np.arange(len(z)));np.savez_compressed(folder/f'VIX_draws/d{day}.npz',**draws[day])
    dump(folder/'VIX_draws/design.json',dict(seeds=seeds,points_per_expiry=1024,dimension=38,proposal='joint defensive Gaussian mixture, exact target/proposal weights, no mass normalization',new_states_and_histories=True))
    mem=histories(k,dated,initial,norm,p,x,draws,folder/'histories',20261029831);H=exposures(k,dated,initial,p,norm,x,draws,mem,folder/'exposures');qv=pd.read_csv(OLDOUT/'VIX_membership.csv');f=pd.read_csv(OLDOUT/'forward_targets.csv');v=VIX(H,draws,x,qv,f);(folder/'VIX').mkdir();scores=[];history_audit=[];diagnostics=[]
    for c in range(5):
        r=v.predict(idx[c],models[c,ACTIVE]);rec=save_vix(v,idx[c],models[c,ACTIVE],f'c{c}',folder/'VIX');rec.update(candidate=c,name=names[c]);tails=[]
        for day,a in r['raw'].items():
            V=a['V'];w=a['weights'];K=float(qv[qv.expiry_days==day].strike.max());pay=w*np.maximum(V-K,0);tot=pay.sum();tails.append(dict(day=int(day),weight_mass=float(w.mean()),right_tail_ESS=float(tot*tot/(pay@pay)) if tot else 0,top10_payoff_share=float(np.sort(pay)[-10:].sum()/tot) if tot else 0))
            levels,floor=physical(models[c]);kind,share,node=cfg[c];h=np.load(folder/f'exposures/d{day}.npz');times=h['times'];tw=h['time_weights'];weight=np.sum(tw*(levels[np.minimum(np.searchsorted(np.array(prior['ends_years']),times,side='right'),15)]-floor)*np.exp(-RATES[node]*(times-day/365)));dm=mem[day][:,node,:,0]-mem[day][:,node,:,1];delta=(1-x[17])*dm[:,0]
            if kind!=1:delta+=x[17]*((1-x[22])*dm[:,1]+x[22]*dm[:,2])
            I=V*V*(30/365)/10000;coarse=100*np.sqrt((I+share*weight*delta)/(30/365));dv=coarse-V;history_audit.append(dict(candidate=c,day=int(day),max_VIX_change=float(abs(dv).max()),weighted_mean_change=float(np.mean(w*dv))))
        rec['tails']=tails;scores.append(rec)
        for stem in ['prices','forwards']:
            path=folder/f'VIX/VIX_{stem}_c{c}.csv';g=pd.read_csv(path)
            if stem=='prices':g['combined_SE_U']=np.sqrt(g.difference_SE_U**2+g.anchor_SE_U**2)
            else:g['combined_SE']=np.sqrt(g.difference_SE**2+g.anchor_SE**2)
            g.to_csv(path,index=False)
    dump(folder/'VIX/scores.json',scores);dump(folder/'history_grid_audit.json',history_audit)
    for day,a in draws.items():
        w=a['weights'];B=(w[:,None,None]*mem[day][:,:,:,1]).reshape(8,-1,3,3).mean(1);diagnostics.append(dict(day=int(day),weighted_memory_means=B.mean(0).tolist(),SE=B.std(0,ddof=1).tolist(),note='SE entries below divided by sqrt8 in mean_SE; no re-normalization',mean_SE=(B.std(0,ddof=1)/np.sqrt(8)).tolist()))
    dump(folder/'memory_mean_diagnostics.json',diagnostics)
    ids={day:np.concatenate([np.flatnonzero(a['nets']==j)[:8] for j in range(8)]) for day,a in draws.items()};sub={day:{key:value[ids[day]] for key,value in a.items()} for day,a in draws.items()};sm={day:mem[day][ids[day]] for day in DAYS};HH=exposures(k,dated,initial,p,norm,x,sub,sm,folder/'quadrature_probe',48,12,64);qa=[]
    for c in [1,2,4]:
        levels,floor=physical(models[c]);j=idx[c]
        for day in DAYS:
            V0=100*np.sqrt((floor*30/365+H[day][j,ids[day]]@(levels-floor))/(30/365));V1=100*np.sqrt((floor*30/365+HH[day][j]@(levels-floor))/(30/365));qa.append(dict(candidate=c,day=int(day),max_VIX_change=float(abs(V1-V0).max()),weighted_mean_change=float(np.mean(sub[day]['weights']*(V1-V0)))))
    dump(folder/'quadrature_audit.json',qa);dump(folder/'run.json',dict(status='complete fresh checking',stock_paths=1024,VIX_terminal_states_per_expiry=1024,new_conditional_history_per_terminal=True,adopted=False));print('MEMORY FRESH CHECK COMPLETE',flush=True)
if __name__=='__main__':run()
