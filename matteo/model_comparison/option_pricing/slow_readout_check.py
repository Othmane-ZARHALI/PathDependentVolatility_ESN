"""Freeze a weak-signal candidate, then check new paired stock and VIX samples."""
from .slow_readout import *
from .asymmetric_readout import OUT as OLDOUT
from .vix_only_curve_fit import ARCHIVE
from .mixed_response_study import ZSHIFTS,ZPROB
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.importance_sampling import DefensiveGaussianMixture
WEAK=HERE/'runs/m2_slow_readout_002'

def run():
    folder=WEAK/'fresh_validation';folder.mkdir(exist_ok=False);cand=json.loads((WEAK/'candidates.json').read_text());ledger=json.loads((WEAK/'VIX_fit_ledger.json').read_text());scores=json.loads((WEAK/'SPX/SPX_scores.json').read_text());feasible=[s for s in scores if s['steps']==256 and s['candidate']>=7 and ledger[s['candidate']-6]['eligible'] and ledger[s['candidate']-6]['loss']<=ledger[0]['loss']+.1];chosen=min(feasible,key=lambda s:s['IV_MAE_pp'])['candidate'];fixed=chosen-6;curve=arr(cand['models'][chosen]);base=arr(cand['models'][0]);strong=json.loads((OUT/'candidates.json').read_text());models=arr([base,base,curve,curve,strong['models'][7]]);alpha=arr([0,cand['alpha'][chosen],cand['alpha'][chosen],0,.5]);idx=int(cand['mode_indices'][chosen]);modes=[idx]*4+[strong['mode_indices'][7]];mi=[0,cand['mean_indices'][chosen],cand['mean_indices'][chosen],0,1];names=['baseline','weak_fixed_curve','weak_refit','curve_only','strong_refit_contrast'];dump(folder/'selection.json',dict(chosen_weak_pilot_candidate=chosen,chosen_name=cand['names'][chosen],criterion='training VIX loss <= alpha0-refit + .1 with four training equalities; smallest 512-path fine-grid SPX IV MAE among weak refits',selection_before_new_stock_and_outer_draws=True,models=models.tolist(),alpha=alpha.tolist(),mode_indices=modes,mean_indices=mi,names=names,adopted=False))
    print('FROZEN SELECTION',chosen,cand['names'][chosen],flush=True)
    prior,dated,norm,_=context();initial=arr(prior['initial_state']);p=arr(prior['selected']);design=json.loads((WEAK/'design.json').read_text());k=Kernel(WEAK/'slow.dylib');d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));oldmom=np.load(WEAK/'component_means.npz')['moments'];strongmom=np.load(OUT/'component_means.npz')['moments'];mom=arr(np.stack([oldmom[:,0],oldmom[:,chosen],oldmom[:,chosen],oldmom[:,0],strongmom[:,7]],axis=1));pieces=[];seeds=[]
    for b in range(8):
        seed=int(np.random.SeedSequence([20261028613,777,b]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed);tick=time.perf_counter();pieces.append(k.stocks(d,64,p,models,alpha,mom,modes,seed));print('NEW STOCK',b,round(time.perf_counter()-tick,2),flush=True)
    stocks=combine(pieces);np.savez_compressed(folder/'stocks.npz',stocks=stocks,seeds=np.array(seeds,dtype=np.uint64));q=pd.read_csv(PARENT/'spx_membership.csv');evaluate_stocks(stocks,np.load(LEVERAGE/'reference.npz')['stats'],models,alpha,q,folder/'SPX');relabel(folder/'SPX')
    # New complete 38-dimensional proposal draws, not new factor dynamics.
    props=json.loads((ARCHIVE/'m2_vix_targeted_check_001/proposals.json').read_text());draws={};(folder/'VIX_draws').mkdir();ds=[]
    for day in DAYS:
        dec=props[str(day)];shifts=arr([np.r_[old,s,f] for old in dec['shifts'] for s in ZSHIFTS for f in ZSHIFTS]);prob=arr([po*ps*pf for po in dec['probabilities'] for ps in ZPROB for pf in ZPROB]);mix=DefensiveGaussianMixture(shifts,prob);sampler=GaussianOuterSampler(dated,day/365,initial);zz=[];ww=[]
        for batch in range(8):
            seed=int(np.random.SeedSequence([20261028613,888,day,batch]).generate_state(1,dtype=np.uint64)[0]);z,w,_=mix.draw(8,seed);zz.append(z);ww.append(w);ds.append(dict(day=day,batch=batch,seed=seed))
        z=np.concatenate(zz);a=dict(states=sampler.mean+z[:,:-2]@sampler.root.T,z_normals=z[:,-2],fast_normals=z[:,-1],weights=np.concatenate(ww),nets=np.repeat(np.arange(8),256),source_indices=np.arange(len(z)));draws[day]=a;np.savez_compressed(folder/f'VIX_draws/d{day}.npz',**a)
    dump(folder/'VIX_draws/design.json',dict(seeds=ds,points_per_expiry=2048,dimension=38,proposal='joint defensive mixture, exact target/proposal weights, no mass normalization',new_draws=True))
    # Reuse exposure engine with a single strength per selected memory and a zero-strength control.
    import matteo.model_comparison.option_pricing.slow_readout as s
    s.LEVELS=arr([alpha[1],.5]);laws=[Law(dated,idx),Law(dated,modes[-1]),Law(dated,idx)]
    M=Means(k,laws,initial,p,norm).build();H=exposures(k,laws,initial,p,norm,M,base,draws,32,8,folder/'VIX_exposures');qv=pd.read_csv(OLDOUT/'VIX_membership.csv');f=pd.read_csv(OLDOUT/'forward_targets.csv');v=VIX(H,draws,base,qv,f);(folder/'VIX').mkdir();out=[];hidx=[0,1,1,0,4]
    for c in range(5):
        rec=save_vix(v,hidx[c],models[c,ACTIVE],f'c{c}',folder/'VIX');rec.update(candidate=c,name=names[c]);r=v.predict(hidx[c],models[c,ACTIVE]);tails=[]
        for day,a in r['raw'].items():
            V=a['V'];w=a['weights'];K=float(qv[qv.expiry_days==day].strike.max());y=w*np.maximum(V-K,0);tot=y.sum();tails.append(dict(day=int(day),weight_mass=float(w.mean()),right_tail_ESS=float(tot*tot/(y@y)) if tot else 0,top10_payoff_share=float(np.sort(y)[-10:].sum()/tot) if tot else 0))
        rec['tails']=tails;out.append(rec)
        for stem in ['prices','forwards']:
            path=folder/f'VIX/VIX_{stem}_c{c}.csv';g=pd.read_csv(path)
            if stem=='prices':g['combined_SE_U']=np.sqrt(g.difference_SE_U**2+g.anchor_SE_U**2)
            else:g['combined_SE']=np.sqrt(g.difference_SE**2+g.anchor_SE**2)
            g.to_csv(path,index=False)
    dump(folder/'VIX/scores.json',out)
    # Matched 64-state refinement of Gaussian/time quadrature.
    ids={day:np.concatenate([np.flatnonzero(a['nets']==j)[:8] for j in range(8)]) for day,a in draws.items()};subset={day:{key:value[ids[day]] for key,value in a.items()} for day,a in draws.items()};fine=exposures(k,laws,initial,p,norm,M,base,subset,48,12,folder/'quadrature_probe');audit=[]
    from .two_speed_joint_core import physical
    for c in [1,2,4]:
        levels,floor=physical(models[c]);j=hidx[c]
        for day in DAYS:
            V0=100*np.sqrt((floor*30/365+H[day][j,ids[day]]@(levels-floor))/(30/365));V1=100*np.sqrt((floor*30/365+fine[day][j]@(levels-floor))/(30/365));audit.append(dict(candidate=c,day=int(day),max_VIX_change=float(abs(V1-V0).max()),weighted_mean_VIX_change=float(np.mean(subset[day]['weights']*(V1-V0)))))
    dump(folder/'quadrature_probe.json',audit);dump(folder/'run.json',dict(status='complete independent paired checking',stock_paths=1024,VIX_states_per_expiry=2048,selection_frozen=True,production_baseline_unchanged=True));print('FRESH VALIDATION COMPLETE',flush=True)
if __name__=='__main__':run()
