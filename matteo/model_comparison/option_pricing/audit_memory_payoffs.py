"""Translate history/quadrature changes into likelihood-weighted native payoffs."""
from .memory_filter import *
from .two_speed_joint_core import physical
from .protected_curve_check import loss

def main():
    folder=OUT/'fresh_validation';sel=json.loads((folder/'selection.json').read_text());models=arr(sel['models']);cfg=sel['config'];idx=sel['exposure_indices'];q=pd.read_csv(OLDOUT/'VIX_membership.csv');prior=json.loads((PRIOR/'run.json').read_text());ends=np.array(prior['ends_years']);x=models[0];records=[];all_diffs={}
    for mode in ['history','future_quadrature']:
        diffs={c:[] for c in range(5)}
        for day in DAYS:
            raw=dict(np.load(folder/f'VIX_draws/d{day}.npz'));hh=dict(np.load(folder/f'exposures/d{day}.npz'));M=np.load(folder/f'histories/d{day}.npz')['memory'];H=hh['H'];weights=raw['weights'];ids=np.arange(len(weights)) if mode=='history' else np.concatenate([np.flatnonzero(raw['nets']==j)[:8] for j in range(8)]);w=weights[ids];g=q[q.expiry_days==day];K=g.strike.to_numpy();sign=np.where(g.option_type.eq('put'),-1,1);baselevels,floor=physical(x);base0=100*np.sqrt((floor*30/365+H[0,ids]@(baselevels-floor))/(30/365));base1=base0.copy()
            if mode!='history':HF=np.load(folder/f'quadrature_probe/d{day}.npz')['H'];base1=100*np.sqrt((floor*30/365+HF[0]@(baselevels-floor))/(30/365))
            for c in range(5):
                levels,floor=physical(models[c]);V0=100*np.sqrt((floor*30/365+H[idx[c],ids]@(levels-floor))/(30/365))
                if mode=='history':
                    kind,share,node=cfg[c];times=hh['times'];tw=hh['time_weights'];coeff=np.sum(tw*(levels[np.minimum(np.searchsorted(ends,times,side='right'),15)]-floor)*np.exp(-RATES[node]*(times-day/365)));dm=M[:,node,:,0]-M[:,node,:,1];delta=(1-x[17])*dm[:,0]
                    if kind!=1:delta+=x[17]*((1-x[22])*dm[:,1]+x[22]*dm[:,2])
                    V1=100*np.sqrt((V0*V0*(30/365)/10000+share*coeff*delta)/(30/365))
                else:V1=100*np.sqrt((floor*30/365+HF[idx[c]]@(levels-floor))/(30/365))
                values=w[:,None]*((np.maximum(sign*(V1[:,None]-K),0)-np.maximum(sign*(base1[:,None]-K),0))-(np.maximum(sign*(V0[:,None]-K),0)-np.maximum(sign*(base0[:,None]-K),0)));deltaP=values.mean(0);B=values.reshape(8,-1,len(g)).mean(1);differences=V1-V0;bound=float(np.mean(w*(abs(differences)+abs(base1-base0))));records.append(dict(mode=mode,candidate=c,day=int(day),points=len(w),weighted_absolute_VIX_change=float(np.mean(w*abs(differences))),native_payoff_Lipschitz_mean_bound=bound,max_absolute_native_price_change=float(abs(deltaP).max()),max_native_price_change_quote_scales=float(np.max(abs(deltaP)/g.price_scale_U.to_numpy())),max_paired_native_price_SE=float((B.std(0,ddof=1)/np.sqrt(8)).max()),max_unweighted_VIX_change=float(abs(differences).max())));diffs[c].extend(deltaP.tolist())
        all_diffs[mode]=diffs
    for mode,diffs in all_diffs.items():
        for c in range(5):
            frame=pd.read_csv(folder/f'VIX/VIX_prices_c{c}.csv');frame[f'{mode}_price_change_U']=diffs[c];frame.to_csv(folder/f'VIX/{mode}_price_changes_c{c}.csv',index=False)
    dump(folder/'weighted_numerical_payoff_audit.json',records);print('WEIGHTED MEMORY PAYOFF AUDIT',max(a['max_native_price_change_quote_scales'] for a in records),'quote scales',flush=True)
if __name__=='__main__':main()
