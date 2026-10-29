"""Gaussian covariance reconstruction and terminal-conditioning rank sensitivity."""
from .memory_filter import *

def main():
    prior,dated,norm,_=context();initial=arr(prior['initial_state']);p=arr(prior['selected']);x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);k=Kernel(OUT/'memory.dylib');r=dated.layout.rates[:24];P=dated.projection[:,:24];S=dated.stationary[:24,:24];draws=load_outer(1,20261029825);records=[]
    for day,a in draws.items():
        b,ba=history_basis(dated,initial,norm,x,day);T=day/365;nt=b['nt'];h=b['dt'];t=np.arange(nt)*h;CT=S*(-np.expm1(-(r[:,None]+r[None,:])*T));ev,Q=np.linalg.eigh(CT);rates,_=stress_parameters(x);zt=np.sqrt(-np.expm1(-2*rates*T))[None,:]*np.column_stack([a['z_normals'],a['fast_normals']]);terminal=a['states'][:,:24]-b['terminalmean'];original=k.history(b,terminal,zt,p,x,20261029825+day)
        lo=dict(b);keep=ev>ev.max()*1e-10;cross=S[None,:,:]*(-np.expm1(-(r[None,:,None]+r[None,None,:])*t[:,None,None]))*np.exp(-r[None,None,:]*(T-t[:,None,None]));beta=(P[None,:,:]@cross@Q[:,keep])/np.sqrt(ev[keep])[None,None,:];BB=np.zeros_like(b['bridge']);BB[:,:,:keep.sum()]=beta;WW=np.zeros((24,24));WW[:keep.sum()]=Q[:,keep].T/np.sqrt(ev[keep])[:,None];lo.update(bridge=arr(BB),whitening=arr(WW));other=k.history(lo,terminal,zt,p,x,20261029825+day)
        # Endpoint projection error is evaluated on actual terminal draws, not arbitrary unattainable state directions.
        keep=ev>ev.max()*1e-11;Bend=np.zeros((3,24));Bend[:,:keep.sum()]=(P@CT@Q[:,keep])/np.sqrt(ev[keep])[None,:];actual=terminal@P.T;reconstructed=(terminal@b['whitening'].T)@Bend.T;endpoint=float(abs(actual-reconstructed).max())
        # Exact covariance of the implemented fixed-step Gaussian recurrence, evaluated in extended precision.
        L=np.longdouble;rr=r.astype(L);AA=rr[:,None]+rr[None,:];QQ=b['root'].astype(L)@b['root'].astype(L).T;CTsim=QQ*(-np.expm1(-AA*L(T)))/(-np.expm1(-AA*L(h)));CTtarget=CT.astype(L);PP=P.astype(L);W=b['whitening'].astype(L);errors=[]
        for j in np.linspace(0,nt-1,9,dtype=int):
            tt=L(t[j]);Cs=QQ*(-np.expm1(-AA*tt))/(-np.expm1(-AA*L(h)));X=Cs*np.exp(-rr[None,:]*(L(T)-tt));K=b['bridge'][j].astype(L)@W;got=PP@Cs@PP.T-PP@X@K.T-K@X.T@PP.T+K@CTsim@K.T+K@CTtarget@K.T;target=PP@(S.astype(L)*(-np.expm1(-AA*tt)))@PP.T;errors.append(float(abs(got-target).max()))
        records.append(dict(day=int(day),actual_terminal_projection_max_error=endpoint,max_unconditional_projected_bridge_covariance_error=max(errors),rank_cutoff_1e10_vs1e11_max_memory_change=float(abs(original-other).max()),rank_cutoff_max_M2_memory_change=float(abs(original[:,:,0]-other[:,:,0]).max()),rank_cutoff_mean_memory_change=np.mean(original-other,axis=0).tolist(),scope='8 matched terminal draws per expiry; rank perturbation changes conditioning, not a uniform pricing bound'))
    dump(OUT/'bridge_precision_audit.json',records);print('MEMORY BRIDGE PRECISION',records,flush=True)
if __name__=='__main__':main()
