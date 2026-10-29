"""Independent scalar evaluation of new readout and both stock time grids."""
import math
from .slow_readout import *
from .audit_leverage_scan import gaussian

def main():
    root=HERE/'runs/m2_slow_readout_002';sel=json.loads((root/'fresh_validation/selection.json').read_text());prior,dated,norm,_=context();p=arr(prior['selected']);models=arr(sel['models']);alpha=arr(sel['alpha']);indices=sel['mode_indices'];d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));Mweak=np.load(root/'component_means.npz')['moments'];Mstrong=np.load(OUT/'component_means.npz')['moments'];chosen=sel['chosen_weak_pilot_candidate'];mom=arr(np.stack([Mweak[:,0],Mweak[:,chosen],Mweak[:,chosen],Mweak[:,0],Mstrong[:,7]],axis=1));n=97;small=dict(d)
    for key in ['time','dt','seg','coarse','reference','expiry','transition','mean','lognorm']:small[key]=d[key][:n].copy()
    small['expiries']=1;small['expiry'][:]=-1;small['expiry'][-1]=0;k=Kernel(root/'slow.dylib');seed=20261028615;got=k.stocks(small,1,p,models,alpha,mom[:n],indices,seed)[:,0];state=np.zeros(24);zs=zf=0;logs=np.zeros((2,5,2));hold=np.zeros((2,5));refv=np.zeros(2);reflog=np.zeros(2);x=models[0];hf=np.exp(np.log(365)*x[21]);hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf]);b=5*x[[16,20]];f=(.05+.07*x[19])**2
    for t in range(n):
        h=d['dt'][t];tt=d['time'][t];tr=d['transition'][t];z=np.array([gaussian(seed,27*t+j) for j in range(27)]);dw=np.sqrt(h)*z[0];B=[]
        for j,a in enumerate(rates):
            q=a*h;res=h*(q*q/12-q**4/120+17*q**6/20160) if q<.01 else h*(1-math.tanh(q/2)/(q/2));B.append(np.sqrt(h-res)*z[25+j]+np.sqrt(res)*gaussian(seed^0x2ad82d187dccb289,2*t+j))
        shock=.6*dw-.8*(np.sqrt(1-x[22])*B[0]+np.sqrt(x[22])*B[1]);proj=d['projection']@state
        for anti,sign in enumerate([1,-1]):
            C,Y,F=d['mean'][t]+sign*proj;eta=p[0]*C+p[1]*(Y-p[2])**2+p[3]*(F-p[4])**2;R=np.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6])-d['lognorm'][t]);AA=np.exp(sign*b*np.array([zs,zf])-.5*b*b*(-np.expm1(-2*rates*tt)));A=(1-x[22])*AA[0]+x[22]*AA[1];H=(1-x[17])*R+x[17]*A
            for c in range(5):
                D=np.exp(alpha[c]*np.logaddexp(0,-sign*state[indices[c]]));MD,MRD=mom[t,c];response=H*D*((1-x[17])+x[17]*A)/((1-x[17])*MRD+x[17]*A*MD);v=f+(.2-f)*models[c,d['seg'][t]]*response
                if d['coarse'][t]:hold[anti,c]=v
                for grid,vv in enumerate([hold[anti,c],v]):logs[anti,c,grid]+=-.5*vv*h+np.sqrt(vv)*sign*shock
            if d['reference'][t]:refv[anti]=f+(.2-f)*x[d['seg'][t]]*((1-x[17])+x[17]*A)
            reflog[anti]+=-.5*refv[anti]*h+np.sqrt(refv[anti])*sign*shock
        state=d['decay'][tr]*state+d['root'][tr]@z[1:25]+d['cross'][tr]*dw;zs,zf=np.exp(-rates*h)*np.array([zs,zf])+np.sqrt(-np.expm1(-2*rates*h))*z[25:27]
    expected=np.concatenate([np.exp(logs),np.repeat(np.exp(reflog)[:,None,None],2,axis=2)],axis=1);err=float(abs(got-expected).max());assert err<1e-12;dump(root/'recurrence_audit.json',dict(maximum_scalar_stock_error=err,steps=n,candidates=5,antithetic_partners=2,grids=[128,256],scope='finite scalar numerical recurrence audit, no continuous-time martingale theorem'));print('SLOW RECURRENCE PASS',err,flush=True)
if __name__=='__main__':main()
