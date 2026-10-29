"""Independent scalar history/filter/stock recurrences and bridge covariance probes."""
import math
from .memory_filter import *
from .audit_leverage_scan import gaussian

def main():
    prior,dated,norm,_=context();initial=arr(prior['initial_state']);p=arr(prior['selected']);x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);k=Kernel(OUT/'memory.dylib');d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));n=97;small=dict(d)
    for key in ['time','dt','seg','coarse','reference','expiry','transition','mean','lognorm']:small[key]=d[key][:n].copy()
    small['expiries']=1;small['expiry'][:]=-1;small['expiry'][-1]=0;cfg=[CONFIG[0],(1,.5,1),(2,.5,1)];models=arr([x]*3);seed=20261029822;got=k.stocks(small,1,p,models,cfg,seed)[:,0];state=np.zeros(24);zs=zf=0;logs=np.zeros((2,3,2));hold=np.zeros((2,3));refv=np.zeros(2);reflog=np.zeros(2);mem=np.ones((2,2,3,3));held=np.zeros((2,3));rates,b=stress_parameters(x);f=(.05+.07*x[19])**2
    for t in range(n):
        h=d['dt'][t];tt=d['time'][t];tr=d['transition'][t];z=np.array([gaussian(seed,27*t+j) for j in range(27)]);dw=np.sqrt(h)*z[0];B=[]
        for j,a in enumerate(rates):
            q=a*h;res=h*(q*q/12-q**4/120+17*q**6/20160) if q<.01 else h*(1-math.tanh(q/2)/(q/2));B.append(np.sqrt(h-res)*z[25+j]+np.sqrt(res)*gaussian(seed^0x2ad82d187dccb289,2*t+j))
        shock=.6*dw-.8*(np.sqrt(1-x[22])*B[0]+np.sqrt(x[22])*B[1]);pr=d['projection']@state
        for anti,sign in enumerate([1,-1]):
            C,Y,F=d['mean'][t]+sign*pr;eta=p[0]*C+p[1]*(Y-p[2])**2+p[3]*(F-p[4])**2;R=np.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6])-d['lognorm'][t]);AA=np.exp(sign*b*np.array([zs,zf])-.5*b*b*(-np.expm1(-2*rates*tt)));A=(1-x[22])*AA[0]+x[22]*AA[1];H=(1-x[17])*R+x[17]*A;inp=np.r_[R,AA]
            if d['coarse'][t]:held[anti]=inp
            for c,(kind,share,j) in enumerate(cfg):
                for grid in range(2):
                    mm=mem[anti,grid,j];filtered=(1-x[17])*mm[0]+x[17]*(A if kind==1 else (1-x[22])*mm[1]+x[22]*mm[2]);response=(1-share)*H+share*filtered;v=f+(.2-f)*x[d['seg'][t]]*response
                    if grid==0 and d['coarse'][t]:hold[anti,c]=v
                    vv=v if grid else hold[anti,c];logs[anti,c,grid]+=-.5*vv*h+np.sqrt(vv)*sign*shock
            for grid in range(2):mem[anti,grid]=np.exp(-RATES*h)[:,None]*mem[anti,grid]+(-np.expm1(-RATES*h))[:,None]*(inp if grid else held[anti])
            if d['reference'][t]:refv[anti]=f+(.2-f)*x[d['seg'][t]]*((1-x[17])+x[17]*A)
            reflog[anti]+=-.5*refv[anti]*h+np.sqrt(refv[anti])*sign*shock
        state=d['decay'][tr]*state+d['root'][tr]@z[1:25]+d['cross'][tr]*dw;zs,zf=np.exp(-rates*h)*np.array([zs,zf])+np.sqrt(-np.expm1(-2*rates*h))*z[25:27]
    expected=np.concatenate([np.exp(logs),np.repeat(np.exp(reflog)[:,None,None],2,axis=2)],axis=1);err=float(abs(got-expected).max());assert err<1e-12
    oldk=OldKernel(OUT/'memory.dylib');a=k.stocks(d,1,p,arr([x]),[CONFIG[0]],seed);bb=oldk.stocks(d,1,p,arr([x]),arr([0]),np.ones((len(d['dt']),1)),1.,seed);replay=float(abs(a-bb).max());assert replay<1e-11
    # Independent conditional Gaussian history recurrence on a short complete interval.
    basis,ba=history_basis(dated,initial,norm,x,1.1,128);terminal=np.zeros((1,24));zt=arr([[.2,-.3]]);hs=20261029823;out=k.history(basis,terminal,zt,p,x,hs)[0];nt=basis['nt'];dt=basis['dt'];S=np.zeros(24);zs=zf=0;paths=[];zz=[]
    for t in range(nt):
        paths.append(S.copy());zz.append([zs,zf]);z=np.array([gaussian(hs,26*t+j) for j in range(26)]);S=basis['decay']*S+basis['root']@z[:24];zs,zf=np.exp(-rates*dt)*np.array([zs,zf])+np.sqrt(-np.expm1(-2*rates*dt))*z[24:26]
    delta=basis['whitening']@(terminal[0]-S);dz=zt[0]-[zs,zf];memory=np.ones((3,3,2));held=np.zeros(3)
    for t in range(nt):
        C,Y,F=basis['mean'][t]+basis['projection']@paths[t]+basis['bridge'][t]@delta;eta=p[0]*C+p[1]*(Y-p[2])**2+p[3]*(F-p[4])**2;R=np.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6])-basis['lognorm'][t]);Z=np.array(zz[t])+basis['stressbridge'][t]*dz;AA=np.exp(b*Z-.5*b*b*basis['stress_q'][t]);inp=np.r_[R,AA]
        if t%2==0:held=inp.copy()
        for j in range(3):
            for grid in range(2):memory[j,:,grid]=np.exp(-RATES[j]*dt)*memory[j,:,grid]+(-np.expm1(-RATES[j]*dt))*(inp if grid else held)
    herr=float(abs(out-memory).max());assert herr<1e-10
    # Positive filters preserve a constant input exactly to roundoff.
    const=np.ones(3)
    for h in [.001,.1,1.]:const=np.exp(-RATES*h)*const-np.expm1(-RATES*h)
    cerr=float(abs(const-1).max());assert cerr<1e-14
    dump(OUT/'implementation_audit.json',dict(scalar_stock_error=err,zero_share_full_stock_replay_error=replay,scalar_conditional_history_error=herr,constant_input_filter_error=cerr,short_history_bridge=ba,scope='finite numerical recurrences, not Lean or a population-forward certificate'));print('MEMORY AUDIT PASS',err,replay,herr,flush=True)
if __name__=='__main__':main()
