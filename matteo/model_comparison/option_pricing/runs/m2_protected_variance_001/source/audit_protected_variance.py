"""Independent scalar stock recurrence, including signed variance corrections."""
import math
from .m2_protected_variance import *
from .audit_leverage_scan import gaussian

def main():
    design=json.loads((OUT/'design.json').read_text());x=arr(design['baseline']);prior,dated,_,_=context();p=arr(prior['selected']);rates=arr(dated.layout.rates[MODES]);d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));n=97;small=dict(d)
    for key in ['time','dt','seg','coarse','reference','expiry','transition','mean','lognorm']:small[key]=d[key][:n].copy()
    # Synthetic clock probes span envelope features; this is a recurrence audit,
    # not a second stochastic model or a pricing estimate.
    small['time']=arr(np.linspace(0,1030/365,n));small['expiries']=1;small['expiry'][:]=-1;small['expiry'][-1]=0
    coeff=arr([[0,0,0],[-1,0,0],[0,1,0],[0,0,-1],[.2,-.3,.5]])
    k=Kernel(OUT/'protected.dylib');seed=20261029099;got=k.stocks(small,1,p,x,coeff,rates,seed)[:,0]
    state=np.zeros(24);zs=zf=0.;logs=np.zeros((2,len(coeff),2));hold=np.zeros((2,len(coeff)));refv=np.zeros(2);reflog=np.zeros(2)
    hf=np.exp(np.log(365)*x[21]);hs=hf*(365/hf)**x[18];sr=np.log(2)*365/np.array([hs,hf]);b=5*x[[16,20]];f=(.05+.07*x[19])**2;minv=1.
    for t in range(n):
        h=small['dt'][t];tt=small['time'][t];tr=small['transition'][t];z=np.array([gaussian(seed,27*t+j) for j in range(27)]);dw=np.sqrt(h)*z[0];B=[]
        for j,a in enumerate(sr):
            q=a*h;res=h*(q*q/12-q**4/120+17*q**6/20160) if q<.01 else h*(1-math.tanh(q/2)/(q/2));B.append(np.sqrt(h-res)*z[25+j]+np.sqrt(res)*gaussian(seed^0x2ad82d187dccb289,2*t+j))
        shock=.6*dw-.8*(np.sqrt(1-x[22])*B[0]+np.sqrt(x[22])*B[1]);pr=small['projection']@state
        for anti,sign in enumerate([1,-1]):
            C,Y,F=small['mean'][t]+sign*pr;eta=p[0]*C+p[1]*(Y-p[2])**2+p[3]*(F-p[4])**2;R=np.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6])-small['lognorm'][t]);AA=np.exp(sign*b*np.array([zs,zf])-.5*b*b*(-np.expm1(-2*sr*tt)));A=(1-x[22])*AA[0]+x[22]*AA[1];H=(1-x[17])*R+x[17]*A
            for c in range(len(coeff)):
                v=f+(.2-f)*x[small['seg'][t]]*H+.5*f*np.dot(coeff[c],correction(tt,sign*state[MODES],rates));minv=min(minv,v)
                if small['coarse'][t]:hold[anti,c]=v
                for grid,vv in enumerate([hold[anti,c],v]):logs[anti,c,grid]+=-.5*vv*h+np.sqrt(vv)*sign*shock
            if small['reference'][t]:refv[anti]=f+(.2-f)*x[small['seg'][t]]*((1-x[17])+x[17]*A)
            reflog[anti]+=-.5*refv[anti]*h+np.sqrt(refv[anti])*sign*shock
        state=small['decay'][tr]*state+small['root'][tr]@z[1:25]+small['cross'][tr]*dw;zs,zf=np.exp(-sr*h)*np.array([zs,zf])+np.sqrt(-np.expm1(-2*sr*h))*z[25:27]
    expected=np.concatenate([np.exp(logs),np.repeat(np.exp(reflog)[:,None,None],2,axis=2)],axis=1);err=float(abs(got-expected).max());assert err<1e-12 and minv>=f/2
    dump(OUT/'recurrence_audit.json',dict(maximum_scalar_stock_error=err,steps=n,candidates=len(coeff),antithetic_partners=2,grids=[128,256],minimum_probe_variance=minv,synthetic_clock=True,scope='finite formula and update-order correspondence; not a continuous-time martingale theorem'));print('PROTECTED RECURRENCE PASS',err,flush=True)
if __name__=='__main__':main()
