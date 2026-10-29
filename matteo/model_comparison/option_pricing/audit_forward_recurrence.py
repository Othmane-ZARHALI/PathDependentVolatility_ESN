"""Scalar stock replay checks the time filtration and both exposure schedules."""
import math,json
import numpy as np
from .forward_leverage import OUT,setup,arr,dump
from .audit_leverage_scan import gaussian

def main():
    prior,dated,norm,_,p,x,d,k=setup(OUT);b=dict(np.load(OUT/'basis.npz'));ss=dict(np.load(OUT/'signals.npz'));n=97;small=dict(d)
    for key in ['time','dt','seg','coarse','reference','expiry','transition','mean','lognorm']:small[key]=d[key][:n].copy()
    small['expiries']=1;small['expiry'][:]=-1;small['expiry'][-1]=0;signal=dict(ss);signal['update']=ss['update'][:n].copy();seed=20261028699;got=k.stocks(small,b,signal,p,x,1,1,seed)[:,0];nd=24;ob=d['diffusion'].shape[1];nb=ob+2;state=np.zeros(nd);zs=zf=0.;logs=np.zeros((2,4,2));held=np.zeros(2);heldinst=np.zeros((2,nb));future=np.zeros((2,2,nb));refv=np.zeros(2);reflog=np.zeros(2);L0=np.log(365);hf=np.exp(L0*x[21]);hs=hf*(365/hf)**x[18];rates=np.log(2)*365/np.array([hs,hf]);bs=5*x[16];bf=5*x[20];floor=(.05+.07*x[19])**2;unit_errors=[];cov_slacks=[]
    def coeff(a,h):
        z=a*h;r=h*(z*z/12-z**4/120+17*z**6/20160) if z<.01 else h*(1-math.tanh(z/2)/(z/2));return math.sqrt(h-r),math.sqrt(r)
    for t in range(n):
        h=d['dt'][t];time=d['time'][t];tr=d['transition'][t];counter=t*27;normal=np.array([gaussian(seed,counter+j) for j in range(25)]);dw=np.sqrt(h)*normal[0];es=gaussian(seed,counter+25);ef=gaussian(seed,counter+26);cs,rs=coeff(rates[0],h);cf,rf=coeff(rates[1],h);dbs=cs*es+rs*gaussian(seed^0x2ad82d187dccb289,2*t);dbf=cf*ef+rf*gaussian(seed^0x2ad82d187dccb289,2*t+1);extra=np.array([gaussian(seed^0x78d904ef6038fb35,t*(ob-1)+j) for j in range(ob-1)]);db=np.r_[dw,d['bmean'][tr]@normal[1:]+d['broot'][tr]@extra,dbs,dbf];fixed=.6*dw-.8*(np.sqrt(1-x[22])*dbs+np.sqrt(x[22])*dbf);projection=d['projection']@state
        for anti,sign in enumerate([1,-1]):
            C,Y,F=d['mean'][t]+sign*projection;eta=p[0]*C+p[1]*(Y-p[2])**2+p[3]*(F-p[4])**2;R=np.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6])-d['lognorm'][t]);As=np.exp(sign*bs*zs-.5*bs*bs*(-np.expm1(-2*rates[0]*time)));Af=np.exp(sign*bf*zf-.5*bf*bf*(-np.expm1(-2*rates[1]*time)));A=(1-x[22])*As+x[22]*Af;H=(1-x[17])*R+x[17]*A;v=floor+(.2-floor)*x[d['seg'][t]]*H;sig=1/(1+np.exp(p[6]*(eta-p[5])));g=2*(1-x[17])*R*sig*np.array([p[0],2*p[1]*(Y-p[2]),2*p[3]*(F-p[4])]);Kinst=np.r_[g@d['diffusion'],x[17]*(1-x[22])*bs*np.sqrt(2*rates[0])*As,x[17]*x[22]*bf*np.sqrt(2*rates[1])*Af];instant=-Kinst/np.linalg.norm(Kinst)
            if d['coarse'][t]:held[anti]=v;heldinst[anti]=instant
            if ss['update'][t]>=0:
                u=ss['update'][t];s={key:ss[key][u] for key in ['futuremean','inverse_norm','hat']};K=k.probe(arr([sign*state]),arr([[sign*zs,sign*zf]]),s,b,p,x,time)[0];a=-K/np.linalg.norm(K);future[anti,1]=a
                if u%2==0:future[anti,0]=a
                unit_errors.append(abs(a@a-1));a0=np.zeros(nb);a0[0]=.6;a0[-2]=-.8*np.sqrt(1-x[22]);a0[-1]=-.8*np.sqrt(x[22]);cov_slacks.append(float(a@K-a0@K));assert a@K<=a0@K+1e-12
            for grid,vv in enumerate([held[anti],v]):
                shocks=[fixed,(instant if grid else heldinst[anti])@db,future[anti,0]@db,future[anti,1]@db];logs[anti,:,grid]+=-.5*vv*h+np.sqrt(vv)*sign*np.array(shocks)
            if d['reference'][t]:refv[anti]=floor+(.2-floor)*x[d['seg'][t]]*((1-x[17])+x[17]*A)
            reflog[anti]+=-.5*refv[anti]*h+np.sqrt(refv[anti])*sign*fixed
        state=d['decay'][tr]*state+d['root'][tr]@normal[1:]+d['cross'][tr]*dw;zs=np.exp(-rates[0]*h)*zs+np.sqrt(-np.expm1(-2*rates[0]*h))*es;zf=np.exp(-rates[1]*h)*zf+np.sqrt(-np.expm1(-2*rates[1]*h))*ef
    expected=np.concatenate([np.exp(logs),np.repeat(np.exp(reflog)[:,None,None],2,axis=2)],axis=1);err=float(np.max(abs(got-expected)));assert err<1e-12;dump(OUT/'recurrence_audit.json',dict(all_checks_passed=True,scalar_stock_recursion_max_absolute_error=err,unit_forward_stock_loading_max_error=max(unit_errors),future_covariance_not_greater_than_baseline_at_update_times=True,scope='97fine steps, both antithetic paths, coarse/fine variance and4/8update schedules; future-loading function independently checked by finite differences in integration_audit.json'));print('FORWARD RECURRENCE AUDIT PASS',err,max(unit_errors),flush=True)

if __name__=='__main__':main()
