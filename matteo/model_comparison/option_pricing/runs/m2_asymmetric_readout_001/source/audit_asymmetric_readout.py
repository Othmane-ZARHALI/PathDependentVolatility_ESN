"""Independent quadrature, stock recursion and forward-fit derivative checks."""
import json,math
import numpy as np
from .asymmetric_readout import OUT,PARENT,PRIOR,Kernel,ALPHA,load_outer,VIX,arr,dump
from .protected_curve_check import context
from .audit_leverage_scan import gaussian
from esn_eq.pricing import _normal_rule
from esn_eq.states import symmetric_square_root

def main():
    k=Kernel(OUT/'asymmetric.dylib');prior,dated,norm,_=context();p=arr(prior['selected']);initial=arr(prior['initial_state']);design=json.loads((OUT/'design.json').read_text());scale=design['scale'];momerr=[]
    for lag in [1e-5,.003,.05]:
        m,cov=dated.gaussian_projection_law(lag,initial);means=np.array([m,m+[0,-2*scale,0],m+[1,scale,2]]);order=32;got=k.moments(means,cov,p,m[1],scale,ALPHA,order);z,w=_normal_rule(order);offset=z@symmetric_square_root(cov).T;expected=np.zeros_like(got)
        for i,mean in enumerate(means):
            C,Y,F=(mean+offset).T;eta=p[0]*C+p[1]*(Y-p[2])**2+p[3]*(F-p[4])**2;G=np.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6]));D=np.exp(np.logaddexp(0,-(Y-m[1])/scale)[:,None]*ALPHA);expected[i,:,0]=w@D;expected[i,:,1]=(w*G)@D
        err=float(np.max(abs(got-expected)/(1+abs(expected))));assert err<2e-8,err;momerr.append(err)
    candidates=json.loads((OUT/'candidates.json').read_text());models=arr(candidates['models']);alpha=arr(candidates['alpha']);d=dict(np.load(PARENT/'fresh_check/grid_and_transitions.npz'));Zfile=np.load(OUT/'normalizers.npz');Z=np.column_stack([np.ones(len(d['time'])),Zfile['Z']]);n=191;small=dict(d)
    for key in ['time','dt','seg','coarse','reference','expiry','transition','mean','lognorm']:small[key]=d[key][:n].copy()
    small['expiries']=1;small['expiry'][:]=-1;small['expiry'][-1]=0;seed=20261028337;got=k.stocks(small,1,p,models,alpha,arr(Z[:n]),scale,seed)[:,0];state=np.zeros(24);zs=zf=0.;logs=np.zeros((2,len(alpha),2));hold=np.zeros((2,len(alpha)));reflog=np.zeros(2);refhold=np.zeros(2);x=models[0];floor=(.05+.07*x[19])**2;af=math.log(2)*365/math.exp(math.log(365)*x[21]);a_s=math.log(2)*365/math.exp(math.log(365)*x[21]+(math.log(365)-math.log(365)*x[21])*x[18]);bs=5*x[16];bf=5*x[20]
    def coeff(a,h):
        z=a*h;rem=h*(z*z/12-z**4/120+17*z**6/20160) if z<.01 else h*(1-math.tanh(z/2)/(z/2));return math.sqrt(h-rem),math.sqrt(rem)
    for t in range(n):
        h=d['dt'][t];time=d['time'][t];counter=t*27;z=np.array([gaussian(seed,counter+j) for j in range(25)]);dw=math.sqrt(h)*z[0];es=gaussian(seed,counter+25);ef=gaussian(seed,counter+26);cs,rs=coeff(a_s,h);cf,rf=coeff(af,h);dbs=cs*es+rs*gaussian(seed^0x2ad82d187dccb289,2*t);dbf=cf*ef+rf*gaussian(seed^0x2ad82d187dccb289,2*t+1);shock=.6*dw-.8*(math.sqrt(1-x[22])*dbs+math.sqrt(x[22])*dbf);proj=d['projection']@state
        for anti,sign in enumerate([1,-1]):
            C,Y,F=d['mean'][t]+sign*proj;eta=p[0]*C+p[1]*(Y-p[2])**2+p[3]*(F-p[4])**2;old=math.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6])-d['lognorm'][t]);stress=(1-x[22])*math.exp(sign*bs*zs-.5*bs*bs*(-math.expm1(-2*a_s*time)))+x[22]*math.exp(sign*bf*zf-.5*bf*bf*(-math.expm1(-2*af*time)));H=(1-x[17])*old+x[17]*stress
            for a in range(len(alpha)):
                response=H*math.exp(alpha[a]*np.logaddexp(0,-sign*proj[1]/scale))/Z[t,a];v=floor+(.2-floor)*models[a,d['seg'][t]]*response
                if d['coarse'][t]:hold[anti,a]=v
                for grid,vv in enumerate([hold[anti,a],v]):logs[anti,a,grid]+=-.5*vv*h+math.sqrt(vv)*sign*shock
            if d['reference'][t]:refhold[anti]=floor+(.2-floor)*x[d['seg'][t]]*((1-x[17])+x[17]*stress)
            reflog[anti]+=-.5*refhold[anti]*h+math.sqrt(refhold[anti])*sign*shock
        tr=d['transition'][t];state=d['decay'][tr]*state+d['root'][tr]@z[1:]+d['cross'][tr]*dw;zs=math.exp(-a_s*h)*zs+math.sqrt(-math.expm1(-2*a_s*h))*es;zf=math.exp(-af*h)*zf+math.sqrt(-math.expm1(-2*af*h))*ef
    expected=np.concatenate([np.exp(logs),np.repeat(np.exp(reflog)[:,None,None],2,axis=2)],axis=1);err=float(np.max(abs(got-expected)));assert err<1e-12
    import pandas as pd
    q=pd.read_csv(OUT/'VIX_membership.csv');f=pd.read_csv(OUT/'forward_targets.csv');draws=load_outer(128,20261028203);H={day:np.load(OUT/f'training_exposures/d{day}.npz')['H'] for day in draws};v=VIX(H,draws,x,q,f);y=np.array([.3,.2,.8,.3,.8,.4,.4]);r=v.predict(2,y);errors=[]
    for j in range(7):
        delta=np.zeros(7);delta[j]=1e-6;fd=(v.predict(2,y+delta)['forwards']-v.predict(2,y-delta)['forwards'])/(2e-6);error=float(np.max(abs(fd-r['fjac'][:,j])));assert error<1e-7,error;errors.append(error)
    # State the statistical budget transfer caused by product normalization.
    share=[]
    for t in [27/365,153/365,1.]:
        m,c=dated.gaussian_projection_law(t,initial);E=k.moments(m,c,p,m[1],scale,ALPHA,48)[0];oldmean=E[:,1]/np.exp(norm(t)[0]);den=(1-x[17])*oldmean+x[17]*E[:,0];share.append(dict(day=365*t,alpha=ALPHA.tolist(),mean_stress_share=(x[17]*E[:,0]/den).tolist()))
    dump(OUT/'implementation_audit.json',dict(all_checks_passed=True,maximum_direct_tensor_quadrature_scaled_error=max(momerr),scalar_stock_recursion_max_error=err,forward_derivative_max_absolute_error=max(errors),stress_mean_budget_share=share,scope='Numerical correspondence checks, not a continuous martingale theorem or finite-wing skew theorem'))
    print('ASYMMETRIC AUDIT PASS',max(momerr),err,max(errors),flush=True)

if __name__=='__main__':main()
