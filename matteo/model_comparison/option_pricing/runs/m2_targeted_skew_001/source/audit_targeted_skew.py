"""Scalar coupled-primitive replay, gradient and conditional-mean audits."""
import json,math
import numpy as np
from .targeted_skew import OUT,Kernel,ALPHA,arr,dump
from .asymmetric_readout import PRIOR
from .audit_leverage_scan import gaussian
from .protected_curve_check import context

def main():
    k=Kernel(OUT/'targeted.dylib');config=json.loads((OUT/'candidates.json').read_text());models=arr(config['models']);alpha=arr(config['alpha']);mode=np.array(config['mode'],np.int32);d=dict(np.load(OUT/'augmented_grid.npz'));data=np.load(OUT/'component_means.npz');M=data['moments'];mom=arr(M[:,config['mean_indices']]);scale=json.loads((OUT/'design.json').read_text())['scale'];p=arr(json.loads((PRIOR/'run.json').read_text())['selected']);na=len(alpha);nb=d['diffusion'].shape[1]+2;oldb=nb-2;n=191;small=dict(d)
    for key in ['time','dt','seg','coarse','reference','expiry','transition','mean','lognorm']:small[key]=d[key][:n].copy()
    small['expiries']=1;small['expiry'][:]=-1;small['expiry'][-1]=0;seed=20261028651;got=k.stocks(small,1,p,models,alpha,mode,arr(mom[:n]),scale,seed)[:,0];s=np.zeros(24);zs=zf=0.;logs=np.zeros((2,na,2));hold=np.zeros((2,na));helddir=np.zeros((2,nb));rl=np.zeros(2);rv=np.zeros(2);x=models[0];floor=(.05+.07*x[19])**2;L=math.log(365);af=math.log(2)*365/math.exp(L*x[21]);a_s=math.log(2)*365/math.exp(L*x[21]+(L-L*x[21])*x[18]);bs=5*x[16];bf=5*x[20];normerrors=[];maxcovslack=0.;gradient_errors=[]
    def coeff(a,h):
        z=a*h;rem=h*(z*z/12-z**4/120+17*z**6/20160) if z<.01 else h*(1-math.tanh(z/2)/(z/2));return math.sqrt(h-rem),math.sqrt(rem)
    def readout(features,t):
        C,Y,F=features;eta=p[0]*C+p[1]*(Y-p[2])**2+p[3]*(F-p[4])**2;return math.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6])-d['lognorm'][t])
    for t in range(n):
        h=d['dt'][t];time=d['time'][t];tr=d['transition'][t];counter=t*27;z=np.array([gaussian(seed,counter+j) for j in range(25)]);dw=math.sqrt(h)*z[0];es=gaussian(seed,counter+25);ef=gaussian(seed,counter+26);cs,rs=coeff(a_s,h);cf,rf=coeff(af,h);dbs=cs*es+rs*gaussian(seed^0x2ad82d187dccb289,2*t);dbf=cf*ef+rf*gaussian(seed^0x2ad82d187dccb289,2*t+1);extra=np.array([gaussian(seed^0x78d904ef6038fb35,t*(oldb-1)+j) for j in range(oldb-1)]);db=np.r_[dw,d['bmean'][tr]@z[1:]+d['broot'][tr]@extra,dbs,dbf];fixed=.6*dw-.8*(math.sqrt(1-x[22])*dbs+math.sqrt(x[22])*dbf);pr=d['projection']@s
        for anti,sign in enumerate([1,-1]):
            features=d['mean'][t]+sign*pr;C,Y,F=features;eta=p[0]*C+p[1]*(Y-p[2])**2+p[3]*(F-p[4])**2;R=readout(features,t);As=math.exp(sign*bs*zs-.5*bs*bs*(-math.expm1(-2*a_s*time)));Af=math.exp(sign*bf*zf-.5*bf*bf*(-math.expm1(-2*af*time)));A=(1-x[22])*As+x[22]*Af;H=(1-x[17])*R+x[17]*A;sig=1/(1+math.exp(p[6]*(eta-p[5])));g=2*R*sig*np.array([p[0],2*p[1]*(Y-p[2]),2*p[3]*(F-p[4])]);kk=np.r_[(1-x[17])*g@d['diffusion'],x[17]*(1-x[22])*bs*math.sqrt(2*a_s)*As,x[17]*x[22]*bf*math.sqrt(2*af)*Af];direction=-kk/np.linalg.norm(kk);normerrors.append(abs(direction@direction-1));a0=np.zeros(nb);a0[0]=.6;a0[-2]=-.8*math.sqrt(1-x[22]);a0[-1]=-.8*math.sqrt(x[22]);maxcovslack=max(maxcovslack,float(direction@kk-a0@kk));assert direction@kk<=a0@kk+1e-12
            if t in [0,90,190]:
                for j in range(3):
                    perturb=np.zeros(3);perturb[j]=1e-5;fd=(readout(features+perturb,t)-readout(features-perturb,t))/(2e-5);gradient_errors.append(abs(fd-g[j])/(1+abs(fd)))
            if d['coarse'][t]:helddir[anti]=direction
            adaptive=np.array([helddir[anti]@db,direction@db]);l=np.logaddexp(0,-sign*pr[1]/scale)
            for a in range(na):
                response=H
                if mode[a]==2:response=H*math.exp(alpha[a]*l)*((1-x[17])+x[17]*A)/((1-x[17])*mom[t,a,1]+x[17]*A*mom[t,a,0])
                v=floor+(.2-floor)*models[a,d['seg'][t]]*response
                if d['coarse'][t]:hold[anti,a]=v
                for grid,vv in enumerate([hold[anti,a],v]):shock=adaptive[grid] if mode[a]==1 else fixed;logs[anti,a,grid]+=-.5*vv*h+math.sqrt(vv)*sign*shock
            if d['reference'][t]:rv[anti]=floor+(.2-floor)*x[d['seg'][t]]*((1-x[17])+x[17]*A)
            rl[anti]+=-.5*rv[anti]*h+math.sqrt(rv[anti])*sign*fixed
        s=d['decay'][tr]*s+d['root'][tr]@z[1:]+d['cross'][tr]*dw;zs=math.exp(-a_s*h)*zs+math.sqrt(-math.expm1(-2*a_s*h))*es;zf=math.exp(-af*h)*zf+math.sqrt(-math.expm1(-2*af*h))*ef
    expected=np.concatenate([np.exp(logs),np.repeat(np.exp(rl)[:,None,None],2,axis=2)],axis=1);error=float(np.max(abs(got-expected)));assert error<1e-12;assert max(gradient_errors)<1e-8
    # Conditional stress-state mean, integrating the complete old-bank product.
    prior,dated,norm,_=context();initial=arr(prior['initial_state']);meanerrs=[];stress_moment_errors=[]
    for t in [27/365,153/365,1.]:
        mean,cov=dated.gaussian_projection_law(t,initial);E=k.moments(mean,cov,p,mean[1],scale,ALPHA,48)[0];E[:,1]/=np.exp(norm(t)[0]);E[0]=1.
        for A in [1e-6,.1,1,10,1e4]:
            numerical=((1-x[17])*E[:,1]+x[17]*A*E[:,0])*((1-x[17])+x[17]*A)/((1-x[17])*E[:,1]+x[17]*A*E[:,0]);desired=(1-x[17])+x[17]*A;meanerrs.append(float(np.max(abs(numerical/desired-1))))
        # Independent vectorized stress GH oracle for arbitrary future states.
        from numpy.polynomial.hermite import hermgauss
        order=24;z,w=hermgauss(order);z=np.sqrt(2)*z;w=w/np.sqrt(np.pi);means=np.array([[-2,1],[0,0],[2,4]]);sd=np.array([.3,.6]);b=np.array([bs,bf]);q=np.array([.8,.9]);got=k.stress(means,sd,b,q,x[22],x[17],E,order);oracle=np.zeros_like(got)
        for j,mm in enumerate(means):
            A=(1-x[22])*np.exp(bs*(mm[0]+sd[0]*z[:,None])-.5*bs*bs*q[0])+x[22]*np.exp(bf*(mm[1]+sd[1]*z[None,:])-.5*bf*bf*q[1]);wt=w[:,None]*w[None,:]
            for a in range(len(ALPHA)):
                ratio=((1-x[17])+x[17]*A)/((1-x[17])*E[a,1]+x[17]*A*E[a,0]);oracle[j,a]=[((1-x[17])*ratio*wt).sum(),(x[17]*A*ratio*wt).sum()]
        stress_moment_errors.append(float(np.max(abs(got-oracle)/(1+abs(oracle)))))
    assert max(meanerrs)<1e-12 and max(stress_moment_errors)<1e-12
    dump(OUT/'implementation_audit.json',dict(all_checks_passed=True,scalar_stock_recursion_max_error=error,unit_stock_loading_max_error=max(normerrors),readout_gradient_max_scaled_error=max(gradient_errors),conditional_stress_mean_max_relative_error=max(meanerrs),stress_quadrature_max_scaled_error=max(stress_moment_errors),instantaneous_covariance_not_greater_than_baseline=True,scope='Numerical primitive/recursion/gradient/conditional-expectation checks; no finite-wing-skew or continuous-stock-martingale theorem'))
    print('TARGETED AUDIT PASS',error,max(gradient_errors),max(stress_moment_errors),flush=True)

if __name__=='__main__':main()
