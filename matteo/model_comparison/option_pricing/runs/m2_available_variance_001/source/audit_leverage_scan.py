"""Independent scalar transition/sufficient-statistic oracle and covariance checks."""
import json
import math
import numpy as np
from scipy.linalg import eigvalsh
from .leverage_scan import Kernel,OUT,PARENT,PRIOR,arr,dump

MASK=(1<<64)-1
def hash64(x):
    x&=MASK;x=((x^(x>>30))*0xbf58476d1ce4e5b9)&MASK;x=((x^(x>>27))*0x94d049bb133111eb)&MASK
    return x^(x>>31)
def gaussian(seed,i):
    a=hash64(seed+2*i*0x9e3779b97f4a7c15);b=hash64(seed+(2*i+1)*0x9e3779b97f4a7c15);u=((a>>11)+.5)*2**-53;v=((b>>11)+.5)*2**-53
    return math.sqrt(-2*math.log(u))*math.cos(2*math.pi*v)

def main():
    d=dict(np.load(PARENT/'fresh_check/grid_and_transitions.npz'));x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);p=arr(json.loads((PRIOR/'run.json').read_text())['selected']);n=191;small=dict(d)
    for key in ['time','dt','seg','coarse','reference','expiry','transition','mean','lognorm']:small[key]=d[key][:n].copy()
    small['expiries']=1;small['expiry'][:]=-1;small['expiry'][-1]=0;seed=202610281931;k=Kernel(OUT/'leverage.dylib');got=k.stats(small,1,p,x,seed)[:,0];state=np.zeros(24);zs=zf=0.;sums=np.zeros((2,3,4));hold=np.zeros((2,3));floor=(.05+.07*x[19])**2;af=math.log(2)*365/math.exp(math.log(365)*x[21]);a_s=math.log(2)*365/math.exp(math.log(365)*x[21]+(math.log(365)-math.log(365)*x[21])*x[18]);bs=5*x[16];bf=5*x[20]
    def coeff(a,h):
        z=a*h;rem=h*(z*z/12-z**4/120+17*z**6/20160) if z<.01 else h*(1-math.tanh(z/2)/(z/2));return math.sqrt(h-rem),math.sqrt(rem)
    maxcov=0.
    for t in range(n):
        h=d['dt'][t];time=d['time'][t];counter=t*27;z=np.array([gaussian(seed,counter+j) for j in range(25)]);dw=math.sqrt(h)*z[0];es=gaussian(seed,counter+25);ef=gaussian(seed,counter+26);cs,rs=coeff(a_s,h);cf,rf=coeff(af,h);dbs=cs*es+rs*gaussian(seed^0x2ad82d187dccb289,2*t);dbf=cf*ef+rf*gaussian(seed^0x2ad82d187dccb289,2*t+1);proj=d['projection']@state
        for anti,sign in enumerate([1,-1]):
            C,Y,F=d['mean'][t]+sign*proj;eta=p[0]*C+p[1]*(Y-p[2])**2+p[3]*(F-p[4])**2;old=math.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6])-d['lognorm'][t]);stress=(1-x[22])*math.exp(sign*bs*zs-.5*bs*bs*(-math.expm1(-2*a_s*time)))+x[22]*math.exp(sign*bf*zf-.5*bf*bf*(-math.expm1(-2*af*time)));v=floor+(.2-floor)*x[d['seg'][t]]*((1-x[17])*old+x[17]*stress)
            if d['coarse'][t]:hold[anti,0]=v
            hold[anti,1]=v
            if d['reference'][t]:hold[anti,2]=floor+(.2-floor)*x[d['seg'][t]]*((1-x[17])+x[17]*stress)
            for g in range(3):sums[anti,g]+=np.r_[hold[anti,g]*h,math.sqrt(hold[anti,g])*sign*np.array([dw,dbs,dbf])]
        tr=d['transition'][t];state=d['decay'][tr]*state+d['root'][tr]@z[1:]+d['cross'][tr]*dw;zs=math.exp(-a_s*h)*zs+math.sqrt(-math.expm1(-2*a_s*h))*es;zf=math.exp(-af*h)*zf+math.sqrt(-math.expm1(-2*af*h))*ef
    err=float(np.max(abs(got-sums)));assert err<1e-12
    # Analytic joint increment covariance, including very small time steps.
    rows=[]
    for a in [a_s,af]:
        for h in np.r_[np.unique(d['dt']),1e-12,1e-9,.1,1.]:
            c,r=coeff(a,h);ns=math.sqrt(-math.expm1(-2*a*h));desired=math.sqrt(2*a)*(-math.expm1(-a*h))/a;error=abs(c*ns-desired)/math.sqrt(h*(-math.expm1(-2*a*h)));maxcov=max(maxcov,error);assert abs((c*c+r*r)/h-1)<1e-14;assert error<1e-12
    # Full retained OU + W0 + stress Brownian covariance PSD, each rotation.
    min_eig=1.
    for theta in [0,.2,.4,.6,.8,.95]:
        for cs,cf in [(math.sqrt(1-x[22]),math.sqrt(x[22])),(1,0),(0,1)]:
            v=np.array([math.sqrt(1-theta*theta),-theta*cs,-theta*cf]);assert abs(v@v-1)<1e-14
            joint=np.eye(4);joint[0,1:]=v;joint[1:,0]=v;min_eig=min(min_eig,float(eigvalsh(joint)[0]));assert eigvalsh(joint)[0]>-1e-12
    # Independent Python reconstruction of cheap reference's exact OU/B increments.
    ref=k.reference(d,1,x,seed);idx=np.flatnonzero(d['reference']);grid=np.r_[d['time'][idx],d['time'][-1]+d['dt'][-1]];out=np.zeros_like(ref);sums=np.zeros((2,4));zs=zf=0.;ends={round(float(d['time'][t]+d['dt'][t]),12):int(e) for t,e in enumerate(d['expiry']) if e>=0}
    for t,h in enumerate(np.diff(grid)):
        time=grid[t];es=gaussian(seed,4*t);ef=gaussian(seed,4*t+1);cs,rs=coeff(a_s,h);cf,rf=coeff(af,h);dbs=cs*es+rs*gaussian(seed,4*t+2);dbf=cf*ef+rf*gaussian(seed,4*t+3)
        for anti,sign in enumerate([1,-1]):
            stress=(1-x[22])*math.exp(sign*bs*zs-.5*bs*bs*(-math.expm1(-2*a_s*time)))+x[22]*math.exp(sign*bf*zf-.5*bf*bf*(-math.expm1(-2*af*time)));v=floor+(.2-floor)*x[d['seg'][idx[t]]]*((1-x[17])+x[17]*stress);sums[anti]+=np.array([v*h,0,math.sqrt(v)*sign*dbs,math.sqrt(v)*sign*dbf]);e=ends.get(round(float(grid[t+1]),12))
            if e is not None:out[anti,e]=sums[anti]
        zs=math.exp(-a_s*h)*zs+math.sqrt(-math.expm1(-2*a_s*h))*es;zf=math.exp(-af*h)*zf+math.sqrt(-math.expm1(-2*af*h))*ef
    erref=float(np.max(abs(out-ref)));assert erref<1e-11
    dump(OUT/'implementation_audit.json',dict(scalar_target_statistics_max_error=err,scalar_reference_statistics_max_error=erref,maximum_standardized_stress_increment_covariance_error=maxcov,rotation_minimum_covariance_eigenvalue=min_eig,positive_spike_shift=float(p[4]),all_checks_passed=True,scope='Independent scalar numerical oracles and analytic finite-step covariance checks; not a proof of continuous-time martingality or sampling convergence'))
    print('AUDIT PASS',err,erref,maxcov,flush=True)

if __name__=='__main__':main()
