"""Continuous Beta measure of normalized geometric M2/stress responses.

R_t(u)=G_t**(1-u)/E[G_t**(1-u)] * exp(b*u*Z_t-.5*(b*u)**2*q(t)).
Shared old factors and one independent OU Z. No mixture of finished prices.
"""
import numpy as np
from scipy.special import roots_jacobi
from scipy.interpolate import CubicSpline
from .mixed_response_core import MixedFit
from .forward_floor_core import LOG_HALF_MAX,unpack,pack
from .vix_only_curve_fit import ACTIVE,DELTA

LOG_K_MIN=np.log(.01)
LOG_K_RANGE=np.log(30/.01)

def measure(x):
    return float(x[8]),float(np.exp(LOG_K_MIN+LOG_K_RANGE*x[11]))

def beta_rule(m,k,n):
    z,w=roots_jacobi(n,(1-m)*k-1,m*k-1)
    w=w/w.sum()
    assert np.all(w>0) and np.all((z>-1)&(z<1))
    return (z+1)/2,w

class ContinuumFit(MixedFit):
    def __init__(self,nodes,*args,measure_order=24,**kwargs):
        super().__init__(nodes,*args,family='additive',**kwargs)
        self.measure_order=measure_order
        for day,d in self.data.items():
            d['spline']=CubicSpline(nodes['u_grid'],nodes[f'd{day}_log_ratios'],axis=0,extrapolate=False)
            d['tw']=nodes[f'd{day}_weights']

    def response(self,day,x):
        d=self.data[day];m,k=measure(x);u,w=beta_rule(m,k,self.measure_order)
        # Differentiate only the small, deterministic quadrature rule numerically.
        # Price/forward derivatives then follow the analytic chain rule.
        rule_der=[]
        for col in [8,11]:
            h=1e-5;xp=x.copy();xm=x.copy();xp[col]+=h;xm[col]-=h
            up,wp=beta_rule(*measure(xp),self.measure_order);um,wm=beta_rule(*measure(xm),self.measure_order)
            rule_der.append(((up-um)/(2*h),(wp-wm)/(2*h)))
        T=day/365;rate=np.log(2)*365/np.exp(LOG_HALF_MAX*x[9]);q=-np.expm1(-2*rate*T)
        lag=d['times']-T;c=np.exp(-rate*lag)[:,None]*np.sqrt(q);eps=d['a']['z_normals'][None,:]
        dh=(rate*lag-rate*T*np.exp(-2*rate*T)/q)[:,None]
        Q=np.zeros_like(c*eps);D=[np.zeros_like(Q) for _ in range(4)]
        sp=d['spline'];lv=sp(u);ld=sp(u,1);b=5*x[7]
        for j,(uj,wj) in enumerate(zip(u,w)):
            z=b*uj*c;v=np.exp(lv[j]+z*eps-.5*z*z);diff=eps-z
            Q+=wj*v;D[0]+=wj*v*diff*uj*c*5
            D[2]+=wj*v*diff*z*dh*LOG_HALF_MAX
            du=ld[j]+diff*b*c
            for l,(un,wn) in enumerate(rule_der):D[[1,3][l]]+=v*(wn[j]+wj*du*un[j])
        return Q,D

    def predict(self,x,detail=False):
        levels,floor,theta=unpack(x,self.p);prices=[];calls=[];pj=[];cj=[];fwd=[];fj=[];raw={};moments=[]
        for day,d in self.data.items():
            Q,D=self.response(day,x)
            W=Q*d['tw'][:,None];coeff=levels[d['seg']]-floor;integral=floor*DELTA+coeff@W
            expo=np.column_stack([W[d['seg']==seg].sum(0) for seg in ACTIVE]);di=np.zeros((len(integral),12));di[:,:7]=(.2-floor)*expo
            for j,col in enumerate([7,8,9,11]):di[:,col]=(coeff*d['tw'])@D[j]
            sigma=.05+.07*x[10];di[:,10]=2*sigma*.07*(DELTA-expo@x[:7])
            V=100*np.sqrt(integral/DELTA);dv=50/np.sqrt(DELTA*integral)[:,None]*di
            a=d['a'];w=a['weights'];g=d['q'];K=g.strike.to_numpy();put=g.option_type.eq('put').to_numpy()
            yf=w*V-d['fc'];yc=w[:,None]*np.maximum(V[:,None]-K,0)-d['correction'];yn=yc-put[None,:]*(yf[:,None]-K)
            call_jac=np.einsum('nq,nj->qj',w[:,None]*(V[:,None]>K),dv)/len(V);f_jac=np.mean(w[:,None]*dv,axis=0);native_jac=call_jac-put[:,None]*f_jac
            prices.extend(yn.mean(0));calls.extend(yc.mean(0));pj.extend(native_jac);cj.extend(call_jac);fwd.append(yf.mean());fj.append(f_jac)
            if detail:
                raw[day]=dict(V=V,weights=w,native=yn,call=yc,forward=yf,nets=a['nets'])
                lengths=np.maximum(0,np.minimum(self.ends,day/365+DELTA)-np.maximum(np.r_[0,self.ends[:-1]],day/365));block=np.array([np.mean((w*V**2)[a['nets']==j]) for j in np.unique(a['nets'])])
                tail=w*np.maximum(V-K.max(),0);total=tail.sum();top=np.sort(tail)[::-1]
                moments.append(dict(day=int(day),exact_second_moment=float(10000/DELTA*(lengths@levels)),sample_second_moment=float(block.mean()),second_moment_se=float(block.std(ddof=1)/np.sqrt(len(block))) if len(block)>1 else None,mass=float(w.mean()),max_VIX=float(V.max()),last_strike=float(K.max()),last_strike_raw_call=float(tail.mean()),last_strike_payoff_ESS=float(total**2/np.sum(tail**2)) if total else 0.,last_strike_top10_fraction=float(top[:10].sum()/total) if total else 0.))
        prices=np.array(prices);pj=np.array(pj);fwd=np.array(fwd);fj=np.array(fj);m=self.market
        dist=np.where(prices<m.price_bid_U,prices-m.price_bid_U,np.where(prices>m.price_ask_U,prices-m.price_ask_U,0.));weights=np.concatenate([np.full(len(g),1/(4*len(g))) for _,g in m.groupby('expiry_days',sort=True)])
        r=np.r_[np.sqrt(weights)*dist/m.price_scale_U,(fwd-self.forwards.forward_vix_points.to_numpy())/(.5*2)];outside=(prices<m.price_bid_U)|(prices>m.price_ask_U)
        jac=np.concatenate([np.sqrt(weights)[:,None]*outside.to_numpy()[:,None]*pj/m.price_scale_U.to_numpy()[:,None],fj/(.5*2)])
        return dict(residuals=r,jacobian=jac,loss=float(r@r),option_loss=float(r[:-4]@r[:-4]),prices=prices,calls=np.array(calls),forwards=fwd,forward_jacobian=fj,native_jacobian=pj,call_jacobian=np.array(cj),raw=raw,moments=moments)
