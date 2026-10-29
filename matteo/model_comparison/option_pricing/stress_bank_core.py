"""Fixed geometric-speed stress bank; normalized positive variance mixture."""
import numpy as np
from .mixed_response_core import MixedFit
from .vix_only_curve_fit import ACTIVE,DELTA

HALF_LIVES=np.array([7.,28.,112.,448.])


def unpack(x,p):
    floor=(.05+.07*x[7])**2
    levels=p[7:-1].copy();levels[ACTIVE]=floor+(.2-floor)*x[:7]
    return levels,floor,None


def weights(x):
    s=x[12:16];alpha=np.array([s[j]*np.prod(1-s[:j]) for j in range(4)])
    jac=np.zeros((4,4))
    for j in range(4):
        jac[j,j]=np.prod(1-s[:j])
        for k in range(j):jac[j,k]=-s[j]*np.prod([1-s[i] for i in range(j) if i!=k])
    return alpha,float(np.prod(1-s)),jac


def stick_units(alpha):
    remaining=1.;s=[]
    for a in alpha:
        s.append(a/remaining if remaining>0 else 0.);remaining-=a
    return np.array(s)


def metadata(x,p):
    levels,floor,_=unpack(x,p);alpha,m2,_=weights(x);total=alpha.sum()
    return dict(half_lives_days=HALF_LIVES.tolist(),amplitudes=(5*x[8:12]).tolist(),
        stress_weights=alpha.tolist(),M2_weight=m2,stress_weight=float(total),
        conditional_stress_weights=(alpha/total).tolist() if total>0 else [0.]*4,
        floor_variance=float(floor),active_levels=levels[ACTIVE].tolist())


def bank_response(d,T,x,reference=False):
    old=d['ratios'];alpha,m2,jac=weights(x);value=m2*old
    derivatives=[];weight_derivatives=[]
    for j,h in enumerate(HALF_LIVES):
        if reference and j<3:
            L=np.ones_like(old);db=np.zeros_like(old)
        else:
            c=d['bank_loadings'][j];b=5*x[8+j];z=b*c[:,None]
            eps=d['a']['bank_normals'][:,j]
            L=np.exp(z*eps[None,:]-.5*z*z)
            db=5*alpha[j]*L*(eps[None,:]-z)*c[:,None]
        value=value+alpha[j]*L;derivatives.append(db);weight_derivatives.append(L-old)
    return value,derivatives+[sum(jac[j,k]*weight_derivatives[j] for j in range(4)) for k in range(4)]


class StressBankFit(MixedFit):
    def __init__(self,*args,reference=False,**kwargs):
        super().__init__(*args,family='additive',**kwargs);self.reference=reference
        for day,d in self.data.items():
            T=day/365;rate=np.log(2)*365/HALF_LIVES
            d['bank_loadings']=np.exp(-rate[:,None]*(d['times']-T)[None,:])*np.sqrt(-np.expm1(-2*rate*T))[:,None]

    def predict(self,x,detail=False):
        levels,floor,theta=unpack(x,self.p);prices=[];calls=[];pj=[];cj=[];fwd=[];fj=[];raw={};moments=[]
        for day,d in self.data.items():
            Q,D=bank_response(d,day/365,x,self.reference)
            W=Q*d['tw'][:,None];coeff=levels[d['seg']]-floor;integral=floor*DELTA+coeff@W
            expo=np.column_stack([W[d['seg']==seg].sum(0) for seg in ACTIVE]);di=np.zeros((len(integral),16));di[:,:7]=(.2-floor)*expo
            for j,col in enumerate(range(8,16)):di[:,col]=(coeff*d['tw'])@D[j]
            sigma=.05+.07*x[7];di[:,7]=2*sigma*.07*(DELTA-expo@x[:7])
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
