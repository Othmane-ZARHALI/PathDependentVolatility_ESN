"""Two positive mean-one variance contributions; analytic conditional pricing.

All mixtures are formed inside the future variance integral on each path.
They are never mixtures of completed option prices.
"""
import numpy as np
from .vix_only_curve_fit import VIXFit, ACTIVE, DELTA


def parameters(family,x):
    if family=='product':
        b1,d,w,lh=x[7:];return dict(b1=float(b1),b2=float(b1+(5-b1)*d),weight=float(w),half_life_days=float(np.exp(lh)))
    b,w,lh=x[7:];return dict(b1=0.,b2=float(b),weight=float(w),half_life_days=float(np.exp(lh)))


def response(family,old,T,times,eps,theta):
    """Conditional E[normalized readout at u|old state_T,Z_T] and derivatives.

    Derivatives are in (b1,d,w,log half-life), or (b,w,log half-life).
    eps is the standard-normal coordinate of Z_T under target Q, not proposal Q.
    """
    if family=='product':b1,d,w,lh=theta;b2=b1+(5-b1)*d
    elif family=='additive':b2,w,lh=theta
    else:raise ValueError(family)
    rate=np.log(2)*365/np.exp(lh);q=-np.expm1(-2*rate*T);lag=times-T
    c=np.exp(-rate*lag)[:,None]*np.sqrt(q);s=(rate*lag-rate*T*np.exp(-2*rate*T)/q)[:,None]
    def comp(b):
        z=b*c;value=np.exp(z*eps[None,:]-.5*z*z);diff=eps[None,:]-z
        return value,value*diff*c,value*diff*z*s
    L2,B2,H2=comp(b2)
    if family=='additive':return (1-w)*old+w*L2,[w*B2,L2-old,w*H2]
    L1,B1,H1=comp(b1)
    value=old*((1-w)*L1+w*L2)
    return value,[old*((1-w)*B1+w*B2*(1-d)),old*w*B2*(5-b1),old*(L2-L1),old*((1-w)*H1+w*H2)]


class MixedFit(VIXFit):
    def __init__(self,*args,family,**kwargs):
        super().__init__(*args,**kwargs);self.family=family
        nodes=args[0] if args else kwargs['nodes']
        for day,d in self.data.items():d['ratios']=nodes[f'd{day}_ratios'];d['tw']=nodes[f'd{day}_weights']

    def predict(self,x,detail=False):
        levels=self.p[7:-1].copy();levels[ACTIVE]=x[:7];npar=len(x)
        prices=[];calls=[];pj=[];fwd=[];fj=[];raw={};moments=[]
        for day,d in self.data.items():
            Q,D=response(self.family,d['ratios'],day/365,d['times'],d['a']['z_normals'],x[7:])
            W=Q*d['tw'][:,None];coeff=levels[d['seg']]-self.p[-1]
            integral=self.p[-1]*DELTA+coeff@W
            di=np.zeros((len(integral),npar))
            for j,seg in enumerate(ACTIVE):di[:,j]=W[d['seg']==seg].sum(0)
            for j,derivative in enumerate(D):di[:,7+j]=(coeff*d['tw'])@derivative
            V=100*np.sqrt(integral/DELTA);dv=50/np.sqrt(DELTA*integral)[:,None]*di
            a=d['a'];w=a['weights'];g=d['q'];K=g.strike.to_numpy();put=g.option_type.eq('put').to_numpy()
            yf=w*V-d['fc'];yc=w[:,None]*np.maximum(V[:,None]-K,0)-d['correction'];yn=yc-put[None,:]*(yf[:,None]-K)
            jac=np.einsum('nq,nj->qj',w[:,None]*(V[:,None]>K),dv)/len(V)
            f_jac=np.mean(w[:,None]*dv,axis=0);jac-=put[:,None]*f_jac
            prices.extend(yn.mean(0));calls.extend(yc.mean(0));pj.extend(jac);fwd.append(yf.mean());fj.append(f_jac)
            if detail:
                raw[day]=dict(V=V,weights=w,native=yn,call=yc,forward=yf,nets=a['nets'])
                lengths=np.maximum(0,np.minimum(self.ends,day/365+DELTA)-np.maximum(np.r_[0,self.ends[:-1]],day/365))
                block=np.array([np.mean((w*V**2)[a['nets']==j]) for j in np.unique(a['nets'])])
                tail=w*np.maximum(V-g.strike.max(),0);total=tail.sum();top=np.sort(tail)[::-1]
                moments.append(dict(day=int(day),exact_second_moment=float(10000/DELTA*(lengths@levels)),sample_second_moment=float(block.mean()),second_moment_se=float(block.std(ddof=1)/np.sqrt(len(block))) if len(block)>1 else None,mass=float(w.mean()),max_VIX=float(V.max()),last_strike=float(g.strike.max()),last_strike_raw_call=float(tail.mean()),last_strike_payoff_ESS=float(total**2/np.sum(tail**2)) if total else 0.,last_strike_top10_fraction=float(top[:10].sum()/total) if total else 0.))
        prices=np.array(prices);pj=np.array(pj);fwd=np.array(fwd);fj=np.array(fj);m=self.market
        dist=np.where(prices<m.price_bid_U,prices-m.price_bid_U,np.where(prices>m.price_ask_U,prices-m.price_ask_U,0.))
        weights=np.concatenate([np.full(len(g),1/(4*len(g))) for _,g in m.groupby('expiry_days',sort=True)])
        r=np.r_[np.sqrt(weights)*dist/m.price_scale_U,(fwd-self.forwards.forward_vix_points.to_numpy())/(.5*2)]
        outside=(prices<m.price_bid_U)|(prices>m.price_ask_U)
        jac=np.concatenate([np.sqrt(weights)[:,None]*outside.to_numpy()[:,None]*pj/m.price_scale_U.to_numpy()[:,None],fj/(.5*2)])
        return dict(residuals=r,jacobian=jac,loss=float(r@r),prices=prices,calls=np.array(calls),forwards=fwd,raw=raw,moments=moments)
