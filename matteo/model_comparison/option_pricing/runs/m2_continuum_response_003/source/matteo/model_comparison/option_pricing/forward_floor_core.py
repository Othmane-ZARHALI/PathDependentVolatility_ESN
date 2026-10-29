"""Positive floor and independently free xi levels, with analytic price/forward Jacobians."""
import numpy as np
from .mixed_response_core import MixedFit,response
from .vix_only_curve_fit import ACTIVE,DELTA

LOG_HALF_MAX=np.log(365.)


def unpack(x,p):
    sigma=.05+.07*x[10];floor=sigma*sigma
    levels=p[7:-1].copy();levels[ACTIVE]=floor+(.2-floor)*x[:7]
    return levels,floor,np.array([5*x[7],x[8],LOG_HALF_MAX*x[9]])


def pack(levels,b,w,half,floor):
    return np.r_[(np.asarray(levels)[ACTIVE]-floor)/(.2-floor),b/5,w,np.log(half)/LOG_HALF_MAX,(np.sqrt(floor)-.05)/.07]


class ForwardFloorFit(MixedFit):
    def __init__(self,*args,**kwargs):super().__init__(*args,family='additive',**kwargs)

    def predict(self,x,detail=False):
        levels,floor,theta=unpack(x,self.p);prices=[];calls=[];pj=[];cj=[];fwd=[];fj=[];raw={};moments=[]
        for day,d in self.data.items():
            Q,D=response('additive',d['ratios'],day/365,d['times'],d['a']['z_normals'],theta)
            W=Q*d['tw'][:,None];coeff=levels[d['seg']]-floor;integral=floor*DELTA+coeff@W
            expo=np.column_stack([W[d['seg']==seg].sum(0) for seg in ACTIVE]);di=np.zeros((len(integral),11));di[:,:7]=(.2-floor)*expo
            for j,scale in enumerate([5.,1.,LOG_HALF_MAX]):di[:,7+j]=scale*((coeff*d['tw'])@D[j])
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
