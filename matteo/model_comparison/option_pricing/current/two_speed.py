"""Conditional VIX pricing for the exploratory two-speed reference.

Extracted from the archived calibration implementation; inputs and controls
are supplied by the caller. See the research record for calibration limits.
"""
import numpy as np

ACTIVE = np.array([1, 2, 3, 4, 7, 8, 10])
DELTA = 30 / 365
LOG_HALF_MAX = np.log(365.)

def canonical_market(quotes):
    q = quotes.copy()
    adjustment = np.where(q.option_type.to_numpy() == 'put',
                          q.forward.to_numpy()-q.strike.to_numpy(), 0.)
    for column in ['price_bid_U', 'price_mid_U', 'price_ask_U']:
        q[column] = q[column].to_numpy()+adjustment
    q['original_option_type'] = q.option_type
    q['option_type'] = 'call'
    return q

def gaussian_controls(z, weights, dimensions=6):
    z = np.asarray(z)[:, :dimensions]; w = np.asarray(weights)
    # z has target N(0,I), even though proposal draws are shifted.
    return np.column_stack([w-1, w[:, None]*z, w[:, None]*(z*z-1)])

class VIXFit:
    def __init__(self, nodes, draws, market, forwards, controls, p, ends):
            self.market=market.reset_index(drop=True);self.calls=canonical_market(self.market)
            self.forwards=forwards;self.p=p;self.ends=ends;self.data={}
            for day,g in self.market.groupby('expiry_days',sort=True):
                a=draws[day];times=nodes[f'd{day}_times'];tw=nodes[f'd{day}_weights']
                seg=np.minimum(np.searchsorted(ends,times,side='right'),len(ends)-1)
                c=gaussian_controls(a['normals'],a['weights']);fixed=controls[str(day)]
                b=np.stack([np.interp(np.log(g.strike),np.log(fixed['anchors']),row) for row in fixed['beta']])
                self.data[day]=dict(a=a, times=times, seg=seg,
                    weighted_ratios=nodes[f'd{day}_ratios']*tw[:,None],
                    correction=c@b, fc=c@np.array(fixed['forward_beta']), q=g)

class MixedFit(VIXFit):
    def __init__(self,*args,family,**kwargs):
            super().__init__(*args,**kwargs);self.family=family
            nodes=args[0] if args else kwargs['nodes']
            for day,d in self.data.items():d['ratios']=nodes[f'd{day}_ratios'];d['tw']=nodes[f'd{day}_weights']

def unpack(x,p):
    sigma=.05+.07*x[10];floor=sigma*sigma
    levels=p[7:-1].copy();levels[ACTIVE]=floor+(.2-floor)*x[:7]
    return levels,floor,np.array([5*x[7],x[8],LOG_HALF_MAX*x[9]])

def pack(levels,b,w,half,floor):
    return np.r_[(np.asarray(levels)[ACTIVE]-floor)/(.2-floor),b/5,w,np.log(half)/LOG_HALF_MAX,(np.sqrt(floor)-.05)/.07]

def exponential(T,times,eps,b,lh):
    a=np.log(2)*365/np.exp(lh);q=-np.expm1(-2*a*T);lag=times-T
    c=(np.exp(-a*lag)*np.sqrt(q))[:,None];z=b*c
    value=np.exp(z*eps[None,:]-.5*z*z)
    db=value*(eps[None,:]-z)*c
    dh=value*(eps[None,:]-z)*z*(a*lag-a*T*np.exp(-2*a*T)/q)[:,None]
    return value,db,dh

def two_speed_response(d,T,x):
    # Slow half-life h_s=h_f*(365/h_f)^x9 ensures h_s>=h_f.
    loghf=LOG_HALF_MAX*x[12]
    loghs=loghf+(LOG_HALF_MAX-loghf)*x[9]
    Ls,Bs,Hs=exponential(T,d['times'],d['a']['z_normals'],5*x[7],loghs)
    Lf,Bf,Hf=exponential(T,d['times'],d['a']['fast_normals'],5*x[11],loghf)
    w=x[8];r=x[13];stress=(1-r)*Ls+r*Lf;old=d['ratios']
    return (1-w)*old+w*stress,[
        5*w*(1-r)*Bs,stress-old,w*(1-r)*Hs*(LOG_HALF_MAX-loghf),
        5*w*r*Bf,LOG_HALF_MAX*w*(r*Hf+(1-r)*(1-x[9])*Hs),w*(Lf-Ls)]

def metadata(x,p):
    levels,floor,_=unpack(x,p);hf=np.exp(LOG_HALF_MAX*x[12]);hs=hf*(365/hf)**x[9]
    return dict(b_slow=float(5*x[7]),b_fast=float(5*x[11]),stress_weight=float(x[8]),
        fast_fraction=float(x[13]),slow_half_days=float(hs),fast_half_days=float(hf),
        floor_variance=float(floor),active_levels=levels[ACTIVE].tolist())

def boundary_x(x):
    y=x[:11].copy();y[9]=x[12]+(1-x[12])*x[9];return y

class TwoSpeedFit(MixedFit):
    def __init__(self,*args,**kwargs):super().__init__(*args,family='additive',**kwargs)

    def predict(self,x,detail=False):
        levels,floor,theta=unpack(x,self.p);prices=[];calls=[];pj=[];cj=[];fwd=[];fj=[];raw={};moments=[]
        for day,d in self.data.items():
            Q,D=two_speed_response(d,day/365,x)
            W=Q*d['tw'][:,None];coeff=levels[d['seg']]-floor;integral=floor*DELTA+coeff@W
            expo=np.column_stack([W[d['seg']==seg].sum(0) for seg in ACTIVE]);di=np.zeros((len(integral),14));di[:,:7]=(.2-floor)*expo
            for j,col in enumerate([7,8,9,11,12,13]):di[:,col]=(coeff*d['tw'])@D[j]
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
