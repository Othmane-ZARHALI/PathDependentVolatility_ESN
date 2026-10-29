"""Four broad VIX-window adjustments with exact frozen off-window curve."""
import numpy as np
from .plot_maturity_curve import iv

WINDOWS=np.array([[27.,57.],[55.,85.],[118.,148.],[153.,183.]])

def window_basis(ends_days):
    ends=np.asarray(ends_days);left=np.r_[0.,ends[:-1]];mid=(left+ends)/2
    for edge in np.unique(WINDOWS):
        if not np.any(ends==edge):raise ValueError('Window edge missing from fixed curve')
    B=((mid[:,None]>WINDOWS[:,0])&(mid[:,None]<WINDOWS[:,1])).astype(float)
    counts=B.sum(1);B[counts>0]/=counts[counts>0,None]
    return B

def adjusted(base,s,B):
    s=np.asarray(s,float)
    if s.shape!=(4,) or not np.isfinite(s).all() or np.max(abs(s))>.15000000001:raise ValueError('Window adjustment outside declared range')
    x=base.copy();x[7:-1]=base[7:-1]+(base[7:-1]-base[-1])*(B@s)
    return x

def curve_jac(base,B):
    D=np.zeros((len(base),4));D[7:-1]=(base[7:-1]-base[-1])[:,None]*B;return D

def prediction(e,p,D):
    v=e.vix_prediction(p);call=v['prices'];J=v['jacobian']@D;forward=v['forwards'];FJ=v['forward_jacobian']@D
    native=call.copy();NJ=J.copy();index=0
    for j,(_,g) in enumerate(e.vix.groupby('expiry_days',sort=True)):
        n=len(g);put=g.option_type.eq('put').to_numpy();native[index:index+n]-=put*(forward[j]-g.strike.to_numpy());NJ[index:index+n]-=put[:,None]*FJ[j];index+=n
    def objective(loss,price,jac):
        q=loss.quotes;dist=np.where(price<q.price_bid_U,price-q.price_bid_U,np.where(price>q.price_ask_U,price-q.price_ask_U,0.));scaled=dist/q.price_scale_U.to_numpy()
        fr=(forward-e.forwards.forward_vix_points.to_numpy())/.5
        val=float(np.sum(loss.weights*scaled**2)+np.mean(fr**2))
        grad=2*np.einsum('q,qj->j',loss.weights*dist/q.price_scale_U.to_numpy()**2,jac)+2*np.mean(fr[:,None]*FJ/.5,axis=0)
        return val,grad
    native_loss,ng=objective(e.otm_loss,native,NJ);common_loss,cg=objective(e.loss,call,J)
    return dict(native_loss=native_loss,native_gradient=ng,common_loss=common_loss,common_gradient=cg,forward=forward,forward_jac=FJ,call=call,native=native)

def shape_metrics(quotes,prices):
    rows=[]
    for day,g in quotes.groupby('expiry_days',sort=True):
        y=iv(g,np.asarray(prices)[g.index]);market=iv(g,g.price_mid_U);xx=np.log(g.K_over_F.to_numpy())
        row=dict(days=int(day),unavailable=int((~np.isfinite(y)).sum()))
        for region,mask in [('ATM',abs(xx)<=.2),('wing',xx>np.log(2))]:
            valid=mask&np.isfinite(y)&np.isfinite(market);row[region+'_quotes']=int(mask.sum());row[region+'_available']=int(valid.sum())
            row[region+'_MAE']=float(np.mean(abs(y[valid]-market[valid]))) if valid.any() else None
            row[region+'_IV_level']=float(y[valid].mean()) if valid.any() else None
            if valid.sum()>=3:
                row[region+'_slope']=float(np.polyfit(xx[valid],y[valid],1)[0]);row[region+'_curvature']=float(2*np.polyfit(xx[valid],y[valid],2)[0])
            else:row[region+'_slope']=row[region+'_curvature']=None
        rows.append(row)
    return dict(rows=rows,ATM_MAE=float(np.mean([r['ATM_MAE'] for r in rows])),wing_MAE=float(np.mean([r['wing_MAE'] for r in rows])),wing_complete=all(r['wing_available']==r['wing_quotes'] for r in rows))
