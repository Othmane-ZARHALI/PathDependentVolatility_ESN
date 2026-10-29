"""23-parameter curve/stress joint mapping and adapted SPX path pricing."""
import ctypes as ct
import numpy as np
from .vix_only_curve_fit import ACTIVE
from .two_speed_core import metadata as vix_metadata
from .parity_adaptive import reference_calls
from .full_readout_fit import arr


def vix_x(x):return np.r_[np.asarray(x)[ACTIVE],x[16:23]]

def vix_jac(j):
    out=np.zeros((len(j),23));out[:,ACTIVE]=j[:,:7];out[:,16:]=j[:,7:];return out


def physical(x):
    floor=(.05+.07*x[19])**2;levels=floor+(.2-floor)*x[:16];return levels,floor


def metadata(x,p):
    d=vix_metadata(vix_x(x),p);levels,floor=physical(x);d['all_curve_levels']=levels.tolist();return d


class StockKernel:
    def __init__(self,path):
        self.lib=ct.CDLL(str(path));P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS');self.lib.two_speed_stock.argtypes=[ct.c_int,ct.c_int,ct.c_int,P,P,P,P,P,I,P,P]
    def stock(self,d,x,derivatives=True):
        n=d['old'].shape[1];out=np.zeros((n,24));self.lib.two_speed_stock(len(d['dt']),n,int(derivatives),d['old'],d['equity'],d['eps'],d['time'],d['dt'],d['seg'],arr(x),out);return out


class SPXFit:
    def __init__(self,data,quotes,kernel,reference_levels):
        self.data=data;self.market=quotes;self.k=kernel;self.reflevels=np.array(reference_levels)
    def predict(self,x,detail=False,derivatives=True):
        prices=[];jac=[];blocks=[];raw={};moments=[]
        for day,g in self.market.groupby('expiry_days',sort=True):
            rows=[];drows=[];weights=[];stocks=[];refs=[]
            for d in self.data[day]:
                out=self.k.stock(d,x,derivatives);S=out[:,0];ds=out[:,1:];ref=np.exp(d['refbrownian']@np.sqrt(self.reflevels)-.5*d['refvariance']);dref=np.zeros_like(ds)
                y,dy=reference_calls(S,ds,ref,dref,d['weights'],g.K_over_F.to_numpy(),d['refvariance'],np.zeros(23))
                F=g.forward.to_numpy();y=y*F;dy=dy*F[None,:,None];puts=g.option_type.eq('put').to_numpy();y=y-puts[None,:]*(F-g.strike.to_numpy())[None,:]
                rows.append(y);drows.append(dy);weights.append(d['weights']);stocks.append(S);refs.append(ref)
            y=np.concatenate(rows);dy=np.concatenate(drows);w=np.concatenate(weights);S=np.concatenate(stocks)
            prices.extend(y.mean(0));jac.extend(dy.mean(0));blocks.append(y.reshape(8,-1,len(g)).mean(1))
            moments.append(dict(day=int(day),mass=float(w.mean()),stock_mean=float(np.mean(w*S)),stock_mean_SE=float(np.std(w*S,ddof=1)/np.sqrt(len(w))),minimum_stock=float(S.min()),maximum_stock=float(S.max())))
            if detail:raw[day]=dict(values=y,stock=S,weights=w,reference=np.concatenate(refs))
        price=np.array(prices);pj=np.array(jac);m=self.market;dist=np.where(price<m.price_bid_U,price-m.price_bid_U,np.where(price>m.price_ask_U,price-m.price_ask_U,0));wt=np.concatenate([np.full(len(g),1/(7*len(g))) for _,g in m.groupby('expiry_days',sort=True)]);r=np.sqrt(wt)*dist/m.price_scale_U
        outside=(price<m.price_bid_U)|(price>m.price_ask_U);J=np.sqrt(wt)[:,None]*outside.to_numpy()[:,None]*pj/m.price_scale_U.to_numpy()[:,None]
        return dict(prices=price,price_jacobian=pj,residuals=np.asarray(r),jacobian=J,loss=float(r@r),blocks=np.concatenate(blocks,axis=1),raw=raw,moments=moments)


class JointFit:
    def __init__(self,spx,vix):self.spx=spx;self.vix=vix
    def predict(self,x,detail=False):
        s=self.spx.predict(x,detail);v=self.vix.predict(vix_x(x),detail);vr=v['residuals'][:-4];vj=vix_jac(v['jacobian'][:-4]);r=np.r_[s['residuals'],vr];j=np.r_[s['jacobian'],vj]
        return dict(loss=float(r@r),spx_loss=s['loss'],vix_loss=v['option_loss'],residuals=r,jacobian=j,forwards=v['forwards'],forward_jacobian=vix_jac(v['forward_jacobian']),spx=s,vix=v)
