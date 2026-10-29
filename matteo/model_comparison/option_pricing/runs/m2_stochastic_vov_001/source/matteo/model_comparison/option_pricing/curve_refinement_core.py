"""Dimension-generic curve pricing; fixed state and primitive sampling law."""
import numpy as np

from .joint_curve_core import path_payoffs, ExploratoryJointLoss


class FlexibleCurvePathBatch:
    def __init__(self,ratio,equity,dt,segment,dimension,floor=1e-4):
        self.ratio,self.equity,self.dt,self.segment=ratio,equity,dt,segment
        self.dimension,self.floor=dimension,floor

    def terminal(self,levels,derivative=True,chunk=2048):
        levels=np.asarray(levels,float)
        if levels.shape!=(self.dimension,) or np.any(levels<=self.floor) or not np.isfinite(levels).all():
            raise ValueError('Invalid curve levels')
        logstock=np.zeros(self.equity.shape[1]);gradient=np.zeros((len(logstock),self.dimension))
        for k in range(0,len(self.dt),chunk):
            a=self.ratio[k:k+chunk];dw=self.equity[k:k+chunk];dt=self.dt[k:k+chunk,None];seg=self.segment[k:k+chunk]
            v=self.floor+(levels[seg,None]-self.floor)*a;sv=np.sqrt(v)
            logstock+=(-v*dt/2+sv*dw).sum(axis=0)
            if derivative:
                terms=a*(-dt/2+dw/(2*sv))
                for j in np.unique(seg):gradient[:,j]+=terms[seg==j].sum(axis=0)
        stock=np.exp(logstock)
        if not np.isfinite(stock).all():raise FloatingPointError('Nonfinite stock')
        return stock,stock[:,None]*gradient


class FlexibleLoss(ExploratoryJointLoss):
    def evaluate(self,prices,forwards,levels,price_jac=None,forward_jac=None):
        result=super().evaluate(prices,forwards,levels)
        if price_jac is not None:
            q=self.quotes;outside=(prices<q.price_bid_U.to_numpy())|(prices>q.price_ask_U.to_numpy())
            pj=self.sqrtweights[:,None]*outside[:,None]*price_jac/q.price_scale_U.to_numpy()[:,None]
            fj=forward_jac/(.5*np.sqrt(len(forwards)))
            cj=np.sqrt(.1)*np.diff(np.eye(len(levels)),axis=0)
            result['jacobian']=np.concatenate([pj,fj,cj])
        return result


class FlexibleFitData:
    def __init__(self,spx,vix,spx_quotes,vix_quotes,forwards,primary,dimension):
        import pandas as pd
        self.spx,self.vix=spx,vix
        self.spx_quotes,self.vix_quotes,self.forwards=spx_quotes,vix_quotes,forwards
        self.all_quotes=pd.concat([spx_quotes,vix_quotes],ignore_index=True)
        self.primary=np.r_[primary,np.ones(len(vix_quotes),bool)]
        self.dimension=dimension;self.loss=FlexibleLoss(self.all_quotes[self.primary],forwards)

    def predict(self,levels,derivative=True):
        n=len(self.all_quotes);d=self.dimension
        prices=np.empty(n);errors=np.empty(n);jac=np.empty((n,d))
        f=np.empty(len(self.forwards));fj=np.empty((len(f),d));fe=np.empty(len(f))
        blocks=np.empty((8,n));fb=np.empty((8,len(f)));archive={}
        for days,data in self.spx.items():
            mean,dv,values=data.price(levels,derivative)
            ix=self.spx_quotes.index[self.spx_quotes.expiry_days==days].to_numpy()
            prices[ix]=mean;jac[ix]=dv;errors[ix]=values.std(axis=0,ddof=1)/np.sqrt(len(values))
            blocks[:,ix]=values.reshape(8,-1,len(ix)).mean(axis=1);archive[f'spx_{days}']=values
        for days,data in self.vix.items():
            mean,dv,b,v=data.price(levels)
            ix=len(self.spx_quotes)+self.vix_quotes.index[self.vix_quotes.expiry_days==days].to_numpy()
            k=int(np.flatnonzero(self.forwards.expiry_days.to_numpy()==days)[0])
            prices[ix]=mean[1:];jac[ix]=dv[1:];errors[ix]=b[:,1:].std(axis=0,ddof=1)/np.sqrt(8)
            f[k]=mean[0];fj[k]=dv[0];fe[k]=b[:,0].std(ddof=1)/np.sqrt(8)
            blocks[:,ix]=b[:,1:];fb[:,k]=b[:,0];archive[f'vix_{days}']=v
        r=self.loss.evaluate(prices[self.primary],f,levels,jac[self.primary] if derivative else None,fj)
        return r|dict(prices=prices,errors=errors,forwards=f,forward_errors=fe,blocks=blocks,forward_blocks=fb,archive=archive)
