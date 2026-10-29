"""Unbiased multilevel estimators; unchanged continuous variance model.

E[P_cont] = E[P_boundary] + E[P_cont-P_boundary]. Both prices use their
own correctly formed VIX variable. This is estimator variance reduction,
not a change to the model nor a mixture of model price laws.
"""
import numpy as np

class MultilevelFit:
    def __init__(self,continuous,boundary_small,boundary_large):
        self.small=continuous;self.bs=boundary_small;self.bl=boundary_large
        for name in ['market','calls','forwards','p','ends','data']:setattr(self,name,getattr(continuous,name))

    def predict(self,x,detail=False):
        c=self.small.predict(x,detail);s=self.bs.predict(x[:11],detail);b=self.bl.predict(x[:11],detail)
        def pad(v):return np.pad(v,[(0,0),(0,1)])
        price=b['prices']+c['prices']-s['prices'];calls=b['calls']+c['calls']-s['calls'];f=b['forwards']+c['forwards']-s['forwards']
        pj=pad(b['native_jacobian'])+c['native_jacobian']-pad(s['native_jacobian']);cj=pad(b['call_jacobian'])+c['call_jacobian']-pad(s['call_jacobian']);fj=pad(b['forward_jacobian'])+c['forward_jacobian']-pad(s['forward_jacobian'])
        m=self.market;distance=np.where(price<m.price_bid_U,price-m.price_bid_U,np.where(price>m.price_ask_U,price-m.price_ask_U,0.));wt=np.concatenate([np.full(len(g),1/(4*len(g))) for _,g in m.groupby('expiry_days',sort=True)]);outside=(price<m.price_bid_U)|(price>m.price_ask_U)
        residual=np.r_[np.sqrt(wt)*distance/m.price_scale_U,(f-self.forwards.forward_vix_points.to_numpy())];jac=np.concatenate([np.sqrt(wt)[:,None]*outside.to_numpy()[:,None]*pj/m.price_scale_U.to_numpy()[:,None],fj])
        out=dict(prices=price,calls=calls,forwards=f,native_jacobian=pj,call_jacobian=cj,forward_jacobian=fj,residuals=residual,jacobian=jac,option_loss=float(residual[:-4]@residual[:-4]),loss=float(residual@residual),raw=c['raw'],moments=c['moments'])
        if detail:
            blocks={}
            for field in ['native','call','forward']:
                terms=[]
                for engine,rr in [(self.bl,b),(self.small,c),(self.bs,s)]:
                    per=[]
                    for day,d in rr['raw'].items():
                        v=d[field];per.append(np.array([v[d['nets']==j].mean(0) for j in range(8)]))
                    terms.append(np.column_stack(per) if field=='forward' else np.concatenate(per,axis=1))
                blocks[field]=terms[0]+terms[1]-terms[2]
            out['blocks']=blocks
        return out
