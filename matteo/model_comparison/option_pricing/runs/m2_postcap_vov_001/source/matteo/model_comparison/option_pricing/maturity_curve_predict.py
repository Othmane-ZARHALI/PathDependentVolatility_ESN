"""Generic full SPX/VIX prediction, preserving the canonical-call estimator.

Reference lognormal control curve is fixed at the previous selected model.
This keeps its known target mean while arbitrary new xi pieces are calibrated.
"""
import numpy as np
from .parity_adaptive import reference_calls, iid_statistics, scrambled_statistics, gaussian_controls
from .parity_hedge_controls import correction_at_strikes

class GenericPredict:
    def predict(self,p,detail=False,derivatives=True):
        if self.controls is None:raise ValueError('Independent pilot controls must be frozen first')
        norm,na=self.normalizer(p);prices=[];errors=[];jac=[];forwards=[];fse=[];fj=[]
        native=[];native_se=[];counts=[];effective=[];raw={};allocations=[];moments=[]
        for days,q0 in self.spx.groupby('expiry_days',sort=True):
            q=q0.reset_index(drop=True);batches,w=self.stock[int(days)];terminal=[];brownian=[];frozen_draws=[]
            grid=np.load(batches[0]/'grid.npy');ng=norm(grid[:-1])
            for folder in batches:
                frozen=None
                if self.frozen_controls:
                    frozen=np.load(self.frozen_controls/'production'/f'd{days}'/f'{folder.name}.npz')
                    frozen_draws.append(frozen['controls'])
                value=self.k.stock(folder,ng,p,derivatives=derivatives)[0]
                terminal.append(value)
                a=np.load(folder/'reference.npz');brownian.append(a['brownian']);duration=a['duration']
            term=np.concatenate(terminal);r=term[:,0];dr=term[:,1:];brownian=np.concatenate(brownian)
            ref=np.exp(brownian@np.sqrt(self.reference_levels)-.5*(duration@self.reference_levels));dref=np.zeros_like(dr)
            da=np.zeros(self.npar)
            y,dy=reference_calls(r,dr,ref,dref,w,q.K_over_F.to_numpy(),float(duration@self.reference_levels),da)
            if frozen_draws:
                hc=self.hedge_coefficients[str(days)]
                correction=correction_at_strikes(np.concatenate(frozen_draws),np.array(hc['beta']),hc['anchors'],q.K_over_F.to_numpy())
                y-=correction
            f=q.forward.to_numpy();y*=f;dy*=f[None,:,None];d=self.allocation_distance(q);s=iid_statistics(y,d)
            prices.extend(s['price']);errors.extend(s['se']);jac.extend(np.einsum('nq,nqj->qj',s['coefficients'],dy))
            # Mean-one SPX is an exact identity of the chosen adapted discrete model.
            adjustment=np.where(q.option_type.eq('put'),f-q.strike.to_numpy(),0.)
            native.extend(s['price']-adjustment);native_se.extend(s['se'])
            counts.extend(s['actual_paths']);effective.extend(s['effective_paths'])
            moments.append(dict(instrument='SPX',days=int(days),mass=float(w.mean()),raw_forward_ratio=float(np.mean(w*r)),
                                raw_forward_ratio_se=float(np.std(w*r,ddof=1)/np.sqrt(len(w)))))
            if detail:
                raw[f'spx_{days}']=y;raw[f'spx_raw_call_{days}']=w[:,None]*np.maximum(r[:,None]-q.K_over_F.to_numpy(),0)*f
                raw[f'stock_{days}']=r;raw[f'weights_spx_{days}']=w;raw[f'reference_stock_{days}']=ref
            allocations.append(dict(instrument='SPX',days=int(days),minimum=int(s['actual_paths'].min()),maximum=int(s['actual_paths'].max())))
        for days,q0 in self.vix.groupby('expiry_days',sort=True):
            q=q0.reset_index(drop=True);v,dv,a=self.vix_values(int(days),p,norm);w=a['weights'];c=gaussian_controls(a['normals'],w)
            record=self.controls[str(days)];anchors=np.array(record['anchors']);beta=np.array(record['beta'])
            b=np.stack([np.interp(np.log(q.strike),np.log(anchors),row) for row in beta])
            y=w[:,None]*np.maximum(v[:,None]-q.strike.to_numpy(),0)-c@b
            s=scrambled_statistics(y.reshape(8,512,len(q)),self.allocation_distance(q));prices.extend(s['price']);errors.extend(s['se'])
            dy=w[:,None,None]*(v[:,None]>q.strike.to_numpy())[:,:,None]*dv[:,None,:]
            jac.extend(np.einsum('nq,nqj->qj',s['coefficients'],dy));counts.extend(s['actual_paths']);effective.extend(s['effective_paths'])
            yf=w*v-c@np.array(record['forward_beta']);bf=yf.reshape(8,512).mean(1);forward=bf.mean()
            forwards.append(forward);fse.append(bf.std(ddof=1)/np.sqrt(8));fj.append(np.mean(w[:,None]*dv,axis=0))
            isput=q.option_type.eq('put').to_numpy();native_blocks=s['blocks']-isput*(bf[:,None]-q.strike.to_numpy())
            native.extend(native_blocks.mean(0));native_se.extend(native_blocks.std(0,ddof=1)/np.sqrt(8))
            moments.append(dict(instrument='VIX',days=int(days),market_forward=float(q.forward.iloc[0]),model_forward=float(forward),mass=float(w.mean()),
                                second_moment=float(np.mean(w*v*v))))
            if detail:
                raw[f'vix_{days}']=v;raw[f'vix_values_{days}']=y;raw[f'vix_forward_values_{days}']=yf
                raw[f'vix_blocks_{days}']=s['blocks'];raw[f'vix_native_blocks_{days}']=native_blocks
                raw[f'weights_vix_{days}']=w
            allocations.append(dict(instrument='VIX',days=int(days),minimum=int(s['actual_paths'].min()),maximum=int(s['actual_paths'].max())))
        prices,errors,jac,forwards,fse,fj=map(np.asarray,[prices,errors,jac,forwards,fse,fj])
        r=self.loss.evaluate(prices,forwards,p[7:-1]);outside=(prices<self.quotes.price_bid_U.to_numpy())|(prices>self.quotes.price_ask_U.to_numpy())
        pj=self.loss.sqrtweights[:,None]*outside[:,None]*jac/self.quotes.price_scale_U.to_numpy()[:,None]
        pen=np.zeros((self.ncurve-1,self.npar));pen[:,7:-1]=np.sqrt(.1)*np.diff(np.eye(self.ncurve),axis=0)
        r['jacobian']=np.concatenate([pj,fj/(.5*np.sqrt(len(forwards))),pen,np.zeros((1,self.npar))]);r['residuals']=np.r_[r['residuals'],np.sqrt(.1/3)]
        r['total']+=.1/3;r['components']['state_penalty']=.1/3
        otm=self.otm_loss.evaluate(np.asarray(native),forwards,p[7:-1]);otm['total']+=.1/3;otm['components']['state_penalty']=.1/3
        r.update(price_jacobian=jac,forward_jacobian=fj,prices=prices,errors=errors,forwards=forwards,forward_errors=fse,raw=raw,normalizer=na,
                 native_prices=np.asarray(native),native_errors=np.asarray(native_se),original_otm_total=otm['total'],original_otm_components=otm['components'],
                 actual_paths=np.asarray(counts),effective_paths=np.asarray(effective),allocations=allocations,moments=moments)
        if not derivatives:r.pop('jacobian',None)
        r['derivatives_computed']=derivatives
        return r
