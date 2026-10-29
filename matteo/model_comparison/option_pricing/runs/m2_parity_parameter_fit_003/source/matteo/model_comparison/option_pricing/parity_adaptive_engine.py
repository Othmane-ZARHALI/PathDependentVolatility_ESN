"""Calibration-ready canonical-call engine on larger reusable OU samples."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .controlled_readout_fast_fit import ControlledEngine
from .all_maturity_core import InSampleLoss
from .full_readout_fit import arr
from .parity_adaptive import (TIERS, canonical_market, sample_coefficients, iid_statistics,
    scrambled_statistics, reference_calls, gaussian_controls, fit_control_coefficients)
from .plot_joint_curve_fit import implied_vol
from .parity_hedge_controls import correction_at_strikes


def distances(q):
    if 'scaled_moneyness' in q and q.scaled_moneyness.notna().all():
        return q.scaled_moneyness.to_numpy(float)
    i=int(np.argmin(abs(np.log(q.K_over_F.to_numpy()))));r=q.iloc[i]
    sigma=implied_vol(float(r.price_mid_U),float(r.forward),float(r.strike),float(r.T_years),r.option_type)/100
    if not np.isfinite(sigma) or sigma<=0:raise ValueError('Invalid frozen market ATM proxy')
    return abs(np.log(q.K_over_F.to_numpy()))/(sigma*np.sqrt(q.T_years.to_numpy()))


class ParityAdaptiveEngine(ControlledEngine):
    def __init__(self, cache, kernel, dated, initial, spx, vix, forwards, controls=None, frozen_controls=None, stock_values=None, allocation='uniform'):
        self.cache=Path(cache);self.role='production';self.folder=self.cache/self.role
        if allocation not in ['uniform','adaptive']:raise ValueError('Unknown path allocation')
        self.allocation=allocation
        self.k,self.dated,self.initial,self.spx,self.vix,self.forwards=kernel,dated,initial,spx,vix,forwards
        self.original_quotes=pd.concat([spx,vix],ignore_index=True)
        self.quotes=canonical_market(self.original_quotes);self.loss=InSampleLoss(self.quotes,forwards)
        self.otm_loss=InSampleLoss(self.original_quotes,forwards)
        self.ends=[.25,.75,2129/365];self.tau=1/max(dated.layout.rates);self.horizon=float(spx.T_years.max())+30/365
        self.pool=ThreadPoolExecutor(max_workers=4);self.normal_laws={};self.vixlaws={}
        self.last_readout=None;self.normal_cached=None;self.exposure_cached={}
        self.stock={};self.controls=controls
        self.frozen_controls=Path(frozen_controls) if frozen_controls else None
        self.stock_values=stock_values
        if self.frozen_controls:
            self.hedge_coefficients=json.loads((self.frozen_controls/'coefficients.json').read_text())
            self.frozen_parameters=np.array(json.loads((self.frozen_controls/'design.json').read_text())['parameters'])
        for days in sorted(spx.expiry_days.unique()):
            batches=sorted((self.folder/f'd{days}').glob('batch*'))
            self.stock[int(days)]=(batches,np.concatenate([np.load(b/'weights.npy') for b in batches]))

    def allocation_distance(self,q):
        return np.full(len(q),3.) if self.allocation=='uniform' else distances(q)

    def vix_values(self,days,p,norm,role='production'):
        key=(days,role)
        if key not in self.vixlaws:
            a=np.load(self.cache/role/f'vix_d{days}.npz');T=days/365
            laws=[self.law(float(t-T),a['states'],32) for t in a['times']]
            self.vixlaws[key]=(a,laws)
        a,laws=self.vixlaws[key];integral=np.zeros(len(a['weights']));di=np.zeros((len(integral),11))
        for t,wt,law in zip(a['times'],a['time_weights'],laws):
            e=self.k.expectation(*law,p);nn=norm(t);m=np.exp(nn[0]);ratio=e[:,0]/m
            s=min(np.searchsorted(self.ends,t,side='right'),2);xi=p[7+s]
            integral+=wt*(p[10]+(xi-p[10])*ratio)
            di[:,:7]+=wt*(xi-p[10])*(e[:,1:]/m-ratio[:,None]*nn[1:])
            di[:,7+s]+=wt*ratio;di[:,10]+=wt*(1-ratio)
        delta=30/365;v=100*np.sqrt(integral/delta);dv=50/np.sqrt(delta*integral)[:,None]*di
        return v,dv,a

    def fit_vix_controls(self,p):
        norm,_=self.normalizer(p);records={}
        for days,q in self.vix.groupby('expiry_days',sort=True):
            v,_,a=self.vix_values(int(days),p,norm,'pilot');w=a['weights'];c=gaussian_controls(a['normals'],w)
            anchors=np.unique(np.r_[np.geomspace(q.strike.min(),q.strike.max(),7),q.forward.iloc[0]])
            y=w[:,None]*np.maximum(v[:,None]-anchors,0)
            beta=fit_control_coefficients(c,y);bf=fit_control_coefficients(c,(w*v)[:,None])[:,0]
            records[str(days)]=dict(anchors=anchors.tolist(),beta=beta.tolist(),forward_beta=bf.tolist(),pilot_points=len(w))
        self.controls=records;return records

    def predict(self,p,detail=False,derivatives=True):
        if self.controls is None:raise ValueError('Independent pilot controls must be frozen first')
        norm,na=self.normalizer(p);prices=[];errors=[];jac=[];forwards=[];fse=[];fj=[]
        native=[];native_se=[];counts=[];effective=[];raw={};allocations=[];moments=[]
        for days,q0 in self.spx.groupby('expiry_days',sort=True):
            q=q0.reset_index(drop=True);batches,w=self.stock[int(days)];terminal=[];brownian=[];frozen_draws=[]
            for folder in batches:
                grid=np.load(folder/'grid.npy')
                frozen=None
                if self.frozen_controls:
                    frozen=np.load(self.frozen_controls/'production'/f'd{days}'/f'{folder.name}.npz')
                    frozen_draws.append(frozen['controls'])
                if not derivatives and frozen is not None and np.array_equal(p,self.frozen_parameters):
                    value=np.zeros((len(frozen['output']),12));value[:,0]=frozen['output'][:,0]
                elif not derivatives and self.stock_values is not None:
                    value=np.zeros((len(np.load(folder/'weights.npy')),12));value[:,0]=self.stock_values(folder,norm(grid[:-1]),p)
                else:value=self.k.stock(folder,norm(grid[:-1]),p)[0]
                terminal.append(value)
                a=np.load(folder/'reference.npz');brownian.append(a['brownian']);duration=a['duration']
            term=np.concatenate(terminal);r=term[:,0];dr=term[:,1:];brownian=np.concatenate(brownian)
            ref=np.exp(brownian@np.sqrt(p[7:10])-.5*(duration@p[7:10]));dref=np.zeros_like(dr)
            dref[:,7:10]=ref[:,None]*(brownian/(2*np.sqrt(p[7:10]))-.5*duration)
            da=np.zeros(11);da[7:10]=duration
            y,dy=reference_calls(r,dr,ref,dref,w,q.K_over_F.to_numpy(),float(duration@p[7:10]),da)
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
        r=self.loss.evaluate(prices,forwards,p[7:10]);outside=(prices<self.quotes.price_bid_U.to_numpy())|(prices>self.quotes.price_ask_U.to_numpy())
        pj=self.loss.sqrtweights[:,None]*outside[:,None]*jac/self.quotes.price_scale_U.to_numpy()[:,None]
        pen=np.zeros((2,11));pen[:,7:10]=np.sqrt(.1)*np.array([[-1,1,0],[0,-1,1]])
        r['jacobian']=np.concatenate([pj,fj/(.5*np.sqrt(len(forwards))),pen,np.zeros((1,11))]);r['residuals']=np.r_[r['residuals'],np.sqrt(.1/3)]
        r['total']+=.1/3;r['components']['state_penalty']=.1/3
        otm=self.otm_loss.evaluate(np.asarray(native),forwards,p[7:10]);otm['total']+=.1/3
        r.update(price_jacobian=jac,forward_jacobian=fj,prices=prices,errors=errors,forwards=forwards,forward_errors=fse,raw=raw,normalizer=na,
                 native_prices=np.asarray(native),native_errors=np.asarray(native_se),original_otm_total=otm['total'],original_otm_components=otm['components'],
                 actual_paths=np.asarray(counts),effective_paths=np.asarray(effective),allocations=allocations,moments=moments)
        if not derivatives:r.pop('jacobian',None)
        r['derivatives_computed']=derivatives
        return r
