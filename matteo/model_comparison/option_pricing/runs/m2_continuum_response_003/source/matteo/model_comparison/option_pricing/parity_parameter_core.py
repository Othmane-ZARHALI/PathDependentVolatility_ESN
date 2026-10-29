"""Full-path joint calibration with checked local SPX approximations.

Only a full joint repricing can accept an outer proposal. Inner SPX prices are
explicitly first-order approximations, never described as actual model prices.
"""
import numpy as np
from .parity_adaptive_engine import ParityAdaptiveEngine
from .parity_adaptive import gaussian_controls


class FitEngine(ParityAdaptiveEngine):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.readout_exposures={};self.exposure_readout_key=None

    def vix_values(self,days,p,norm,role='production'):
        readout=tuple(p[:7])
        if readout!=self.exposure_readout_key:
            self.readout_exposures={};self.exposure_readout_key=readout
        key=(days,role)
        if key not in self.vixlaws:
            a=np.load(self.cache/role/f'vix_d{days}.npz');T=days/365
            laws=[self.law(float(t-T),a['states'],32) for t in a['times']]
            self.vixlaws[key]=(a,laws)
        a,laws=self.vixlaws[key]
        if key not in self.readout_exposures:
            expo=np.zeros((len(a['weights']),3));de=np.zeros((len(expo),3,7))
            for t,wt,law in zip(a['times'],a['time_weights'],laws):
                e=self.k.expectation(*law,p);nn=norm(t);m=np.exp(nn[0]);ratio=e[:,0]/m
                s=min(np.searchsorted(self.ends,t,side='right'),2)
                expo[:,s]+=wt*ratio;de[:,s]+=wt*(e[:,1:]/m-ratio[:,None]*nn[1:])
            self.readout_exposures[key]=(expo,de)
        expo,de=self.readout_exposures[key];delta=30/365
        integral=p[10]*delta+expo@(p[7:10]-p[10]);di=np.zeros((len(expo),11))
        di[:,:7]=np.einsum('nsj,s->nj',de,p[7:10]-p[10]);di[:,7:10]=expo;di[:,10]=delta-expo.sum(1)
        v=100*np.sqrt(integral/delta);dv=50/np.sqrt(delta*integral)[:,None]*di
        return v,dv,a

    def predict(self,p,detail=False,derivatives=True):
        r=super().predict(p,detail,derivatives)
        # Complete the component ledger; scalar formula is unchanged.
        r['original_otm_components']['state_penalty']=.1/3
        return r

    def vix_prediction(self,p):
        norm,na=self.normalizer(p);prices=[];jac=[];forwards=[];fj=[]
        for days,q in self.vix.groupby('expiry_days',sort=True):
            v,dv,a=self.vix_values(int(days),p,norm);w=a['weights'];c=gaussian_controls(a['normals'],w)
            fixed=self.controls[str(days)];anchors=np.array(fixed['anchors']);beta=np.array(fixed['beta'])
            b=np.stack([np.interp(np.log(q.strike),np.log(anchors),row) for row in beta])
            y=w[:,None]*np.maximum(v[:,None]-q.strike.to_numpy(),0)-c@b
            prices.extend(y.mean(0));jac.extend(np.einsum('nq,nj->qj',w[:,None]*(v[:,None]>q.strike.to_numpy()),dv)/len(v))
            forwards.append(np.mean(w*v-c@np.array(fixed['forward_beta'])));fj.append(np.mean(w[:,None]*dv,axis=0))
        return dict(prices=np.asarray(prices),jacobian=np.asarray(jac),forwards=np.asarray(forwards),
                    forward_jacobian=np.asarray(fj),normalizer=na)


def local_objective(engine,p,anchor_p,anchor_result):
    """Nonlinear VIX + local SPX at the exact same price/loss scales."""
    n=len(engine.spx);vix=engine.vix_prediction(p)
    spx=anchor_result['prices'][:n]+anchor_result['price_jacobian'][:n]@(p-anchor_p)
    prices=np.r_[spx,vix['prices']];jac=np.concatenate([anchor_result['price_jacobian'][:n],vix['jacobian']])
    r=engine.loss.evaluate(prices,vix['forwards'],p[7:10]);q=engine.quotes
    outside=(prices<q.price_bid_U.to_numpy())|(prices>q.price_ask_U.to_numpy())
    pj=engine.loss.sqrtweights[:,None]*outside[:,None]*jac/q.price_scale_U.to_numpy()[:,None]
    pen=np.zeros((2,11));pen[:,7:10]=np.sqrt(.1)*np.array([[-1,1,0],[0,-1,1]])
    r['jacobian']=np.concatenate([pj,vix['forward_jacobian']/(.5*np.sqrt(len(vix['forwards']))),pen,np.zeros((1,11))])
    r['residuals']=np.r_[r['residuals'],np.sqrt(.1/3)];r['total']+=.1/3;r['components']['state_penalty']=.1/3
    r['normalizer']=vix['normalizer'];return r


def accept_proposal(anchor,trial,predicted_reduction):
    """True joint and VIX improvement in BOTH price conventions; no check data."""
    actual=anchor['total']-trial['total']
    cv=lambda r:r['components']['VIX']+r['components']['forwards']
    nv=lambda r:r['original_otm_components']['VIX']+r['original_otm_components']['forwards']
    reasons=[]
    if predicted_reduction<=0:reasons.append('no predicted joint decrease')
    if actual<=1e-8:reasons.append('no true canonical joint decrease')
    if cv(trial)>=cv(anchor)-1e-8:reasons.append('no true canonical VIX-plus-forward decrease')
    if nv(trial)>=nv(anchor)-1e-8:reasons.append('no native OTM VIX-plus-forward decrease')
    if trial['original_otm_total']>=anchor['original_otm_total']-1e-8:reasons.append('no native OTM joint decrease')
    return dict(accepted=not reasons,reasons=reasons,actual_reduction=float(actual),
                predicted_reduction=float(predicted_reduction),ratio=float(actual/predicted_reduction) if predicted_reduction>0 else None,
                canonical_vix_change=float(cv(trial)-cv(anchor)),native_vix_change=float(nv(trial)-nv(anchor)))
