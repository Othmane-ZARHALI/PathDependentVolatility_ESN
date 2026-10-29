"""VIX-focused local search constrained in both economic price conventions."""
import numpy as np
from .parity_parameter_core import FitEngine,local_objective


class ConstrainedFitEngine(FitEngine):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);self.partial_key=None;self.partial_value=None
    def vix_prediction(self,p):
        key=tuple(p)
        if key!=self.partial_key:
            self.partial_value=super().vix_prediction(p);self.partial_key=key
        return self.partial_value


def constrained_local(engine,p,anchor_p,anchor):
    local=local_objective(engine,p,anchor_p,anchor);n=len(engine.spx);q=engine.vix
    v=engine.vix_prediction(p)  # Cached at the same readout by FitEngine.
    # Cache the complete partial prediction at driver level to avoid redoing
    # this call: engine.vix_prediction has an optional per-vector memo wrapper.
    native=v['prices'].copy();native_jac=v['jacobian'].copy()
    index=0
    for j,(_,g) in enumerate(q.groupby('expiry_days',sort=True)):
        put=g.option_type.eq('put').to_numpy();m=len(g);native[index:index+m]-=put*(v['forwards'][j]-g.strike.to_numpy())
        native_jac[index:index+m]-=put[:,None]*v['forward_jacobian'][j]
        index+=m
    bid=q.price_bid_U.to_numpy();ask=q.price_ask_U.to_numpy();scale=q.price_scale_U.to_numpy()
    residual=np.where(native<bid,native-bid,np.where(native>ask,native-ask,0))/scale
    residual*=engine.otm_loss.sqrtweights[n:]
    outside=(native<bid)|(native>ask)
    jac=engine.otm_loss.sqrtweights[n:,None]*outside[:,None]*native_jac/scale[:,None]
    native_loss=float(residual@residual);native_gradient=2*residual@jac
    nv=len(q);nf=len(v['forwards']);res=local['residuals'];jres=local['jacobian']
    common_loss=local['components']['VIX'];common_gradient=2*res[n:n+nv]@jres[n:n+nv]
    f_loss=local['components']['forwards'];f_gradient=2*res[n+nv:n+nv+nf]@jres[n+nv:n+nv+nf]
    joint_gradient=2*res@jres
    native_joint=local['total']-common_loss+native_loss
    native_joint_gradient=joint_gradient-common_gradient+native_gradient
    criteria=np.array([local['total'],native_joint,common_loss+f_loss,native_loss+f_loss])
    gradients=np.stack([joint_gradient,native_joint_gradient,common_gradient+f_gradient,native_gradient+f_gradient])
    merit=.5*(common_loss+native_loss)+f_loss
    merit_gradient=.5*(common_gradient+native_gradient)+f_gradient
    return dict(merit=float(merit),merit_gradient=merit_gradient,criteria=criteria,criteria_gradients=gradients,
                canonical_joint=local['total'],native_joint=float(native_joint),canonical_vix=float(common_loss+f_loss),native_vix=float(native_loss+f_loss),
                normalizer=local['normalizer'])


def anchor_criteria(anchor):
    return np.array([anchor['total'],anchor['original_otm_total'],anchor['components']['VIX']+anchor['components']['forwards'],
                     anchor['original_otm_components']['VIX']+anchor['original_otm_components']['forwards']])
