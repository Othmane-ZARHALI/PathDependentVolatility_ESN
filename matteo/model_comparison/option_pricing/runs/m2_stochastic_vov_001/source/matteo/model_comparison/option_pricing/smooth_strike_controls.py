"""C2 deterministic interpolation of fixed known-zero control draws."""
import numpy as np
from scipy.interpolate import CubicSpline
from .continuous_strike_controls import interpolation_weights


def smooth_weights(strikes,anchors):
    strikes,anchors=np.asarray(strikes,float),np.asarray(anchors,float)
    interpolation_weights(strikes,anchors)
    if len(anchors)<3:return interpolation_weights(strikes,anchors)
    spline=CubicSpline(anchors,np.eye(len(anchors)),axis=0,bc_type='natural',extrapolate=False)
    out=spline(np.clip(strikes,anchors[0],anchors[-1]))
    for bound,mask in [(anchors[0],strikes<anchors[0]),(anchors[-1],strikes>anchors[-1])]:
        out[mask]=spline(bound)+(strikes[mask]-bound)[:,None]*spline(bound,1)
    return out


def smooth_controls(quotes,frozen_draws,spx=True):
    frozen_draws=np.asarray(frozen_draws,float);out=np.empty_like(frozen_draws);anchors=[]
    for side,sub in quotes.groupby('option_type'):
        if spx:ix=sorted(set(int((sub.scaled_moneyness-z).abs().idxmin()) for z in [.5,1.5,2.5]),key=lambda j:quotes.strike.iloc[j])
        else:
            indices=sub.sort_values('strike').index.to_numpy();ix=sorted(set(indices[np.round(np.linspace(0,len(indices)-1,3)).astype(int)]))
        k=quotes.strike.iloc[ix].to_numpy();out[:,sub.index]=frozen_draws[:,ix]@smooth_weights(sub.strike.to_numpy(),k).T
        anchors.extend(dict(side=side,strike=float(quotes.strike.iloc[j])) for j in ix)
    return out,anchors
