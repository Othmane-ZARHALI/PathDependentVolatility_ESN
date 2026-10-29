"""Explicit all-observation feasibility scope and fixed-control transfer.

Historical source splits are retained in original_split. The new split is an
explicit in-sample label accepted by the existing production objective API.
No processed market artifact or production train/test default is modified.
"""
import numpy as np
import pandas as pd

from matteo.preprocessing.options import black76_implied_vol
from esn_eq.calibration_loss import JointPriceObjective
from .joint_curve_core import ExploratoryJointLoss
from .joint_curve_fit import FitData

SPX_DAYS=(29,57,119,211,302,666,1030)
EXCLUDED_DAYS=(1,1401,1765,2129)
VIX_DAYS=(27,55,118,153)


def populations(quotes,cutoff=3.):
    if cutoff!=3.:raise ValueError('The authorized cutoff is 3')
    q=quotes.copy(deep=True)
    q['original_split']=q['split'];q['original_maturity_role']=q['maturity_role']
    q['split']='in_sample_fit';q['maturity_role']='in_sample_fit'
    eligible=q[(q.instrument=='SPX')&q.expiry_days.isin(SPX_DAYS)].copy()
    if set(eligible.expiry_days)!=set(SPX_DAYS):raise ValueError('Missing retained SPX maturity')
    if eligible.duplicated(['instrument','expiry','strike','option_type']).any():raise ValueError('Duplicate contract')
    pieces=[];anchors=[]
    for days,g in eligible.groupby('expiry_days',sort=True):
        g=g.copy();distance=np.abs(np.log(g.K_over_F))
        a=g.assign(distance=distance).sort_values(['distance','strike','option_type']).iloc[0]
        sigma=black76_implied_vol(float(a.forward),float(a.strike),float(a.T_years),float(a.price_mid_U),a.option_type)
        if not np.isfinite(sigma) or sigma<=0:raise ValueError('Invalid ATM proxy')
        g['atm_proxy']=sigma;g['scaled_moneyness']=distance/(sigma*np.sqrt(g.T_years))
        g['priority_selected']=g.scaled_moneyness<=cutoff
        pieces.append(g)
        anchors.append(dict(expiry_days=int(days),T_years=float(a.T_years),forward=float(a.forward),
            strike=float(a.strike),option_type=a.option_type,atm_proxy=sigma,source_original_split=a.original_split,
            eligible=len(g),selected=int(g.priority_selected.sum())))
    membership=pd.concat(pieces).sort_values(['expiry','strike','option_type']).reset_index(drop=True)
    spx=membership[membership.priority_selected].reset_index(drop=True)
    vix=q[(q.instrument=='VIX')&q.expiry_days.isin(VIX_DAYS)].sort_values(['expiry','strike','option_type']).reset_index(drop=True)
    if set(vix.expiry_days)!=set(VIX_DAYS):raise ValueError('Missing VIX maturity')
    return spx,vix,membership,pd.DataFrame(anchors)


def nearest_same_side(target,donors):
    """Market-key-only mapping; controls and their coefficients remain frozen."""
    indices=[]
    for r in target.itertuples():
        available=np.flatnonzero(donors.option_type.to_numpy()==r.option_type)
        if not len(available):raise ValueError('Missing same-side donor')
        distance=np.abs(np.log(donors.K_over_F.to_numpy()[available]/r.K_over_F))
        indices.append(available[np.argmin(distance)])
    return np.array(indices,dtype=int)


class InSampleLoss(ExploratoryJointLoss):
    def __init__(self,quotes,forwards):
        JointPriceObjective(quotes,forwards,split='in_sample_fit')
        self.quotes=quotes.copy().reset_index(drop=True);self.forwards=forwards.copy().reset_index(drop=True)
        self.weights=np.zeros(len(quotes))
        for instrument in ['SPX','VIX']:
            q=self.quotes[self.quotes.instrument==instrument]
            for _,g in q.groupby('expiry'):self.weights[g.index]=1/(q.expiry.nunique()*len(g))
        self.sqrtweights=np.sqrt(self.weights)


class AllMaturityData(FitData):
    def __init__(self,spx,vix,spx_quotes,vix_quotes,forwards):
        self.spx,self.vix=spx,vix;self.spx_quotes,self.vix_quotes=spx_quotes,vix_quotes;self.forwards=forwards
        self.all_quotes=pd.concat([spx_quotes,vix_quotes],ignore_index=True)
        self.primary=np.ones(len(self.all_quotes),dtype=bool)
        self.loss=InSampleLoss(self.all_quotes,forwards)
