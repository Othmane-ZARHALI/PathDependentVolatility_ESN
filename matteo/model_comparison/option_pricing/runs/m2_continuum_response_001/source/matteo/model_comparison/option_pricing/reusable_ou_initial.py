"""Exact starting-state translations under FIXED exogenous OU Q dynamics.

No simulation or pricing occurs here. Full initial OU states are required;
three readout projections alone generally do not identify the full state.
Use dated.initial_state if the project's chosen projection-to-state map is intended.
"""
import numpy as np


def state_translation(dated,times,old_initial,new_initial):
    t=np.asarray(times,float)
    if not np.isfinite(t).all() or np.any(t<0):raise ValueError('Nonnegative finite times required')
    old=dated._state(old_initial);new=dated._state(new_initial)
    return np.exp(-t[...,None]*dated.layout.rates)*(new-old)


def projection_translation(dated,times,old_initial,new_initial):
    return state_translation(dated,times,old_initial,new_initial)@dated.projection.T
