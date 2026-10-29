"""Exact OU conditional innovation, predictor derivative and payoff oracles."""
from dataclasses import replace
import unittest
import numpy as np
from esn_eq.innovation_controls import InnovationRecorder,spx_innovation_controls,readout_derivatives
from esn_eq.hedging_controls import spx_hedging_controls
from esn_eq.pricing import ForwardVarianceCurve,coarsen_paths
from esn_eq.states import joint_transition_covariance
from test_pricing import small_pricer,constant_readout


class InnovationControlTests(unittest.TestCase):
    def assert_mc(self,x,target):
        error=np.abs(x.mean(axis=0)-target);se=x.std(axis=0,ddof=1)/np.sqrt(len(x))
        self.assertTrue(np.all(error<6*se+1e-12),(error,se))

    def test_exact_full_state_centering_covariance_and_read_only_observer(self):
        p=small_pricer();times=np.array([0.,.03,.1]);initial=np.linspace(-.1,.2,len(p.layout.rates))
        observer=InnovationRecorder(p,times,dict(fine=times,coarse=times[[0,2]]),80000)
        cache=p.prepare(times,initial,80000,63,state_observer=observer)
        for name in ('fine','coarse'):
            recorded=observer.result(name);coarse=coarsen_paths(cache,recorded.times)
            for step,dt in enumerate(np.diff(recorded.times)):
                eps=recorded.values[step];q=joint_transition_covariance(p.layout,dt)
                expected=p.projection@q[1:,1:]@p.projection.T
                np.testing.assert_allclose(recorded.covariance[step],expected,atol=1e-15)
                self.assert_mc(eps,np.zeros(3))
                self.assert_mc(eps[:,:,None]*eps[:,None,:],expected)
                self.assert_mc(eps*coarse.equity_increments[step,:,None],p.projection@q[1:,0])
        # Identical projected initial state can have different full-state next means.
        _,_,v=np.linalg.svd(p.projection,full_matrices=True);null=v[-1]
        np.testing.assert_allclose(p.projection@null,0.,atol=1e-14)
        self.assertGreater(np.linalg.norm(p.projection@(np.exp(-p.layout.rates*.1)*null)),1e-4)
        recordings=[]
        for initial in (np.zeros(len(null)),null):
            obs=InnovationRecorder(p,times,dict(coarse=times[[0,2]]),8)
            observed=p.prepare(times,initial,8,29,state_observer=obs)
            plain=p.prepare(times,initial,8,29)
            np.testing.assert_array_equal(plain.projections,observed.projections)
            np.testing.assert_array_equal(plain.equity_increments,observed.equity_increments)
            recordings.append(obs.result('coarse').values)
        np.testing.assert_allclose(recordings[0],recordings[1],atol=1e-14)

    def test_smooth_readout_derivatives_against_finite_differences(self):
        p=small_pricer();parameters=replace(constant_readout(),cluster_loading=.2,feedback_curvature=.1,spike_curvature=.05)
        x=np.array([[.1,.3,-.2],[-.2,-.4,.5]])
        first,second=readout_derivatives(x,parameters);g=p.multiplier(x,parameters);h=1e-4
        for axis in range(3):
            dx=np.zeros_like(x);dx[:,axis]=h;plus=p.multiplier(x+dx,parameters);minus=p.multiplier(x-dx,parameters)
            np.testing.assert_allclose((np.log(plus)-np.log(minus))/(2*h),first[:,axis],rtol=1e-6,atol=1e-8)
            np.testing.assert_allclose((plus-2*g+minus)/(2*h*h*g),second[:,axis],rtol=1e-5,atol=1e-7)

    def test_martingale_controls_and_exact_existing_stock_scheme(self):
        p=small_pricer();times=np.array([0.,.02,.05,.1]);initial=np.zeros(len(p.layout.rates));paths=60000
        obs=InnovationRecorder(p,times,dict(fine=times),paths);cache=p.prepare(times,initial,paths,81,state_observer=obs)
        params=replace(constant_readout(),cluster_loading=.15,feedback_curvature=.02,spike_curvature=.03)
        curve=ForwardVarianceCurve((.1,),(.04,));maturities=[.1,.1];strikes=[.95,1.05];kinds=['put','call']
        new=spx_innovation_controls(p,cache,obs.result('fine'),params,curve,maturities,strikes,kinds,order=12)
        from unittest.mock import patch
        with patch('esn_eq.innovation_controls.readout_derivatives',side_effect=AssertionError('quadratic path evaluated')):
            linear=spx_innovation_controls(p,cache,obs.result('fine'),params,curve,maturities,strikes,kinds,order=12,linear_only=True)
        np.testing.assert_array_equal(linear.volatility_hedges,new.volatility_hedges[:,:,:3])
        np.testing.assert_array_equal(linear.stock_hedges,new.stock_hedges)
        np.testing.assert_array_equal(linear.terminal_ratios[.1],new.terminal_ratios[.1])
        self.assert_mc(linear.volatility_hedges,np.zeros((2,3)))
        old=spx_hedging_controls(p,cache,params,curve,maturities,strikes,kinds,order=12)
        np.testing.assert_array_equal(new.terminal_ratios[.1],old.terminal_ratios[.1])
        np.testing.assert_array_equal(new.stock_hedges,old.hedges)
        self.assert_mc(new.volatility_hedges,np.zeros((2,9)))
        constant=spx_innovation_controls(p,cache,obs.result('fine'),constant_readout(),curve,maturities,strikes,kinds,order=4)
        np.testing.assert_array_equal(constant.volatility_hedges,0.)

    def test_weighted_proposal_martingales_and_independent_black_oracle(self):
        from scipy.special import ndtr
        from esn_eq.importance_sampling import DefensiveGaussianMixture
        from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths
        p=small_pricer();times=np.array([0.,.02,.05,.1]);initial=np.zeros(len(p.layout.rates));paths=80000
        basis=volatility_path_basis(p,times,.05);shifts=np.zeros((4,4));shifts[1,0]=2.;shifts[2,2]=-2.;shifts[3,3]=2.
        mixture=DefensiveGaussianMixture(shifts,[.5,.15,.2,.15])
        observer=InnovationRecorder(p,times,dict(fine=times),paths)
        sampled=prepare_volatility_paths(p,times,initial,paths,93,basis,mixture,state_observer=observer)
        recorded=observer.result('fine');w=sampled.weights
        for step in range(3):
            self.assert_mc(w[:,None]*recorded.values[step],np.zeros(3))
            self.assert_mc(w[:,None]*recorded.values[step]**2,np.diag(recorded.covariance[step]))
            self.assert_mc(w[:,None]*recorded.values[step]*sampled.cache.equity_increments[step,:,None],recorded.equity_covariance[step])
        parameters=replace(constant_readout(),cluster_loading=.15,feedback_curvature=.02,spike_curvature=.03)
        curve=ForwardVarianceCurve((.1,),(.04,));contracts=([.1,.1],[.95,1.05],['put','call'])
        result=spx_innovation_controls(p,sampled.cache,recorded,parameters,curve,*contracts,order=12)
        self.assert_mc(w[:,None,None]*result.volatility_hedges,np.zeros((2,9)))
        constant=spx_innovation_controls(p,sampled.cache,recorded,constant_readout(),curve,*contracts,order=4)
        stock=constant.terminal_ratios[.1];sigma=np.sqrt(.04*.1)
        for k,kind in zip(contracts[1],contracts[2]):
            d1=(-np.log(k)+sigma**2/2)/sigma;d2=d1-sigma
            exact=ndtr(d1)-k*ndtr(d2) if kind=='call' else k*ndtr(-d2)-ndtr(-d1)
            self.assert_mc(w*np.maximum((1 if kind=='call' else -1)*(stock-k),0.),exact)
