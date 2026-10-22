"""Reflection algebra, represented Gaussian laws and full pricing oracles."""
import unittest
from dataclasses import replace
import numpy as np
from scipy.special import ndtr
from esn_eq.bridge_coupling import joint_root,endpoint_reflector,EndpointBridge
from esn_eq.importance_sampling import DefensiveGaussianMixture
from esn_eq.volatility_path_importance import volatility_path_basis,prepare_volatility_paths
from esn_eq.pricing import ForwardVarianceCurve,coarsen_paths
from esn_eq.states import joint_transition_covariance
from test_pricing import small_pricer,constant_readout


class BridgeCouplingTests(unittest.TestCase):
    def assert_mc(self,x,target):
        error=np.abs(x.mean(axis=0)-target);se=x.std(axis=0,ddof=1)/np.sqrt(len(x))
        self.assertTrue(np.all(error<6*se+2e-12),(error,se))

    def test_reflection_endpoints_and_gaussian_covariance_at_unequal_steps(self):
        p=small_pricer();l1=joint_root(p,.013);l2=joint_root(p,.027);decay=np.exp(-p.layout.rates*.027)
        r,diagnostics=endpoint_reflector(l1,l2,decay)
        np.testing.assert_allclose(r@r.T,np.eye(len(r)),atol=1e-13)
        np.testing.assert_allclose(r@r,np.eye(len(r)),atol=1e-13)
        m=np.c_[np.r_[1.,decay][:,None]*l1,l2]
        np.testing.assert_allclose(m@r,m,atol=1e-13)
        np.testing.assert_allclose(m@m.T,joint_transition_covariance(p.layout,.04),atol=1e-13)
        # Singular roots also require no inversion and preserve represented endpoints.
        zero=np.zeros_like(l1);r,_=endpoint_reflector(zero,zero,decay)
        np.testing.assert_allclose(r@r.T,np.eye(len(r)),atol=1e-13)

    def test_joint_gaussian_black_and_original_sampler_correspondence(self):
        p=small_pricer();times=np.array([0.,.02,.05,.09,.15]);initial=np.linspace(-.1,.2,len(p.layout.rates))
        basis=volatility_path_basis(p,times,.05);shifts=np.zeros((3,4));shifts[1,0]=2.;shifts[2,2]=-2.
        mixture=DefensiveGaussianMixture(shifts,[.5,.2,.3]);bridge=EndpointBridge(p,times,basis)
        x=bridge.draw(initial,100000,612,mixture,checkpoints=[.05,.15]);w=x.original.weights
        self.assertLess(x.endpoint_error,1e-11)
        self.assert_mc(w,1.);self.assert_mc(w[:,None]*x.reflected_scores,np.zeros(4))
        self.assert_mc(w[:,None,None]*x.reflected_scores[:,:,None]*x.reflected_scores[:,None,:],np.eye(4))
        for t in (.05,.15):
            np.testing.assert_allclose(x.original.cache.checkpoint_states[t],x.reflected.checkpoint_states[t],atol=1e-11)
        curve=ForwardVarianceCurve((.05,.15),(.04,.07));params=constant_readout()
        original=p.spx(x.original.cache,params,curve,[.05,.15],order=4).terminal_ratios
        anti=p.spx(x.reflected,params,curve,[.05,.15],order=4).terminal_ratios
        coarse=p.spx(coarsen_paths(x.original.cache,times[::2]),params,curve,[.05,.15],order=4).terminal_ratios
        for t in (.05,.15):
            np.testing.assert_allclose(original[t],anti[t],atol=1e-12)
            np.testing.assert_allclose(original[t],coarse[t],atol=1e-12)
            self.assert_mc(w*anti[t],1.)
            k=1.1;sigma=np.sqrt(curve.integral(0.,t));d1=(-np.log(k)+sigma*sigma/2)/sigma
            self.assert_mc(w*np.maximum(anti[t]-k,0),ndtr(d1)-k*ndtr(d1-sigma))
        for step in (1,3):
            mean,cov=p.gaussian_projection_law(times[step],initial)
            centered=x.reflected.projections[step]-mean
            self.assert_mc(w[:,None]*centered,np.zeros(3))
            self.assert_mc(w[:,None,None]*centered[:,:,None]*centered[:,None,:],cov)
        stochastic=replace(params,cluster_loading=.12,feedback_curvature=.015,spike_curvature=.02)
        for cache in (x.original.cache,x.reflected):
            stock=p.spx(cache,stochastic,curve,[.15],order=12).terminal_ratios[.15]
            self.assert_mc(w*stock,1.)
        # The streaming original reproduces the previous full-path sampler.
        small=bridge.draw(initial,16,28,mixture)
        reference=prepare_volatility_paths(p,times,initial,16,28,basis,mixture)
        np.testing.assert_allclose(small.original.cache.projections,reference.cache.projections,atol=1e-12)
        np.testing.assert_allclose(small.original.weights,reference.weights,atol=1e-13)
        np.testing.assert_allclose(small.original.cache.equity_increments,reference.cache.equity_increments,atol=1e-13)
