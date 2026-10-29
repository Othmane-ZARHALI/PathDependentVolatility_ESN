import tempfile
from pathlib import Path
import subprocess
from types import SimpleNamespace
import unittest
import numpy as np
import pandas as pd
from .m2_spx import frozen_stage3
from .full_readout_fit import HERE,Kernel,NAMES
from .parity_adaptive_engine import ParityAdaptiveEngine
from .parity_parameter_core import FitEngine,local_objective,accept_proposal
from esn_eq.pricing import DatedM2Pricer

class ParameterFitTests(unittest.TestCase):
    def test_cached_vix_matches_original_value_and_all_derivatives(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);lib=folder/'conv.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'controlled_readout_convolution_kernel.cpp'),'-o',str(lib)],check=True)
            a,params,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial,_=dated.initial_state([1,0,0])
            p=np.array([getattr(params,n) for n in NAMES[:7]]+[.05,.06,.055,.001]);p[4]=.3
            state=dated.prepare([0,.1],initial,16,19,checkpoints=[.1]).checkpoint_states[.1]
            (folder/'production').mkdir();T=55/365;np.savez(folder/'production/vix_d55.npz',states=state,weights=np.ones(16),times=[T+.01,T+.06],time_weights=[.03,30/365-.03])
            e=object.__new__(FitEngine);e.cache=folder;e.k=Kernel(lib);e.dated=dated;e.vixlaws={};e.ends=[.25,.75,2129/365];e.readout_exposures={};e.exposure_readout_key=None
            norm=lambda t:np.zeros(8)
            v,dv,_=e.vix_values(55,p,norm);old,dold,_=ParityAdaptiveEngine.vix_values(e,55,p,norm)
            np.testing.assert_allclose(v,old,rtol=3e-14);np.testing.assert_allclose(dv,dold,rtol=5e-13,atol=3e-10)
            p2=p.copy();p2[7:10]+=[.002,.003,-.001]
            reference,dreference,_=ParityAdaptiveEngine.vix_values(e,55,p2,norm)
            e.k.expectation=lambda *a:(_ for _ in ()).throw(AssertionError('Unnecessary readout quadrature on curve-only move'))
            got,dgot,_=e.vix_values(55,p2,norm);np.testing.assert_allclose(got,reference,rtol=3e-14);np.testing.assert_allclose(dgot,dreference,rtol=5e-13,atol=3e-10)

    def test_local_objective_jacobian_and_no_hidden_SPX_repricing(self):
        q=pd.DataFrame(dict(price_bid_U=[0.]*4,price_ask_U=[0.]*4,price_scale_U=[1.]*4))
        class Loss:
            sqrtweights=np.ones(4)
            def evaluate(self,prices,forwards,levels):
                residual=np.r_[prices,(forwards-10)/.5,np.sqrt(.1)*np.diff(levels)]
                return dict(total=float(residual@residual),residuals=residual,
                            components=dict(SPX=float(np.sum(prices[:2]**2)),VIX=float(np.sum(prices[2:]**2)),forwards=float(np.sum(((forwards-10)/.5)**2)),state_penalty=0.))
        def vix(p):
            values=np.array([p[0]**2+2*p[1],np.exp(p[2])+.2*p[3]]);jac=np.zeros((2,11));jac[0,0]=2*p[0];jac[0,1]=2;jac[1,2]=np.exp(p[2]);jac[1,3]=.2
            fj=np.zeros((1,11));fj[0,4]=1
            return dict(prices=values,jacobian=jac,forwards=np.array([10+p[4]]),forward_jacobian=fj,normalizer={})
        e=SimpleNamespace(spx=[0,1],quotes=q,loss=Loss(),vix_prediction=vix)
        rng=np.random.default_rng(18);anchor_p=np.array([.1,.02,.3,.15,.2,5.,2.,.05,.06,.055,.001]);anchor=dict(prices=np.array([1.,2.,3.,4.]),price_jacobian=rng.normal(size=(4,11)))
        p=anchor_p+.001;result=local_objective(e,p,anchor_p,anchor);h=1e-6
        for j in range(11):
            x=p.copy();y=p.copy();x[j]+=h;y[j]-=h
            fd=(local_objective(e,x,anchor_p,anchor)['residuals']-local_objective(e,y,anchor_p,anchor)['residuals'])/(2*h)
            np.testing.assert_allclose(fd,result['jacobian'][:,j],rtol=2e-7,atol=2e-8)

    def test_acceptance_requires_true_joint_and_both_VIX_conventions(self):
        def r(total,native,vix,native_vix):
            return dict(total=total,original_otm_total=native,components=dict(VIX=vix,forwards=1.),original_otm_components=dict(VIX=native_vix,forwards=1.))
        a=r(100,90,40,35)
        self.assertTrue(accept_proposal(a,r(80,75,30,25),25)['accepted'])
        self.assertFalse(accept_proposal(a,r(110,75,30,25),25)['accepted'])
        self.assertFalse(accept_proposal(a,r(80,75,30,40),25)['accepted'])
        self.assertFalse(accept_proposal(a,r(80,95,30,25),25)['accepted'])
        self.assertFalse(accept_proposal(a,r(80,75,30,25),-1)['accepted'])

if __name__=='__main__':unittest.main()
