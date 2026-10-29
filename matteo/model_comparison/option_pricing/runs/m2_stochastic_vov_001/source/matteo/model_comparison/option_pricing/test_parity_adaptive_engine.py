"""End-to-end price/Jacobian correspondence on a deliberately small time grid.

This is an implementation test, not production time-grid convergence evidence.
"""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
import pandas as pd
from .m2_spx import frozen_stage3
from .full_readout_fit import HERE,Kernel
from .parity_adaptive_engine import ParityAdaptiveEngine
from esn_eq.pricing import DatedM2Pricer
from esn_eq.outer_sampling import GaussianOuterSampler

class EngineTests(unittest.TestCase):
    def test_all_parameter_jacobian_with_frozen_controls_and_budgets(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);lib=root/'k.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'controlled_readout_convolution_kernel.cpp'),'-o',str(lib)],check=True)
            k=Kernel(lib);a,_,premium,_=frozen_stage3();p=DatedM2Pricer(a,premium);initial,_=p.initial_state([1,0,0])
            selected=np.array(json.loads((HERE/'runs/m2_controlled_readout_fit_002/selected.json').read_text())['parameters'])
            parent=HERE/'runs/m2_controlled_readout_cache_001'
            spx=pd.read_csv(parent/'spx_contracts.csv');spx=spx[spx.expiry_days==57].iloc[[0,80,140,180,240,297]].reset_index(drop=True)
            vix=pd.read_csv(parent/'vix_contracts.csv');vix=vix[vix.expiry_days==55].iloc[[0,10,20,30,40,54]].reset_index(drop=True)
            forwards=pd.read_csv(parent/'forward_targets.csv');forwards=forwards[forwards.expiry_days==55].reset_index(drop=True)
            grid=np.linspace(0,57/365,12);cache=p.prepare(grid,initial,4096,913);folder=root/'production/d57/batch000';folder.mkdir(parents=True)
            for name,value in [('features',cache.projections),('equity',cache.equity_increments),('grid',grid),('weights',np.ones(4096)),('events',np.array([len(grid)-1],dtype='int32'))]:np.save(folder/f'{name}.npy',value)
            np.savez(folder/'reference.npz',brownian=np.column_stack([cache.equity_increments.sum(0),np.zeros((4096,2))]),duration=[57/365,0,0])
            sampler=GaussianOuterSampler(p,55/365,initial);rng=np.random.default_rng(314)
            for role,n in [('production',4096),('pilot',1024)]:
                d=root/role;d.mkdir(exist_ok=True);z=rng.normal(size=(n,len(initial)));states=sampler.mean+z@sampler.root.T
                np.savez(d/'vix_d55.npz',states=states,normals=z,weights=np.ones(n),nets=np.repeat(np.arange(8),n//8),times=np.array([55/365+15/365]),time_weights=np.array([30/365]))
            # A fixed nonzero zero-mean-control draw correction must be treated
            # as CONSTANT when differentiating subsequent candidate parameters.
            frozen=root/'frozen';dest=frozen/'production/d57';dest.mkdir(parents=True)
            anchors=[float(spx.K_over_F.min()),1.,float(spx.K_over_F.max())]
            (frozen/'coefficients.json').write_text(json.dumps({'57':{'anchors':anchors,'beta':[[.1],[.2],[.15]]}}))
            (frozen/'design.json').write_text(json.dumps({'parameters':selected.tolist()}))
            controls=rng.normal(size=(4096,3,1))
            np.savez(dest/'batch000.npz',controls=controls,output=np.zeros((4096,1)))
            e=ParityAdaptiveEngine(root,k,p,initial,spx,vix,forwards,frozen_controls=frozen)
            try:
                e.fit_vix_controls(selected);result=e.predict(selected)
                direction=np.array([.01,.002,.03,.004,.02,.04,.02,.001,-.001,.002,.00001]);h=1e-4
                plus=e.predict(selected+h*direction);minus=e.predict(selected-h*direction)
                fd=(plus['prices']-minus['prices'])/(2*h)
                # Price Jacobian is recovered from the residual Jacobian away
                # from the bid/ask interval; check the whole residual vector.
                fdres=(plus['residuals']-minus['residuals'])/(2*h)
                relative=np.linalg.norm(fdres-result['jacobian']@direction)/max(1.,np.linalg.norm(fdres))
                np.testing.assert_allclose(fd,result['price_jacobian']@direction,rtol=2e-5,atol=2e-7)
                self.assertLess(relative,.003)
                np.testing.assert_allclose(result['prices'][:len(spx)]-result['native_prices'][:len(spx)],
                    np.where(spx.option_type.eq('put'),spx.forward-spx.strike,0.),atol=1e-12)
                q=vix;canonical=result['prices'][len(spx):];native=result['native_prices'][len(spx):]
                np.testing.assert_allclose(canonical-native,np.where(q.option_type.eq('put'),result['forwards'][0]-q.strike,0),atol=1e-12)
            finally:e.close()

if __name__=='__main__':unittest.main()
