import ctypes as ct
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
import pandas as pd
from .full_readout_fit import HERE,Kernel,arr
from .maturity_curve import ENDS_DAYS,embed_three_piece,window_nodes,bounds,validate,CurveKernel,CurveEngine,vix_only_loss
from .maturity_curve_local import CurvePriceLocal


class MaturityCurveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.old=HERE/'runs/m2_joint_local_fit_001';cls.p0=np.array(json.loads((cls.old/'run.json').read_text())['selected']);cls.p=embed_three_piece(cls.p0)
        cls.tmp=tempfile.TemporaryDirectory();cls.lib=Path(cls.tmp.name)/'curve.dylib'
        subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'maturity_curve_stock.cpp'),'-o',str(cls.lib)],check=True)
        cls.k=CurveKernel(cls.old/'convolution.dylib',cls.lib,ENDS_DAYS/365)
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()

    def test_nested_curve_and_hard_bounds(self):
        validate(self.p);self.assertEqual(len(self.p),24)
        lo,hi,_=bounds(16);self.assertTrue(np.all(lo<hi))
        t=np.linspace(0,1030,10001);a=np.minimum(np.searchsorted(ENDS_DAYS,t,side='right'),15);b=np.searchsorted([91.25,273.75],t,side='right')
        np.testing.assert_array_equal(self.p[7:-1][a],self.p0[7:10][b])
        p=self.p.copy();p[9]=p[-1];validate(p);p[9]=np.nextafter(p[-1],0)
        with self.assertRaises(ValueError):validate(p)

    def test_window_integration_splits_every_curve_boundary(self):
        for d in [27,55,118,153]:
            t,w=window_nodes(d/365,ENDS_DAYS/365,1e-4)
            j=np.minimum(np.searchsorted(ENDS_DAYS/365,t,side='right'),15)
            measured=np.bincount(j,weights=w,minlength=16)
            left=np.r_[0,ENDS_DAYS[:-1]];exact=np.maximum(0,np.minimum(ENDS_DAYS,d+30)-np.maximum(left,d))/365
            np.testing.assert_allclose(measured,exact,atol=1e-15)

    def test_stock_kernel_and_all_24_derivatives(self):
        rng=np.random.default_rng(89);grid=np.sort(np.r_[0.,ENDS_DAYS/365]);nt=len(grid)-1;npth=12
        x=arr(rng.normal(0,.2,(nt,npth,3)));dw=arr(rng.normal(0,.1,(nt,npth)));dt=arr(np.diff(grid));seg=np.arange(16,dtype='int32')
        norm=np.zeros((nt,8));norm[:,0]=1.4;norm[:,1:]=rng.normal(0,.02,(nt,7));norm=arr(norm)
        def calc(p):
            nn=norm.copy();nn[:,0]+=norm[:,1:]@(p[:7]-self.p[:7]);out=np.zeros((npth,len(p)+1))
            self.k.curve_lib.stock_curve(nt,npth,16,1,x,dw,dt,seg,nn,arr(p),out);return out
        a=calc(self.p);_,_,sc=bounds(16)
        for j in range(24):
            h=1e-6*sc[j];p=self.p.copy();m=p.copy();p[j]+=h;m[j]-=h
            np.testing.assert_allclose((calc(p)[:,0]-calc(m)[:,0])/(2*h),a[:,j+1],rtol=1e-5,atol=1e-7)

    def test_three_piece_kernel_correspondence_on_actual_cached_paths(self):
        folder=HERE/'runs/m2_parity_adaptive_cache_001/production/d29/batch000';grid=np.load(folder/'grid.npy')
        norm=np.zeros((len(grid)-1,8));norm[:,0]=1.4
        old=Kernel(self.old/'convolution.dylib');new=CurveKernel(self.old/'convolution.dylib',self.lib,[.25,.75,2129/365])
        np.testing.assert_allclose(new.stock(folder,norm,self.p0),old.stock(folder,norm,self.p0),rtol=1e-12,atol=1e-12)
        expanded=self.k.stock(folder,norm,self.p)[0];previous=old.stock(folder,norm,self.p0)[0]
        np.testing.assert_allclose(expanded[:,:8],previous[:,:8],rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(expanded[:,-1],previous[:,-1],rtol=1e-12,atol=1e-12)
        mids=(np.r_[0,ENDS_DAYS[:-1]]+ENDS_DAYS)/2;group=np.searchsorted([91.25,273.75],mids,side='right')
        for j in range(3):np.testing.assert_allclose(expanded[:,8:24][:,group==j].sum(1),previous[:,8+j],rtol=1e-12,atol=1e-12)

    def test_generic_local_anchor_and_all_24_criterion_columns(self):
        q=pd.read_csv(self.old/'prices_canonical_selected.csv');o=pd.read_csv(self.old/'prices_otm_selected.csv');f=pd.read_csv(self.old/'forwards_selected.csv')
        rng=np.random.default_rng(91);a=dict(prices=q.model_price_U.to_numpy(),forwards=f.model_forward.to_numpy(),price_jacobian=rng.normal(size=(len(q),24)),forward_jacobian=rng.normal(size=(4,24)))
        model=CurvePriceLocal(self.p,a,q,o,f);r=model.evaluate(self.p)
        s=json.loads((self.old/'score_selected.json').read_text())
        np.testing.assert_allclose(r['criteria'][:2],[s['canonical_loss'],s['original_otm_loss']],rtol=1e-12)
        _,_,sc=bounds(16)
        for j in range(24):
            h=1e-6*sc[j];p=self.p.copy();m=p.copy();p[j]+=h;m[j]-=h
            fd=(model.evaluate(p)['criteria']-model.evaluate(m)['criteria'])/(2*h)
            np.testing.assert_allclose(fd,r['criteria_gradients'][:,j],rtol=1e-4,atol=1e-4)

    def test_full_prediction_assembly_has_24_columns_and_correct_parity(self):
        from .m2_spx import frozen_stage3
        from esn_eq.pricing import DatedM2Pricer
        cache=HERE/'runs/m2_parity_adaptive_cache_001';ctrl=HERE/'runs/m2_parity_adaptive_controls_001'
        spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv')
        a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial,_=dated.initial_state([1,0,0])
        class StubStock:
            def stock(self,folder,norm,p,derivatives=True):
                out=np.zeros((1,len(np.load(folder/'weights.npy')),len(p)+1));out[:,:,0]=1.;return out
        vc=json.loads((self.old/'vix_controls.json').read_text())
        e=CurveEngine(cache,StubStock(),dated,initial,spx,vix,f,controls=vc,frozen_controls=ctrl,ends=ENDS_DAYS/365,reference_levels=self.p0[7:10])
        e.normalizer=lambda p:(lambda t:np.zeros((len(np.atleast_1d(t)),8)),{})
        def values(day,p,norm):
            with np.load(cache/'production'/f'vix_d{day}.npz') as saved:data={k:saved[k] for k in saved.files}
            return np.full(4096,21.),np.zeros((4096,24)),data
        e.vix_values=values
        try:r=e.predict(self.p,True,True)
        finally:e.close()
        self.assertEqual(r['price_jacobian'].shape,(1511,24));self.assertEqual(r['forward_jacobian'].shape,(4,24))
        self.assertAlmostEqual(sum(r['original_otm_components'].values()),r['original_otm_total'])
        index=len(spx)
        for j,(_,g) in enumerate(vix.groupby('expiry_days',sort=True)):
            expected=np.where(g.option_type.eq('put'),r['forwards'][j]-g.strike.to_numpy(),0)
            np.testing.assert_allclose(r['prices'][index:index+len(g)]-r['native_prices'][index:index+len(g)],expected,atol=1e-11)
            index+=len(g)
        e.spx=spx.iloc[:0];e.original_quotes=vix.reset_index(drop=True)
        from .parity_adaptive import canonical_market
        e.quotes=canonical_market(vix);e.loss=vix_only_loss(e.quotes,f);e.otm_loss=vix_only_loss(vix,f)
        r=e.predict(self.p,True,True)
        self.assertEqual(r['price_jacobian'].shape,(217,24));self.assertEqual(r['components']['SPX'],0)


if __name__=='__main__':unittest.main()
